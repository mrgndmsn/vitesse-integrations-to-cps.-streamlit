"""Vitesse reader, conservative rebinning, and v7-compatible CSV writer.

Binary layout/scaling reference: https://github.com/djdt/pewlib/blob/master/src/pewlib/io/nu.py
No isotope deconvolution or background subtraction is performed.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
import math
import re
import shutil
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from zoneinfo import ZoneInfo

import numpy as np
import pandas as pd


class DataError(ValueError):
    pass


def safe_relative(name: str) -> Path:
    name = name.replace('\\', '/')
    p = PurePosixPath(name)
    if p.is_absolute() or '..' in p.parts or not p.parts or ':' in p.parts[0]:
        raise DataError(f'Unsafe upload path: {name}')
    return Path(*p.parts)


def unpack_zip(source, destination: Path, limit=30 * 1024**3):
    """Extract only acquisition data, never executable content or links."""
    allowed = {'.info', '.index', '.integ', '.pulse', '.autob', '.dat', '.method', '.tuning', '.raw'}
    total = 0
    seen = set()
    with zipfile.ZipFile(source) as z:
        for item in z.infolist():
            path = safe_relative(item.filename)
            if item.is_dir() or '__MACOSX' in path.parts or path.name.startswith('.'):
                continue
            if (item.external_attr >> 16) & 0o170000 == 0o120000:
                raise DataError('ZIP symlinks are not supported.')
            if path.suffix.lower() not in allowed:
                continue
            if path in seen:
                raise DataError(f'Duplicate ZIP member: {path}')
            seen.add(path)
            total += item.file_size
            if total > limit:
                raise DataError('Uncompressed upload exceeds the 30 GB limit. Use a local folder instead.')
            target = destination / path
            target.parent.mkdir(parents=True, exist_ok=True)
            with z.open(item) as src, target.open('wb') as dst:
                shutil.copyfileobj(src, dst)


def discover(root: Path) -> dict[str, Path]:
    root = root.expanduser().resolve()
    if not root.is_dir():
        raise DataError('Choose an existing acquisition folder.')
    paths = [root] if (root/'run.info').is_file() else sorted({p.parent for p in root.rglob('run.info')})
    result = {str(p.relative_to(root)) if p != root else root.name: p for p in paths if (p/'integrated.index').is_file()}
    if not result:
        raise DataError('No acquisition found. Include run.info, integrated.index and the .integ files. A .raw waveform file alone is not enough for this exporter.')
    return result


@dataclass
class Acquisition:
    info: dict
    start: np.ndarray
    end: np.ndarray
    counts: np.ndarray
    mz: np.ndarray
    pulses: np.ndarray
    dwell: float
    candidates: list[dict]
    transit: float
    quality: dict


def read_acquisition(root: Path) -> Acquisition:
    info = json.loads((root/'run.info').read_text())
    segments = info['SegmentInfo']
    if len(segments) != 1:
        raise DataError('This version supports one acquisition segment; multi-segment data needs additional validation.')
    seg = segments[0]
    accum = int(info['NumAccumulations1'] * info['NumAccumulations2'])
    period = float(seg['AcquisitionPeriodNs']) * 1e-9
    area = float(info['AverageSingleIonArea'])
    dwell = accum * period
    if not all(math.isfinite(x) and x > 0 for x in [dwell, area, info['SampleRateGsSec']]):
        raise DataError('Invalid timing or single-ion calibration.')
    index = json.loads((root/'integrated.index').read_text())
    chunks, counters, first_centers = [], [], None
    flags, auxiliary, last = 0, 0, None
    result_type = np.dtype([('center','<f4'), ('signal','<f4'), ('aux','<f4'), ('flag','u1')])
    for item in index:
        p = root/f"{int(item['FileNum'])}.integ"
        if not p.is_file():
            raise DataError(f'Missing {p.name}; export stopped rather than joining across missing data.')
        with p.open('rb') as f:
            compressed = f.read(2) == b'\x1f\x8b'
        if compressed:
            with gzip.open(p,'rb') as f:
                payload=f.read()
            head=np.frombuffer(payload[:16],'<u4')
        else:
            with p.open('rb') as f:
                head=np.frombuffer(f.read(16),'<u4')
        if len(head)!=4 or not 1 <= head[3] <= 10000:
            raise DataError(f'Invalid integration header: {p.name}')
        if tuple(head[:3]) != tuple(item[k] for k in ('FirstCycNum','FirstSegNum','FirstAcqNum')):
            raise DataError(f'Index/header mismatch: {p.name}')
        dtype=np.dtype([('cycle','<u4'),('segment','<u4'),('acq','<u4'),('n','<u4'),('results',result_type,int(head[3]))])
        size=len(payload) if compressed else p.stat().st_size
        if size % dtype.itemsize:
            raise DataError(f'Truncated integration file: {p.name}')
        # Bulk reads avoid retaining mappings into a removable drive while UI reruns.
        records=np.frombuffer(payload,dtype=dtype) if compressed else np.fromfile(p,dtype=dtype)
        if not np.all(records['cycle']==1) or not np.all(records['segment']==seg['Num']) or not np.all(records['n']==head[3]):
            raise DataError('Changing cycles/segments/result counts are not supported in this version.')
        acq=records['acq'].astype('int64')
        if np.any(np.diff(acq)!=accum) or (last is not None and acq[0]-last!=accum):
            raise DataError('Missing, duplicated, or out-of-order integrations. Nothing was exported.')
        last=int(acq[-1])
        centers=records['results']['center']
        if first_centers is None:
            first_centers=centers[0].copy()
        if centers.shape[1]!=len(first_centers) or not np.all(centers==first_centers):
            raise DataError('Mass-channel centers change during acquisition; mapping needs validation.')
        flags+=int(np.count_nonzero(records['results']['flag']))
        auxiliary+=int(np.count_nonzero(records['results']['aux']))
        chunks.append(records['results']['signal'].copy())
        counters.append(acq)
    if not chunks:
        raise DataError('The integration index is empty.')
    adc=np.concatenate(chunks)
    if len(adc)!=info['TotalAcquisitions']:
        raise DataError('Integration count does not match run.info.')
    if not np.isfinite(adc).all():
        raise DataError('Non-finite stored signals found; this dataset needs review.')
    end=np.concatenate(counters)*period
    start=end-dwell
    counts=adc.astype('float64')/area
    a,b=info['MassCalCoefficients']
    mz=(a+b*(first_centers.astype(float)/info['SampleRateGsSec']+seg['AcquisitionTriggerDelayNs']))**2
    pulse_list=[]
    if (root/'pulse.index').exists():
        dtype=np.dtype([('cycle','<u4'),('segment','<u4'),('acq','<u4'),('overflow','u1')])
        for item in json.loads((root/'pulse.index').read_text()):
            path=root/f"{int(item['FileNum'])}.pulse"
            if path.stat().st_size%dtype.itemsize:
                raise DataError('Truncated pulse file.')
            p=np.fromfile(path,dtype=dtype)
            if np.any(p['overflow']) or np.any(p['cycle']!=1) or np.any(p['segment']!=seg['Num']):
                raise DataError('Pulse counter overflow or unsupported pulse segment.')
            pulse_list.extend(p['acq'].astype(float)*period)
    candidates=[]
    laser_path=root.parent/'laser.info'
    if laser_path.exists():
        laser=json.loads(laser_path.read_text())
        first=info['FirstLaserLineNumber']
        candidates=[x for x in laser['LaserLineInfo'] if first<=x['ln']<first+laser['AcquisitionLineGroupSize']]
    transit=0.0
    correction=root.parent/'TriggerCorrections.dat'
    if correction.exists():
        c=json.loads(correction.read_text())
        if c['CorrectionMode']!=0:
            raise DataError('Only constant trigger-transit correction is supported.')
        transit=c['Transit1Time']*.001
    autob_files=list(root.glob('*.autob'))
    autob_events=any(p.stat().st_size>25 for p in autob_files)
    return Acquisition(info,start,end,counts,mz,np.array(pulse_list),dwell,candidates,transit,
                       {'flagged_values':flags,'nonzero_auxiliary_values':auxiliary,'autoblank_events_present':autob_events})


def template_lines(text: str) -> tuple[list[str], list[str]]:
    lines=text.lstrip('\ufeff').splitlines()
    i=next((i for i,line in enumerate(lines) if line.strip().startswith('Cycle time (ms)')),None)
    if i is None:
        raise DataError('Template is missing Cycle time (ms).')
    columns=next(csv.reader([lines[i]]))
    if [x.strip() for x in columns[:3]]!=['Cycle time (ms)','x [um]','y [um]']:
        raise DataError('Template must use the v7 time, x and y columns.')
    required=['Timestamp','Laser line number','Laser line name','Laser image name','Starting X','Starting Y','Starting Z','Spot size','Spot spacing','Number of shots','Laser rep rate','Laser power','Laser fluence','Direction of ablation']
    found=[next(csv.reader([line]))[0].strip().rstrip(':').strip() for line in lines[:i]]
    if found!=required or len(columns)<4 or len(set(columns))!=len(columns):
        raise DataError('Template metadata order/columns differ from the supplied v7 layout.')
    return lines[:i+1],columns


def propose_mapping(columns: list[str], mz: np.ndarray) -> pd.DataFrame:
    rows=[]
    for label in columns[3:]:
        match=re.match(r'(\d+)',label.strip())
        mass=float(match[1]) if match else math.nan
        k=int(np.argmin(abs(mz-mass))) if math.isfinite(mass) else -1
        source=f'ch_{k:03d}' if k>=0 and abs(mz[k]-mass)<.25 else 'unmapped'
        rows.append({'CSV column':label,'Source channel':source})
    return pd.DataFrame(rows)


def mapping_indices(mapping: pd.DataFrame, columns: list[str], channels: int) -> np.ndarray:
    if list(mapping['CSV column'])!=columns[3:]:
        raise DataError('Channel mapping must preserve every template column in its original order.')
    result=[]
    for source in mapping['Source channel']:
        if source=='unmapped':
            result.append(-1)
        elif re.fullmatch(r'ch_\d+',str(source)) and 0<=int(source[3:])<channels:
            result.append(int(source[3:]))
        else:
            raise DataError(f'Invalid source channel: {source}')
    return np.array(result,dtype=int)


WINDOW_COLUMNS=['acquisition','name','kind','start_s','end_s','x_um','y_um','z_um']


def empty_windows():
    return pd.DataFrame({k:pd.Series(dtype='float64' if k.endswith(('_s','_um')) else 'str') for k in WINDOW_COLUMNS})


def trigger_windows(acq, key, offset_s=0.0):
    """Associate ordered pulse trains with stationary spots; never invent absent triggers.

    Name suffixes provide an explicitly reported fallback for this acquisition's
    zero shot-count metadata. Ambiguous associations retain data as unassigned.
    """
    import itertools
    if not math.isfinite(offset_s):
        raise DataError('Trigger offset must be finite.')
    spots=acq.candidates
    if not spots or not len(acq.pulses):
        return empty_windows(), ['No usable spot metadata or recorded triggers; all signal retained as unassigned.']
    if any(s.get('lt')!=4 for s in spots):
        raise DataError('Automatic coordinate assignment currently supports stationary spots only.')
    rates=[float(s.get('Metadata',{}).get('RepRate',0)) for s in spots]
    if any(not math.isfinite(r) or r<=0 for r in rates):
        raise DataError('Recorded spot repetition rates are missing or invalid.')
    pulses=acq.pulses
    trains=np.split(pulses,np.flatnonzero(np.diff(pulses)>2.5/min(rates))+1)
    notes=[]
    if len(trains)==len(spots):
        assignment=tuple(range(len(spots)))
    else:
        expected=[]
        for s in spots:
            n=int(s.get('ns',0))
            match=re.search(r'-(\d+)(?:\(\d+\))?$',s['na'])
            expected.append(n if n>0 else int(match[1]) if match else None)
        possible=[ids for ids in itertools.combinations(range(len(spots)),len(trains))
                  if all(expected[i] is not None and len(t) in (expected[i],expected[i]-1)
                         for t,i in zip(trains,ids))]
        if len(possible)!=1:
            return empty_windows(), ['Trigger trains cannot be uniquely linked to recorded spot names; signal retained as unassigned.']
        assignment=possible[0]
        notes.append('Spot association inferred from ordered trigger counts and shot-count name suffixes because raw shot counts are zero; verify labels against the preview.')
    rows=[]
    for train,i in zip(trains,assignment):
        s=spots[i];period=1/rates[i]
        # Cover the recorded train through one repetition period after its last trigger.
        start=float(train[0]+acq.transit+offset_s)
        end=float(train[-1]+period+acq.transit+offset_s)
        lo=max(start,float(acq.start[0]));hi=min(end,float(acq.end[-1]))
        if hi<=lo:
            raise DataError(f"Offset moves {s['na']} outside the acquisition.")
        if lo!=start or hi!=end:
            notes.append(f"{s['na']}: trigger interval clipped to recorded acquisition limits.")
        rows.append([key,s['na'],'signal',lo,hi,s['sx'],s['sy'],s['sz']])
    missing=[s['na'] for i,s in enumerate(spots) if i not in assignment]
    if missing:
        notes.append('No recorded trigger for: '+', '.join(missing)+'. Their signal remains in unassigned/background data; no coordinates were guessed.')
    return pd.DataFrame(rows,columns=WINDOW_COLUMNS),notes


@dataclass
class Binned:
    start: np.ndarray
    end: np.ndarray
    counts: np.ndarray
    n: np.ndarray
    names: list[str]
    kinds: list[str]
    xyz: np.ndarray
    segment_ids: np.ndarray

    @property
    def exposure(self):
        return self.end-self.start

    @property
    def cps(self):
        return self.counts/self.exposure[:,None]


def rebin(acq: Acquisition, windows: pd.DataFrame, target_ms: float, pixel_sum: bool=False) -> Binned:
    if not math.isfinite(target_ms) or target_ms<=0:
        raise DataError('Integration interval must be positive and finite.')
    step=max(1,int(math.floor(target_ms/1000/acq.dwell+.5)))
    middle=(acq.start+acq.end)/2
    segments=[]
    for row in windows.to_dict('records'):
        s,e=float(row['start_s']),float(row['end_s'])
        if not math.isfinite(s+e) or s<acq.start[0]-1e-8 or e>acq.end[-1]+1e-8 or e<=s:
            raise DataError(f"Invalid interval for {row['name']}; use times within this acquisition.")
        if row['kind'] not in ('signal','background') or not str(row['name']).strip():
            raise DataError('Each interval needs a name and a signal/background type.')
        xyz=np.array([row[k] for k in ['x_um','y_um','z_um']],float)
        if not np.isfinite(xyz).all():
            raise DataError('Mapped intervals require finite X, Y and Z coordinates.')
        lo,hi=np.searchsorted(middle,[s,e])
        if lo==hi:
            raise DataError(f"Interval {row['name']} is shorter than one integration.")
        segments.append((int(lo),int(hi),row))
    segments.sort(key=lambda x:x[0])
    filled=[];cursor=0
    for lo,hi,row in segments:
        if lo<cursor:
            raise DataError('Intervals overlap after snapping to native integrations.')
        if lo>cursor:
            filled.append((cursor,lo,None))
        filled.append((lo,hi,row));cursor=hi
    if cursor<len(acq.start):
        filled.append((cursor,len(acq.start),None))
    starts=[];ends=[];counts=[];sizes=[];names=[];kinds=[];positions=[];segment_ids=[]
    for segment_id,(lo,hi,row) in enumerate(filled):
        kind=row['kind'] if row else 'unassigned'
        width=hi-lo if pixel_sum and kind=='signal' else step
        idx=np.arange(lo,hi,width)
        stop=np.minimum(idx+width,hi)
        block=np.add.reduceat(acq.counts[lo:hi],idx-lo,axis=0)
        starts.extend(acq.start[idx]);ends.extend(acq.end[stop-1]);counts.extend(block)
        sizes.extend(stop-idx)
        names.extend([str(row['name']) if row else 'Unassigned / background retained']*len(idx))
        kinds.extend([kind]*len(idx))
        positions.extend([[row['x_um'],row['y_um'],row['z_um']] if row else [np.nan]*3]*len(idx))
        segment_ids.extend([segment_id]*len(idx))
    result=Binned(np.array(starts),np.array(ends),np.array(counts),np.array(sizes),names,kinds,np.array(positions,float),np.array(segment_ids))
    if result.n.sum()!=len(acq.start) or not np.allclose(result.counts.sum(0),acq.counts.sum(0),rtol=1e-10,atol=1e-6):
        raise DataError('Internal count-preservation check failed.')
    return result


def v7_csv(template: str, acq: Acquisition, binned: Binned, indices: np.ndarray,
           line_number: int, image_name: str, timezone='America/Los_Angeles',
           sample_name: str|None=None, endpoint_marker: bool=False) -> str:
    lines,columns=template_lines(template)
    if len(indices)!=len(columns)-3:
        raise DataError('Mapping length differs from the template.')
    stamp=datetime.fromisoformat(acq.info['AnalysisDateTime'])
    if stamp.tzinfo is not None:
        stamp=stamp.astimezone(ZoneInfo(timezone))
    stamp=stamp.replace(tzinfo=None).isoformat(timespec='microseconds')
    candidates=acq.candidates
    first=next((s for s in candidates if s['na']==sample_name),candidates[0] if candidates else {})
    # Unknown coordinates remain explicitly NaN; never synthesize a raster path.
    values={'Timestamp':stamp,'Laser line number':str(line_number),
            'Laser line name':sample_name if sample_name is not None else f"{acq.info.get('SampleName','Acquisition')} group {line_number}",
            'Laser image name':image_name,'Starting X':str(binned.xyz[0,0]),
            'Starting Y':str(binned.xyz[0,1]),'Starting Z':str(binned.xyz[0,2]),
            'Spot size':str(first.get('ss','nan')),'Spot spacing':str(first.get('sp','nan')),
            'Number of shots':str(first.get('ns',0)),
            'Laser rep rate':str(first.get('Metadata',{}).get('RepRate','nan')),
            'Laser power':'nan','Laser fluence':str(first.get('Metadata',{}).get('Fluence','nan'))}
    output=[]
    for line in lines[:-1]:
        prefix,old=line.split(',',1)
        key=prefix.strip().rstrip(':').strip()
        if key in values:
            value=values[key]
            if any(x in value for x in [',','\n','\r']):
                value=value.replace(',',';').replace('\n',' ').replace('\r',' ')
            # Keep original key spacing and surrounding value whitespace.
            left=old[:len(old)-len(old.lstrip())]
            right=old[len(old.rstrip()):]
            output.append(prefix+','+left+value+right)
        else:
            output.append(line)  # Direction field stays exactly as in the supplied v7 template.
    output.append(lines[-1])
    f=io.StringIO();f.write('\n'.join(output)+'\n')
    writer=csv.writer(f,lineterminator='\n')
    cps=binned.cps
    valid=indices>=0
    for i in range(len(binned.start)):
        signals=np.full(len(indices),np.nan)
        signals[valid]=cps[i,indices[valid]]
        row=[binned.start[i]*1000,binned.xyz[i,0],binned.xyz[i,1],*signals]
        writer.writerow(['nan' if not np.isfinite(v) else format(float(v),'.12g') for v in row])
    if endpoint_marker:
        # v7 has no duration column. A missing-value boundary establishes the true
        # sample end even for a single summed pixel, without duplicating counts.
        writer.writerow([format(float(binned.end[-1]*1000),'.12g')]+['nan']*(len(columns)-1))
    return f.getvalue()


def named_parts(bins: Binned, acquisition_key: str):
    boundaries=np.r_[0,np.flatnonzero(np.diff(bins.segment_ids))+1,len(bins.start)]
    for part,(lo,hi) in enumerate(zip(boundaries[:-1],boundaries[1:])):
        kind=bins.kinds[lo]
        name=bins.names[lo]
        if kind=='unassigned':
            name=f'Unassigned background - {acquisition_key} - {part+1}'
        elif kind=='background':
            name=f'Background - {name}'
        piece=Binned(bins.start[lo:hi],bins.end[lo:hi],bins.counts[lo:hi],bins.n[lo:hi],
            bins.names[lo:hi],bins.kinds[lo:hi],bins.xyz[lo:hi],bins.segment_ids[lo:hi])
        yield name,kind,piece


def diagnostics_frame(acq_id: str, bins: Binned) -> pd.DataFrame:
    result=pd.DataFrame({'acquisition':acq_id,'start_s':bins.start,'end_s':bins.end,
        'exposure_s':bins.exposure,'native_integrations':bins.n,'interval':bins.names,'kind':bins.kinds,'segment_id':bins.segment_ids,
        'x_um':bins.xyz[:,0],'y_um':bins.xyz[:,1],'z_um':bins.xyz[:,2]})
    return pd.concat([result,pd.DataFrame(bins.counts,columns=[f'ch_{i:03d}_counts' for i in range(bins.counts.shape[1])])],axis=1)


def build_archives(groups, selected, read_fn, template, mapping, windows, target_ms, pixel_sum,
                   image_name, timezone, output: Path, progress=None, window_fn=None):
    output.mkdir(parents=True,exist_ok=True)
    vit_path=output/'vitesse_export.vit'
    diag_path=output/'vitesse_diagnostics.zip'
    _,columns=template_lines(template)
    manifest={'requested_interval_ms':target_ms,'pixel_sum':pixel_sum,'background_subtracted':False,
              'sample_endpoint_markers':'One final all-NaN signal row per CSV at the exact interval end; not a measurement.',
              'estimated_cps':True,'timezone':timezone,'groups':[],'limitations':[
              'Direct mass-channel estimates; NuQuant isotope processing is not reproduced.',
              'Unmapped isotope columns and unassigned coordinates are NaN.',
              'Review isotope mapping; sharing a channel does not resolve isobaric interference.',
              'Pixel mode sums counts but exports counts/exposure as CPS. Exact counts and exposure are in diagnostics.',
              'Use recorded X/Y on import; direction metadata is retained for v7 compatibility.']}
    try:
        file_counter=0
        with zipfile.ZipFile(vit_path,'w',zipfile.ZIP_DEFLATED) as vz, zipfile.ZipFile(diag_path,'w',zipfile.ZIP_DEFLATED) as dz:
            dz.writestr('channel_mapping.csv',mapping.to_csv(index=False))
            dz.writestr('intervals.csv',windows.to_csv(index=False))
            for i,key in enumerate(selected):
                acq=read_fn(groups[key])
                indices=mapping_indices(mapping,columns,len(acq.mz))
                # Mapping is by stable channel position; reject changed masses across groups.
                if i==0:
                    reference_mz=acq.mz
                elif acq.mz.shape!=reference_mz.shape or not np.allclose(acq.mz,reference_mz,atol=.05,rtol=0):
                    raise DataError(f'Channel masses differ in {key}. Export separately with a reviewed mapping.')
                w,notes=window_fn(acq,key) if window_fn else (windows[windows['acquisition']==key],[])
                dz.writestr(f'group_{i}_intervals.csv',w.to_csv(index=False))
                b=rebin(acq,w,target_ms,pixel_sum)
                files=[]
                for sample_name,kind,piece in named_parts(b,key):
                    content=v7_csv(template,acq,piece,indices,file_counter,image_name,timezone,
                                   sample_name=sample_name,endpoint_marker=True)
                    stamp=re.sub(r'\D','',next(csv.reader([content.splitlines()[0]]))[1])
                    name=f'line_{file_counter}_{stamp}.csv'
                    vz.writestr(name,content)
                    files.append({'file':name,'sample_name':sample_name,'kind':kind,
                        'start_s':float(piece.start[0]),'end_s':float(piece.end[-1]),
                        'measurement_rows':len(piece.start),'endpoint_rows':1})
                    file_counter+=1
                dz.writestr(f'group_{i}_counts.csv',diagnostics_frame(key,b).to_csv(index=False))
                dz.writestr(f'group_{i}_masses.csv',pd.DataFrame({'channel':[f'ch_{j:03d}' for j in range(len(acq.mz))],'calibrated_mz':acq.mz}).to_csv(index=False))
                manifest['groups'].append({'acquisition':key,'alignment_notes':notes,'csv_files':len(files),'samples':files,'native_rows':len(acq.start),'exported_rows':len(b.start),
                    'effective_regular_interval_ms':max(1,int(math.floor(target_ms/1000/acq.dwell+.5)))*acq.dwell*1000,
                    'unassigned_rows':b.kinds.count('unassigned'),'source_quality':acq.quality})
                if progress:
                    progress((i+1)/len(selected))
                del acq,b
            dz.writestr('manifest.json',json.dumps(manifest,indent=2))
        with zipfile.ZipFile(vit_path) as z:
            if z.testzip() is not None or len(z.namelist())!=file_counter:
                raise DataError('Archive verification failed.')
        return vit_path,diag_path,manifest
    except Exception:
        vit_path.unlink(missing_ok=True);diag_path.unlink(missing_ok=True)
        raise

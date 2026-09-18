"""Spatial grid summaries with an explicit summary clock and source provenance."""
from dataclasses import replace
from pathlib import Path
import json
import math
import re
import zipfile
import numpy as np
import pandas as pd
from core import Binned, DataError, mapping_indices, template_lines, v7_csv


class SpatialGrid:
    def __init__(self, width, height, origin_x=0., origin_y=0.):
        if not all(math.isfinite(v) for v in (width,height,origin_x,origin_y)) or min(width,height)<=0:
            raise DataError('Pixel dimensions must be finite and positive.')
        self.width,self.height,self.ox,self.oy=width,height,origin_x,origin_y
        self.cells={};self.sources=[];self.reference=None;self.assigned_counts=None

    def add(self, counts, exposures, xyz, names, acquisition, starts, ends, native_n=None):
        counts=np.asarray(counts,float);exposures=np.asarray(exposures,float);xyz=np.asarray(xyz,float)
        if counts.ndim!=2 or xyz.shape!=(len(counts),3) or exposures.shape!=(len(counts),):
            raise DataError('Spatial contributions have inconsistent shapes.')
        if not np.isfinite(xyz).all() or not np.isfinite(counts).all() or not np.isfinite(exposures).all() or np.any(exposures<=0):
            raise DataError('Spatial contributions require known coordinates and positive exposure.')
        native_n=np.ones(len(counts),int) if native_n is None else native_n
        for j in range(len(counts)):
            ix=math.floor((xyz[j,0]-self.ox)/self.width)
            iy=math.floor((xyz[j,1]-self.oy)/self.height)
            key=(ix,iy)
            if key not in self.cells:
                self.cells[key]={'counts':np.zeros(counts.shape[1]),'exposure':0.,'n':0,'names':[], 'z':set()}
            cell=self.cells[key];cell['counts']+=counts[j];cell['exposure']+=float(exposures[j]);cell['n']+=int(native_n[j]);cell['z'].add(float(xyz[j,2]))
            if names[j] not in cell['names']:cell['names'].append(names[j])
            self.sources.append({'pixel_x':ix,'pixel_y':iy,'acquisition':acquisition,'spot':names[j],
                'source_start_s':float(starts[j]),'source_end_s':float(ends[j]),'exposure_s':float(exposures[j]),
                'native_integrations':int(native_n[j]),'source_x_um':float(xyz[j,0]),'source_y_um':float(xyz[j,1]),'source_z_um':float(xyz[j,2])})
        total=counts.sum(0)
        self.assigned_counts=total if self.assigned_counts is None else self.assigned_counts+total

    def add_acquisition(self, acq, key, windows):
        if self.reference is None:
            self.reference=replace(acq,counts=np.empty((0,len(acq.mz))),start=np.empty(0),end=np.empty(0),candidates=list(acq.candidates))
        elif not np.allclose(acq.mz,self.reference.mz,atol=.05,rtol=0):
            raise DataError('Mass channels changed between spatial contributions.')
        middle=(acq.start+acq.end)/2
        for row in windows.to_dict('records'):
            if row['kind']!='signal':continue
            lo,hi=np.searchsorted(middle,[row['start_s'],row['end_s']])
            if hi<=lo:continue
            # These raw files supply one stationary position per interval. Grid.add
            # also accepts individually positioned measurements when available.
            self.add(acq.counts[lo:hi].sum(0)[None,:],[(hi-lo)*acq.dwell],
                [[row['x_um'],row['y_um'],row['z_um']]],[row['name']],key,
                [acq.start[lo]],[acq.end[hi-1]],[hi-lo])
        self.reference.candidates.extend(s for s in acq.candidates if s not in self.reference.candidates)

    def table(self, channel=0):
        return pd.DataFrame([{'pixel_x':ix,'pixel_y':iy,'x_um':self.ox+(ix+.5)*self.width,
            'y_um':self.oy+(iy+.5)*self.height,'counts':c['counts'][channel],
            'exposure_s':c['exposure'],'CPS':c['counts'][channel]/c['exposure'],
            'spots':'; '.join(c['names'])} for (ix,iy),c in sorted(self.cells.items())])

    def export(self, template, mapping, timezone, output):
        if not self.cells:raise DataError('No measurements have linked X/Y coordinates; no spatial pixels can be exported.')
        summed=np.sum([c['counts'] for c in self.cells.values()],axis=0)
        if not np.allclose(summed,self.assigned_counts,rtol=1e-10,atol=1e-6):
            raise DataError('Spatial count conservation failed.')
        output=Path(output);output.mkdir(parents=True,exist_ok=True)
        path=output/'vitesse_spatial.vit';diag=output/'vitesse_spatial_details.zip'
        _,columns=template_lines(template)
        indices=mapping_indices(mapping,columns,len(self.reference.mz))
        clock=0.;records=[]
        with zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED) as z:
            for line,((ix,iy),cell) in enumerate(sorted(self.cells.items())):
                x=self.ox+(ix+.5)*self.width;y=self.oy+(iy+.5)*self.height
                duration=cell['exposure'];sample='; '.join(cell['names'])
                zpos=next(iter(cell['z'])) if len(cell['z'])==1 else np.nan
                bins=Binned(np.array([clock]),np.array([clock+duration]),cell['counts'][None,:],
                    np.array([cell['n']]),[sample],['signal'],np.array([[x,y,zpos]]),np.array([0]))
                matching=[s for s in self.reference.candidates if s['na'] in cell['names']]
                meta={'na':sample}
                for field in ['ss','sp','ns']:
                    vals={s.get(field) for s in matching}
                    meta[field]=next(iter(vals)) if len(vals)==1 else 'nan'
                meta['Metadata']={}
                for field in ['RepRate','Fluence']:
                    vals={s.get('Metadata',{}).get(field) for s in matching}
                    meta['Metadata'][field]=next(iter(vals)) if len(vals)==1 else 'nan'
                ref=replace(self.reference,candidates=[meta])
                content=v7_csv(template,ref,bins,indices,line,'Spatial grid summary',timezone,sample_name=sample,endpoint_marker=True)
                stamp=re.sub(r'\D','',content.splitlines()[0].split(',',1)[1])
                filename=f'line_{line}_{stamp}.csv';z.writestr(filename,content)
                records.append({'file':filename,'pixel_x':ix,'pixel_y':iy,'x_um':x,'y_um':y,
                    'sample_name':sample,'summary_start_s':clock,'exposure_s':duration,'native_integrations':cell['n'],
                    **{f'ch_{k:03d}_counts':v for k,v in enumerate(cell['counts'])}})
                clock+=duration
        manifest={'width_um':self.width,'height_um':self.height,'origin_x_um':self.ox,'origin_y_um':self.oy,
            'pixels':len(records),'background_subtracted':False,'units_in_v7':'CPS = summed counts / summed exposure',
            'time_axis':'Synthetic cumulative exposure clock for spatial summaries, NOT original acquisition time.',
            'coordinate_rule':'Half-open grid cells; rows use cell centers. Counts are assigned by recorded position, not beam-footprint overlap.',
            'scope':'All linked signal contributions across selected acquisitions; unassigned/background times remain in companion timeline .vit.',
            'warning':'Import spatial and timeline packages into separate sessions to avoid double-counting.'}
        with zipfile.ZipFile(diag,'w',zipfile.ZIP_DEFLATED) as z:
            z.writestr('pixels_counts.csv',pd.DataFrame(records).to_csv(index=False))
            z.writestr('pixel_contributions.csv',pd.DataFrame(self.sources).to_csv(index=False))
            z.writestr('manifest.json',json.dumps(manifest,indent=2))
        return path,diag,manifest

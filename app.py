from __future__ import annotations

import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
import uuid

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from core import (DataError, WINDOW_COLUMNS, build_archives, discover, empty_windows,
                  mapping_indices, propose_mapping, read_acquisition, rebin,
                  safe_relative, template_lines, unpack_zip)

HERE=Path(__file__).parent
st.set_page_config(page_title='Vitesse · Time & Pixels',page_icon='🔬',layout='wide')
st.markdown('''<style>
.block-container{max-width:1450px;padding-top:2rem}
h1{letter-spacing:-.045em} div[data-testid="stMetric"]{background:white;padding:16px;border-radius:12px;border:1px solid #dfe7ec}
</style>''',unsafe_allow_html=True)


@st.cache_resource(max_entries=1,show_spinner='Reading continuous integrations…')
def cached_read(path,signature):
    return read_acquisition(Path(path))


def load_group(path):
    files=[path/'run.info',path/'integrated.index',*path.glob('*.integ'),*path.glob('*.pulse')]
    signature=tuple((str(p),p.stat().st_size,p.stat().st_mtime_ns) for p in files)
    return cached_read(str(path),signature)


def activate(path):
    groups=discover(Path(path))
    st.session_state.groups=groups
    st.session_state.windows=empty_windows()
    st.session_state.pop('mapping',None)
    st.session_state.pop('exports',None)
    st.session_state.dataset_id=uuid.uuid4().hex
    st.session_state.windows_revision=0
    st.session_state.mapping_revision=0


def extrema(x,y,limit=18000):
    if len(x)<=limit:
        return x,y
    width=math.ceil(len(x)/(limit//2))
    keep=[]
    for lo in range(0,len(x),width):
        block=y[lo:lo+width]
        keep.extend(sorted({lo+int(np.argmin(block)),lo+int(np.argmax(block))}))
    return x[keep],y[keep]


if 'temp_root' not in st.session_state:
    st.session_state.temp_root=tempfile.mkdtemp(prefix='vitesse-app-')
    st.session_state.dataset_id='empty'
    st.session_state.windows_revision=0
    st.session_state.mapping_revision=0

# Keep in-progress sessions usable after a source update.
st.session_state.setdefault('windows_revision',0)
st.session_state.setdefault('mapping_revision',0)

st.title('Vitesse · Time & Pixels')
st.caption('Keep the complete acquisition. Choose the time resolution. Export in your NICE v7 layout.')

with st.sidebar:
    st.header('1 · Load acquisition')
    source=st.radio('Source',['Local folder','Upload folder','Upload ZIP'])
    try:
        if source=='Local folder':
            path=st.text_input('Acquisition or results folder',value=os.environ.get('VITESSE_DEFAULT_DATA',''))
            st.caption('Best for large datasets. Files stay on this computer and are read in place.')
            if st.button('Load folder',type='primary'):
                activate(path)
        elif source=='Upload folder':
            uploaded=st.file_uploader('Choose the original results folder',accept_multiple_files='directory',
                type=['info','index','integ','pulse','autob','dat','method','tuning','raw'])
            if st.button('Load uploaded folder',disabled=not uploaded):
                dest=Path(st.session_state.temp_root)/uuid.uuid4().hex
                dest.mkdir()
                seen=set()
                for f in uploaded:
                    rel=safe_relative(f.name)
                    if rel in seen:
                        raise DataError(f'Duplicate uploaded path: {rel}')
                    seen.add(rel)
                    target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True)
                    with target.open('wb') as out:
                        f.seek(0);shutil.copyfileobj(f,out)
                activate(dest)
        else:
            uploaded=st.file_uploader('ZIP of the original results folder',type=['zip'])
            if st.button('Load ZIP',disabled=uploaded is None):
                dest=Path(st.session_state.temp_root)/uuid.uuid4().hex
                dest.mkdir()
                with st.spinner('Unpacking acquisition…'):
                    uploaded.seek(0);unpack_zip(uploaded,dest)
                    activate(dest)
    except (DataError,OSError,ValueError,KeyError) as exc:
        st.error(str(exc))
    st.divider()
    st.caption('Required: run.info, integrated.index and .integ files. Include laser.info and .pulse files for timing and coordinates. No background subtraction is applied.')

if 'groups' not in st.session_state:
    st.info('Load a results folder or ZIP to start. Your original files are never modified.')
    st.markdown('**Workflow:** inspect the signal → assign time intervals and coordinates → adjust integration → download `.vit`.')
    st.stop()

groups=st.session_state.groups
with st.sidebar:
    key=st.selectbox('Preview acquisition',list(groups),key=f'group_{st.session_state.dataset_id}')
    st.header('2 · Output resolution')
    target_ms=st.number_input('Target integration interval (ms)',min_value=.001,value=10.,step=1.,format='%.4f')
    mode=st.radio('Signal output',['Time-resolved bins','Sum each mapped signal interval into one pixel'])
    pixel_sum=mode.startswith('Sum')
    st.caption('Background and unassigned intervals always retain time-resolved bins. Pixel CPS = summed counts ÷ total exposure.')

try:
    acq=load_group(groups[key])
except Exception as exc:
    st.error(f'Could not read this acquisition: {exc}')
    st.stop()

effective=max(1,int(math.floor(target_ms/1000/acq.dwell+.5)))*acq.dwell*1000
a,b,c,d=st.columns(4)
a.metric('Recorded integrations',f'{len(acq.start):,}')
b.metric('Recorded duration',f'{acq.end[-1]-acq.start[0]:.3f} s')
c.metric('Native interval',f'{acq.dwell*1000:.4f} ms')
d.metric('Regular output interval',f'{effective:.4f} ms')
st.caption('Intervals use whole native integrations. Short bins at interval boundaries and the end keep their actual exposure.')
st.warning('Signals are estimated CPS from stored mass channels. NuQuant isotope corrections are not yet reproduced; review channel mapping before export.',icon='⚠️')
if any(acq.quality.values()):
    st.warning(f'Recorded quality information needs review: {acq.quality}. Stored signals are retained without applying new blanking corrections.')

template=st.session_state.get('template',(HERE/'v7_template.csv').read_text())
try:
    header_lines,columns=template_lines(template)
except DataError as exc:
    st.error(str(exc));st.stop()
mapping_key=hashlib.sha256((template+str(len(acq.mz))).encode()).hexdigest()
if st.session_state.get('mapping_key')!=mapping_key:
    st.session_state.mapping=propose_mapping(columns,acq.mz)
    st.session_state.mapping_key=mapping_key
    st.session_state.mapping_revision+=1

view,intervals,channels,export=st.tabs(['Signal preview','Intervals & coordinates','v7 columns & mapping','Export'])

with intervals:
    st.subheader('Assign time intervals to fixed coordinates')
    st.write('Use seconds from this acquisition’s start. Unassigned time is preserved automatically. Overlapping intervals are rejected; intervals snap to native integration boundaries.')
    with st.expander('Recorded spot metadata and pulse times',expanded=False):
        if acq.candidates:
            st.dataframe(pd.DataFrame([{'name':x['na'],'x_um':x['sx'],'y_um':x['sy'],'z_um':x['sz'],
              'line':x['ln'],'stored_shot_count':x['ns']} for x in acq.candidates]),hide_index=True)
        if len(acq.pulses):
            pulse_groups=np.split(acq.pulses,np.flatnonzero(np.diff(acq.pulses)>.3)+1)
            st.dataframe(pd.DataFrame([{'pulse_train':i+1,'recorded_pulses':len(g),'first_trigger_s':g[0],
               'last_trigger_s':g[-1],'configured_transit_s':acq.transit} for i,g in enumerate(pulse_groups)]),hide_index=True)
        st.caption('Pulse trains are not automatically assigned to spot names. Trigger counts can be incomplete; confirm boundaries from the trace or laser log.')
    with st.form('add_interval'):
        opts=['Enter coordinates']+[x['na'] for x in acq.candidates]
        coord_source=st.selectbox('Copy coordinates from recorded spot',opts)
        name=st.text_input('Spot / interval name',value='',help='Leave blank to use the recorded spot name when copying its coordinates.')
        c1,c2,c3=st.columns(3)
        kind=c1.selectbox('Interval type',['signal','background'])
        start_s=c2.number_input('Start (s)',value=max(0.,float(acq.start[0])),format='%.6f')
        end_s=c3.number_input('End (s)',value=float(acq.end[-1]),format='%.6f')
        x,y,z=st.columns(3)
        xval=x.number_input('X (µm)',value=0.,format='%.3f')
        yval=y.number_input('Y (µm)',value=0.,format='%.3f')
        zval=z.number_input('Z (µm)',value=0.,format='%.3f')
        st.caption('Choosing a recorded spot uses its X/Y/Z instead of the manual fields.')
        if st.form_submit_button('Add interval'):
            try:
                if coord_source!='Enter coordinates':
                    spot=next(x for x in acq.candidates if x['na']==coord_source)
                    xval,yval,zval=spot['sx'],spot['sy'],spot['sz']
                    name=name.strip() or spot['na']
                row={'acquisition':key,'name':name,'kind':kind,'start_s':start_s,'end_s':end_s,'x_um':xval,'y_um':yval,'z_um':zval}
                proposed=pd.concat([st.session_state.windows,pd.DataFrame([row])],ignore_index=True)
                rebin(acq,proposed[proposed.acquisition==key],target_ms,pixel_sum)
                st.session_state.windows=proposed
                st.session_state.windows_revision+=1
                st.rerun()
            except (ValueError,DataError) as exc:
                st.error(str(exc))
    st.caption('Edit or delete existing intervals below, then apply. Background coordinates can be assigned explicitly too.')
    edited=st.data_editor(st.session_state.windows,num_rows='dynamic',hide_index=True,key=f'windows_editor_{st.session_state.dataset_id}_{st.session_state.windows_revision}',
       column_config={'acquisition':st.column_config.SelectboxColumn(options=list(groups),required=True),
                      'kind':st.column_config.SelectboxColumn(options=['signal','background'],required=True)},width='stretch')
    if st.button('Apply interval edits'):
        try:
            if not set(edited.acquisition).issubset(groups):
                raise DataError('An interval refers to an unknown acquisition.')
            for group_id in edited.acquisition.unique():
                source_acq=acq if group_id==key else load_group(groups[group_id])
                rebin(source_acq,edited[edited.acquisition==group_id],target_ms,pixel_sum)
            st.session_state.windows=edited.copy();st.session_state.windows_revision+=1;st.rerun()
        except Exception as exc:
            st.error(str(exc))
    csv_upload=st.file_uploader('Import interval / coordinate table',type='csv')
    if st.button('Apply uploaded intervals',disabled=csv_upload is None):
        try:
            new=pd.read_csv(csv_upload)
            if list(new.columns)!=WINDOW_COLUMNS or not set(new.acquisition).issubset(groups):
                raise DataError('Use the downloaded table’s column names and acquisition IDs.')
            for group_id in new.acquisition.unique():
                source_acq=acq if group_id==key else load_group(groups[group_id])
                rebin(source_acq,new[new.acquisition==group_id],target_ms,pixel_sum)
            st.session_state.windows=new;st.session_state.windows_revision+=1;st.rerun()
        except Exception as exc:
            st.error(str(exc))
    st.download_button('Download interval table',st.session_state.windows.to_csv(index=False),'intervals.csv','text/csv')

windows=st.session_state.windows
try:
    bins=rebin(acq,windows[windows.acquisition==key],target_ms,pixel_sum)
except Exception as exc:
    st.error(f'Check intervals: {exc}');st.stop()

with view:
    st.subheader('Recorded signal and proposed export')
    choice=st.multiselect('Mass channels to plot',range(len(acq.mz)),default=[int(np.argmin(abs(acq.mz-177.92)))],
       format_func=lambda i:f'ch_{i:03d} · m/z {acq.mz[i]:.4f}',key=f'plot_channels_{len(acq.mz)}')
    show_native=st.checkbox('Overlay native integrations',value=False)
    show_pulses=st.checkbox('Show recorded laser pulses',value=True)
    fig=go.Figure()
    mid=(bins.start+bins.end)/2
    for k in choice:
        px,py=extrema(mid,bins.cps[:,k])
        fig.add_trace(go.Scattergl(x=px,y=py,name=f'ch_{k:03d} output',mode='lines+markers' if pixel_sum else 'lines',line=dict(width=1.4)))
        if show_native:
            nx,ny=extrema((acq.start+acq.end)/2,acq.counts[:,k]/acq.dwell)
            fig.add_trace(go.Scattergl(x=nx,y=ny,name=f'ch_{k:03d} native',opacity=.3,line=dict(width=.7)))
    for row in windows[windows.acquisition==key].to_dict('records'):
        fig.add_vrect(x0=row['start_s'],x1=row['end_s'],fillcolor='#10b981' if row['kind']=='signal' else '#94a3b8',
                     opacity=.1,line_width=0,annotation_text=row['name'])
    if show_pulses and len(acq.pulses):
        fig.add_trace(go.Scattergl(x=acq.pulses+acq.transit,y=np.full(len(acq.pulses),.97),yaxis='y2',
              mode='markers',marker=dict(symbol='line-ns',size=8,color='#c89524'),name='Triggers + configured transit'))
    fig.update_layout(height=480,template='plotly_white',margin=dict(l=20,r=20,t=20,b=20),hovermode='x unified',
       xaxis_title='Time from acquisition origin (s)',yaxis_title='Estimated CPS',
       yaxis2=dict(overlaying='y',range=[0,1],visible=False),legend=dict(orientation='h'),uirevision=key)
    st.plotly_chart(fig,width='stretch',config={'scrollZoom':True})
    st.caption('Drag to zoom. Display uses peak-preserving downsampling for long traces; export always includes all integrations. Yellow ticks use the saved transit delay, not a fitted alignment.')
    st.write(f'**{len(bins.start):,} output rows** · {bins.n.sum():,} native integrations retained · {bins.kinds.count("unassigned"):,} rows with unassigned coordinates')
    if choice:
        sample=pd.DataFrame({'start_s':bins.start,'end_s':bins.end,'exposure_s':bins.exposure,'type':bins.kinds,
                            'X':bins.xyz[:,0],'Y':bins.xyz[:,1]})
        for k in choice:
            sample[f'ch_{k:03d} counts']=bins.counts[:,k]
            sample[f'ch_{k:03d} CPS']=bins.cps[:,k]
        st.dataframe(sample.head(100),hide_index=True,width='stretch')

with channels:
    st.subheader('Exact v7 layout, reviewed channel mapping')
    st.write(f'The supplied template has {len(columns)-3} isotope columns. This acquisition has {len(acq.mz)} stored mass channels. The original header order is preserved, including unmapped columns.')
    template_upload=st.file_uploader('Optional replacement v7 CSV template',type='csv',key='template_upload')
    if st.button('Use uploaded v7 template',disabled=template_upload is None):
        try:
            text=template_upload.getvalue().decode('utf-8-sig');template_lines(text)
            st.session_state.template=text;st.rerun()
        except Exception as exc:
            st.error(str(exc))
    st.caption('Initial suggestions use nominal mass proximity only. They do not separate isotopes at the same mass or apply NuQuant corrections. Unmapped values export as nan, never zero.')
    with st.expander('Original metadata and header'):
        st.code('\n'.join(header_lines),language=None)
    chmap=st.data_editor(st.session_state.mapping,hide_index=True,disabled=['CSV column'],
       column_config={'Source channel':st.column_config.SelectboxColumn(options=['unmapped']+[f'ch_{i:03d}' for i in range(len(acq.mz))],required=True)},
       key=f'mapping_editor_{mapping_key}_{st.session_state.mapping_revision}',height=400,width='stretch')
    if st.button('Apply channel mapping'):
        try:
            mapping_indices(chmap,columns,len(acq.mz));st.session_state.mapping=chmap.copy();st.session_state.mapping_revision+=1;st.rerun()
        except DataError as exc:
            st.error(str(exc))
    with st.expander('Stored channel masses'):
        st.dataframe(pd.DataFrame({'Source channel':[f'ch_{i:03d}' for i in range(len(acq.mz))],'calibrated_mz':acq.mz}),hide_index=True)
    st.download_button('Download channel mapping',st.session_state.mapping.to_csv(index=False),'channel_mapping.csv','text/csv')
    mapfile=st.file_uploader('Restore channel mapping CSV',type='csv',key='map_upload')
    if st.button('Apply uploaded mapping',disabled=mapfile is None):
        try:
            new=pd.read_csv(mapfile);mapping_indices(new,columns,len(acq.mz))
            st.session_state.mapping=new;st.session_state.mapping_key=mapping_key
            st.session_state.mapping_revision+=1;st.rerun()
        except Exception as exc:
            st.error(str(exc))

with export:
    st.subheader('Build the iolite package')
    selected=st.multiselect('Acquisitions to export',list(groups),default=[key])
    image_name=st.text_input('Laser image name',value=acq.info.get('SampleName','Vitesse export'))
    timezone=st.text_input('Timestamp timezone',value='America/Los_Angeles')
    st.write('The `.vit` ZIP contains one v7-layout CSV per named spot interval, plus separate CSVs for background and unassigned time. Spot names become iolite Samples via the Laser line name field. A separate diagnostics ZIP contains counts, exposure, intervals and processing settings.')
    named_count=int(((windows.kind=='signal') & windows.acquisition.isin(selected)).sum())
    st.caption(f'{named_count} named spot intervals assigned for this export. Set their names and time boundaries in Intervals & coordinates; gaps remain separate background samples.')
    if not named_count:
        st.warning('No named spot intervals are assigned yet. Add intervals to get individual spot Samples in iolite; this export would contain only unassigned background samples.')
    st.info('No additional background subtraction. Every recorded integration contributes exactly once. Unknown coordinates and unmapped isotope columns remain nan.')
    with st.expander('Iolite import notes',expanded=False):
        st.write('Tested in iolite 4.10.12: regular and pixel-summed exports imported with signal values, time coordinates and X/Y intact. Iolite appends a NaN separator row.')
        st.warning('Your existing Vitesse “Treat channels as already background subtracted” default was left unchanged. It can label these raw-background signals as background-subtracted. Review that setting before baseline-dependent reduction; CSV contents cannot override it.')
        st.write('Iolite may ask for channel dwell times. Exact exposures are in diagnostics. Variable-duration pixels cannot be described by one global dwell value.')
    if pixel_sum:
        st.caption('Pixel mode has one measurement row per mapped signal interval, with background time bins retained. Each CSV also has a final NaN boundary row so iolite knows the full sample duration. Counts are not duplicated; exact counts and exposures are in diagnostics.')
    review=st.checkbox('I have reviewed the channel mapping and accept estimated CPS, including any shared channels or unmapped values.')
    coordinates=st.checkbox('I have reviewed interval coordinates; any unassigned X/Y values may remain nan.')
    settings={'groups':selected,'template':template,'mapping':st.session_state.mapping.to_csv(index=False),
       'windows':windows.to_csv(index=False),'ms':target_ms,'pixel_sum':pixel_sum,'image_name':image_name,'timezone':timezone,
       'dataset':st.session_state.dataset_id,
       'source_versions':[(k,[(p.name,p.stat().st_size,p.stat().st_mtime_ns) for p in sorted(groups[k].iterdir()) if p.is_file()]) for k in selected]}
    signature=hashlib.sha256(json.dumps(settings,sort_keys=True).encode()).hexdigest()
    if st.button('Build .vit and diagnostics',type='primary',disabled=not(selected and review and coordinates)):
        try:
            from zoneinfo import ZoneInfo
            ZoneInfo(timezone)
            progress=st.progress(0.,text='Exporting complete timelines…')
            out=Path(st.session_state.temp_root)/'exports'/uuid.uuid4().hex
            with st.spinner('Writing CSVs and verifying counts…'):
                result=build_archives(groups,selected,read_acquisition,template,st.session_state.mapping,windows,
                    target_ms,pixel_sum,image_name,timezone,out,lambda value:progress.progress(value))
            st.session_state.exports=(signature,result)
            progress.empty()
        except Exception as exc:
            st.error(f'Export stopped: {exc}')
    saved=st.session_state.get('exports')
    if saved and saved[0]==signature:
        vit,diag,manifest=saved[1]
        st.success(f'Package ready: {sum(g["csv_files"] for g in manifest["groups"])} CSV files. All count-preservation checks passed.')
        one,two,three=st.columns(3)
        with vit.open('rb') as f:
            one.download_button('Download .vit',f,'vitesse_export.vit','application/zip')
        with vit.open('rb') as f:
            two.download_button('Download identical .zip',f,'vitesse_export.zip','application/zip')
        with diag.open('rb') as f:
            three.download_button('Download diagnostics',f,'vitesse_diagnostics.zip','application/zip')
        st.dataframe(pd.DataFrame([s for g in manifest['groups'] for s in g['samples']]),hide_index=True)
    elif saved:
        st.caption('Settings changed. Rebuild to download an export with the current settings.')
    st.caption('Use iolite’s Vitesse text importer. Review its background-subtraction flag and dwell-time settings separately; the app never changes those defaults. Full import checks and scientific limits are in the README.')

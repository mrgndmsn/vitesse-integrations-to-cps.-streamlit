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
                  safe_relative, template_lines, unpack_zip, trigger_windows)

from spatial import SpatialGrid

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
    source=st.radio('Source',['Upload ZIP','Upload folder','Local folder (server only)'])
    try:
        if source=='Local folder (server only)':
            path=st.text_input('Acquisition or results folder',value=os.environ.get('VITESSE_DEFAULT_DATA',''))
            st.caption('This path must exist on the server running Streamlit. A cloud app cannot read your Mac or external drive; use Upload ZIP instead.')
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
    st.caption('Cloud uploads leave your computer and are processed on the hosting server. Large runs may exceed the host’s memory limits.')
    st.caption('Required: run.info, integrated.index and .integ files. Include laser.info and .pulse files for timing and coordinates. No background subtraction is applied.')

if 'groups' not in st.session_state:
    st.info('Load a results folder or ZIP to start. Your original files are never modified.')
    st.markdown('**Workflow:** preview → align triggers → choose integration → export `.vit`.')
    st.stop()

groups=st.session_state.groups
key=st.selectbox('Preview acquisition',list(groups),key=f'group_{st.session_state.dataset_id}')
try:
    acq=load_group(groups[key])
except Exception as exc:
    st.error(f'Could not read this acquisition: {exc}');st.stop()

left,right=st.columns(2)
offset=left.number_input('Trigger time offset (s) — all acquisitions',value=0.0,step=.01,format='%.4f',
    help='Added to the saved transit correction. Negative moves triggers earlier; positive moves them later. Signal timestamps and coordinates stay unchanged.')
target_ms=right.number_input('Output integration (ms)',min_value=.001,value=10.,step=1.,format='%.4f')
st.caption(f'Negative offset ← earlier · positive offset → later. Saved transit: {acq.transit:.4f} s; total shift here: {acq.transit+offset:+.4f} s. The same added offset applies to every exported acquisition.')
choice=st.multiselect('Signals to preview',range(len(acq.mz)),default=[int(np.argmin(abs(acq.mz-177.92)))],
    format_func=lambda i:f'm/z {acq.mz[i]:.4f}',key=f'plot_channels_{len(acq.mz)}')
spatial_mode=st.checkbox('Bin counts into an X/Y pixel grid',value=False)
width=height=5.0
origin_x=origin_y=0.0
if spatial_mode:
    c1,c2=st.columns(2)
    width=c1.number_input('Pixel width (µm)',min_value=.001,value=5.,step=1.)
    height=c2.number_input('Pixel height (µm)',min_value=.001,value=5.,step=1.)
    with st.expander('Grid origin',expanded=False):
        origin_x=st.number_input('Grid origin X (µm)',value=0.,step=1.)
        origin_y=st.number_input('Grid origin Y (µm)',value=0.,step=1.)
    st.caption('All linked measurements in the same X/Y cell are combined, including revisits and different spots. Counts are summed before dividing by summed exposure for the v7 CPS columns.')
try:
    windows,notes=trigger_windows(acq,key,offset)
    bins=rebin(acq,windows,target_ms,False)
except (DataError,ValueError) as exc:
    st.error(str(exc));st.stop()
fig=go.Figure()
for k in choice:
    px,py=extrema((bins.start+bins.end)/2,bins.cps[:,k])
    fig.add_trace(go.Scattergl(x=px,y=py,name=f'm/z {acq.mz[k]:.4f}',mode='lines'))
for row in windows.to_dict('records'):
    fig.add_vrect(x0=row['start_s'],x1=row['end_s'],fillcolor='#10b981',opacity=.10,line_width=0,annotation_text=row['name'])
if len(acq.pulses):
    fig.add_trace(go.Scattergl(x=acq.pulses+acq.transit+offset,y=np.full(len(acq.pulses),.97),yaxis='y2',
        mode='markers',marker=dict(symbol='line-ns',size=10,color='#c89524'),name='Shifted laser triggers'))
fig.update_layout(height=510,template='plotly_white',margin=dict(l=20,r=20,t=35,b=20),hovermode='x unified',
    xaxis_title='Time from acquisition origin (s)',yaxis_title='Estimated CPS',
    yaxis2=dict(overlaying='y',range=[0,1],visible=False),legend=dict(orientation='h'),uirevision=key)
st.plotly_chart(fig,width='stretch',config={'scrollZoom':True})
st.caption(f'{len(acq.start):,} recorded integrations → {len(bins.start):,} output measurements. Background is retained. Drag to zoom; changing offset keeps your zoom.')
for note in notes:
    st.warning(note)
if spatial_mode:
    preview_grid=SpatialGrid(width,height,origin_x,origin_y)
    preview_grid.add_acquisition(acq,key,windows)
    st.subheader('Spatial preview · current acquisition')
    if preview_grid.cells and choice:
        table=preview_grid.table(choice[0])
        spatial_fig=go.Figure()
        for row in table.to_dict('records'):
            spatial_fig.add_shape(type='rect',x0=row['x_um']-width/2,x1=row['x_um']+width/2,
                y0=row['y_um']-height/2,y1=row['y_um']+height/2,line=dict(color='#cbd5e1',width=1))
        spatial_fig.add_trace(go.Scatter(x=table.x_um,y=table.y_um,mode='markers',text=table.spots,
            customdata=table[['counts','exposure_s']].to_numpy(),
            marker=dict(symbol='square',size=14,color=table.counts,colorscale='Viridis',showscale=True,colorbar=dict(title='Counts')),
            hovertemplate='%{text}<br>X %{x} µm · Y %{y} µm<br>Counts %{customdata[0]}<br>Exposure %{customdata[1]} s<extra></extra>'))
        spatial_fig.update_layout(height=330,template='plotly_white',xaxis_title='X (µm)',yaxis_title='Y (µm)',yaxis=dict(scaleanchor='x',scaleratio=1))
        st.plotly_chart(spatial_fig,width='stretch')
        st.caption(f'{len(table)} occupied cells at {width:g} × {height:g} µm. Export merges contributions across all selected acquisitions. Blank cells contain no assigned measurements.')
    else:
        st.info('No located measurements to show in this acquisition.')
    st.caption('Stationary holes retain one recorded position. Smaller grid cells do not create additional spatial resolution within a hole.')
    st.info('Spatial .vit uses a summary exposure clock, not acquisition time. A separate complete timeline .vit preserves all background and original times. Import these into separate iolite sessions.')


# The supplied v7 template stays built in. Routine users do not edit its columns.
template=(HERE/'v7_template.csv').read_text()
_,columns=template_lines(template)
mapping_key=hashlib.sha256((template+str(len(acq.mz))).encode()).hexdigest()
if st.session_state.get('mapping_key')!=mapping_key:
    st.session_state.mapping=propose_mapping(columns,acq.mz)
    st.session_state.mapping_key=mapping_key
with st.expander('Details & advanced settings',expanded=False):
    st.caption('The supplied v7 metadata layout and column order are built in. Spot names, X/Y/Z, spot size, repetition rate and fluence come from the raw metadata. Time without a reliable spot association retains NaN coordinates.')
    st.dataframe(windows[['name','start_s','end_s','x_um','y_um','z_um']],hide_index=True)
    st.caption('Trigger windows extend from the first recorded trigger to one repetition period after the last. Missing first-shot triggers are not reconstructed. Check the shaded regions against the signal, especially washout tails.')
    st.warning('CPS are estimates from stored mass channels; NuQuant isotope corrections are not reproduced. Unmapped isotope columns remain NaN. Channel assignments below are nominal-mass suggestions.')
    edited=st.data_editor(st.session_state.mapping,hide_index=True,disabled=['CSV column'],
       column_config={'Source channel':st.column_config.SelectboxColumn(options=['unmapped']+[f'ch_{i:03d}' for i in range(len(acq.mz))],required=True)},
       key=f'mapping_{mapping_key}',height=260)
    mapping_indices(edited,columns,len(acq.mz))
    st.session_state.mapping=edited.copy()
    timezone=st.text_input('Timestamp timezone',value='America/Los_Angeles')
    st.caption('Your iolite background-subtracted default remains unchanged. Exact exposures are in diagnostics; a single dwell value cannot represent variable-duration summed spots.')

st.divider()
all_groups=st.checkbox('Export all loaded acquisitions',value=True)
selected=list(groups) if all_groups else [key]
st.caption(f'{len(selected)} acquisition(s). One CSV per linked spot, plus separate background/unassigned CSVs, packaged as .vit. Recorded coordinates are fixed within each stationary spot.')
st.caption('Estimated CPS · no additional background subtraction · original v7 layout')
signature=hashlib.sha256(json.dumps({'dataset':st.session_state.dataset_id,'groups':selected,'offset':offset,
    'ms':target_ms,'grid':[spatial_mode,width,height,origin_x,origin_y],'mapping':edited.to_csv(index=False),'timezone':timezone},sort_keys=True).encode()).hexdigest()
if st.button('Export .vit',type='primary'):
    try:
        from zoneinfo import ZoneInfo
        ZoneInfo(timezone)
        progress=st.progress(0.,text='Exporting acquisitions…')
        output=Path(st.session_state.temp_root)/'exports'/uuid.uuid4().hex
        grid=SpatialGrid(width,height,origin_x,origin_y) if spatial_mode else None
        def prepare_windows(a,k):
            w,notes=trigger_windows(a,k,offset)
            if grid is not None:
                grid.add_acquisition(a,k,w)
            return w,notes
        result=build_archives(groups,selected,read_acquisition,template,edited,empty_windows(),target_ms,False,
            acq.info.get('LaserImageName',''),timezone,output,lambda v:progress.progress(v),
            window_fn=prepare_windows)
        vit,diag,manifest=result
        manifest['additional_trigger_offset_s']=offset
        # Export settings accompany the diagnostic archive, never the v7 CSV headers.
        import zipfile
        with zipfile.ZipFile(diag,'a',zipfile.ZIP_DEFLATED) as z:
            z.writestr('alignment_settings.json',json.dumps({'additional_trigger_offset_s':offset,'uses_saved_transit':True}))
        spatial_result=grid.export(template,edited,timezone,output) if grid is not None else None
        st.session_state.exports=(signature,result,spatial_result)
        progress.empty()
    except Exception as exc:
        st.error(f'Export stopped: {exc}')
saved=st.session_state.get('exports')
if saved and saved[0]==signature:
    vit,diag,manifest=saved[1]
    spatial_result=saved[2] if len(saved)>2 else None
    if spatial_result:
        spatial_vit,spatial_diag,spatial_manifest=spatial_result
        st.success(f'{spatial_manifest["pixels"]} spatial pixels ready; counts conserved across all assigned measurements.')
        with spatial_vit.open('rb') as f:
            st.download_button('Download spatial .vit',f,'vitesse_spatial.vit','application/zip',type='primary')
        with spatial_diag.open('rb') as f:
            st.download_button('Download pixel counts & contributions',f,'vitesse_spatial_details.zip','application/zip')
    st.success(f'Ready: {sum(g["csv_files"] for g in manifest["groups"])} CSVs. All recorded integrations retained.')
    with vit.open('rb') as f:
        st.download_button('Download complete timeline .vit' if spatial_mode else 'Download .vit',f,'vitesse_export.vit','application/zip',type='primary')
    with st.expander('Export details'):
        with diag.open('rb') as f:
            st.download_button('Download counts & diagnostics',f,'vitesse_diagnostics.zip','application/zip')
        for g in manifest['groups']:
            for note in g['alignment_notes']:
                st.write(f'{g["acquisition"]}: {note}')
elif saved:
    st.caption('Settings changed. Export again to update the file.')

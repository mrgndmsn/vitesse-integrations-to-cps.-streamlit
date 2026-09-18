import csv
import io
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
import pytest

from core import (Acquisition, DataError, WINDOW_COLUMNS, build_archives, empty_windows,
                  mapping_indices, propose_mapping, read_acquisition, rebin,
                  template_lines, unpack_zip, v7_csv)

HERE=Path(__file__).parent


def sample(n=11):
    end=np.arange(1,n+1)*.001
    return Acquisition({'AnalysisDateTime':'2026-07-02T17:16:19.286342-06:00','SampleName':'Test'},
       end-.001,end,np.arange(n*2,dtype=float).reshape(n,2)+1,np.array([27.9724,177.921]),
       np.array([.004]),.001,[],.03,{})


def windows(*rows):
    return pd.DataFrame(rows,columns=WINDOW_COLUMNS)


def test_partial_bins_and_background_preserve_all_counts():
    a=sample()
    w=windows(['g','pixel','signal',.003,.007,12.,13.,14.])
    b=rebin(a,w,2,True)
    np.testing.assert_array_equal(b.n,[2,1,4,2,2])
    np.testing.assert_allclose(b.counts.sum(0),a.counts.sum(0))
    np.testing.assert_allclose((b.cps*b.exposure[:,None]).sum(0),a.counts.sum(0))
    assert b.kinds==['unassigned','unassigned','signal','unassigned','unassigned']
    np.testing.assert_array_equal(b.xyz[2],[12,13,14])
    assert np.isnan(b.xyz[0]).all()


def test_disjoint_same_pixel_does_not_swallow_background():
    a=sample()
    w=windows(['g','same','signal',0,.003,1,2,3],['g','same','signal',.007,.011,1,2,3])
    b=rebin(a,w,2,True)
    assert b.kinds==['signal','unassigned','unassigned','signal']
    np.testing.assert_array_equal(b.n,[3,2,2,4])


@pytest.mark.parametrize('w',[
    windows(['g','a','signal',0,.006,1,2,3],['g','b','signal',.005,.008,1,2,3]),
    windows(['g','a','signal',0,.02,1,2,3]),
    windows(['g','a','signal',.008,.006,1,2,3]),
    windows(['g','a','signal',0,.006,float('nan'),2,3]),
])
def test_invalid_windows_rejected(w):
    with pytest.raises(DataError):
        rebin(sample(),w,2)


def test_native_lower_bound_and_short_final_bin():
    a=sample()
    assert len(rebin(a,empty_windows(),.01).start)==11
    b=rebin(a,empty_windows(),3)
    np.testing.assert_array_equal(b.n,[3,3,3,2])
    np.testing.assert_allclose(b.exposure,[.003,.003,.003,.002])


def test_v7_header_and_column_order_exact():
    template=(HERE/'v7_template.csv').read_text()
    lines,cols=template_lines(template)
    a=sample();b=rebin(a,empty_windows(),2)
    m=propose_mapping(cols,a.mz);idx=mapping_indices(m,cols,2)
    exported=v7_csv(template,a,b,idx,0,'image')
    output=exported.splitlines()
    assert output[14]==lines[14]
    assert [x.split(',')[0] for x in output[:14]]==[x.split(',')[0] for x in lines[:14]]
    assert output[0]=='Timestamp:,2026-07-02T16:16:19.286342'
    assert output[13]==lines[13]
    rows=list(csv.reader(io.StringIO(exported)))
    assert all(len(r)==len(cols) for r in rows[15:])
    assert len(cols)==275
    assert rows[15][1:3]==['nan','nan']
    assert rows[15][cols.index('238U cps')]=='nan'
    assert float(rows[15][cols.index('28Si cps')])==2000


def test_zip_rejects_traversal(tmp_path):
    source=io.BytesIO()
    with zipfile.ZipFile(source,'w') as z:
        z.writestr('../run.info','{}')
    source.seek(0)
    with pytest.raises(DataError):
        unpack_zip(source,tmp_path)
    assert not (tmp_path.parent/'run.info').exists()


def test_two_group_package_has_only_v7_csvs(tmp_path):
    a=sample();template=(HERE/'v7_template.csv').read_text();_,cols=template_lines(template)
    mapping=propose_mapping(cols,a.mz)
    vit,diag,manifest=build_archives({'g1':Path('g1'),'g2':Path('g2')},['g1','g2'],lambda _:a,
         template,mapping,empty_windows(),2,False,'test','America/Los_Angeles',tmp_path)
    with zipfile.ZipFile(vit) as z:
        assert len(z.namelist())==2
        assert all('/' not in name and name.endswith('.csv') for name in z.namelist())
        assert all(template_lines(z.read(name).decode())[1]==cols for name in z.namelist())
    with zipfile.ZipFile(diag) as z:
        frame=pd.read_csv(z.open('group_0_counts.csv'))
        assert frame.native_integrations.sum()==11
        assert json.loads(z.read('manifest.json'))['background_subtracted'] is False


def test_named_spots_and_background_get_separate_samples(tmp_path):
    a=sample();template=(HERE/'v7_template.csv').read_text();_,cols=template_lines(template)
    mapping=propose_mapping(cols,a.mz)
    w=windows(['g','915001-1','signal',.003,.007,11937.,43284.,1281.])
    vit,diag,manifest=build_archives({'g':Path('g')},['g'],lambda _:a,template,mapping,w,2,True,'test',
                                    'America/Los_Angeles',tmp_path)
    with zipfile.ZipFile(vit) as z:
        assert len(z.namelist())==3
        parsed=[list(csv.reader(io.StringIO(z.read(name).decode()))) for name in z.namelist()]
        names=[rows[2][1].strip() for rows in parsed]
        assert names[1]=='915001-1'
        assert names[0].startswith('Unassigned background') and names[2].startswith('Unassigned background')
        signal=parsed[1]
        assert len(signal[15:])==2  # one real pixel plus one missing-value endpoint
        assert float(signal[15][0])==3 and float(signal[16][0])==7
        assert all(value=='nan' for value in signal[16][1:])
        assert float(signal[15][cols.index('28Si cps')])==10000
        assert all(rows[14]==cols for rows in parsed)
    assert manifest['groups'][0]['exported_rows']==5


def test_reader_rejects_missing_chunk(tmp_path):
    (tmp_path/'run.info').write_text(json.dumps({'SegmentInfo':[{'AcquisitionPeriodNs':1000,'Num':1}],
        'NumAccumulations1':1,'NumAccumulations2':1,'AverageSingleIonArea':1,'SampleRateGsSec':2}))
    (tmp_path/'integrated.index').write_text('[{"FileNum":1}]')
    with pytest.raises(DataError,match='Missing'):
        read_acquisition(tmp_path)


@pytest.mark.parametrize('compressed',[False,True])
def test_uploaded_nested_zip_decodes_and_preserves_native_counts(tmp_path,compressed):
    import gzip
    import struct
    source=io.BytesIO()
    info={'SegmentInfo':[{'AcquisitionPeriodNs':1000,'Num':1,'AcquisitionTriggerDelayNs':0}],
      'NumAccumulations1':2,'NumAccumulations2':1,'AverageSingleIonArea':10,
      'SampleRateGsSec':2,'TotalAcquisitions':3,'MassCalCoefficients':[0,1],
      'AnalysisDateTime':'2026-07-02T17:00:00-06:00','SampleName':'Test'}
    binary=b''.join(struct.pack('<IIIIfffB',1,1,acq,1,10,signal,0,0) for acq,signal in [(2,10),(4,20),(6,30)])
    if compressed:
        binary=gzip.compress(binary)
    with zipfile.ZipFile(source,'w') as z:
        z.writestr('Image001/00001/run.info',json.dumps(info))
        z.writestr('Image001/00001/integrated.index','[{"FileNum":1,"FirstCycNum":1,"FirstSegNum":1,"FirstAcqNum":2}]')
        z.writestr('Image001/00001/1.integ',binary)
    source.seek(0);unpack_zip(source,tmp_path)
    a=read_acquisition(tmp_path/'Image001/00001')
    np.testing.assert_allclose(a.counts[:,0],[1,2,3])
    np.testing.assert_allclose(a.end,[.000002,.000004,.000006])
    np.testing.assert_allclose(rebin(a,empty_windows(),1).counts[:,0],[6])

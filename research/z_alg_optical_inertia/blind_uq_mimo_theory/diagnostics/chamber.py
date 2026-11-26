"""Measured empty-repeat lower scale for square subsampled MIMO layouts.
No antenna impedances or unrecorded reflection entries are inferred.
"""
import csv,json,os,sys,time
from pathlib import Path
from zipfile import ZipFile
import numpy as np
P=Path(__file__).resolve().parents[1]; historical=P.parent/'validation/manitoba'
sys.path.insert(0,str(historical if historical.exists() else P/'vendor/chamber'))
import physical as chamber
out=P/'results'; start=time.perf_counter()
with ZipFile(os.environ.get('MANITOBA_TOMOGRAPHIC_ZIP',chamber.RAW)) as z: data={key:chamber.load_file(z,val)[0] for key,val in chamber.FILE_MAP.items()}
rows=[]
for q in [6,8,12,16]:
    tx=np.floor(np.arange(q)*24/q).astype(int);rx=(tx+1)%24
    for k,fr in enumerate(chamber.FREQUENCIES):
        ij=np.ix_(rx,tx)
        e1=data['empty_first'][k][ij];e2=data['empty_last'][k][ij]
        ok=np.isfinite(e1)&np.isfinite(e2)
        a=(e1[ok]+e2[ok])/2; d=np.linalg.norm(e2[ok]-e1[ok])/2
        for case in chamber.TARGETS:
            m=data[case][k][ij][ok];norm=np.linalg.norm(m-a)
            rows.append(dict(q=q,case=case,freq_GHz=fr,recorded_channels=int(ok.sum()),empty_radius_floor=d,
                             specimen_signal=norm,min_radius_as_fraction_of_contrast=d/norm))
with (out/'manitoba_empty_repeat.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
res={}
for q in [6,8,12,16]:
    rr=[i for i in rows if i['q']==q];v=np.array([i['min_radius_as_fraction_of_contrast'] for i in rr])
    res[q]=dict(materials=4,conditions=len(v),channels=rr[0]['recorded_channels'],quantiles=np.quantile(v,[0,.5,1]).tolist(),over_3pct=int(sum(v>.03)))
(out/'manitoba_empty_repeat.json').write_text(json.dumps(dict(scope='24-port ring, fixed positions subsampled. Empty-repeat minimum common additive radius, normalized to measured contrast signal. No feed calibration or loaded error upper.',scales=res,wall_s=time.perf_counter()-start),indent=2))
print(res)

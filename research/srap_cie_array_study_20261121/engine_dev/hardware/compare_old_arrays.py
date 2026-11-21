"""Compare every numerical old training/held-out array directly, including missing masks.
The remote copy is compressed and streamed in memory. No remote file is created.
"""
import base64
import io
import json
import os
from pathlib import Path
import subprocess
import sys
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[k]='1'
import numpy as np
here=Path(__file__).resolve().parent
hist=here.parents[2]/'limited_aperture_routed_dr_20261120'
remote_script="""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
import base64,io,json,sys
from pathlib import Path
import numpy as np
names={};unique={}
for folder in ['public','heldout']:
    for file in sorted((Path('.runs')/folder).glob('*.npz')):
        with np.load(file,allow_pickle=False) as z:
            for k in z.files:
                a=z[k]
                if a.dtype.kind not in 'iufcb':
                    continue
                nun=next((n for n,prev in unique.items() if prev.shape==a.shape and prev.dtype==a.dtype and np.array_equal(prev,a,equal_nan=True)),None)
                if nun is None:
                    nun='a'+str(len(unique))
                    unique[nun]=a
                names[folder+'/'+file.stem+'/'+k]=nun
buf=io.BytesIO(); np.savez_compressed(buf,**unique,manifest=np.array(json.dumps(names)))
print(base64.b64encode(buf.getvalue()).decode())
"""
cmd=r'cd /d D:\AI\LIMITED_APERTURE_20261120 && D:\python\python.exe -B -'
ssh=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8']
if os.getenv('SRAP_SSH_KEY'):
    ssh += ['-i',os.environ['SRAP_SSH_KEY']]
ssh.append(os.environ['SRAP_SSH_HOST'])
p=subprocess.run(ssh+[cmd],input=remote_script.encode(),capture_output=True)
if p.returncode:raise RuntimeError(p.stderr.decode()[-1000:])
with np.load(io.BytesIO(base64.b64decode(p.stdout))) as rec:
    names=json.loads(str(rec['manifest']))
    unique={k:rec[k] for k in set(names.values())}
    rows=[]
    for folder in ('public','heldout'):
        for f in sorted((hist/'.runs'/folder).glob('*.npz')):
            prefix=folder+'/'+f.stem+'/'
            with np.load(f,allow_pickle=False) as z:
                keys=sorted(set([k for k in z.files if z[k].dtype.kind in 'iufcb']+[k[len(prefix):] for k in names if k.startswith(prefix)]))
                for k in keys:
                    row=dict(folder=folder,case=f.stem,key=k)
                    if prefix+k not in names or k not in z:
                        row.update(status='missing_array',in_mac=k in z,in_xinan=prefix+k in names)
                    else:
                        a=z[k]
                        b=unique[names[prefix+k]]
                        row.update(mac_shape=list(a.shape),xinan_shape=list(b.shape))
                        if a.shape!=b.shape:row['status']='shape_mismatch'
                        elif not np.array_equal(np.isfinite(a),np.isfinite(b)):row['status']='finite_mask_mismatch'
                        else:
                            ix=np.isfinite(a)
                            ad,bd=a[ix],b[ix]
                            # Work on real/complex floating arrays for differences in unsigned/bool masks.
                            exact=np.array_equal(ad,bd)
                            er=0. if exact else float(np.linalg.norm(ad.astype(complex)-bd.astype(complex))/max(np.linalg.norm(ad),1.e-300))
                            row.update(status='equal' if exact else 'numerical_difference',
                                finite_elements=int(ix.sum()),relative_F=er,maximum_absolute=0. if exact else float(np.max(abs(ad.astype(complex)-bd.astype(complex)))))
                    rows.append(row)
    (here/'all_old_arrays.json').write_text(json.dumps(rows,indent=2))
    biggest=sorted((r for r in rows if r.get('relative_F',0)>0),key=lambda r:r['relative_F'],reverse=True)
    print('rows',len(rows),'largest',[(r['case'],r['key'],r['relative_F']) for r in biggest[:6]])

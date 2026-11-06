"""Saved-array postprocessing: coverage, soft-filter stability, no Maxwell calls."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import csv,json,time
import numpy as np
from threadpoolctl import threadpool_limits
HERE=Path(__file__).resolve().parents[1];OUT=HERE/'results'

def readcsv(p):return list(csv.DictReader(p.open()))
def writecsv(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)));w.writeheader();w.writerows(rows)

def main():
    start=time.perf_counter();cases=readcsv(OUT/'e1/cases.csv');ref=np.load(OUT/'e1/public_reference_geometry.npz')
    vref=ref['J_V'];gref=np.linalg.qr(np.array(ref['leakage_vectors'],order='C'))[0]
    rows=[];soft=[];cums=[];ids=[]
    for c in cases:
        tid=c['task_id'];a=np.load(OUT/'e1'/f'{tid}.npz');err=a['error_OFFLINE'];en=float(err@err);y=float(np.linalg.norm(a['Y']))
        bases={'anchor_J':a['J_V'][:,::-1],'anchor_Green':np.linalg.qr(np.array(a['leakage_vectors'],order='C'))[0],
               'reference_J':vref[:,::-1],'reference_Green':gref}
        for name,b in bases.items():
            energy=(b.T@err)**2;cum=np.cumsum(energy)/en
            row={k:c[k]for k in ('task_id','observation_id','shape','amplitude','noise','method')}
            row.update(basis=name,r90=int(np.searchsorted(cum,.90)+1),r95=int(np.searchsorted(cum,.95)+1),
                       r99=int(np.searchsorted(cum,.99)+1),coverage_32=float(cum[31]),coverage_64=float(cum[63]),
                       coverage_128=float(cum[127]),coverage_256=float(cum[255]),coverage_512=float(cum[511]))
            rows.append(row);cums.append(cum);ids.append(tid+'|'+name)
        va=a['J_V'];dot2=(vref.T@va)**2
        eps=1e-4 if c['noise']=='clean' else 1e-3;r=.03;lam=(eps/(r*24))**2
        wr=lam/((ref['J_s']/y)**2+lam);wa=lam/((a['J_s']/y)**2+lam)
        difference2=max(0,float(wr@wr+wa@wa-2*(wr[:,None]*dot2*wa[None,:]).sum()))
        weak=a['J_s']/y*r*24<=eps;coef=va.T@err
        soft.append(dict(task_id=tid,shape=c['shape'],amplitude=c['amplitude'],noise=c['noise'],method=c['method'],
             epsilon=eps,RMS=r,soft_filter_relative_frobenius_difference=float(np.sqrt(difference2/(wr@wr))),
             data_tolerance_weak_dimension=int(weak.sum()),data_tolerance_weak_GT_energy=float((coef[weak]@coef[weak])/en),
             generalized_onepercent_dimension=int(a['approx_basis'].shape[1])))
    writecsv(OUT/'coverage_dimension.csv',rows);writecsv(OUT/'fixed_filter_stability.csv',soft)
    np.savez_compressed(OUT/'coverage_cumulative.npz',key=np.array(ids),fraction=np.array(cums))
    low=[r for r in rows if float(r['amplitude'])==.4];slow=[r for r in soft if float(r['amplitude'])==.4]
    summary={}
    for name in bases:
        rr=[r for r in low if r['basis']==name]
        summary[name]={k:float(np.median([r[k]for r in rr]))for k in ('r90','r95','r99','coverage_32','coverage_64','coverage_128','coverage_256','coverage_512')}
    summary['soft_filter']={k:float(np.median([r[k]for r in slow]))for k in ('soft_filter_relative_frobenius_difference','data_tolerance_weak_dimension','data_tolerance_weak_GT_energy')}
    summary['runtime']=dict(wall_seconds=time.perf_counter()-start,new_physics=0,source='saved E1 arrays only')
    (OUT/'GEOMETRY_SUMMARY.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary))

if __name__=='__main__':
    with threadpool_limits(limits=1):main()

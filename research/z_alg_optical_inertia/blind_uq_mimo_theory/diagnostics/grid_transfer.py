"""Paired resolution challenge after all primary readouts and inverse predictions were fixed."""
import csv,json,os,sys,time
from pathlib import Path
for k in ('VECLIB_MAXIMUM_THREADS','OPENBLAS_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'imaging'))
from shared import Model
OUT=P/'results'; TARGETS=[f'n32_c{i:02}' for i in range(1,9)]+['n96_c01','n96_c02']
def reps():
    # All labels are used only for the separate post-freeze forward comparison.
    for c in TARGETS:
        with np.load(P/'ports/private'/f'{c}.npz') as z:yield c,'scoring_material',z['material']
    for c in ['n32_c01','n32_c02','n32_c03','n32_c04','n32_c05','n32_c06','n32_c08']:
        f=OUT/'extended'/f'{c}.json'
        if not f.exists():continue
        meta=json.loads(f.read_text())
        with np.load((P/meta['zero_file']).with_suffix('.npz')) as a:yield c,'zero',a['material'].reshape(32,32)
        with np.load(f.with_suffix('.npz')) as a:yield c,'extended',a['material'].reshape(32,32)
def run():
    start=time.perf_counter(); cpu=time.process_time();models={s:Model.read('n%d_c01'%s) for s in [32,96]};rows=[];fields={}
    for c,tag,x in reps():
        wall=time.perf_counter(); cpu0=time.process_time();n=x.shape[0];other=96 if n==32 else 32
        x2=np.repeat(np.repeat(x,3,0),3,1) if n==32 else x.reshape(32,3,32,3).mean((1,3))
        a=models[n].forward(x.ravel(),fresh=True);b=models[other].forward(x2.ravel(),fresh=True)
        with np.load(P/'ports/public/cases'/f'{c}.npz') as obs:Y=obs['Y'];r=float(obs['radius'])
        f,g=a['pred'],b['pred'];fields[c+'_'+tag+'_native']=f;fields[c+'_'+tag+'_changed']=g
        rows.append(dict(case=c,role=tag,from_side=n,to_side=other,contract='piecewise constant exact prolongation' if n==32 else 'conservative cell average restriction',
                         relative_readout_change=float(np.linalg.norm(f-g)/np.linalg.norm(f)),difference_over_r=float(np.linalg.norm(f-g)/r),
                         native_residual_over_r=float(np.linalg.norm(f-Y)/r),other_residual_over_r=float(np.linalg.norm(g-Y)/r),
                         elapsed_wall=time.perf_counter()-wall,elapsed_cpu=time.process_time()-cpu0))
    with (OUT/'grid_transfer.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    np.savez_compressed(OUT/'grid_transfer.npz',**fields)
    (OUT/'grid_transfer_cost.json').write_text(json.dumps(dict(materials=10,wall=time.perf_counter()-start,cpu=time.process_time()-cpu,solves={s:dict(m.cost) for s,m in models.items()},note='Original grid remains the only numerical guarantee. Intergrid differences are samples.'),indent=2))
if __name__=='__main__':
    with threadpool_limits(1):run()

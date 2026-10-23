"""ORACLE DIAGNOSTIC only. Not imported by online_inference."""
import gc,time,json,sys
import numpy as np
from .accounting import ROOT,Book,write
from wq.fem import FEM,PMLConfig
from wq.meshes import wall_loops,interface_mesh,circle_loop
from p0.geometry import wall_values,in_query,antennas,probes

def coefficient(seed,f,probe=None):
    def eps(pts):
        re,sigma=wall_values(pts,seed)
        re[in_query(pts)]=1.; sigma[in_query(pts)]=0.
        if probe is not None:
            center,radius,contrast=probe; mask=np.linalg.norm(pts-np.asarray(center),axis=1)<=radius
            re[mask]=1+contrast; sigma[mask]=.02
        return re+1j*sigma/f
    return eps

def reference(seed,band,label='cal'):
    old=json.loads((ROOT/'configs/AUTHORIZATION.json').read_text())['historical_array_root']
    with np.load(__import__('pathlib').Path(old)/'wall-pml/ARRAYS.npz') as z:
        return z[f'seed{seed}_{band}_{label}_volume'].copy()

def persist(path,p,s,**extras):
    path.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(path,mesh_p=p.mesh.p,mesh_t=p.mesh.t,points=s.material_points,
       weights=s.material_weights,epsilon=s.epsilon,field=s.field,current=s.current,
       scattered_dofs=s.scattered_dofs,data=s.data,probe_data=s.probe_data,direct=s.direct,
       residuals=s.residuals,**extras)

def main(mode='r0'):
    book=Book('R'); tx,rx=antennas(); cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text())
    out=ROOT/'results/R0';out.mkdir(parents=True,exist_ok=True); rows=[]
    frequencies=sorted(set(sum(cfg['frequencies'].values(),[]))); probe=probes()[8]
    for seed in cfg['truth_seeds']:
        pm=PMLConfig(strength=6)
        mesh=interface_mesh(h=.145,loops=wall_loops(seed,5e-5,2)+[circle_loop(probe[0],probe[1],256)])
        outputs={}; probeouts={}
        for f in frequencies:
            p=FEM(f,tx,rx,book,h=.145,order=3,mesh=mesh,pml=pm,receiver_mode='volume',phase='R0_truth_closure')
            for label,pr in [('cal',None),('task',probe)]:
                t=time.perf_counter(); s=p.solve(coefficient(seed,f,pr)); outputs[f,label]=s.data.copy();probeouts[f,label]=s.probe_data.copy()
                persist(out/f'{seed}_{f:.8f}_{label}.npz',p,s)
                row=dict(seed=seed,frequency=f,label=label,ndofs=p.ndofs,residual=s.residual,wall_seconds=time.perf_counter()-t,metadata=s.metadata,identity='ORACLE_DIAGNOSTIC')
                rows.append(row);write(out/'STATES.json',rows)
                if len(rows)==1: write(ROOT/'BENCHMARK.json',dict(source='first original truth endpoint, charged to R',**row,forecast='16 RHS per frequency for forward+receiver Full-J; 80 per full five-frequency iterate',peak_memory_policy='one frequency factor at a time; old factors released'))
                del s
            del p;gc.collect()
        for band,fs in cfg['frequencies'].items():
            rng=np.random.default_rng(2026102307+seed+(100 if band=='wide' else 0))
            for label in ['cal','task']:
                pred=np.stack([outputs[f,label] for f in fs]);old=reference(seed,band,label)
                noise=.01*(rng.normal(size=pred.shape)+1j*rng.normal(size=pred.shape))
                row=dict(seed=seed,band=band,label=label,d_mean=float(np.linalg.norm(pred-old)/.01),d_observed_truth=float(np.linalg.norm(noise)/.01),volume_probe=float(np.linalg.norm(pred-np.stack([probeouts[f,label] for f in fs]))/.01),same_contract_pass=bool(np.linalg.norm(pred-old)/.01<1e-7),identity='ORACLE_DIAGNOSTIC')
                np.savez_compressed(out/f'OBSERVATION_{seed}_{band}_{label}.npz',mean=old,replay=pred,observed=old+noise,noise=noise,frequencies=fs)
                rows.append(row);write(out/'STATES.json',rows); print(json.dumps(row),flush=True)
    write(out/'RESULTS.json',dict(status='PASS' if all(r.get('same_contract_pass',True) for r in rows) else 'FAIL',rows=rows,totals=book.totals(),new_noisy_observations='diagnostic only; original inverse search remains mean-data target',old_cross_solver_qualification='reused .00653–.03322, not rerun entire P0-Q'))

if __name__=='__main__':main(*sys.argv[1:])

"""Substantive new-interface audit: native parity, replay and M1 pullback.

This does not repeat a frozen Maxwell identity or score material truth.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import argparse,json,time,sys
from pathlib import Path
from collections import Counter


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--input',type=Path,required=True)
    a=p.parse_args();cpu=time.process_time();wall=time.perf_counter()
    sys.path[:0]=[str(a.root/'src'),str(a.legacy_repo),str(a.legacy_repo/'vendor')]
    import numpy as np
    from research.research_fbg_r3_experiments.src.physics import MaxwellModel
    from research.research_max_struct_r8.src.cie_fft import _run_cie
    from research.research_max_struct_r7.src.aplus_readout import mobius_material_readout
    from run_experiment import build_config,dump
    from dualmetric_core import StatefulDualMetric
    out=a.root/'audit/wiring_v1';out.mkdir(parents=True,exist_ok=False)
    with np.load(a.input,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
    Y=public['Y'];model=MaxwellModel(public);config,dec=build_config(model,public,out)
    targets={32,501,1001};refs={};n=[0]
    class StopNative(BaseException):pass
    def observe(row):
        n[0]+=1
        if n[0] in targets:refs[n[0]]={k:row[k].copy() for k in ('m','J','coefficients')}
        if n[0]==max(targets):raise StopNative()
    inherited=dict(config,current_transform='scalar',material_transform='onehop',_counter=Counter(),
                   material_readout=mobius_material_readout,iteration_callback=observe)
    try:_run_cie(model,Y,public['m_init'],inherited)
    except StopNative:pass
    solver=StatefulDualMetric(MaxwellModel(public),Y,public['m_init'],config,pair='01')
    parity=[];replay=[];pullback=[]
    for k in range(1,max(targets)+1):
        row=solver.step()
        if k in targets:
            for name in ('m','J','coefficients'):
                v=row[name];r=refs[k][name]
                parity.append(dict(iteration=k,array=name,exact=bool(np.array_equal(v,r)),relative_error=float(np.linalg.norm(v-r)/max(np.linalg.norm(r),1e-300))))
        if k in (32,499,999):
            path=out/f'checkpoint_{k}.npz';solver.checkpoint(path)
            with np.load(path,allow_pickle=False) as z:keys=z.files
            restored=StatefulDualMetric.from_checkpoint(MaxwellModel(public),Y,config,path)
            twin=StatefulDualMetric.from_state_dict(MaxwellModel(public),Y,config,solver.state_dict())
            for step in range(3):
                r=restored.step();t=twin.step()
                replay.append(dict(checkpoint_iteration=k,steps=step+1,exact=bool(all(np.array_equal(r[x],t[x]) for x in ('m','J','Q','coefficients','old_grad','direction'))),stage=int(r['stage']),readable_npz_keys=keys))
    rng=np.random.default_rng(14014)
    for pair in ('00','01','10','11'):
        clone=StatefulDualMetric.from_state_dict(model,Y,config,solver.state_dict(),pair=pair)
        x=rng.normal(size=clone.A.shape)+1j*rng.normal(size=clone.A.shape)
        z=rng.normal(size=clone.J.shape)+1j*rng.normal(size=clone.J.shape)
        lhs=np.vdot(clone._H(x,clone.m),z);rhs=np.vdot(x,clone._HH(z,clone.m))
        pullback.append(dict(pair=pair,relative_adjoint_defect=float(abs(lhs-rhs)/max(abs(lhs),abs(rhs),1e-300))))
    result=dict(status='COMPLETE',truth_access='NONE',purpose='new current-metric wiring and full-state re-entry',parity=parity,
                replay=replay,pullback=pullback,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                physical_cost=dict(model.cost),dec=dec)
    dump(out/'RECEIPT.json',result)
    print(json.dumps(dict(max_native_relative_error=max(x['relative_error'] for x in parity),all_replay_exact=all(x['exact'] for x in replay),max_pullback_defect=max(x['relative_adjoint_defect'] for x in pullback),CPU_seconds=result['CPU_seconds'])))


if __name__=='__main__':main()

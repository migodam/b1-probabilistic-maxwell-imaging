"""Separate controlled-prior simulation process. Never imported by inference.

The only oracle labels in P1 are generated here AFTER public prior freezing.
Online files contain observations, parent and noise identities only.
"""
import sys,json
import numpy as np
from .accounting import ROOT,Book,write
from .probability_family import Family

def main(stage):
    if stage not in ['P1A','P1B']:raise ValueError(stage)
    protocol=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text())
    family=Family.load();book=Book(stage)
    out=ROOT/'results'/stage/'data';out.mkdir(parents=True,exist_ok=True)
    identity=out/'PUBLIC_MODEL_SNAPSHOT.json'
    if identity.exists():
        if json.loads(identity.read_text())!=family.config:raise RuntimeError('data resume public family differs')
    else:write(identity,family.config)
    wall_copy=out/'PUBLIC_WALLS.npz'
    if wall_copy.exists():
        z=np.load(wall_copy)
        if any(not np.array_equal(z[f'q{i}'],q) for i,q in enumerate(family.walls)):raise RuntimeError('data resume wall fields differ')
    else:np.savez_compressed(wall_copy,q0=family.walls[0],q1=family.walls[1])
    seed=protocol['data_seeds'][stage];count=protocol['independent_parents'];replicas=protocol['noise_replicas']
    rng=np.random.default_rng(seed);plan=[]
    # Draw ALL labels/noise before physical evaluation. No residual-based case
    # selection, adaptive stopping or use of labels in the inference code.
    for index in range(count):
        b=int(rng.integers(2));u=float(rng.uniform(-.01,.01));c=float(rng.uniform(2.94,3.06))
        a=float(rng.uniform(.99,1.01));p=float(rng.uniform(-.01,.01))
        noise=rng.normal(size=(replicas,2,5,4,12))*.01
        plan.append(dict(parent=index,branch=b,u=u,contrast=c,amplitude=a,phase=p,noise=noise))
    plan_path=out/'OFFLINE_GENERATION_PLAN.npz'
    if plan_path.exists():
        old=np.load(plan_path);expected=np.array([[r[k] for k in ['branch','u','contrast','amplitude','phase']] for r in plan])
        if not np.array_equal(old['parameters'],expected) or not np.array_equal(old['noise'],np.array([r['noise'] for r in plan])) or int(old['seed'])!=seed:raise RuntimeError('resume simulation plan mismatch')
    else:np.savez_compressed(plan_path,parameters=np.array([[r[k] for k in ['branch','u','contrast','amplitude','phase']] for r in plan]),noise=np.array([r['noise'] for r in plan]),seed=seed)
    ys=[];truths=[];params=[];new_points=0;reused_points=0
    for r in plan:
        dest=out/f'parent_{r["parent"]:03d}';complete=dest/'COMPLETE.npz'
        if complete.exists():
            z=np.load(complete);f=z['F'];y=z['y'];reused_points+=1
            if any(float(z[k])!=float(r[k]) for k in ['branch','u','contrast','amplitude','phase']):raise RuntimeError('saved parent labels differ from frozen plan')
            if f.shape!=(5,4,12) or y.shape!=(replicas,5,4,12) or not np.all(np.isfinite(f)) or not np.all(np.isfinite(y)):raise RuntimeError('saved parent response/observation contract mismatch')
            expected=r['amplitude']*np.exp(1j*r['phase'])*f[None]+r['noise'][:,0]+1j*r['noise'][:,1]
            if not np.array_equal(y,expected):raise RuntimeError('saved noisy data do not match frozen generation plan')
        else:
            # Every incomplete attempt stays on disk and in the same ledger.
            dest.mkdir(exist_ok=True)
            attempt=len(list(dest.glob('physical_attempt_*')))
            book.event(event='new_controlled_truth_simulation',parent=r['parent'],evidence_mode='DATA_GENERATION_NOT_ONLINE_INFERENCE',stage_data_seed=seed)
            f=family.point(r['branch'],r['u'],r['contrast'],book,str((dest/f'physical_attempt_{attempt:03d}').relative_to(ROOT/'results')))
            new_points+=1
            mean=r['amplitude']*np.exp(1j*r['phase'])*f
            y=mean[None]+r['noise'][:,0]+1j*r['noise'][:,1]
            np.savez_compressed(complete,F=f,y=y,branch=r['branch'],u=r['u'],contrast=r['contrast'],amplitude=r['amplitude'],phase=r['phase'])
        ys.extend(y);truths.extend([f]*replicas);params.extend([[r['branch'],r['u'],r['contrast'],r['amplitude'],r['phase']]]*replicas)
        write(out/'PROGRESS.json',dict(completed_parents=r['parent']+1,total_parents=count,stage=stage,seed=seed))
    # Public online file has no latent labels or noiseless means.
    np.savez_compressed(out/'OBSERVATIONS.npz',y=np.array(ys),parent_ids=np.repeat(np.arange(count),replicas),noise_ids=np.tile(np.arange(replicas),count))
    np.savez_compressed(out/'OFFLINE_TRUTH.npz',parameters=np.array(params),F=np.array(truths),parent_ids=np.repeat(np.arange(count),replicas))
    write(out/'COMPLETE.json',dict(stage=stage,independent_parents=count,observations=len(ys),new_full_vectors_this_invocation=new_points,reused_complete_parent_vectors=reused_points,total_parent_vectors=count,identity='INDEPENDENT_CONTROLLED_PRIOR_SIMULATION_AFTER_PRIOR_FREEZE',online_truth_access=False))

if __name__=='__main__':main(sys.argv[1])

"""Measurement-only proposal bank and separately stored prior-generated data."""
from pathlib import Path
import json,gc,time
import numpy as np
from wslam.accounting import BudgetBook,json_default
from wslam.finite_bayes import build_coarse_grid
from wslam.forward import Evaluator
from wslam.physics import acquisition
from wslam.world import World
root=Path(__file__).resolve().parents[1]
if json.loads((root/'results/e0/gate.json').read_text())['decision']!='PASS':raise RuntimeError('E0 failed')
out=root/'results/a0/preparation';out.mkdir(parents=True,exist_ok=False)
online=root/'inputs/online';offline=root/'inputs/offline';online.mkdir(parents=True,exist_ok=True);offline.mkdir(parents=True,exist_ok=True)
cfg=json.loads((root/'EXPERIMENT_CONTRACT.json').read_text());book=BudgetBook();tx,rx=acquisition(outer=True)
grid=build_coarse_grid();N=grid['state_count'];data=np.empty((N,4,16),complex);resids=[]
np.savez_compressed(online/'COARSE_GRID.npz',**{k:v for k,v in grid.items() if isinstance(v,np.ndarray)})
solver=Evaluator(book,'A0',tx=tx,rx=rx);started=time.perf_counter()
for i in range(N):
    world=World(grid['support'][i].reshape(2,2),grid['epsinf'][i],grid['ell'][i],grid['pose'][i])
    ev=solver.evaluate(world);data[i]=ev.data;resids.append(ev.state.residual)
    if (i+1)%25==0 or i+1==N:
        np.savez_compressed(online/'COARSE_BANK.npz',data=data[:i+1],residuals=np.asarray(resids),completed=i+1,total=N,tx=tx,rx=rx)
        print(json.dumps(dict(case='coarse_bank',completed=i+1,total=N,counts=book.totals('A0'))),flush=True)
    del ev
solver.clear();gc.collect()
# New worlds are drawn once from the frozen proper prior; every support is
# unknown to inference. The source and optimizer never load this offline file.
rng=np.random.default_rng(cfg['seeds']['truth']);truths=[];means=[];observations=[]
support_prior=grid['support_prior'];legal=grid['legal_supports']
def truncated(sd):
    while True:
        x=rng.normal(0,sd)
        if abs(x)<=.01:return float(x)
for parent in range(cfg['a0']['development_parents']):
    sid=int(rng.choice(len(legal),p=support_prior))
    world=World(legal[sid],float(rng.uniform(1.5,6)),float(rng.uniform(0,.25)),[rng.uniform(-1,1),rng.uniform(1.6,2.6)])
    amp,phase=truncated(.005),truncated(.005);gamma=(1+amp)*np.exp(1j*phase)
    ev=solver.evaluate(world);mean=gamma*ev.data;means.append(mean)
    row=dict(parent_id=parent,support_id=sid,**world.metadata(),instrument=dict(amplitude_offset=amp,phase=phase))
    truths.append(row)
    for repeat in range(cfg['a0']['noise_repeats']):
        nrng=np.random.default_rng(cfg['seeds']['noise']+parent*37+repeat)
        noise=.01*(nrng.normal(size=mean.shape)+1j*nrng.normal(size=mean.shape))
        observations.append(mean+noise)
    if parent in (0,1):
        np.savez_compressed(offline/f'PARENT_{parent:02d}_FIELD.npz',data=mean,points=ev.state.material_points,epsilon=ev.state.epsilon,field=ev.state.field,mesh_p=ev.backend.mesh.p,mesh_t=ev.backend.mesh.t)
    (offline/'CONTINUOUS_TRUTHS.json').write_text(json.dumps(truths,indent=2,default=json_default))
    print(json.dumps(dict(case='prior_parent',parent=parent,counts=book.totals('A0'))),flush=True)
    del ev
solver.clear()
np.savez_compressed(offline/'CONTINUOUS_CLEAN_MEANS.npz',data=np.stack(means))
np.savez_compressed(online/'A0_OBSERVATIONS.npz',base=np.stack(observations)[:,:2],parent_id=np.repeat(np.arange(len(truths)),2),noise_repeat=np.tile([0,1],len(truths)),tx=tx[:2],rx=rx)
# Holdout is separate and will be read only by the post-optimization scorer.
np.savez_compressed(offline/'A0_HOLDOUT.npz',data=np.stack(observations)[:,2:],tx=tx[2:],rx=rx)
cal=np.random.default_rng(cfg['seeds']['calibration']).normal(size=(4096,64))
upper=float(np.quantile(np.sum(cal*cal,axis=1),.99))
(out/'DISCREPANCY_RULE.json').write_text(json.dumps(dict(upper_noise_energy=upper,dimension=64,quantile=.99,
    source='4096 frozen conditional prior-predictive noise draws; F(theta) cancels at calibration parameter',
    warning='one-sided compatibility envelope, not fitted residual degrees of freedom or calibrated posterior',fit_and_holdout_separate=True),indent=2))
(out/'COMPLETE.json').write_text(json.dumps(dict(coarse_states=N,independent_parents=len(truths),observations=len(observations),newtruthlabels=len(truths),counts=book.totals('A0'),wall_seconds=time.perf_counter()-started,online_loader_forbidden_fields=['support truth','true material','true pose','holdout observations']),indent=2))
print('Preparation complete',flush=True)

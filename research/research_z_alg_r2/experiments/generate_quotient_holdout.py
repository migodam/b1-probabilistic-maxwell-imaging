"""Isolated fresh-data generator; never imported by an inverse solver."""
import os
for v in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[v]='1'
import pathlib,sys,json,secrets,time
import numpy as np
P=pathlib.Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model
from validate import CountingBook
out=P/'experiments/quotient_holdout_public';private=P/'private/quotient_holdout'
if out.exists():raise RuntimeError('Fresh holdout already generated')
c=json.loads((P/'experiments/HIGH4_CONTRACT.json').read_text())
contract={'status':'FROZEN_BEFORE_GENERATION','frozen_unix':time.time(),'geometry':c['geometry'],'box':[.5,30.],'targets':4,'law':'independent uniform real contrasts, new private seed; no selection','methods':['quotient','six','ReducedFullJ8','ReducedFullJ64'],'parameters':'same quotient and continuation source and reduced-TRF settings as frozen QUOTIENT_BENCHMARK_CONTRACT; no tuning','root_validation':{'relative_response':1e-8,'state':1e-10,'sigma_min':1e-6,'merge_RMS':1e-7},'blindness':'separate generator, solver, scorer processes; solver public input only; not OS access isolation','source_snapshot':[]}
for p in [P/'src/quotient_solver.py',P/'src/branch_enumerator.py',P/'src/physics.py',P/'src/algebraic_solver.py']:
 st=p.stat();contract['source_snapshot'].append({'path':str(p.relative_to(P)),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns})
(P/'experiments/QUOTIENT_HOLDOUT_CONTRACT.json').write_text(json.dumps(contract,indent=2))
out.mkdir();private.mkdir();seed=secrets.randbits(63);rng=np.random.default_rng(seed);b=CountingBook();model=Model(b,**c['geometry']);manifest=[]
for i in range(4):
 m=rng.uniform(.5,30,4);case=f'QH4_{i:02d}';y=model.forward(m)['y']
 np.savez(out/f'{case}.npz',y=y)
 np.savez(private/f'{case}.npz',m=m)
 manifest.append({'case':case,'observations':f'{case}.npz'})
(out/'manifest.json').write_text(json.dumps(manifest,indent=2))
(private/'generator_state.json').write_text(json.dumps({'seed':seed,'distribution':'uniform','box':[.5,30.]}))
(P/'results/quotient_holdout_generation.json').write_text(json.dumps({'generated_targets':4,'seed_exposed_to_solver':False,'accounting':b.summary(model.ns,model.nr)},indent=2))
print('4 fresh targets generated. Only observations and geometry public; GT values and seed withheld.')

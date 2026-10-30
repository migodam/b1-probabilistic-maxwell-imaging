"""Separate generation process. Inverse modules never import this file."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name]='1'
import argparse, json, sys, time, uuid
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
sys.path.insert(0,str(REPO))
from research.research_z_alg_r3_digits.src.physics import Model
from research.research_z_alg_r3_digits.experiments.generate_data import digit
HISTORY=ROOT.parent/'research_z_alg_r3_digits'

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('split',choices=['factor','fresh'])
    args=parser.parse_args()
    begun,cpu=time.perf_counter(),time.process_time()
    operators=dict(np.load(HISTORY/'public_data/operators.npz',allow_pickle=False))
    model=Model(operators)
    public=ROOT/'public_data'/args.split
    truth=ROOT/('truth_development' if args.split=='factor' else 'private_truth/fresh')
    public.mkdir(parents=True,exist_ok=True);truth.mkdir(parents=True,exist_ok=True)
    rows=[]
    if args.split=='factor':
        entries=[]
        for label,split,case in [('DEV7','dev_v2','dev_2a5c221387b1'),('holdout7','holdout','obs_c3604cf3555d')]:
            source=dict(np.load(HISTORY/'truth_after_freeze'/split/(case+'.npz'),allow_pickle=False))
            if 'alpha' in source:
                alpha=source['alpha'].ravel()
            else:
                mapping=json.loads((HISTORY/'truth_after_freeze'/split/'mapping_seedless.json').read_text())
                params=next(r['parameters'] for r in mapping['records'] if r['observation_id']==case)
                alpha=(source['m'].ravel()-params['background'])/(params['foreground']-params['background'])
            for amplitude in (1.3,6.,15.):
                entries.append((7,1.2+amplitude*alpha,alpha,dict(shape=label,background=1.2,amplitude=amplitude,scope='EXPOSED_TRUTH_FACTOR_DEVELOPMENT')))
    else:
        contract=json.loads((ROOT/'configs/FRESH_CONTRACT.json').read_text())
        rng=np.random.default_rng(int.from_bytes(os.urandom(16),'big'))
        entries=[]
        for label,font in zip((2,7,8),('Helvetica Bold.ttf','Georgia Bold.ttf','Helvetica Bold.ttf')):
            old,mask,params=digit(label,font,rng)
            alpha=(old-params['background'])/(params['foreground']-params['background'])
            background=float(rng.uniform(1.,1.4));foreground=float(rng.uniform(13.,18.))
            material=background+(foreground-background)*alpha
            params.update(background=background,foreground=foreground,scope='FRESH_AFTER_METHOD_FREEZE')
            entries.append((label,material,alpha,params))
        rng.shuffle(entries)
    with threadpool_limits(limits=1):
        for label,material,alpha,params in entries:
            case=args.split+'_'+uuid.uuid4().hex[:12]
            state=model.forward(material,label='generation')
            np.savez_compressed(public/(case+'.npz'),Y=state['y'],grid_shape=np.array([24,24]),case_id=case)
            np.savez_compressed(truth/(case+'.npz'),m=material,alpha=alpha,mask=alpha>=.5)
            rows.append(dict(case_id=case,digit=label,parameters=params))
    (truth/'mapping.json').write_text(json.dumps({'records':rows,'seeds_saved':False},indent=2)+'\n')
    receipt=dict(split=args.split,generated_utc=datetime.now(timezone.utc).isoformat(),case_ids=[r['case_id'] for r in rows],
       model_counts=dict(model.counts),wall_seconds=time.perf_counter()-begun,cpu_seconds=time.process_time()-cpu,
       physics='unchanged fixed R3-DIGITS operators; no new Green setup',truth_input_to_inverse=False)
    (ROOT/'audit'/('GENERATION_'+args.split.upper()+'.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()

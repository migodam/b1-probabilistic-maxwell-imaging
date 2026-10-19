"""Recover raw, unpaid D draws from fixed seeds/checkpoints; zero physics."""
from pathlib import Path
from contextlib import nullcontext
from collections import Counter
import sys,json,time
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from material_fiber.core import OUT,Target,write,save
from material_fiber.learn import build_model,feature
from material_fiber.explore import geometry,direction
from p1_material.proposals import frozen_latent_step

class ZeroBook:
    stage='D';counts=Counter()
    def require(self,*a,**kw):raise RuntimeError('Replay forbids paid physics')
    def reserve_state(self,*a,**kw):raise RuntimeError('Replay forbids new states')
    def reserve_rhs(self,*a,**kw):raise RuntimeError('Replay forbids new RHS')
    def span(self,*a,**kw):return nullcontext()
    def physical(self,*a,**kw):return nullcontext()
    def append(self,*a,**kw):pass

def main():
    start=time.perf_counter();torch.set_num_threads(2);t=Target(ZeroBook());aa=json.loads((OUT/'A/ANCHORS.json').read_text());rr=json.loads((OUT/'D/SUMMARY.json').read_text())['runs'];models={};contexts={};n=0;failed=0
    for r in rr:
        if 'method' not in r:continue
        parent,kind,seed,j,s=r['parent'],r['method'],r['seed'],r['trial'],r['radius']
        if parent not in contexts:contexts[parent]=t.load(aa[parent])
        p=contexts[parent];y=t.observation(parent);g=geometry(t,p,s);rng=np.random.default_rng(np.random.SeedSequence([seed,j,int(s*10000)]));torch.manual_seed(seed+j);raw=p.q.copy()
        if kind=='nonlearned':
            v,_=direction(t,p,g,'soft_dynamic',rng)
            if v is not None:raw[:64]+=s*v[:64]
        else:
            if (kind,seed) not in models:
                cp=torch.load(OUT/f'D/models/{kind}_{seed}.pt',weights_only=False);model=build_model(kind,cp['nfeatures']);model.load_state_dict(cp['state_dict']);model.eval();models[kind,seed]=model
            with torch.no_grad():raw[:64]+=models[kind,seed].sample(torch.tensor(feature(t,p,y,s,g))[None])[0].numpy().astype(float)
        prng=np.random.default_rng(np.random.SeedSequence([seed,j,int(s*10000),991]));u=np.mod(raw[64:]+.04*prng.normal(size=2),2);raw[64:]=np.where(u<=1,u,2-u)
        chi,_=t.prior.decode(raw);save(OUT/f'D/raw_proposal_replay/{parent}_{kind}_{seed}_{j}.npz',raw_q=raw,raw_chi=chi,original_status=np.array(r.get('status','PHYSICS_EVALUATED')),radius=s,seed=seed,trial=j)
        n+=1;failed+=int('status' in r)
    write(OUT/'D/RAW_PROPOSAL_REPLAY.json',dict(draws=n,original_decoder_failures=failed,provenance='deterministic regeneration of raw pre-retraction draws from frozen seed/checkpoint/context; not newly discovered candidates',new_states=0,new_RHS=0,new_LU=0,wall_seconds=time.perf_counter()-start,postfit_materials_for_evaluated_candidates='original endpoint arrays',failed_postfit_materials='not falsely reconstructed as measured endpoints'))
    print(json.dumps(dict(draws=n,failed=failed,new_physics=0)))
if __name__=='__main__':main()

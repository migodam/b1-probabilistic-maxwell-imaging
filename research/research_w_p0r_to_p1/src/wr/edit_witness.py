"""Direct joint compensation and finite cavity birth: ORACLE existence only."""
import sys,json
import numpy as np
from p0.geometry import probes
from .physical_witness import ModulatedWorld,mesh_for
from .physics import Forward,pack,packj
from .accounting import ROOT,Book,ScopedBook,write,BudgetStop
from .optimization import trust_region

def main(seed,band='narrow',kind='joint',iterations=14):
    book=ScopedBook('B2',f'{kind}_{seed}_{band}');cfg=json.loads((ROOT/'configs/ORIGINAL_P0.json').read_text());freq=cfg['frequencies'][band]
    probe=probes()[8];zero=np.zeros(48);folder=f'B2/{kind}_direct_{seed}_{band}'
    base_model=ModulatedWorld(seed,probe if kind=='joint' else None)
    base=Forward(base_model,mesh_for(seed,probe),freq,book,folder+'/base',phase='B2_direct_reference',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6)
    d0,_=base.evaluate(zero)
    if kind=='joint':
        candidate_probe=((probe[0][0]+.1,probe[0][1]),probe[1],probe[2]);cavities=()
    else:candidate_probe=None;cavities=(((-.85,.85),.06),)
    model=ModulatedWorld(seed,candidate_probe,cavities=cavities)
    mesh=mesh_for(seed,candidate_probe,cavities=cavities)
    fw=Forward(model,mesh,freq,book,folder+'/edit',phase='B2_direct_nonlinear_compensation',identity='PHYSICAL_EXISTENCE_SEARCH',h=.145,pml_strength=6,retain_current=True)
    out=ROOT/'results'/folder;out.mkdir(parents=True,exist_ok=True);rows=[]
    raw,_=fw.evaluate(zero)
    for method in ['FULL_J','RANDOM_MULTISTART']:
        rng=np.random.default_rng(2026102321+seed);x0=zero.copy() if method=='FULL_J' else rng.normal(size=48)*.1
        history=[]
        def cb(row):history.append(row);write(out/f'{method}_ITERATIONS.json',history)
        fit=trust_region(fw.residual(d0),x0,np.full(48,-2.),np.full(48,2.),max_iter=iterations,max_calls=2*iterations+6,radius=.2,max_radius=.7,l2_radius=False,callback=cb,objective_floor=.2,fatal_exceptions=(BudgetStop,))
        d,_=fw.evaluate(fit.x)
        target_mass_change=0.0 if kind=='joint' else None
        row=dict(seed=seed,band=band,kind=kind,method=method,raw_distance=float(np.linalg.norm(raw-d0)/.01),final_distance=float(np.linalg.norm(d-d0)/.01),target_shift=.1 if kind=='joint' else 0.,target_mass_relative_change=target_mass_change,edit='target translation with identical circle mass' if kind=='joint' else 'new strict-air circular cavity of diameter .12',stop=str(fit.message),iterations=fit.nit,gradient=float(fit.projected_gradient),identity='ORACLE_CONSTRUCTED_PHYSICAL_WITNESS',status='NEEDS_INDEPENDENT_NUMERICAL_VALIDATION',same_source_noise_prior_bounds=True)
        np.savez_compressed(out/f'{method}_FINAL.npz',x=fit.x,data=d,reference_data=d0,raw_data=raw,probe=np.asarray(candidate_probe[0]) if candidate_probe else np.empty(0))
        rows.append(row);write(out/'RESULTS.json',rows);print(json.dumps(row),flush=True)
    write(out/'COMPLETE.json',dict(totals=book.totals(),probabilities=None))
if __name__=='__main__':main(int(sys.argv[1]),*(sys.argv[2:4]),**({'iterations':int(sys.argv[4])} if len(sys.argv)>4 else {}))

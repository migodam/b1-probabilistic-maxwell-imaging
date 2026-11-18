"""Noisy old data alone determine the bandwidth, spectrum, and candidate initializations."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from scipy.stats import qmc
from freq_model import Fit, FrequencyTM
from freq_fold import candidate_basis
from spectral_moments import coefficients, predict, ROOT, PATTERNS

sys.path.insert(0,str(ROOT.parent/'research_z_alg_r10/src'))
from homotopy_control import run as find_roots

SIGMA=3e-7
GAIN_SD=4e-4
LIMIT=.02
LAMBDA=(SIGMA/GAIN_SD)**2
PATTERNS_USED=['opposite_edges','two_positive','two_negative']


def nominal_gap(model,roots,band):
    if len(roots)<2:return 0.
    p=np.sqrt([1.,.5,.5])
    vec=[np.concatenate([w*model.response(q,x*band)[0] for w,x in zip(p,[0.,-1.,1.])]+
                        [np.array([np.sqrt(LAMBDA)])]) for q in roots]
    def directed(v,u):
        return np.linalg.norm(v-u*np.vdot(u,v)/np.vdot(u,u))
    return float(min(min(directed(u,v),directed(v,u)) for u,v in combinations(vec,2)))


def choose(model,old):
    start=perf_counter();cost0=model.factorizations
    scan=find_roots(model.geometry,old)
    roots=np.asarray(scan['validated_roots']).reshape(-1,6)
    roots=sorted(roots,key=lambda q:tuple(q))
    pairs=sorted(combinations(range(len(roots)),2),key=lambda ij:np.linalg.norm(roots[ij[0]]-roots[ij[1]]))
    b_dyn=0.; pair_chosen=None
    for i,j in pairs:
        pair=np.asarray([roots[i],roots[j]])
        basis=candidate_basis(model,pair)
        if not basis[-1]:continue
        c=coefficients(model,basis)
        eta=LAMBDA/(LAMBDA+2*c[-1])
        b_dyn=1.25*predict('opposite_edges',c,basis[3],1,eta)[0]
        pair_chosen=[i,j]
        break
    gap_grid=[]
    for b in np.geomspace(1.e-6,LIMIT,29) if len(roots)>=2 else []:
        gap=nominal_gap(model,roots,b)
        gap_grid.append((float(b),gap))
        if gap>=4*SIGMA:break
    uses_gap=bool(gap_grid and gap_grid[-1][1]>=4*SIGMA)
    gap_b=gap_grid[-1][0] if uses_gap else LIMIT
    b_choice=float(min(LIMIT,max(gap_b,b_dyn))) if pair_chosen else LIMIT
    return roots,b_choice,dict(central_roots=[x.tolist() for x in roots],
        path_classifications=scan['path_classifications'],paths_attempted=scan['paths_attempted'],
        path_jacobian_evaluations=scan['path_jacobian_evaluations'],root_wall_s=scan['total_wall_seconds'],
        selected_pair=pair_chosen,predicted_dynamics_band=b_dyn,minimum_gap_band=gap_b,nominal_gaps=gap_grid,
        design_factorizations=model.factorizations-cost0,design_wall_s=perf_counter()-start)


def compare_reads(model,old,extra,band,pattern,starts,roots):
    offsets,powers=PATTERNS[pattern]; amp=np.sqrt(powers)
    f=Fit(model,[x*band for x in offsets],[old,*extra],weights=amp,gain_calibration=(LAMBDA,1.))
    def candidates(seq):
        fit=[]
        for z in seq:
            y,work=f.solve(z,130)
            res=f.evaluate(y)[0]
            fit.append((float(np.linalg.norm(res)**2/SIGMA**2),y,work))
        best=min(fit,key=lambda k:k[0])
        return {'prediction':best[1].tolist(),'loss_over_sigma2':best[0],
                'factorizations':int(sum(v[2] for v in fit)),'starts':len(seq)}
    # Both alternatives have three starts. Branch ranking and all final selection use received data only.
    ordered=sorted(roots,key=lambda z:np.linalg.norm(f.evaluate(z)[0]))
    chosen=ordered[:3]+starts[:max(0,3-len(ordered))]
    return {'ordinary_three':candidates(starts),'old_roots_three':candidates(chosen)}


def study(repeats):
    archive=json.loads((ROOT/'results/spectral_new_layouts.json').read_text())
    record=[]
    for i,row in enumerate(archive):
        if not row['included']:continue
        pair=np.array(row['pair']); model=FrequencyTM(row['geometry'])
        for t,truth in enumerate(pair):
            for rep in range(repeats):
                rng=np.random.default_rng(667029+rep*131+i*13+t)
                n=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))*SIGMA/np.sqrt(2)
                g=1+.65*GAIN_SD*np.exp(1j*(rep*2.13+.41))
                old=g*model.response(truth,0.)[0]+n[0]
                roots,b,info=choose(model,old)
                starts=[np.full(6,12.),*list(qmc.Sobol(6,scramble=True,seed=90341+rep).random_base2(1)*29.5+.5)]
                arms={}
                for label,pattern,band in [('adaptive_split','opposite_edges',b),('adaptive_plus','two_positive',b),
                    ('adaptive_minus','two_negative',b),('fixed_wide_split','opposite_edges',LIMIT),
                    ('fixed_wide_plus','two_positive',LIMIT),('repeat','opposite_edges',0.)]:
                    offsets=PATTERNS[pattern][0]
                    extra=[g*model.response(truth,di*band)[0]+np.sqrt(2)*n[k] for k,di in enumerate(offsets[1:],1)]
                    a=compare_reads(model,old,extra,band,pattern,starts,roots)
                    for result in a.values():
                        q=np.array(result['prediction'])
                        result['RMS']=float(np.linalg.norm(q-truth)/np.sqrt(6))
                        result['closer_to_generation_than_other_old_material']=bool(
                            np.linalg.norm(q-truth)<np.linalg.norm(q-pair[1-t]))
                    arms[label]=dict(band=band,extra_ri=np.stack((np.real(extra),np.imag(extra)),axis=-1).tolist(),fits=a)
                record.append(dict(label=row['label'],truth_side=t,rep=rep,common_gain=[g.real,g.imag],old_ri=np.c_[old.real,old.imag].tolist(),
                                   sigma=SIGMA,gain_sd=GAIN_SD,design=info,arms=arms))
            print(row['label'],t,'completed',len(record),flush=True)
            (ROOT/'results/spectral_online.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--repeats',type=int,default=3)
    study(p.parse_args().repeats)

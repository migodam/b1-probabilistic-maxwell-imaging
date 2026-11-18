"""Near-fold spectral moment prediction versus full multi-frequency Maxwell inversion."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import root

from freq_model import Fit, FrequencyTM, ri
from freq_fold import candidate_basis, exact_profile, fold_case

ROOT=Path(__file__).resolve().parents[1]
R10=ROOT.parent/'research_z_alg_r10/results'

# Weights here are transmitted energies. All new tones together cost one old block.
PATTERNS={
    'two_positive':([0.,1.,1.],[1.,.5,.5]),
    'two_negative':([0.,-1.,-1.],[1.,.5,.5]),
    'opposite_edges':([0.,-1.,1.],[1.,.5,.5]),
    'inner_symmetric':([0.,-2**-.5,-2**-.5,2**-.5,2**-.5],[1.,.25,.25,.25,.25]),
    'edges_and_repeats':([0.,-1.,0.,0.,1.],[1.,.25,.25,.25,.25])
}


def coefficients(model,basis):
    center,v,_,_,_=basis
    y,j,h=model.response(center,hessian=True)
    *_,dj=model.response(center,derivative=True)
    left=np.linalg.svd(ri(j))[0][:,-1]
    curv=float(left@ri(np.einsum('rij,i,j->r',h,v,v)))/2
    derivative=ri(dj@v)
    gain_normal=left@np.column_stack((ri(y),ri(1j*y)))
    return curv,float(left@derivative),float(np.linalg.norm(derivative)),float(gain_normal@gain_normal)


def energy_moments(pattern):
    x,p=(np.asarray(v) for v in PATTERNS[pattern])
    mu=float(np.average(x,weights=p))
    return mu,float(np.average((x-mu)**2,weights=p)),float(np.average(x*x,weights=p))


def predict(pattern,coef,a,side,eta):
    c,dw,dn,_=coef
    mean,var,_=energy_moments(pattern)
    scale=np.sqrt(2*var)*dn-np.sqrt(eta)*np.sign(c*side)*mean*dw/2
    if scale<=0:return None,None
    b=np.sqrt(eta)*abs(c)*a/scale
    x=-side/2-3*mean*b*dw/(4*c*a)
    return b,x


def investigate(model,pair,label,effective_gains,patterns,must_have_center):
    wall=perf_counter()
    basis=candidate_basis(model,pair)
    center,v,_,a,valid=basis
    if not valid and must_have_center:
        return {'label':label,'included':False,'reason':'no regular candidate fold located inside half-separation'}
    c=coefficients(model,basis)
    result={'label':label,'included':True,'candidate_center_solved':bool(valid),'a':a,
            'coefficients_c_dw_dnorm_gainnorm2':list(c),'old_RMS_separation':float(2*a/np.sqrt(len(v))),
            'old_data_difference':float(np.linalg.norm(model.response(pair[1])[0]-model.response(pair[0])[0])),
            'moment_rules':{k:energy_moments(k) for k in patterns},'pair':pair.tolist(),'runs':[]}
    for eta in effective_gains:
        lam=None if eta==1 else 2*c[-1]*eta/(1-eta)
        nuisance=None if lam is None else (lam,1.)
        best=predict('opposite_edges',c,a,1,eta)[0]
        worst=max(predict('two_positive',c,a,s,eta)[0] for s in (-1,1))
        decision_band=np.sqrt(best*worst)
        for side,truth,start in ((-1,pair[0],pair[1]),(1,pair[1],pair[0])):
            @lru_cache(maxsize=256)
            def measure(delta):
                return model.response(truth,delta)[0]

            for pattern in patterns:
                offsets,power=PATTERNS[pattern]
                amp=np.sqrt(power)
                pred_b,pred_x=predict(pattern,c,a,side,eta)
                def full_fit(b):
                    ds=tuple(float(x*b) for x in offsets)
                    return Fit(model,ds,[measure(q) for q in ds],weights=amp,gain_calibration=nuisance)
                scale=max(1e-23,eta*(abs(c[0])*a*a)**2)
                def two_derivatives(z):
                    _,g,h,_,_=exact_profile(full_fit(z[1]*pred_b),basis,z[0])
                    return np.array([g,h])/scale
                sol=root(two_derivatives,[pred_x,1.],tol=2.e-8)
                if not np.all(np.isfinite(sol.x)):
                    raise RuntimeError((label,pattern,sol.message))
                exact_b=pred_b*sol.x[1]
                exact_x=sol.x[0]
                xf,g,h,strong_g,strong_h=exact_profile(full_fit(exact_b),basis,exact_x)
                ok=bool(np.linalg.norm(two_derivatives(sol.x))<1e-4 and b_positive_ratio(exact_b,pred_b)
                        and strong_h>0 and min(xf)>.5 and max(xf)<30 and abs(exact_x)<1.8)
                direction=np.eye(2)*.0003
                ja=np.column_stack([(two_derivatives(sol.x+e)-two_derivatives(sol.x-e))/.0006 for e in direction])
                # A simple double critical point needs nonzero d_b gradient and d_x Hessian.
                ok=bool(ok and abs(np.linalg.det(ja))>1e-5)
                joint=full_fit(decision_band)
                rec,count=joint.solve(start,limit=300)
                classification=pair[np.argmin([np.linalg.norm(joint.evaluate(q)[0]) for q in pair])]
                other_measurements=np.concatenate([w*(model.response(pair[1],x*decision_band)[0]-
                                  model.response(pair[0],x*decision_band)[0]) for x,w in zip(offsets,amp)])
                decision={'band':decision_band,'conditional_TRF_error_RMS':float(np.linalg.norm(rec-truth)/np.sqrt(len(v))),
                          'conditional_TRF_correct_branch':bool(np.linalg.norm(rec-truth)<np.linalg.norm(rec-start)),
                          'TRF_factorizations':count,'old_pair_score_error_RMS':float(np.linalg.norm(classification-truth)/np.sqrt(len(v))),
                          'fixed_pair_signal_separation':float(np.linalg.norm(other_measurements))}
                result['runs'].append(dict(eta=eta,gain_penalty=lam,side=side,pattern=pattern,
                       predicted_band=pred_b,predicted_x=pred_x,numerical_band=exact_b,collision_valid=ok,
                       scaled_equation_norm=float(np.linalg.norm(two_derivatives(sol.x))),
                       transverse_slope=float(ja[0,1]),third_derivative_scaled=float(ja[1,0]),
                       strong_hessian_min=strong_h,strong_gradient=strong_g,numerical_x=exact_x,
                       interior_margin=float(min(np.min(xf-.5),np.min(30-xf))),decision=decision))
    result['cost']={'wall_s':perf_counter()-wall,'factorizations_including_synthesis':model.factorizations}
    return result


def b_positive_ratio(exact,pred):
    return np.isfinite(exact) and .1 < exact/pred < 8


def main(omit_natural=False):
    rows=[]
    out=ROOT/'results/spectral_moments.json'
    def record(result):
        rows.append(result)
        out.write_text(json.dumps(rows,indent=2,allow_nan=False)+'\n')
        print(result['label'],result.get('included'),[(round(q['numerical_band'],8),q['collision_valid']) for q in result.get('runs',[])[:5]],flush=True)
    for a in (.18,.4,.72):
        model,_,pair=fold_case(a,1)
        record(investigate(model,pair,'fold_'+str(a),(1.,.25) if a==.4 else (1.,),list(PATTERNS),True))
    if not omit_natural:
        cases=json.loads((R10/'branch_holdout.json').read_text())['cases']
        with np.load(R10/'natural_inputs.npz') as data:
            for tag in data['case_ids'].tolist():
                candidates=data[tag+'/clean_witnesses']
                if len(candidates)<2:
                    record({'label':tag,'included':False,'reason':'fewer than two old numerical branches'})
                    continue
                model=FrequencyTM(cases[tag]['geometry'])
                record(investigate(model,candidates[:2],tag,(1.,),list(PATTERNS)[:3],True))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--omit-natural',action='store_true')
    main(parser.parse_args().omit_natural)

"""Prove that no choice of two global scalar multipliers can rescue a
single-source energy bound at the stated cap, on two historical full rings.
The lossy auxiliary current is only a witness against the relaxation.
"""
import json
from pathlib import Path
import sys
import time
import numpy as np
from scipy.linalg import eigh, solve
from flint import acb, arb, ctx
from threadpoolctl import threadpool_limits

from hidden_radiation import STUDY, SOURCE, load_operators, herm, star, ball, gershgorin, lower, upper, frobenius_squared

CASES={'r10_twod_0001': ('history__r12_dr__m',5.),'austria_eps4_400MHz':('history__cie_dual__m',1.)}

def witness(case, field, u, threshold):
    t0=time.perf_counter()
    a=load_operators(case)
    G,B,S=a['G'],a['E0'],a['S']
    with np.load(SOURCE/(case+'.npz')) as data:m,y=data[field],data['Y_archive']
    rbank=json.loads((SOURCE.parent/'observation_radii.json').read_text())[case]
    # This fixed small loss produces strict optical slack; only the entries of X
    # are used by the certificate. No material physics is assumed of this trial.
    loss=1e-6
    material=m.ravel()+1j*loss
    X=solve(np.eye(len(m))-material[:,None]*G,material[:,None]*B,assume_a='gen',check_finite=False)
    E=B+G@X; w=X.conj().T@E
    real=herm(X.conj().T@X)-u*herm(w); T=-herm(w/1j)
    r0=eigh(real,T,eigvals_only=True)[-1]
    t=max(1.,float(np.ceil(1.1*r0+1)))
    A,Bb,Sb,Yb,Xb=map(ball,(G,B,S,y,X))
    fieldb=Bb+A*Xb
    we=star(Xb)*fieldb
    Tb=(we-star(we))*acb(0,.5)
    realb=star(Xb)*Xb-arb(u)*(we+star(we))/2
    residual=upper(frobenius_squared(Sb*Xb-Yb).sqrt())
    results={}
    for name,candidate,actual in [('T',T,Tb),('tT_minus_real',t*T-real,arb(t)*Tb-realb)]:
        V=ball(eigh(herm(candidate))[1])
        results[name]=lower(gershgorin(star(V)*actual*V,upper=False))
        results[name+'_basis']=lower(gershgorin(star(V)*V,upper=False))
    positive=min(results.values())>0
    enabled=threshold>t and positive
    return dict(cap=u,prior_public_endpoint=field,loss_added_only_for_a_relaxation_witness=loss,
                residual_upper=residual,radii=rbank,ratio_high=t,dark_negative_for_ratio_below=threshold,
                arb_proof=results,status={name:('single_source_global_multipliers_exhausted' if enabled and residual<r else 'unresolved')
                                        for name,r in rbank.items()},elapsed_wall_seconds=time.perf_counter()-t0)

def run():
    ctx.prec=128;ctx.threads=1
    depths=json.loads((STUDY/'results/theory/hidden_radiation.json').read_text())['cases']
    report={c:witness(c,field,u,depths[c]['hidden_threshold']['minimum_ratio_lambda_over_mu_to_remove_this_vector'])
            for c,(field,u) in CASES.items()}
    (STUDY/'results/theory/all_scalar_multipliers.json').write_text(json.dumps(report,indent=2)+'\n')
    for k,v in report.items():print(k,v['residual_upper'],v['ratio_high'],v['status'],flush=True)
if __name__=='__main__':
    with threadpool_limits(1):run()

"""Independent numerical checks on stored TM operators and measurements."""

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "results"


def cofactor_check(H, R, B, host, m, ids, q, saved_clean=None):
    d = len(ids)-1
    left = np.linalg.svd(B, full_matrices=True)[0][:, -d:].conj().T
    U = np.eye(len(m))[:, ids] - H[:, ids]*q
    C = left @ U
    weights = np.array([(-1)**j*np.linalg.det(np.delete(C,j,axis=1)) for j in range(d+1)])
    K = U @ weights
    alpha = np.linalg.lstsq(B,K,rcond=None)[0]
    c = host @ alpha + (q*R[ids]) @ weights
    F = host + R @ np.linalg.solve(np.diag(1/m)-H,B)
    v = np.linalg.solve((np.eye(len(m))-m[:,None]*H).T,R)
    lhs = F @ alpha-c
    terms = weights*(m[ids]-q)*v[ids]
    rhs = np.sum(terms)
    scale = abs(F @ alpha)+abs(c)+np.sum(abs(terms))
    return dict(d=d, ids=ids.tolist(), q=q.tolist(), condition_B=np.linalg.cond(B),
                lhs_real_imag=[lhs.real,lhs.imag], rhs_real_imag=[rhs.real,rhs.imag],
                absolute_error=abs(lhs-rhs), cancellation_scale=scale,
                error_over_scale=abs(lhs-rhs)/scale,
                relative_nonzero_error=abs(lhs-rhs)/abs(rhs) if abs(rhs)>1e-7*scale else None,
                signed_null_relative=np.linalg.norm(C@weights)/np.linalg.norm(abs(C)@abs(weights)),
                synthesis_relative=np.linalg.norm(B@alpha-K)/np.linalg.norm(K),
                saved_clean_max_error=np.max(abs(F-saved_clean)) if saved_clean is not None else None)


def gain_check(H,R,B,host,m,y,eps,gamma):
    A = np.eye(len(m))-H*m
    F = host + (R*m) @ np.linalg.solve(A,B)
    M = np.column_stack((gamma*F,eps*np.eye(len(y))))
    r = y-F
    normalized_noise = np.linalg.lstsq(M,r,rcond=None)[0]
    ls = np.linalg.norm(normalized_noise)**2
    gram = np.vdot(r,np.linalg.solve(M@M.conj().T,r)).real
    D = np.linalg.det(A)
    Q = np.array([np.linalg.det(np.block([[A,b[:,None]],[-(R*m)[None,:],np.array([[h]])]]))
                  for b,h in zip(B.T,host)])
    c = (eps/gamma)**2
    wedge = np.triu(y[:,None]*Q[None,:]-Q[:,None]*y[None,:],1)
    psi = np.linalg.norm(wedge)**2+c*np.linalg.norm(D*y-Q)**2-eps**2*(np.linalg.norm(Q)**2+c*abs(D)**2)
    normalizer = eps**2*(np.linalg.norm(Q)**2+c*abs(D)**2)
    polynomial = 1+psi/normalizer
    return dict(ls=ls,gram=gram,polynomial=polynomial,psi=psi,normalizer=normalizer,
                error_over_score_scale=max(abs(ls-gram),abs(ls-polynomial))/max(1.,ls),
                feasible=[bool(ls<=1),bool(gram<=1),bool(psi<=0)])


natural = np.load(ROOT/"natural_inputs.npz")
designed = np.load(ROOT/"lift_k24.npz")
name = "fresh6-00"
H,R,m = (natural[f"{name}/{key}"] for key in ("H","R","truth"))
B2,host2,y = (natural[f"{name}/active/{key}"] for key in ("B","y_host","y"))
B3 = natural[f"{name}/B0"]
clean = (y-natural[f"{name}/additive_noise"])/(1+natural[f"{name}/gain_error"])
clean3 = clean[:3]
cases = [(designed["H"],designed["R"],designed["B_old"],designed["y_host"]@designed["p_old"],
          designed["pair_materials"][0],np.array([19,12]),
          designed["pair_old_data"][0]+designed["y_host"]@designed["p_old"]),
         (H,R,B2,host2,m,np.array([5,0,2]),clean),
         (H,R,B3,np.zeros(3,dtype=complex),m,np.array([5,0,2,4]),clean3)]
out = dict(cofactors=[],gain_elimination={})
for Hd,Rd,Bd,host,md,ids,saved in cases:
    for shift in (0.,1e-9,1.):
        q = md[ids]+shift*np.array([.4,-.3,.2,-.1])[:len(ids)]
        out["cofactors"].append(cofactor_check(Hd,Rd,Bd,host,md,ids,q,saved))
for key,z in (("truth",m),("old_alternate",natural[f"{name}/clean_witnesses"][1]),
              ("perturbed_truth",m+np.linspace(-.1,.1,len(m))),("box_center",np.full(len(m),15.25))):
    out["gain_elimination"][key] = gain_check(H,R,B2,host2,z,y,2e-4,.01)
assert max(t[key] for t in out["cofactors"]
           for key in ("error_over_scale","signed_null_relative","synthesis_relative","saved_clean_max_error")) < 1e-10
assert all(len(set(t["feasible"]))==1 and t["error_over_score_scale"]<1e-10
           for t in out["gain_elimination"].values())
(ROOT/"formula_checks.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")

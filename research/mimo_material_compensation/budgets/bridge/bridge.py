"""A fixed-mean ROI with an independent within-ROI contrast tolerance."""
from pathlib import Path
import sys
from fractions import Fraction
import json

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from bound import (CFG, Response, audited_state, block, down, elapsed, explicit_difference, fro,
                   load_case, now, np, optimize_direction, optimize_witness, rd, real, support, up, write_json)
from shape_mean import best_task_multiplier
from budget_proof import VerifiedWave
from scipy.optimize import minimize


def bound(state,ref,roi,k,j,ej,u,e,d,ed,lam,beta,alpha):
    ext=~roi
    lo=-state.x[ext]
    hi=np.nextafter(4.5-state.x[ext],np.inf)
    hi[state.x[ext]==4.5]=0.
    cap=np.minimum(np.maximum(-lo,hi),beta)
    dE=min(beta,float(np.max(cap)))
    r0,r1=np.sqrt(alpha),np.sqrt(dE)
    scales=np.array([r0,r1])
    # The material diagonal factors out on each side of K. The weighted 2x2
    # block-norm matrix bounds the norm on the full ROI plus its exterior.
    contraction=fro(k*np.outer(scales,scales))
    if contraction>=1:
        return dict(lower=0.,contraction_upper=contraction)
    scale=np.where(roi,r0,r1)
    wu,we=u*scale,e*scale[:,None]
    err_u=up(max(scales)*(state.u_error+rd(4)*fro(u)))
    err_e=up(max(scales)*(state.e_error+rd(4)*fro(e)))
    p=lam/(up(max(1.,fro(lam)))*(1.+rd(4)))
    if fro(p)>1:
        raise ArithmeticError("Direction exceeds unit norm")
    m=p[:64].reshape(8,8)+1j*p[64:].reshape(8,8)
    elu=up(fro(m.conj().T@wu)+err_u+rd(8)*fro(wu)*fro(m))
    ele=up(fro(we@m.conj().T)+err_e+rd(8)*fro(we)*fro(m))
    remainder=up(min(elu*(fro(we)+err_e),ele*(fro(wu)+err_u))*contraction/(1.-contraction))
    q=j.T@p
    lower_phase,theta=best_task_multiplier(q[roi], np.full(roi.sum(),-min(alpha,float(state.x[roi][0]))),
                            np.full(roi.sum(),min(alpha,4.5-float(state.x[roi][0]))),0.,False)
    s,_,_,_=support(-q[ext],lo,hi,beta,"max")
    d2=up(np.hypot(fro(cap),alpha*np.sqrt(roi.sum())))
    dq=up(d2*(ej+rd(len(p))*fro(j)*fro(p)))
    center=down(float(p@d)-ed-rd(len(p))*fro(p)*fro(d))
    lower=max(0.,down(center+lower_phase-s-dq-remainder-
                     rd(8)*(abs(center)+abs(lower_phase)+s+dq+remainder)))
    return dict(lower=lower,contraction_upper=contraction,roi_support_upper=-lower_phase,
                exterior_support_upper=s,directional_remainder_upper=remainder,
                derivative_error_upper=dq,signal_projection_lower=center,
                task_multiplier=theta,lambda_norm_upper=fro(p))


def exact_roi(v,t,alpha):
    """Repair an SLSQP proposal to a binary64 material of exactly mean t."""
    # A dyadic contraction leaves room for the rounding-size mean repair.
    v=t+(v-t)*(1.-1e-8)
    target=len(v)*Fraction(t)
    val=sum(Fraction(x) for x in v)
    gap=target-val
    for i in np.argsort(v):
        exact=Fraction(v[i])+gap
        candidate=float(exact)
        if Fraction(candidate)==exact and 0<=candidate<=4.5 and abs(candidate-t)<alpha:
            v[i]=candidate
            if sum(Fraction(x) for x in v)==target:
                return v
    raise ArithmeticError("Could not encode exact registered mean")


def propose_shape(physical,roi,c,target,r,t,alpha):
    lay=CFG["layout"]["B"]
    def objective(p):
        x=c.copy()
        x[roi]=p
        z=physical.full(x,derivatives=True)
        d=(block(z["p"],"B")-target)/(2*r)
        u=z["u"][lay["rx"]][:,roi]
        e=z["e"][roi][:,lay["tx"]]
        g=np.einsum("ij,ik,kj->k",d.conj(),u,e,optimize=True).real/(2*r)
        return .5*np.linalg.norm(d)**2,g
    lo,hi=max(0.,t-alpha),min(4.5,t+alpha)
    s=minimize(objective,np.full(roi.sum(),t),jac=True,method="SLSQP",bounds=[(lo,hi)]*int(roi.sum()),
         constraints=[dict(type="eq",fun=lambda x:np.sum(x)-t*len(x),jac=lambda x:np.ones(len(x)))],
         options=dict(maxiter=50,ftol=1e-11))
    x=c.copy()
    x[roi]=exact_roi(s.x,t,alpha)
    return x,int(s.nfev)


def main():
    start=now()
    t=0.27
    side,base,exc,rec,roi,c,r,obs,observed=load_case("joint32")
    rec=rec[:16]
    ext=~roi
    physical=Response(side,base,exc,rec)
    wave=VerifiedWave(dict(base,S=rec,E0=exc))
    ref=audited_state(wave,physical,rec,exc,c)
    state=audited_state(wave,physical,rec,exc,c+t*roi,derivatives=True)
    j,ej,u,e=state.jacobian(CFG["layout"]["B"])
    d,ed=explicit_difference(block(state.p,"B"),block(ref.p,"B"))
    ed=up(ed+state.response_error+ref.response_error)
    ks=np.array([[state.k_upper(a,b,dense=True) for b in (roi,ext)] for a in (roi,ext)])
    focus=json.loads((HERE.parent/"budget_focus_joint32.json").read_text())
    budgets=[0.,0.001,0.002,0.0025,0.003]
    alphas=[0.,0.0001,0.0003,0.0005,0.001,0.002,0.003,0.005,0.0075,0.01,0.015,0.02,0.03,0.04,
            0.05,0.07,0.10,0.13,0.16,0.20,0.23,0.27,0.32,0.40,0.55,0.70,0.85,1.]
    candidate_alphas=[.0075,.015,.03,.05,.07,.10,.13,.16,.20,.23,.27,.32,.40,.55,.70,.85]
    results=[]
    arrays=dict(reference=c,roi=roi,p_reference=ref.p,anchor=state.x,anchor_p=state.p,radius=r)
    ref_block=block(ref.p,"B")
    proposals={}
    for bi,beta in enumerate(budgets):
        query=next(rec for rec in focus["records"] if
                   (rec["t"],rec["layout"],rec["beta"])==(t,"B",beta))
        with np.load(HERE.parent/"raw_focus_joint32.npz",allow_pickle=False) as raw:
            lam=raw[query["dual_key"]]
        certs=[]
        for alpha in alphas:
            info=bound(state,ref,roi,ks,j,ej,u,e,d,ed,lam,beta,alpha)
            info.update(alpha=alpha,certified_lower_over_two_r=down(info["lower"]/(2*r)) if info["lower"]>0 else 0.)
            certs.append(info)
        last=max((q for q in certs if q["certified_lower_over_two_r"]>1),key=lambda q:q["alpha"])
        # Using the same feasible direction and monotonicity, refine the
        # largest strictly certified contrast tolerance around the grid bracket.
        upper_alpha=min(a for a in alphas if a>last["alpha"])
        lower_alpha=last["alpha"]
        for step in range(22):
            mid=(lower_alpha+upper_alpha)*.5
            info=bound(state,ref,roi,ks,j,ej,u,e,d,ed,lam,beta,mid)
            if info["lower"]>2*r:
                lower_alpha=mid
            else:
                upper_alpha=mid
        # Round the usable interval endpoint down by 10 micro-contrast.
        published_alpha=max(0.,np.floor(lower_alpha*1e5)/1e5-1e-5)
        cert=bound(state,ref,roi,ks,j,ej,u,e,d,ed,lam,beta,published_alpha)
        arrays[f"b{bi:02d}_lambda"]=lam
        best=None
        attempts=[]
        for alpha in candidate_alphas:
            if alpha<=published_alpha:
                continue
            if alpha not in proposals:
                proposals[alpha]=propose_shape(physical,roi,c,ref_block,r,t,alpha)
            x0,nfev=proposals[alpha]
            x=x0
            if beta>0:
                p0=physical.full(x0,derivatives=True)
                dz=real(block(p0["p"],"B")-ref_block)
                layout=CFG["layout"]["B"]
                u0=p0["u"][layout["rx"]][:,ext]
                e0=p0["e"][ext][:,layout["tx"]]
                jp=(u0[:,None,:]*e0.T[None,:,:]).reshape(64,-1)
                jr=np.concatenate((jp.real,jp.imag),axis=0)
                lo,hi=-x0[ext],4.5-x0[ext]
                _,z=optimize_direction(dz,jr,lo,hi,beta,"max")
                x=optimize_witness(physical,x0,ref_block,"B",z,ext,lo,hi,beta,"max",2*r,maxiter=8)
            pred=physical.full(x)["p"]
            proposal=float(np.linalg.norm(block(pred,"B")-ref_block)/(2*r))
            attempts.append(dict(alpha=alpha,proposal=proposal,roi_objective_evaluations=nfev))
            if proposal>=.9998:
                continue
            assert sum(Fraction(y) for y in x[roi])==int(roi.sum())*Fraction(t)
            da=up(np.max(abs(x[roi]-t)))
            db=up(np.max(abs((x-c)[ext]))) if beta else 0.
            assert da<=alpha and db<=beta and np.all((x>=0)&(x<=4.5))
            w=audited_state(wave,physical,rec,exc,x)
            a,ar=explicit_difference(block(w.p,"B"),ref_block)
            wu=up((fro(a)+ar+w.response_error+ref.response_error)/(2*r))
            if wu<1:
                ampl=f"b{bi:02d}_witness"
                arrays[ampl+"_x"]=x
                arrays[ampl+"_p"]=w.p
                arrays[ampl+"_j"]=w.j
                best=dict(alpha_upper=da,beta_upper=db,strict_W_over_two_r=wu,raw_prefix=ampl,
                          roi_values=x[roi].tolist(),exact_mean_as_ratio=t.as_integer_ratio(),
                          weighted_L1=float(np.sum(abs((x-c)[ext]))*(2/side)**2),
                          weighted_L2=float(np.linalg.norm((x-c)[ext])*(2/side)))
                break
        results.append(dict(beta=beta,certified_alpha_through=published_alpha,
                            fixed_direction_bisection=[lower_alpha,upper_alpha],
                            B_lower_over_two_r=down(cert["lower"]/(2*r)),certified_endpoint=cert,
                            witness=best,attempts=attempts,grid=certs))
        print(beta,published_alpha,results[-1]["B_lower_over_two_r"],best,flush=True)
    with (HERE/"raw.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    write_json(HERE/"receipt.json",dict(t=t,case=obs,r=r,cell_count=int(roi.sum()),
               k_block_upper=ks.tolist(),region="fixed mean plus a pixelwise contrast tolerance",
               records=results,cost=dict(physical.model.cost),**elapsed(start)))


if __name__=="__main__":
    main()

"""Cold-start convex QCQP duals for all preregistered public structural nodes."""
import argparse
import json
from pathlib import Path
import signal
import sys
from time import perf_counter
import warnings

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from common import PUBLIC, RESULTS, RATIOS, digest, herm, load, save, star, write_json
from bounds import (Acquisition, arb, ball, bt, ctx, high, low, midpoint,
                    propose_frame, select, validate_frame)
from baseline_psd import strict_schur
from known_zero_control import check as inertia_known_zero
from replay import original
from run import ArrayStore
import cvxpy as cp
import numpy as np
from scipy.linalg import eigh, hadamard, solve

ctx.prec=192
ctx.threads=1
RATIO=3.
BUDGET=9.4  # reserve 0.6 s for signal latency and saved-proof replay
STAGE_ORDER=("spatial_psd","signed_local_qcqp","source_diagonal_qcqp","source_hadamard_matrix_qcqp")
SHRINK=(1.,.999,.99,.9,.7,.5,.3,.1,0.)


class NodeDeadline(Exception):
    pass


def time_limit(_signum,_frame):
    raise NodeDeadline("10 second node budget")


def select_public():
    arrays=load(PUBLIC)
    selection=[]
    for k,pair in enumerate(arrays["pair_ids"]):
        if arrays["families"][k]!="structure":
            continue
        frozen=json.loads((RESULTS/"raw"/f"pair_{k+1:02d}.json").read_text())
        for row in frozen["rows"]:
            target,bg,k_width,slot=row["box_positions"]
            if bg==0 and RATIOS[k_width]==RATIO:
                selection.append((str(pair),str(row["target_id"]),slot,k,row))
    return arrays,sorted(selection,key=lambda entry:entry[:3])


def forms(G,S,B,Y,ell,upp):
    n,s=B.shape
    j=np.hstack((np.eye(n),np.zeros((n,s))))
    e=np.hstack((G,B))
    o=np.hstack((S,-Y))
    a0=herm(star(o)@o)
    m=(ell+upp)/2
    h=(upp-ell)/2
    diag,optic=[],[]
    for i in range(n):
        q=j[i]-m[i]*e[i]
        diag.append(np.outer(q.conj(),q)-h[i]**2*np.outer(e[i].conj(),e[i]))
        cross=-.5j*np.outer(j[i].conj(),e[i])
        optic.append(cross+star(cross))
    return a0,np.array(diag),np.array(optic)


def cone_expression(basis,W):
    n,q,_=basis.shape
    return cp.reshape(cp.Constant(basis.reshape(n,q*q).T)@W,(q,q),order="C")


def solve_cone(acq,ell,upp,variant,expiry):
    start=perf_counter()
    if variant=="source_diagonal_qcqp":
        U=np.eye(acq.s)
    elif variant=="source_hadamard_matrix_qcqp":
        U=hadamard(acq.s)/np.sqrt(acq.s)
    else:
        U=None
    matrices=[]
    if U is None:
        matrices.append(forms(acq.G,acq.S,acq.B,acq.Y,ell,upp))
    else:
        for q in range(acq.s):
            u=U[:,q:q+1]
            matrices.append(forms(acq.G,acq.S,acq.B@u,acq.Y@u,ell,upp))
    scales=[]
    variables=[]
    objective=[]
    constraints=[]
    for a0,D,L in matrices:
        s=a0.shape[0]-acq.n
        unit=max(1.e-4,np.linalg.norm(a0,ord=2))
        scales.append(float(unit))
        w=cp.Variable(acq.n,nonneg=True)
        opt=None if variant=="spatial_psd" else cp.Variable(acq.n)
        coeff = a0/unit+cone_expression(D,w)
        if opt is not None:
            coeff+=cone_expression(L,opt)
        theta=cp.Variable((s,s),hermitian=True)
        zero=np.zeros((acq.n,acq.n))
        Z_s=np.zeros((acq.n,s))
        offset=cp.bmat([[zero,Z_s],[Z_s.T,theta]])
        # Bound search coefficients uniformly and reserve a small state gap.
        # This is a stricter proposal cone. Replays use the unconstrained
        # mathematical dual and only require original K to be positive.
        guard = np.diag(np.r_[np.ones(acq.n)*1.e-7,np.zeros(s)])
        constraints.extend((coeff-offset-guard >> 0, w <= 64.))
        if opt is not None:
            constraints.extend((opt <= 64., opt >= -64.))
        objective.append(unit*cp.real(cp.trace(theta)))
        variables.append((w,opt))
    problem=cp.Problem(cp.Maximize(sum(objective)),constraints)
    build=perf_counter()-start
    deadline=expiry-perf_counter()
    if deadline<=0.:
        raise NodeDeadline("graph preparation reached node budget")
    solver_start=perf_counter()
    proposal_history=[]
    proposals = [[] for _ in variables]
    with warnings.catch_warnings(record=True) as issued:
        warnings.simplefilter("always")
        for method in ("CLARABEL","SCS"):
            left = expiry-perf_counter()
            if left<=.05:
                break
            try:
                if method=="CLARABEL":
                    optimum=problem.solve(solver=method,warm_start=False,time_limit=max(.03,.5*left),
                                          max_iter=100,tol_gap_abs=1.e-7,tol_gap_rel=1.e-6,
                                          tol_feas=1.e-7)
                else:
                    optimum=problem.solve(solver=method,warm_start=False,time_limit_secs=max(.02,left-.04),
                                          max_iters=15000,eps=5.e-5)
                receipt=dict(solver=method,status=str(problem.status),
                             engine_seconds=problem.solver_stats.solve_time,
                             iterations=problem.solver_stats.num_iters,
                             objective=None if optimum is None or not np.isfinite(optimum) else float(optimum))
                proposal_history.append(receipt)
                for options,(w,o),scale in zip(proposals,variables,scales,strict=True):
                    if w.value is not None and (o is None or o.value is not None):
                        weights=np.maximum(np.asarray(w.value)*scale,0.)
                        optical=np.zeros(acq.n) if o is None else np.asarray(o.value)*scale
                        options.append((weights,optical,method))
                if problem.status=="optimal":
                    break
            except cp.SolverError as exc:
                proposal_history.append(dict(solver=method,status="optimizer_error",message=str(exc)))
    return dict(status=proposal_history[-1]["status"] if proposal_history else "no solver allowance",
                build_seconds=build,solver_wall_seconds=perf_counter()-solver_start,
                solver_receipts=proposal_history,
                warnings=sorted(set(str(w.message) for w in issued))),proposals,U


def project(H,idx,n_sources,basis):
    if idx is None:
        return sum((H[i,i].real for i in range(n_sources)),arb(0))
    source=np.eye(n_sources) if basis=="I" else hadamard(n_sources)
    u=ball(source[:,idx:idx+1])
    return (bt(u)*H*u)[0,0].real/(n_sources if basis=="H" else 1)


def selected_modes(acq,ell,upp,variant,candidates,store,expiry):
    lower=arb(0)
    terms=[]
    basis="H" if variant=="source_hadamard_matrix_qcqp" else "I"
    for k,parameters in enumerate(candidates):
        if perf_counter()>=expiry:
            break
        mode=k if variant.startswith("source_") else None
        mains=[]
        A_base,D,L=forms(acq.G,acq.S,acq.B,acq.Y,ell,upp)
        for w,o,method in parameters:
            for shrink in SHRINK:
                wr,lr=w*shrink,o*shrink
                kernel=A_base+np.einsum("k,kij->ij",wr,D)+np.einsum("k,kij->ij",lr,L)
                K=herm(kernel[:acq.n,:acq.n])
                ev=eigh(K,eigvals_only=True,check_finite=False)
                if ev[0]<=1.e-13:
                    continue
                ab=kernel[:acq.n,acq.n:]
                H=herm(kernel[acq.n:,acq.n:]-star(ab)@solve(K,ab,assume_a="pos",check_finite=False))
                if mode is None:
                    hint=float(np.trace(H).real)
                else:
                    ve=np.eye(acq.s)[:,mode] if basis=="I" else hadamard(acq.s)[:,mode]/np.sqrt(acq.s)
                    hint=float(np.vdot(ve,H@ve).real)
                mains.append((hint,shrink,wr,lr,K,method))
        for hint,shrink,wr,lr,K,method in sorted(mains,key=lambda q:q[0],reverse=True)[:3]:
            if perf_counter()>=expiry:
                break
            frame,edge=propose_frame(K,acq.n)
            if frame is None:
                continue
            try:
                H,cert=strict_schur(acq,ell,upp,wr,lr,frame)
            except (ZeroDivisionError,ValueError):
                continue
            if H is None:
                continue
            value=project(H,mode,acq.s,basis)
            lower+=max(arb(0),value.lower())
            terms.append(dict(weights=wr.tolist(),optical=lr.tolist(),S_space=basis,source_mode=mode,
                              strict_value_lower=low(value),strict_value_upper=high(value),
                              shrink=shrink,proposal_solver=method,frame=store.add(frame),positive_state_proof=cert,
                              K_smallest_proposed_eigenvalue=edge))
            break
    return dict(strict_squared_residual_lower=max(0.,low(lower)),
                strict_residual_over_radius_lower=max(0.,low(lower.sqrt()/arb(acq.radius))),
                excluded=bool(lower>arb(acq.radius)**2),proof_terms=terms)


def run_node(B,ell,upp,frames):
    start=perf_counter()
    expiry=start+BUDGET
    signal.signal(signal.SIGALRM,time_limit)
    rows={}
    for variant in STAGE_ORDER:
        left=expiry-perf_counter()
        if left<=0.:
            rows[variant]=dict(status="unattempted_wall_cap",excluded=False)
            continue
        # Split the remaining allowance to give each announced baseline a try.
        budget=left/(4-len(rows))
        phase_expiry=perf_counter()+max(.01,budget-.10)
        signal.setitimer(signal.ITIMER_REAL,left)
        try:
            proposal,params,_=solve_cone(B,ell,upp,variant,phase_expiry)
            t0=perf_counter()
            strict=selected_modes(B,ell,upp,variant,params,frames,expiry)
            rows[variant]=dict(proposal=proposal,certificate_seconds=perf_counter()-t0,**strict)
        except NodeDeadline as exc:
            rows[variant]=dict(status="wall_cap",reason=str(exc),excluded=False)
        finally:
            signal.setitimer(signal.ITIMER_REAL,0.)
    return rows,perf_counter()-start


def independent(acq,ell,upp,variants,frames):
    """Reassemble raw local 2x2 blocks, then verify the saved modal positive frames."""
    replays={}
    for variant,report in variants.items():
        lower=arb(0)
        t0=perf_counter()
        for term in report.get("proof_terms",[]):
            nu=acq.n
            W,L=np.array(term["weights"]),np.array(term["optical"])
            if np.any(W<0):
                raise ArithmeticError("nonnegative PSD weight requirement failed replay")
            full=original(acq,ell,upp,W,L)
            idx,src=range(nu),range(nu,nu+acq.s)
            K=select(full,idx,idx)
            frame=frames[term["frame"]]
            if frame.shape!=(nu,nu) or not validate_frame(K,frame)[0]>0:
                raise ArithmeticError("proposed positive PSD-dual state frame failed replay")
            couple=select(full,idx,src)
            H=select(full,src,src)-bt(couple)*K.solve(couple)
            value=project(H,term["source_mode"],acq.s,term["S_space"])
            lower+=max(arb(0),value.lower())
        excludes=bool(lower>arb(acq.radius)**2)
        if report["excluded"]!=excludes:
            raise ArithmeticError("frozen PSD decision failed independent reassembly")
        replays[variant]=dict(excluded=excludes,strict_squared_residual_lower=max(0.,low(lower)),
                              elapsed_seconds=perf_counter()-t0)
    return replays


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--begin",type=int,default=0)
    parser.add_argument("--end",type=int,default=32)
    args = parser.parse_args()
    acq,targets=select_public()
    assert len(targets)==32
    out=HERE/"nodes"
    out.mkdir(exist_ok=True)
    for pos,(pair,target,slot,k,row) in enumerate(targets):
        if not args.begin<=pos<args.end:
            continue
        location=out/f"{pos:02d}.json"
        if location.exists():
            continue
        before=perf_counter()
        t,b,width,_=row["box_positions"]
        ell,upp=acq["boxes"][k,b,width,slot]
        keep=np.flatnonzero(upp>0.)
        shared=Acquisition(acq["G"][np.ix_(keep,keep)],acq["S"][:,keep],acq["B"][keep],
                           acq["Y"][k,t,b],float(acq["radii"][k,t,b]))
        setup=perf_counter()-before
        st=ArrayStore()
        s=perf_counter()
        inertia=inertia_known_zero(acq | dict(Y_block=shared.Y,radius=shared.radius),(ell,upp),row,st)
        inertia_online_and_replay=perf_counter()-s
        variants,wall=run_node(shared,ell[keep],upp[keep],st)
        replays=independent(shared,ell[keep],upp[keep],variants,st.arrays)
        eliminated=np.flatnonzero(upp==0.).tolist()
        witness=row["witness"]["audit"]["proved"]
        exclusion=any(v["excluded"] for v in variants.values())
        if exclusion and witness:
            raise ArithmeticError(f"independent PSD search excluded a certified known witness at {pos}")
        status="EXCLUDED" if exclusion else "WITNESS_FOUND" if witness else "UNRESOLVED"
        payload=dict(pair_id=pair,target_id=target,box_slot=slot,remaining_state_rows=keep.tolist(),
                     eliminated_known_zero_rows=eliminated,state_dimension=len(keep),
                     condition_id=row["condition_id"],box_positions=row["box_positions"],
                     C_known_zero=inertia,old_B_status=row["status"]["B"],B_independent_state=status,
                     shared_setup_seconds=setup,B_end_to_end_online_seconds=wall,
                     C_total_proposal_and_replay_seconds=inertia_online_and_replay,
                     B_independent_replay=replays,B_variants=variants)
        save(out/f"{pos:02d}_frames.npz",**st.arrays)
        write_json(location,payload)
        print(pos,pair,target,slot,status,inertia["C_known_zero"],f"{wall:.3f}s",
              [v for v in variants if variants[v]["excluded"]],flush=True)


if __name__=="__main__":
    main()

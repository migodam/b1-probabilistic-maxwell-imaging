"""Verify two reference material neighborhoods using the nominal native arrays."""
from pathlib import Path
import argparse
import json
import os
import sys

for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
    os.environ[variable]="1"
sys.dont_write_bytecode=True

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"budgets"))
from bound import explicit_difference, directional, down, elapsed, fro, now, np, rd, up, write_json
from budget_proof import VerifiedState, VerifiedWave, opnorm, norm_lower
from passive_replay import load_cover, singular_upper
from scipy.ndimage import label
from scipy.optimize import minimize
from scipy.linalg import cholesky
from quick import STUDY

GRID=(0.,0.0001,0.001,0.003,0.01,0.03)
NATIVE_REFINEMENT=(.011,.012,.013)
HARDWARE=dict(shots=8,complex_reads=64,available_W=1.,source_and_termination_ohms=50,
              hold_ports=list(range(16,24)), tx=list(range(8)), rx=list(range(8,16)))
LAYOUT=dict(tx=list(range(8)),rx=list(range(8)))


class AuditedState(VerifiedState):
    """Reuse the interval state audit, with a tight dense norm on the small grid."""
    def __init__(self,wave,rec,exc,x,require_k=True):
        self.require_k=require_k
        super().__init__(wave,rec,exc,x,full_k=len(x)==1024)

    def k_upper(self,rows,cols,dense):
        if not self.require_k:
            return 0.
        left,right=self.cross_left[rows],self.cross_right[:,cols]
        el,er,ew=self.cross_left_error,self.cross_right_error,self.weighted_error
        if not dense:
            return self.compressed_native_norm(left,right,el,er,ew)
        direct,ed=self.green_block(np.flatnonzero(rows),np.flatnonzero(cols))
        d=len(self.active)
        wright=self.weighted@right
        ewr=up(ew*(opnorm(right)+er)+opnorm(self.weighted)*er+rd(d)*fro(self.weighted)*fro(right))
        correction=left@wright
        ec=up(el*(opnorm(wright)+ewr)+opnorm(left)*ewr+rd(d)*fro(left)*fro(wright))
        self.k_matrix=direct+correction
        self.k_matrix_error=up(ed+ec+rd(2)*(fro(direct)+fro(correction)))
        return up(singular_upper(self.k_matrix)+self.k_matrix_error)

    def compressed_native_norm(self,a,b,ea,eb,ew):
        # Numerically compressed correction K-G = A W B on all native pixels.
        # Hermitian Gram majorization verifies ||A R^-1||<=1. Reciprocity
        # residual is bounded explicitly before reusing R on the right.
        n,d=a.shape
        a_norm=opnorm(a)
        b_norm=opnorm(b)
        w=self.weighted
        w_norm=opnorm(w)
        ah_a=a.conj().T@a
        gram=(ah_a+ah_a.conj().T)*.5
        round_gram=up(rd(n)*(fro(a)**2+fro(ah_a))+rd(4)*fro(ah_a))
        shift=up(4.*(round_gram+rd(d)*(fro(gram)+1.)))
        majorant=gram.copy()
        majorant.flat[::d+1]+=shift
        r=cholesky(majorant,lower=False,check_finite=False)
        defect=up(fro(r.conj().T@r-majorant)+rd(d+3)*(fro(r)**2+fro(majorant))+
                  rd(2)*(fro(gram)+shift*np.sqrt(d)))
        if not round_gram+defect<shift:
            raise ArithmeticError("Failed Gram domination check")
        wr=w@r.T
        central=r@wr
        central_err=up(rd(d)*(fro(r)*fro(wr)+opnorm(r)*fro(w)*fro(r)))
        central_upper=up(min(opnorm(central),singular_upper(central))+central_err)
        reciprocity_error=up(fro(b-a.T)+rd(3)*(fro(a)+fro(b)))
        correction_error=up(ea*(w_norm+ew)*(b_norm+eb)+a_norm*ew*(b_norm+eb)+
                            a_norm*w_norm*eb+a_norm*w_norm*reciprocity_error)
        self.norm_audit=dict(gram_majorant_shift=shift,gram_roundoff_upper=round_gram,
                             cholesky_defect_upper=defect,reciprocity_difference_upper=reciprocity_error,
                             compressed_correction_norm_upper=central_upper,correction_input_error_upper=correction_error,
                             bare_plus_load_green_norm_upper=self.wave.gnorm)
        return up(self.wave.gnorm+central_upper+correction_error)


def material(side):
    t=2.0 if side==32 else 1.0
    source=HERE/f"phase_{side}_two_shells.npz"
    with np.load(source,allow_pickle=False) as raw:
        ix=raw["ix"]
        group=raw["group"]
        plus=np.zeros(side**2)
        minus=np.zeros(side**2)
        roi=np.zeros(side**2,bool)
        plus[ix]=raw[f"plus_{t}"]
        minus[ix]=raw[f"minus_{t}"]
        roi[ix[group==0]]=True
        radius=float(raw["r"])
    with np.load(STUDY/"ports/public/cases"/f"n{side}_c01.npz",allow_pickle=False) as raw:
        if radius!=float(raw["radius"]):
            raise ValueError("Phase packet disagrees with nominated radius")
        obs_id=str(raw["observation_id"].item())
    if not np.all(plus[roi]==t) or not np.all(minus[roi]==0.) or np.min(minus)<0 or np.max(plus)>4.5:
        raise ValueError("Unexpected material state")
    if label(roi.reshape(side,side))[1]!=1:
        raise ValueError("ROI should be connected")
    return t,source.name,radius,obs_id,ix,roi,plus,minus


def released_operator(side, frequency):
    name=f"n{side}_f{frequency}_a08"+("" if frequency==400 else "_c048")
    with np.load(STUDY/"ports/public/operators"/(name+".npz"),allow_pickle=False) as f:
        return name,{key:f[key] for key in f.files if key!="G"}


def bounds(x,beta):
    lo=np.maximum(-beta,-x)
    hi=np.minimum(beta,np.nextafter(4.5-x,np.inf))
    hi[x==4.5]=0.
    return lo,hi


def choose_direction(d,j,first,second,beta):
    lo0,hi0=bounds(first,beta)
    lo1,hi1=bounds(second,beta)
    # A realified center difference and an approximate least-squares residual
    # supply arbitrary directions. Only the verifier decides their validity.
    u=d/np.linalg.norm(d)
    if beta==0:
        return u,0
    j0,j1=j
    s=np.linalg.norm(d)
    n=len(first)
    def evaluate(v):
        p=(d+j0@v[:n]-j1@v[n:])/s
        adj=np.r_[j0.T@p,-j1.T@p]
        return .5*(p@p),adj/s
    opt=minimize(evaluate,np.zeros(2*n),jac=True,method="L-BFGS-B",
                 bounds=list(zip(np.r_[lo0,lo1],np.r_[hi0,hi1])),options=dict(maxiter=25,ftol=1e-11,gtol=1e-9))
    resid=d+j0@opt.x[:n]-j1@opt.x[n:]
    direction=resid/np.linalg.norm(resid) if np.linalg.norm(resid)>0 else u
    return direction,int(opt.nfev)


def check_direction(d,ed,lam,states,jacs,beta):
    scale=up(max(1.,fro(lam)))
    q=lam/(scale*(1.+rd(4)))
    if fro(q)>1:
        raise ArithmeticError("Unverified direction norm")
    # directional() applies precisely the same normalization to lam once.
    terms=[]
    for sign,state,(j,ej,u,e) in zip((1.,-1.),states,jacs):
        x=state.x
        lo=-x
        hi=np.nextafter(4.5-x,np.inf)
        hi[x==4.5]=0.
        info,checked=directional(np.zeros_like(d),lam,sign*j,ej,u,e,state.u_error,
                                 state.e_error,state.k_total,beta,"max",lo,hi,0.)
        if not np.array_equal(checked,q):
            raise RuntimeError("Mismatched verification directions")
        terms.append({k:info[k] for k in ("safe_support","derivative_term_error","nonlinear","k_times_max_edit")})
    projected=down(float(q@d)-ed-rd(len(q))*fro(q)*fro(d))
    if any(term["nonlinear"] is None for term in terms):
        return 0.,q,projected,terms
    penalty=up(sum(term["safe_support"]+term["derivative_term_error"]+term["nonlinear"] for term in terms))
    value=max(0.,down(projected-penalty-rd(8)*(abs(projected)+penalty)))
    return value,q,projected,terms


def parity_to_open(side,t,band,pair,states):
    prior_file=HERE/f"verified_responses_{side}_{t}_phase_{side}_two_shells_{band}.npz"
    with np.load(prior_file,allow_pickle=False) as prior:
        values=[(prior["s1"],float(prior["e1"])),(prior["s0"],float(prior["e0"]))]
    results=[]
    for (wave,e),state in zip(values,states):
        open_value,open_error,_=load_cover(wave,e,1.,0.,0.)
        numerical=open_value[8:16,:8]/np.sqrt(8)
        independent_error=up(open_error/np.sqrt(8)+rd(4)*fro(numerical))
        # Differences of total responses are background-invariant. For each
        # center the released stored-state readout adds the same empty total.
        results.append((numerical,independent_error))
    empty=pair["M_nominal"]
    empty_scatter=np.eye(len(empty))-2*float(pair["reference_ohms"])*np.linalg.solve(empty,np.eye(len(empty)))
    empty_scaled=empty_scatter[8:16,:8]/np.sqrt(8)
    residuals=[]
    for (reference,e),state in zip(results,states):
        delta=reference-(state.p+empty_scaled)
        observed=fro(delta)
        slack=up(e+state.response_error+rd(32)*(fro(reference)+fro(state.p)+fro(empty_scaled)))
        residuals.append(dict(open_direct_difference_norm=observed, enclosure_compatibility_slack=slack,
                              within=observed<=slack))
    open_delta=results[0][0]-results[1][0]
    open_err=up(results[0][1]+results[1][1]+rd(3)*(fro(results[0][0])+fro(results[1][0])))
    return residuals,open_delta,open_err,empty_scaled


def verify_midpoint(a,b,e0,e1,r):
    midpoint=.5*(a+b)
    guards=[]
    for p,e in ((a,e0),(b,e1)):
        inc=up(fro(p-midpoint)+e+rd(3)*(fro(p)+fro(midpoint)))
        guards.append(up(inc/r))
    return midpoint,guards


def run(side,bands):
    start=now()
    t,name,r,obs_id,ix,roi,plus,minus=material(side)
    total=dict(side=side, reference_packet=name, radius=r, radius_source_id=obs_id,
               observation_access="Only observation_id and radius read. No object measurement or material truth.",
               contract="Both candidate neighborhoods have independent all-pixel maximum budget, intersect [0,4.5]. Single fixed A layout. Shared chosen frequency and real non-dispersive contrast assumed; equal hard Frobenius readout radius.",
               hardware=HARDWARE,core_cells=int(roi.sum()),core_area_m2=float(roi.sum())*(2/side)**2,
               nominal_core_mean_gap=t, frequency={})
    arrays=dict(plus=plus,minus=minus,core=roi,active_union=ix,radius=r)
    for band in bands:
        timer=now()
        oid,op=released_operator(side,band)
        rec,exc=op["S"],op["E0"]
        wave=VerifiedWave(op)
        a=AuditedState(wave,rec,exc,plus,require_k=band==460)
        b=AuditedState(wave,rec,exc,minus,require_k=band==460)
        d,rounding=explicit_difference(a.p,b.p)
        ed=up(rounding+a.response_error+b.response_error)
        d_lower=max(0.,down(norm_lower(d)-ed))
        d_upper=up(fro(d)+ed)
        parity,open_d,open_e,empty=parity_to_open(side,t,band,op,(a,b))
        if not all(q["within"] for q in parity):
            raise ArithmeticError(f"Released operator and independent open replay disagree: {parity}")
        result=dict(operator=oid, lower_center_distance_over_two_r=down(d_lower/(2*r)),
                    upper_center_distance_over_two_r=up(d_upper/(2*r)),
                    independent_open_band_center_distance_over_two_r=float(np.linalg.norm(open_d)/(2*r)),
                    independent_open_error_over_two_r=up(open_e/(2*r)), open_vs_frozen=parity,
                    plus=dict(a.log,native_norm_audit=getattr(a,"norm_audit",None)),
                    minus=dict(b.log,native_norm_audit=getattr(b,"norm_audit",None)))
        arrays.update({f"p_{band}_plus":a.p,f"p_{band}_minus":b.p,f"p_{band}_empty_total":empty,
                       f"currents_{band}_plus":a.j,f"currents_{band}_minus":b.j,
                       f"delta_{band}_independent_open":open_d})
        if band==400:
            mid,rel=verify_midpoint(a.p,b.p,a.response_error,b.response_error,r)
            assert max(rel)<1.
            arrays[f"midpoint_{band}"]=mid
            result.update(shared_midpoint_verified=True,midpoint_residual_upper_over_r=rel)
        else:
            jacs=[state.jacobian(LAYOUT) for state in (a,b)]
            initial=d/np.linalg.norm(d)
            result["budget"]=[]
            budgets=tuple(sorted(GRID+NATIVE_REFINEMENT)) if side==96 else GRID
            for i,beta in enumerate(budgets):
                physical_j=[x[0] for x in jacs]
                new,num=choose_direction(d,physical_j,plus,minus,beta)
                alternatives=[check_direction(d,ed,l,(a,b),jacs,beta) for l in (initial,new)]
                score,chosen,projection,terms=max(alternatives,key=lambda z:z[0])
                status="verified_excluded" if score>2*r else "unresolved"
                key=f"lambda_{band}_{i:02d}"
                arrays[key]=chosen
                outcome=dict(beta=beta,certified_lower_over_two_r=max(0.,down(score/(2*r))),
                             empirical_lower_over_two_r=max(0.,down(alternatives[0][0]/(2*r))),
                             status=status,optimized_direction_nfev=num,
                             projection_lower=projection,endpoint_terms=terms,lambda_key=key,
                             core_mean_gap_lower=down(t-2*beta))
                result["budget"].append(outcome)
            result["largest_tested_excluded_budget"]=max((q["beta"] for q in result["budget"]
                                                            if q["status"]=="verified_excluded"),default=None)
        result["cost"]=dict(**elapsed(timer),verified_material_LU=2,verified_inverse_rhs=len(a.active)+len(b.active),
                            material_current_rhs=16,reciprocal_adjoint_rhs_represented_by_inverse=16,
                            verified_gram_factorizations=2 if side==96 and band==460 else 0,
                            verified_svd_factorizations=2 if band==460 else 0)
        total["frequency"][str(band)]=result
        print(side,band, result["lower_center_distance_over_two_r"],
              result.get("midpoint_residual_upper_over_r"),
              [(b["beta"],b["certified_lower_over_two_r"]) for b in result.get("budget",[])],
              result["cost"],flush=True)
    total.update(**elapsed(start))
    suffix="" if set(bands)=={400,460} else "_"+str(bands[0])
    write_json(HERE/f"frequency_classes_{side}{suffix}.json",total)
    with (HERE/f"frequency_classes_{side}{suffix}.npz").open("wb") as f:
        np.savez_compressed(f,**arrays)
    return total


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("side",type=int,choices=[32,96],nargs="+")
    p.add_argument("--bands",nargs="+",choices=[400,460],type=int,default=[400,460])
    args=p.parse_args()
    for side in args.side:
        run(side,args.bands)

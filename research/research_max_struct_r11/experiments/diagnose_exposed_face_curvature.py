"""Exact original-Maxwell face Hessian at three previously selected exposed endpoints."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.linalg import lu_solve
ROOT=Path(__file__).absolute().parents[1]
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_fbg_r3_experiments.src.physics import MaxwellModel

def load(path):
    with np.load(path,allow_pickle=False) as z:return {k:z[k].copy() for k in z.files}

def main():
    started=time.perf_counter()
    out=ROOT/"RESULTS/exposed_face_curvature";out.mkdir(parents=True,exist_ok=False)
    cfg=dict(status="FROZEN_BEFORE_EVALUATION",cases=["r10_twod_0001","r10_twod_0004","r10_twod_0022"],
        state="t=0 of existing exposed_endpoint_paths: larger original-Phi endpoint, no GT",
        lower=.02,upper=20.,active_distance=1e-8,
        active_scope="Numerical face diagnostic; near-bound and exact-bound counts both retained. Not a rigorous KKT/local-minimum certificate.",
        hessian="Exact data Hessian including residual-weighted second derivative, plus exact smoothed-TV/squared-difference prior Hessian",
        negative_direction="Only if smallest face eigenvalue is negative: retain direction and feasible distance; no automatic new inversion",
        wiring="One centered analytic-gradient JVP on first case, fixed seeded free-face unit vector, h=min(1e-4,quarter legal distance); cost charged",
        truth_access="NONE",seed=202611071117)
    (out/"FROZEN_CONFIG.json").write_text(json.dumps(cfg,indent=2)+"\n")
    inv=json.loads((ROOT/"delegated/diagnostics/STATE_INVENTORY.json").read_text())
    rows=[];costs=[];rng=np.random.default_rng(cfg["seed"])
    for ci,cid in enumerate(cfg["cases"]):
        obj=next(x for x in inv["objects"] if x["case_id"]==cid)
        public={}
        for path in obj["portable_public_files"]:public.update(load(path))
        model=MaxwellModel(public);Y=public["Y"]
        m=load(ROOT/"RESULTS/exposed_endpoint_paths"/(cid+".npz"))["m"][0]
        phi,g=model.data_value_gradient(m,Y);st=model.state(m)
        lo,hi=cfg["lower"],cfg["upper"];atlo=m<=lo+cfg["active_distance"];athi=m>=hi-cfg["active_distance"]
        free=np.flatnonzero(~(atlo|athi));nf=len(free)
        directions=np.eye(model.n)[free]
        dj=model.tangent(m,directions).transpose(1,0,2)
        flat=dj.reshape(model.n,-1)
        dq=(model.G@flat).reshape(model.n,nf,model.ns)
        dy=(model.S@flat).reshape(model.nr,nf,model.ns).transpose(0,2,1).reshape(-1,nf)
        scale=float(np.vdot(Y,Y).real);res=st["Y"]-Y
        with model._physics.operation("FACE_HESSIAN_ADJOINT",rhs=model.ns,adjoint_RHS=model.ns):
            adj=lu_solve(st["lu"],model.S.conj().T@res/scale,trans=2,check_finite=False)
        gn=(dy.conj().T@dy).real/scale
        one=np.einsum("is,ijs->ij",adj[free].conj(),dq[free]).real
        nonlin=one+one.T
        D=model.reg.D;d=D@m;q=np.sqrt(d*d+model.reg.epsilon**2)
        weights=(2*model.reg.smooth+model.reg.tv*model.reg.epsilon**2/q**3)/model.reg.edges
        DF=D[:,free];prior=(DF.T@DF.multiply(weights[:,None])).toarray()
        H=gn+nonlin+prior;eig,vec=np.linalg.eigh(H)
        ray=np.zeros(model.n);ray[free]=vec[:,0]
        feasible=[]
        for sign in [1,-1]:
            v=sign*ray
            caps=np.r_[(hi-m[v>0])/v[v>0],(lo-m[v<0])/v[v<0]]
            feasible.append(float(np.min(caps)))
        row=dict(case_id=cid,original_Phi=phi,n_free=nf,near_lower=int(atlo.sum()),near_upper=int(athi.sum()),
            exact_lower=int(np.sum(m==lo)),exact_upper=int(np.sum(m==hi)),
            free_gradient_inf=float(np.max(abs(g[free]))),
            projected_gradient_inf=float(np.max(abs(m-np.clip(m-g,lo,hi)))),
            lower_multiplier_min=float(np.min(g[atlo])) if atlo.any() else None,
            lower_multiplier_max=float(np.max(g[atlo])) if atlo.any() else None,
            wrong_sign_lower_multipliers=int(np.sum(g[atlo]<0)),
            wrong_sign_upper_multipliers=int(np.sum(g[athi]>0)),
            Hessian_min=float(eig[0]),Hessian_max=float(eig[-1]),negative_eigenvalues=int(np.sum(eig<0)),
            GN_prior_min=float(np.linalg.eigvalsh(gn+prior)[0]),
            nonlinear_Hessian_norm=float(np.linalg.norm(nonlin,2)),
            ray_gradient=float(g@ray),feasible_plus_minus=feasible,
            statement="Face curvature only; active-bound critical-cone directions and approximate stationarity require separate interpretation.")
        if ci==0:
            hdir=np.zeros(model.n);hdir[free]=rng.standard_normal(nf);hdir/=np.linalg.norm(hdir)
            nz=abs(hdir)>0
            eps=min(1e-4,.25*float(np.min(np.minimum(m[nz]-lo,hi-m[nz])/abs(hdir[nz]))))
            gp=model.data_value_gradient(m+eps*hdir,Y)[1]
            gm=model.data_value_gradient(m-eps*hdir,Y)[1]
            numeric=(gp-gm)[free]/(2*eps);analytic=H@hdir[free]
            row["gradient_JVP_step"]=eps
            row["gradient_JVP_relative_error"]=float(np.linalg.norm(numeric-analytic)/max(np.linalg.norm(numeric),np.linalg.norm(analytic)))
            row["gradient_JVP_material_box_legal"]=bool(np.all(m+eps*hdir>=lo)&np.all(m+eps*hdir<=hi)&np.all(m-eps*hdir>=lo)&np.all(m-eps*hdir<=hi))
        np.savez_compressed(out/(cid+".npz"),m=m,g=g,free=free,H=H,GN=gn,nonlinear=nonlin,prior=prior,eigenvalues=eig,weakest_face_direction=ray)
        rows.append(row);costs.append(dict(case_id=cid,operation_counts=dict(model.cost),extra_dense_Green_actions=nf*model.ns,extra_receiver_actions=nf*model.ns))
    (out/"SUMMARY.json").write_text(json.dumps(rows,indent=2)+"\n")
    (out/"COSTS.json").write_text(json.dumps(costs,indent=2)+"\n")
    (out/"RECEIPT.json").write_text(json.dumps(dict(status="COMPLETE",truth_access="NONE",new_inverse_runs=0,
        wall_seconds=time.perf_counter()-started,CPU_seconds_near_exit=time.process_time(),cost_scope="External whole-child CPU authoritative"),indent=2)+"\n")
    print(json.dumps(rows,indent=2))
if __name__=="__main__":main()

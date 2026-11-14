"""Authorized post-replay native-state analysis; never an online selector.

No NumPy/model/material import/read precedes parent comparison authorization,
terminal diagnostic queue and reviewed diagnostic seal validation. No GT NPZ,
optimization, new projection, parameter sweep or identity-verification test.
"""
import os,time
START_CPU=time.process_time();START_WALL=time.perf_counter()
for _k in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS","BLIS_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[_k]="1"
import argparse,csv,importlib.util,json,math,platform,sys
from collections import Counter,defaultdict
from pathlib import Path

TERMINAL={"COMPLETED","COMPLETED_WITH_ERRORS","COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY"}


def read(p):return json.loads(Path(p).read_text(encoding="utf-8"))


def path(root,value):
    p=Path(str(value).replace("\\","/"))
    return p if p.is_absolute() else root/p


def dump(p,value):
    def encode(x):
        if hasattr(x,"tolist"):return x.tolist()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(p).write_text(json.dumps(value,indent=2,default=encode,allow_nan=False)+"\n",encoding="utf-8")


def table(out,name,rows):
    dump(out/(name+".json"),rows)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with (out/(name+".csv")).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields,lineterminator="\n");w.writeheader()
        for r in rows:
            w.writerow({k:json.dumps(r[k],separators=(",",":")) if isinstance(r.get(k),(dict,list)) else r.get(k) for k in fields})


def authorize(root,sup,diagroot,manifest_file,auth_file):
    # compare48 and seal_results import stdlib only.
    auth=read(auth_file)
    if auth.get("authorized_diagnostic_analysis") is not True:
        raise PermissionError("Parent authorized_diagnostic_analysis=true required")
    spec=importlib.util.spec_from_file_location("native_state_compare_guard",sup/"src/compare48.py")
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    _,seals,locations,cases=module.authorize(root,sup,auth_file)
    manifest=read(manifest_file);campaign=manifest["campaign"]
    if auth.get("diagnostic_campaign")!=campaign:
        raise PermissionError("Diagnostic authorization campaign differs")
    if not auth.get("diagnostic_seal"):raise PermissionError("Reviewed diagnostic_seal required")
    seal_path=path(sup,auth["diagnostic_seal"]).resolve();seal=read(seal_path)
    sys.path.insert(0,str(root/"src"));from seal_results import validate_seal
    validate_seal(seal,require_authorization=True)
    if Path(seal["root"]).resolve()!=diagroot or seal["campaign"]!=campaign:
        raise ValueError("Diagnostic seal root/campaign mismatch")
    if seal["manifest_snapshot"]!=manifest:
        raise ValueError("Diagnostic manifest differs from sealed snapshot")
    queue_path=diagroot/"queues"/campaign/"RECEIPT.json";queue=read(queue_path)
    if queue.get("status") not in TERMINAL or queue.get("pending") or queue.get("running"):
        raise RuntimeError("Diagnostic queue must be terminal with no pending/running tasks")
    ids=[t["task_id"]for t in manifest["tasks"]]
    if len(ids)!=len(set(ids)):raise ValueError("Duplicate diagnostic task ID")
    allowed={}
    for name,s in seals.items():
        base=sup if name=="aplus48_extension32" else root
        for t in s["manifest_snapshot"]["tasks"]:
            allowed.setdefault((t["case_id"],t["method"]),set()).add((base/"RESULTS"/name/t["task_id"]).resolve())
    original={t["case_id"]:t for t in seals["confirmation_core48"]["manifest_snapshot"]["tasks"]if t["method"]=="g19"}
    for t in manifest["tasks"]:
        if t["case_id"] not in cases or t.get("algorithm",t["method"]) not in ("cie","dual_cie","g19","a_plus","dr"):
            raise ValueError("Diagnostic task outside declared same48/method scope")
        expected=allowed.get((t["case_id"],t["method"]))
        if not expected or path(root,t["primary_reference_folder"]).resolve()not in expected:
            raise ValueError("Primary reference is not the sealed original/extension run")
        if t["case_id"] not in original:raise ValueError("Missing original public operator reference")
    return auth,seals,locations,manifest,seal,seal_path,queue,original


def load_public(np,root,original_task):
    # Resolve from original sealed public cohort, never a Windows execution path.
    p=path(root,original_task["input_file"]).resolve()
    with np.load(p,allow_pickle=False) as z:
        public={k:z[k].copy()for k in ("G","S","Y","grid_shape")}
        public["E0"]=z["E0"if"E0"in z else"E"].copy()
        for k in ("lower","upper","delta","noise_radius","sigma"):
            if k in z:public[k]=float(z[k])
        if "components"in z and int(z["components"])!=1:raise ValueError("Scalar native-state analysis only")
    return public,p


def counts_delta(after,before):
    return {k:after.get(k,0)-before.get(k,0)for k in set(after)|set(before)}


def jsonl(p):
    if not p.is_file():return []
    with p.open(encoding="utf-8")as f:return [json.loads(line)for line in f if line.strip()]


class CasePhysics:
    """One public case, exact-equal materials reuse physics; LU not banked."""
    def __init__(self,np,model,Y,lo,hi,cid,records,gradient_bank):
        self.np,self.model,self.Y=np,model,Y
        self.lo,self.hi,self.cid=lo,hi,cid
        self.entries=[];self.records=records;self.bank=gradient_bank
        self.actions=Counter()

    def norm(self,x):
        self.actions["explicit_Frobenius_norm_calls"]+=1
        self.actions["explicit_Frobenius_norm_elements"]+=x.size
        return float(self.np.linalg.norm(x))

    def cosine(self,a,b):
        na,nb=self.norm(a),self.norm(b)
        self.actions["real_inner_products"]+=1
        return float(self.np.real(self.np.vdot(a,b))/(na*nb))if na>0 and nb>0 else None

    def physical(self,m,reference):
        np=self.np
        old=next((e for e in self.entries if np.array_equal(e["m"],m)),None)
        if old is not None:
            old["references"].append(reference);return old
        eid=self.cid+"/m%04d"%len(self.entries)
        cpu,wall=time.process_time(),time.perf_counter();before=dict(self.model.cost)
        e=dict(id=eid,m=m.copy(),references=[reference],status="FAILED")
        self.entries.append(e)
        try:
            value,g=self.model.data_value_gradient(m,self.Y) # original beta=None
            state=self.model.state(m)                       # retained same LU
            R,gR,_=self.model.reg.evaluate(m)
            if not np.isfinite(value)or not np.isfinite(g).all()or not np.isfinite(gR).all():
                raise FloatingPointError("Nonfinite original physical value/gradient")
            residual=state["Y"]-self.Y;rr=self.norm(residual);yn=self.norm(self.Y)
            pg=m-np.clip(m-g,self.lo,self.hi)
            metrics=dict(Phi=float(value),D=.5*(rr/yn)**2,R=float(R),
                physical_residual_norm=rr,physical_relative_residual=rr/yn,
                gradient_L2=self.norm(g),gradient_inf=float(np.linalg.norm(g,np.inf)),
                unit_step_box_pg_inf=float(np.linalg.norm(pg,np.inf)),
                unit_step_box_pg_L2=self.norm(pg),gradient_step_scale=1.,
                lower_active_exact=int(np.count_nonzero(m==self.lo)),
                upper_active_exact=int(np.count_nonzero(m==self.hi)),
                outside_box_cells=int(np.count_nonzero((m<self.lo)|(m>self.hi))))
            e.update(status="AVAILABLE",metrics=metrics,gPhi=g.copy(),gR=gR.copy(),
                     Jphys=state["j"].copy(),Yphys=state["Y"].copy())
            prefix=eid.replace("/","__")
            self.bank[prefix+"__m"]=m.copy();self.bank[prefix+"__gPhi"]=g.copy()
            self.bank[prefix+"__gPrior"]=gR.copy()
        except Exception as ex:e["error"]=repr(ex)
        cost=counts_delta(self.model.cost,before)
        record=dict(evaluation_id=eid,case_id=self.cid,status=e["status"],
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            physical_cost_delta=cost,physical_LU=cost.get("LU",0),
            complex_RHS=cost.get("actual_complex_RHS",0),
            adjoint_RHS=cost.get("adjoint_RHS",0),source_RHS=cost.get("source_RHS",0),
            references=e["references"],error=e.get("error"),metrics=e.get("metrics"))
        self.records.append(record);return e

    def field_context(self,J,Q):
        # Exactly one diagnostic G and S action for this fixed field block;
        # reused before/after material. Physical model's internal actions are
        # already covered by physical evaluation costs/whole-process CPU.
        self.actions["diagnostic_G_actions"]+=1
        self.actions["diagnostic_G_complex_RHS"]+=J.shape[1]
        GJ=self.model.G@J
        self.actions["diagnostic_S_actions"]+=1
        self.actions["diagnostic_S_complex_RHS"]+=J.shape[1]
        SJ=self.model.S@J
        return dict(J=J,Q=Q,v=Q-self.model.E-GJ,SJ=SJ)

    def defects(self,f,m,e,delta):
        np=self.np;J,Q,v=f["J"],f["Q"],f["v"]
        u=J-m[:,None]*Q
        forcing=u+m[:,None]*v
        row=dict(u_norm=self.norm(u),v_norm=self.norm(v),
                 constitutive_plus_propagation_forcing_norm=self.norm(forcing),
                 auxiliary_data_residual_over_delta=self.norm(f["SJ"]-self.Y)/delta)
        if e["status"]=="AVAILABLE":
            dj=J-e["Jphys"];nj=self.norm(dj);nf=row["constitutive_plus_propagation_forcing_norm"]
            row.update(current_to_physical_norm=nj,
                port_current_to_physical_over_delta=self.norm(f["SJ"]-e["Yphys"])/delta,
                physical_residual_over_delta=e["metrics"]["physical_residual_norm"]/delta,
                directional_amplification=nj/nf if nf>0 else None,
                amplification_undefined_reason=None if nf>0 else"Zero forcing norm",
                amplification_scope="Observed directional response, not resolvent operator norm")
        else:
            for k in ("current_to_physical_norm","port_current_to_physical_over_delta",
                      "physical_residual_over_delta","directional_amplification"):row[k]=None
        return row


def counterfactual(np,C,blocks,J,Q,m,beta,rho,Os,onehop_fn,cost):
    """Fixed J/Q gradients; explicitly NOT the G19 native gradient."""
    cpu,wall=time.process_time(),time.perf_counter()
    R,gR,_=C.model.reg.evaluate(m)
    den=1+beta*m[:,None]
    residual=beta*(J-m[:,None]*Q)/den
    derivative=-beta*(Q+beta*J)/(den**2)
    weighted=residual@Os
    cost["ordinary_CIE_source_metric_products"]+=1
    ordinary=.5*rho*float(np.vdot(residual,weighted).real)+float(R)
    gOrd=rho*np.real(np.sum(derivative.conj()*weighted,axis=1))+gR
    onehop,gOne=onehop_fn(blocks,J,Q,m,beta,Os,rho,
        lambda z:C.model.reg.evaluate(z)[:2],cost)
    if not np.isfinite(ordinary)or not np.isfinite(onehop)or not np.isfinite(gOrd).all()or not np.isfinite(gOne).all():
        raise FloatingPointError("Nonfinite fixed-current counterfactual value/gradient")
    return dict(ordinary_value=ordinary,onehop_value=float(onehop),
        gOrd=gOrd,gOne=np.array(gOne,copy=True),
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)


def projected(np,m,g,lo,hi):return m-np.clip(m-g,lo,hi)


def current_force_transport(np,C,U,J,Q,SJ,m0,m,beta,rho,Od,Os,native_grad,direction,cost):
    """Parent-specified current-gradient transport at fixed postcurrent J/Q/A.

    This is not a numerical line-orthogonality identity test, a new PR step,
    or the material derivative. No Maxwell solve or coefficient update occurs.
    """
    cpu,wall=time.process_time(),time.perf_counter()
    cost["current_transport_data_metric_products"]+=1
    cost["current_transport_S_adjoint_calls"]+=1
    cost["current_transport_S_adjoint_complex_RHS"]+=J.shape[1]
    data_back=C.model.S.conj().T@((SJ-C.Y)@Od)
    GH=C.model.G.conj().T
    def force(material):
        p=beta/(1+beta*material)
        R=p[:,None]*(J-material[:,None]*Q)
        weighted=R@Os
        cost["current_transport_state_metric_products"]+=1
        cost["current_transport_G_adjoint_calls"]+=1
        cost["current_transport_G_adjoint_complex_RHS"]+=J.shape[1]
        back=p[:,None]*weighted-GH@((material*p)[:,None]*weighted)
        return U.adjoint(data_back+rho*back)
    before,after=force(m0),force(m)
    if before.shape!=direction.shape or before.shape!=native_grad.shape:
        raise ValueError("Saved native coefficient gradient/direction shape differs from saved basis")
    if not np.isfinite(before).all()or not np.isfinite(after).all():
        raise FloatingPointError("Nonfinite recomputed current force")
    tau=after-before
    nb,na,nd=C.norm(before),C.norm(after),C.norm(direction)
    cost["current_transport_real_inner_products"]+=2
    db=float(np.vdot(before,direction).real);da=float(np.vdot(after,direction).real)
    metrics=dict(beta_current=beta,current_force_transport_norm=C.norm(tau),
        postcurrent_prereadout_gradient_norm=nb,postreadout_current_gradient_norm=na,
        native_precurrent_gradient_norm=C.norm(native_grad),native_current_direction_norm=nd,
        cosine_transport_native_direction=C.cosine(tau,direction),
        cosine_postreadout_prereadout_current_gradients=C.cosine(after,before),
        postcurrent_prereadout_gradient_dot_native_direction=db,
        postreadout_gradient_dot_native_direction=da,
        normalized_prereadout_gradient_dot_native_direction=db/(nb*nd)if nb>0 and nd>0 else None,
        normalized_postreadout_gradient_dot_native_direction=da/(na*nd)if na>0 and nd>0 else None,
        gradient_transport_relative_to_prereadout=C.norm(tau)/nb if nb>0 else None,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        scope="Fixed postcurrent auxiliary J/Q/A; beta_current; saved native grad is PRE-current, not recomputed prereadout force",
        operator_scope="Original public dense G/S adjoints and native Fourier/DPCS basis; no new SVD/Maxwell solve/update",
        interpretation="Material-readout-induced change of current force/PR context; no chaos, basin or causal-prediction certificate",
        undefined_reason="Normalized values/cosines null when a participating norm is zero")
    return metrics,dict(gA_postcurrent_prereadout=before,gA_postreadout=after,tau_current_force=tau)


def native_rows(task,folder):
    dense=[]
    wanted=("stage","iteration","global_iteration","beta_current","beta_material","mf",
        "data_term","state_term_current","state_term_material","regularizer",
        "data_before_current","state_before_current","data_after_current","state_after_current",
        "current_step_objective_change","readout_material_change","PR_plus","line_alpha",
        "direction_restarted","material_old_block_loss","material_scalar_seed_block_loss",
        "material_final_block_loss","material_fallback","aplus_sweeps","event","status","wall")
    for r in jsonl(folder/"DIAGNOSTIC_NATIVE_HISTORY.jsonl"):
        dense.append(dict(case_id=task["case_id"],task_id=task["task_id"],method=task["method"],
            stream="native_CIE",**{k:r.get(k)for k in wanted},
            CPU_seconds=None,CPU_reason="Native scalar history stores wall, not per-update process CPU"))
    for r in jsonl(folder/"DIAGNOSTIC_DR_STEPS.jsonl"):
        dense.append(dict(case_id=task["case_id"],task_id=task["task_id"],method=task["method"],
            stream="native_DR",iteration=r["iteration"],CPU_seconds=r.get("CPU_seconds"),
            shadow_gap=r.get("shadow_gap"),multiplier=r.get("A_info",{}).get("multiplier"),
            gap_scope="Existing native audit-time gap only; unsampled gaps null"))
    return dense


def snapshot_base(task,entry,observer):
    scalar=entry.get("scalar",{})
    return dict(case_id=task["case_id"],task_id=task["task_id"],method=task["method"],
        snapshot_id=entry["snapshot_id"],kind=entry["kind"],stage=entry["stage"],
        iteration=entry["iteration"],global_iteration=entry["global_iteration"],
        snapshot_CPU_seconds=scalar.get("CPU_seconds"),phase=scalar.get("phase"),
        retained_slots=entry.get("slots",[]),
        requested_target_iterations=observer.get("requested_target_iterations"),
        reached_target_iterations=observer.get("reached_target_iterations"),
        unreached_target_iterations=observer.get("unreached_target_iterations"))


def target_checks(np,task,folder,entries,z,observer,bank):
    requests=task.get("diagnostic_target_iterations",[])
    rawfolder=Path(str(task["primary_reference_folder"]).replace("\\","/"))
    rawmeta=read(rawfolder/"raw.json")
    with np.load(rawfolder/"raw.npz",allow_pickle=False)as a:raw=a["m"].copy()
    aux=rawmeta.get("auxiliary")or{}
    original_it=aux.get("global_iteration",aux.get("iteration"))
    rows=[]
    for target in requests:
        selected=[e for e in entries if e["global_iteration"]==target]
        row=dict(case_id=task["case_id"],task_id=task["task_id"],method=task["method"],
            requested_iteration=target,original_raw_iteration=original_it,
            request_matches_original_raw_iteration=target==original_it,
            original_raw_candidate=rawmeta.get("candidate"),
            original_raw_path=str(rawfolder/"raw.npz"),
            observer_declares_reached=target in observer.get("reached_target_iterations",[]),
            target_snapshot_count=len(selected),status="UNREACHED",
            comparison_scope="All material cells only; auxiliary J/Q/PR equality is not established",
            material_array_equal=None,material_max_abs_difference=None,
            material_RMS_difference=None,qualification="No replay equality asserted for absent/different state")
        if selected:
            entry=selected[0];key=entry["array_prefix"]+("cm"if entry["kind"]=="dr"else"m")
            m=z[key].copy()
            if m.shape!=raw.shape:
                row["status"]="SHAPE_MISMATCH"
            else:
                d=m-raw;equal=bool(np.array_equal(m,raw))
                row.update(status="EXACT_RECORDED_MATERIAL_MATCH"if equal else"DIFFERENT_RECORDED_MATERIAL",
                    material_array_equal=equal,material_max_abs_difference=float(np.max(np.abs(d))),
                    material_RMS_difference=float(np.linalg.norm(d)/np.sqrt(d.size)),
                    snapshot_id=entry["snapshot_id"])
                bank[task["task_id"]+"__target_%s__minus_primary_raw"%target]=d
        rows.append(row)
    return rows


def plot_case(np,plt,out,cid,dense,snapshots,transitions,fieldrows,audits,transports):
    figures=out/"FIGURES";figures.mkdir(exist_ok=True)
    def save(fig,name):
        fig.tight_layout()
        for ext in ("png","pdf"):fig.savefig(figures/(name+"."+ext),dpi=160)
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8))
    for method in sorted({r["method"]for r in dense}):
        rr=[r for r in dense if r["method"]==method and r.get("global_iteration")is not None]
        for ax,key in zip(axes.flat,("state_term_current","state_term_material","readout_material_change","PR_plus")):
            ss=[r for r in rr if isinstance(r.get(key),(int,float))]
            if ss:ax.plot([r["global_iteration"]for r in ss],[r[key]for r in ss],label=method,lw=.8)
            ax.set_xlabel("Completed native update");ax.set_ylabel(key);ax.set_yscale("symlog",linthresh=1e-10)
        rest=[r for r in rr if r.get("direction_restarted")]
        if rest:axes[1,1].scatter([r["global_iteration"]for r in rest],[r.get("PR_plus",0)for r in rest],marker="x",s=10)
    for ax in axes.flat:
        if ax.lines:ax.legend(fontsize=7)
    fig.suptitle(cid+" — dense native history; auxiliary losses are not physical Phi")
    save(fig,cid+"__native_history")
    fig,axes=plt.subplots(1,3,figsize=(15,4))
    for method in sorted({r["method"]for r in transitions}):
        rr=[r for r in transitions if r["method"]==method and r.get("status")=="AVAILABLE"]
        x=[r["global_iteration"]for r in rr]
        for ax,key in zip(axes,("actual_delta_Phi","cosine_step_negative_physical_gradient",
                               "cosine_step_negative_onehop_gradient")):
            yy=[r.get(key)for r in rr]
            ax.plot(x,yy,"o-",ms=3,label=method)
            ax.set_xlabel("Saved native update");ax.set_ylabel(key);ax.axhline(0,color="gray",lw=.5)
    for ax in axes:
        if ax.lines:ax.legend(fontsize=7)
    fig.suptitle(cid+" — readout movement; onehop direction is counterfactual for G19")
    save(fig,cid+"__alignment")
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for method in sorted({r["method"]for r in fieldrows}):
        rr=[r for r in fieldrows if r["method"]==method and r.get("port_current_to_physical_over_delta")is not None]
        if rr:
            axes[0].scatter([r["u_norm"]for r in rr],[r["port_current_to_physical_over_delta"]for r in rr],s=14,label=method)
            axes[1].scatter([r["v_norm"]for r in rr],[r["directional_amplification"]for r in rr],s=14,label=method)
    for ax in axes:
        ax.set_xscale("symlog",linthresh=1e-10);ax.set_yscale("symlog",linthresh=1e-10)
        if ax.collections:ax.legend(fontsize=7)
    axes[0].set_xlabel("Constitutive defect ||u||");axes[0].set_ylabel("||S(J-Jphys)|| / delta")
    axes[1].set_xlabel("Propagation defect ||v||");axes[1].set_ylabel("Observed directional amplification")
    fig.suptitle(cid+" — saved field blocks; no resolvent operator-norm certificate")
    save(fig,cid+"__observable_defects")
    fig,ax=plt.subplots(figsize=(9,4))
    for method in sorted({r["method"]for r in audits}):
        rr=[r for r in audits if r["method"]==method]
        ax.scatter([r["CPU_seconds"]for r in rr],[r["Phi"]for r in rr],s=3,alpha=.4,label=method)
    ax.set_xlabel("Actual replay process CPU at paid physical audit")
    ax.set_ylabel("Original physical Phi (all paid trials)");ax.set_yscale("symlog",linthresh=1e-10)
    if ax.collections:ax.legend(fontsize=7)
    fig.suptitle(cid+" — actual physical audit candidates; diagnostic replay only")
    save(fig,cid+"__paid_physical_audits")
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for method in sorted({r["method"]for r in transports}):
        rr=[r for r in transports if r["method"]==method and r.get("status")=="AVAILABLE"]
        for ax,key in zip(axes,("gradient_transport_relative_to_prereadout","normalized_postreadout_gradient_dot_native_direction")):
            ss=[r for r in rr if r.get(key)is not None]
            if ss:ax.plot([r["global_iteration"]for r in ss],[r[key]for r in ss],"o-",ms=3,label=method)
            ax.set_xlabel("Saved native update");ax.set_ylabel(key);ax.axhline(0,color="gray",lw=.5)
    for ax in axes:
        if ax.lines:ax.legend(fontsize=7)
    axes[0].set_yscale("symlog",linthresh=1e-8)
    fig.suptitle(cid+" — current-force transport at fixed postcurrent fields; posthoc diagnostic")
    save(fig,cid+"__current_force_transport")
    dr=[r for r in dense if r["stream"]=="native_DR"]
    if dr:
        fig,axes=plt.subplots(1,2,figsize=(10,4))
        axes[0].plot([r["iteration"]for r in dr],[r.get("CPU_seconds")for r in dr]);axes[0].set_ylabel("Process CPU at native update")
        ss=[r for r in dr if r.get("shadow_gap")is not None]
        axes[1].plot([r["iteration"]for r in ss],[r["shadow_gap"]for r in ss],"o-")
        axes[1].set_ylabel("Native audit-time reflector gap")
        for ax in axes:ax.set_xlabel("Completed DR update")
        fig.suptitle(cid+" — DR scalar history, gaps sparse as originally computed")
        save(fig,cid+"__DR_history")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root","--parent-root",dest="root",type=Path,required=True)
    p.add_argument("--supplement-root","--sup",dest="supplement_root",type=Path)
    p.add_argument("--diag-root",type=Path,required=True)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--compareauth",type=Path,required=True)
    p.add_argument("--legacy-repo",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--no-figures",action="store_true")
    a=p.parse_args();root=a.root.resolve();sup=(a.supplement_root or root/"supplement_aplus48").resolve()
    diagroot=a.diag_root.resolve()
    auth,seals,locations,manifest,seal,sealpath,queue,original=authorize(root,sup,diagroot,a.manifest.resolve(),a.compareauth.resolve())
    out=a.out.resolve()
    banned=[base/name for base in (root,sup,diagroot)for name in
            ("RESULTS","queues","sealed_truth","data_public","src","runtime","configs")]
    if out in (root,sup,diagroot)or any(out.is_relative_to(x)for x in banned):
        raise ValueError("Analysis output must be a new directory outside predictions/data/source")
    for s in list(seals.values())+[seal]:
        if any(out.is_relative_to(Path(d["path"]).resolve())for d in s["directories"]):
            raise ValueError("Analysis output overlaps a sealed directory")
    out.mkdir(parents=True,exist_ok=False)
    # Authorization, all seals, terminal queue and task scope precede these imports.
    import numpy as np
    sys.path.insert(0,str(a.legacy_repo.resolve()))
    import research.research_fbg_r3_experiments.src.physics as physics
    if not Path(physics.__file__).resolve().is_relative_to(a.legacy_repo.resolve()):
        raise RuntimeError("Imported model is outside requested inherited repository")
    from research.research_max_struct_r4.src.cie_onehop import OneHopAbsorption,onehop_material_value_gradient
    from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
    plt=None
    if not a.no_figures:
        import matplotlib;matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    scores={}
    score_file=sup/"analysis/comparison48/ALL48_SCORES.csv"
    comp_receipt=sup/"analysis/comparison48/COMPARISON_RECEIPT.json"
    if score_file.is_file():
        comp=read(comp_receipt)
        if comp.get("seals")!=locations:raise ValueError("Existing comparison labels have different seals")
        with score_file.open(newline="",encoding="utf-8")as f:
            for r in csv.DictReader(f):
                if r.get("view")in("raw","primary"):
                    scores[(r["case_id"],r["method"],r["view"])]={
                        k:r.get(k)for k in ("NRMSE","IoU","relativeL2","mean_amplitude_bias")}
    tasks_by_case=defaultdict(list)
    for task in manifest["tasks"]:tasks_by_case[task["case_id"]].append(task)
    coverage=[];snaprows=[];transitions=[];fieldrows=[];targets=[];evaluations=[];costs=[];issues=[];all_dense=[];audits=[];transports=[]
    gradient_dir=out/"GRADIENT_BANKS";gradient_dir.mkdir()
    protocol=manifest["protocol"];prior=protocol["prior"]
    for cid,tasks in tasks_by_case.items():
        case_cpu,case_wall=time.process_time(),time.perf_counter();bank={};model=None;C=None;blockcost=Counter();transportcost=Counter()
        try:
            public,public_file=load_public(np,root,original[cid])
            model=physics.MaxwellModel(public)
            model.reg=physics.Regularizer(model.shape,**{k:float(prior[k])for k in("smooth","tv","epsilon")})
            lo=float(public.get("lower",protocol["material_box"][0]));hi=float(public.get("upper",protocol["material_box"][1]))
            C=CasePhysics(np,model,public["Y"],lo,hi,cid,evaluations,bank)
            blocks=OneHopAbsorption(model.G,model.shape,blockcost)
            case_dense=[];case_trans=[];case_fields=[];case_audits=[];case_transports=[]
            for task in tasks:
                folder=diagroot/"RESULTS"/manifest["campaign"]/task["task_id"]
                cov=dict(case_id=cid,task_id=task["task_id"],method=task["method"],status="UNAVAILABLE",
                    requested_target_iterations=task.get("diagnostic_target_iterations",[]),
                    reached_target_iterations=None,unreached_target_iterations=None,saved_snapshots=None,
                    public_input_path=str(public_file),issues=[])
                coverage.append(cov)
                try:
                    observer=read(folder/"DIAGNOSTIC_RECEIPT.json")
                    if observer.get("role")!="POST_SCORE_OFFLINE_DIAGNOSTIC"or observer.get("primary_selection_changed")is not False:
                        raise ValueError("Wrong diagnostic role/selection declaration")
                    if observer.get("metadata",{}).get("task_id")!=task["task_id"]or observer.get("method")!=task.get("algorithm",task["method"]):
                        raise ValueError("Observer task/method differs from sealed manifest")
                    idx=read(folder/"DIAGNOSTIC_STATE_INDEX.json");entries=idx["states"]
                    cov.update(status="AVAILABLE",saved_snapshots=len(entries),
                        observer_status=observer.get("observer_status"),
                        reached_target_iterations=observer.get("reached_target_iterations"),
                        unreached_target_iterations=observer.get("unreached_target_iterations"),
                        replay_termination=observer.get("final_context"),peak_full_states=observer.get("peak_full_states"))
                    dense=native_rows(task,folder);case_dense.extend(dense);all_dense.extend(dense)
                    dr_cpu={r["iteration"]:r["CPU_seconds"]for r in dense if r["stream"]=="native_DR"}
                    for h in read(folder/"HISTORY.json"):
                        ar=dict(case_id=cid,task_id=task["task_id"],method=task["method"],
                            candidate=h["candidate"],phase=h["phase"],CPU_seconds=h["CPU_seconds"],
                            Phi=h["Phi"],D=h["D"],R=h["R"],eligible=h["eligible"],
                            global_iteration=(h.get("auxiliary")or{}).get("global_iteration"),
                            iteration=(h.get("auxiliary")or{}).get("iteration"))
                        audits.append(ar);case_audits.append(ar)
                    delta=float(read(folder/"RECEIPT.json")["delta"])
                    stages={r["stage"]:r for r in idx.get("stages",[])}
                    bases={}
                    with np.load(folder/"DIAGNOSTIC_STATES.npz",allow_pickle=False)as z:
                        targets.extend(target_checks(np,task,folder,entries,z,observer,bank))
                        Os=z["static__Os"].copy()if"static__Os"in z else None
                        Od=z["static__Od"].copy()if"static__Od"in z else None
                        Ls=z["static__Ls"].copy()if"static__Ls"in z else None
                        Vd=z["static__Vd"].copy()if"static__Vd"in z else None
                        JD=z["static__JD"].copy()if"static__JD"in z else None
                        blockmatch=None
                        if"static__block_goff"in z:
                            blockmatch=bool(np.array_equal(z["static__block_goff"],blocks.goff)and np.array_equal(z["static__block_diag_g"],blocks.diag_g)and np.array_equal(z["static__block_indices"],blocks.indices))
                        previous_dr=None
                        for entry in sorted(entries,key=lambda e:e["global_iteration"]):
                            base=snapshot_base(task,entry,observer);prefix=entry["array_prefix"]
                            if base["snapshot_CPU_seconds"]is None and entry["kind"]=="dr":
                                base["snapshot_CPU_seconds"]=dr_cpu.get(entry["global_iteration"])
                            def arr(key):return z[prefix+key].copy()
                            row=dict(base,status="UNAVAILABLE",native_primary_label=scores.get((cid,task["method"],"primary")),
                                     saved_block_context_matches_public_G=blockmatch)
                            snaprows.append(row)
                            try:
                                if not math.isfinite(delta)or delta<=0:raise ValueError("Positive original data-ball delta required")
                                if entry["kind"]=="cie":
                                    m0,m=arr("m_before_readout"),arr("m");J,Q=arr("J"),arr("Q")
                                    stage=stages[entry["stage"]];rho=float(stage["rho"]);beta=float(stage["beta_material"])
                                    if Os is None or Ls is None:raise ValueError("Saved source metrics missing; no default invented")
                                    row.update(rho=rho,beta_material=beta,beta_current=stage["beta_current"],
                                        sourcewise=stage.get("sourcewise"),receiver_rank=stage.get("receiver_rank"),mf=stage.get("mf"),
                                        Ls_context_retained=True,counterfactual_scope="Ordinary scalar CIE and exact onehop at fixed auxiliary J/Q, NOT G19 native gradient")
                                    before=C.physical(m0,dict(base,phase="before_readout"))
                                    after=C.physical(m,dict(base,phase="after_readout"))
                                    physical_ok=before["status"]=="AVAILABLE"and after["status"]=="AVAILABLE"
                                    row["physical_evaluation_status"]={"before":before["status"],"after":after["status"]}
                                    f=C.field_context(J,Q)
                                    for phase,material,e in(("before_readout",m0,before),("after_readout",m,after)):
                                        defects=dict(base,field_block="fixed_current_"+phase,physical_evaluation_id=e["id"],**C.defects(f,material,e,delta))
                                        fieldrows.append(defects);case_fields.append(defects)
                                        if e["status"]=="AVAILABLE":
                                            for k,v in e["metrics"].items():row[phase+"_"+k]=v
                                    if Vd is not None and JD is not None:
                                        C.actions["DPCS_adjoint_actions"]+=1;C.actions["DPCS_complex_RHS"]+=J.shape[1]
                                        row["DPCS_locked_coordinate_norm"]=C.norm(Vd.conj().T@(J-JD))
                                    cb=counterfactual(np,C,blocks,J,Q,m0,beta,rho,Os,onehop_material_value_gradient,blockcost)
                                    ca=counterfactual(np,C,blocks,J,Q,m,beta,rho,Os,onehop_material_value_gradient,blockcost)
                                    for phase,material,counter in(("before",m0,cb),("after",m,ca)):
                                        for key in("ordinary_value","onehop_value","CPU_seconds","wall_seconds"):row["counterfactual_"+phase+"_"+key]=counter[key]
                                        for key in("gOrd","gOne"):
                                            pg=projected(np,material,counter[key],lo,hi)
                                            row["counterfactual_"+phase+"_"+key+"_pg_inf"]=float(np.linalg.norm(pg,np.inf))
                                            bank[task["task_id"]+"__s%s__%s__%s"%(entry["snapshot_id"],phase,key)]=counter[key]
                                    step=m-m0;gdiff=cb["gOne"]-cb["gOrd"]
                                    pgOrd=projected(np,m0,cb["gOrd"],lo,hi)
                                    pgOne=projected(np,m0,cb["gOne"],lo,hi)
                                    C.actions["real_inner_products"]+=2
                                    tr=dict(base,status="UNAVAILABLE",before_evaluation_id=before["id"],after_evaluation_id=after["id"],material_step_norm=C.norm(step),
                                        cosine_step_negative_ordinary_gradient=C.cosine(step,-cb["gOrd"]),
                                        cosine_step_negative_onehop_gradient=C.cosine(step,-cb["gOne"]),
                                        cosine_step_negative_gradient_difference=C.cosine(step,-gdiff),
                                        cosine_ordinary_onehop_gradients=C.cosine(cb["gOrd"],cb["gOne"]),
                                        counterfactual_gradient_difference_norm=C.norm(gdiff),
                                        cosine_step_negative_ordinary_pg=C.cosine(step,-pgOrd),
                                        cosine_step_negative_onehop_pg=C.cosine(step,-pgOne),
                                        cosine_step_negative_pg_difference=C.cosine(step,-(pgOne-pgOrd)),
                                        counterfactual_pg_difference_norm=C.norm(pgOne-pgOrd),
                                        ordinary_first_order_step=float(np.dot(cb["gOrd"],step)),
                                        onehop_first_order_step=float(np.dot(cb["gOne"],step)),
                                        counterfactual_ordinary_loss_change=ca["ordinary_value"]-cb["ordinary_value"],
                                        counterfactual_onehop_loss_change=ca["onehop_value"]-cb["onehop_value"])
                                    if before["status"]=="AVAILABLE"and after["status"]=="AVAILABLE":
                                        C.actions["real_inner_products"]+=1
                                        tr.update(status="AVAILABLE",actual_delta_Phi=after["metrics"]["Phi"]-before["metrics"]["Phi"],
                                            physical_first_order_step=float(np.dot(before["gPhi"],step)),
                                            cosine_step_negative_physical_gradient=C.cosine(step,-before["gPhi"]),
                                            cosine_step_negative_physical_pg=C.cosine(step,-projected(np,m0,before["gPhi"],lo,hi)),
                                            cosine_physical_ordinary_gradients=C.cosine(before["gPhi"],cb["gOrd"]),
                                            cosine_physical_onehop_gradients=C.cosine(before["gPhi"],cb["gOne"]))
                                    transitions.append(tr);case_trans.append(tr)
                                    transport=dict(base,status="UNAVAILABLE")
                                    transports.append(transport);case_transports.append(transport)
                                    try:
                                        if Vd is None or Od is None:raise ValueError("Saved Vd/Od unavailable for current-force transport")
                                        if entry["stage"]not in bases:
                                            U=NativeFourierAPCS(model.shape,stage["mf"],Vd,transportcost)
                                            transportcost["current_transport_basis_constructors"]+=1
                                            if not np.array_equal(U.frequencies,np.asarray(stage["frequencies"])):
                                                raise ValueError("Native basis frequency order differs from saved stage")
                                            bases[entry["stage"]]=U
                                        tm,tb=current_force_transport(np,C,bases[entry["stage"]],J,Q,f["SJ"],m0,m,
                                            float(stage["beta_current"]),rho,Od,Os,arr("grad"),arr("direction"),transportcost)
                                        transport.update(status="AVAILABLE",**tm)
                                        for key,value in tb.items():bank[task["task_id"]+"__s%s__"%entry["snapshot_id"]+key]=value
                                    except Exception as ex:
                                        transport["error"]=repr(ex);issues.append(dict(base,error=repr(ex),scope="current_force_transport"))
                                else:
                                    m=arr("cm");scale=float(protocol["scale"])
                                    e=C.physical(m,dict(base,phase="DR_constitutive_shadow"))
                                    physical_ok=e["status"]=="AVAILABLE"
                                    row["physical_evaluation_status"]=e["status"]
                                    for block,jkey,qkey in(("constitutive","cX","cQ"),("reflected_A","reflected_A_X","reflected_A_Q"),("pre_iterator","pre_X","pre_Q")):
                                        f=C.field_context(scale*arr(jkey),arr(qkey))
                                        defects=dict(base,field_block=block,physical_evaluation_id=e["id"],**C.defects(f,m,e,delta))
                                        fieldrows.append(defects);case_fields.append(defects)
                                    def pairnorm(x,q):return float(np.hypot(C.norm(x),C.norm(q)))
                                    row.update(full_state_memory_norm=pairnorm(arr("pre_X")-arr("cX"),arr("pre_Q")-arr("cQ")),
                                        reflector_shadow_gap_norm=pairnorm(arr("reflected_A_X")-arr("cX"),arr("reflected_A_Q")-arr("cQ")),
                                        iterator_step_norm=pairnorm(arr("post_X")-arr("pre_X"),arr("post_Q")-arr("pre_Q")),
                                        A_info=entry.get("scalar",{}).get("A_info"),
                                        physical_evaluation_id=e["id"],physical_metrics=e.get("metrics"),
                                        post_iterator_material_projection=None,
                                        post_projection_reason="No extra C or A projection; next shadow comparison uses retained checkpoint only")
                                    if previous_dr is not None:
                                        row.update(previous_saved_shadow_iteration=previous_dr[0],
                                            checkpoint_shadow_material_change_norm=C.norm(m-previous_dr[1]),
                                            checkpoint_change_scope="Between saved checkpoints, not a native one-step change")
                                    previous_dr=(entry["global_iteration"],m.copy())
                                row["status"]="AVAILABLE"if physical_ok else"PARTIAL_PHYSICS_UNAVAILABLE"
                            except Exception as ex:row["error"]=repr(ex);issues.append(dict(base,error=repr(ex),scope="snapshot"))
                except Exception as ex:
                    cov["status"]="UNAVAILABLE";cov["issues"].append(repr(ex))
                    issues.append(dict(case_id=cid,task_id=task["task_id"],error=repr(ex),scope="task"))
                    for target in task.get("diagnostic_target_iterations",[]):
                        if not any(r["task_id"]==task["task_id"]and r["requested_iteration"]==target for r in targets):
                            targets.append(dict(case_id=cid,task_id=task["task_id"],method=task["method"],
                                requested_iteration=target,status="UNAVAILABLE",error=repr(ex)))
            np.savez(gradient_dir/(cid+".npz"),**bank)
            if plt is not None:plot_case(np,plt,out,cid,case_dense,snaprows,case_trans,case_fields,case_audits,case_transports)
            costs.append(dict(case_id=cid,status="PROCESSED",CPU_seconds=time.process_time()-case_cpu,
                wall_seconds=time.perf_counter()-case_wall,physical_cost=dict(model.cost),
                explicit_diagnostic_actions=dict(C.actions),onehop_and_counterfactual_cost=dict(blockcost),
                current_force_transport_cost=dict(transportcost),
                unique_physical_materials=len(C.entries),gradient_bank_arrays=len(bank),
                cost_scope="Whole case includes model/setup/cache lookup/fields/gradients/export/figures; per-evaluation/counterfactual timers are nested, not additive"))
        except Exception as ex:
            issues.append(dict(case_id=cid,error=repr(ex),scope="case"))
            costs.append(dict(case_id=cid,status="FAILED",CPU_seconds=time.process_time()-case_cpu,
                wall_seconds=time.perf_counter()-case_wall,error=repr(ex),
                physical_cost=dict(model.cost)if model is not None else None,
                explicit_diagnostic_actions=dict(C.actions)if C is not None else None,
                onehop_and_counterfactual_cost=dict(blockcost),current_force_transport_cost=dict(transportcost)))
            for task in tasks:
                if not any(r["task_id"]==task["task_id"]for r in coverage):
                    coverage.append(dict(case_id=cid,task_id=task["task_id"],method=task["method"],status="UNAVAILABLE",issues=[repr(ex)]))
    for task in manifest["tasks"]:
        for target in task.get("diagnostic_target_iterations",[]):
            if not any(r["task_id"]==task["task_id"]and r["requested_iteration"]==target for r in targets):
                targets.append(dict(case_id=task["case_id"],task_id=task["task_id"],method=task["method"],
                    requested_iteration=target,status="UNAVAILABLE",error="No retained target evidence after task/case failure"))
    for name,rows in(("TASK_COVERAGE",coverage),("SNAPSHOTS",snaprows),("READOUT_TRANSITIONS",transitions),
        ("FIELD_DEFECTS",fieldrows),("RAW_TARGET_CHECKS",targets),("PHYSICAL_EVALUATIONS",evaluations),
        ("CURRENT_FORCE_TRANSPORT",transports),
        ("DENSE_NATIVE_HISTORY",all_dense),("PAID_PHYSICAL_AUDITS",audits),("CASE_COSTS",costs),("ISSUES",issues)):
        table(out,name,rows)
    receipt=dict(status="POST_AUTHORIZED_NATIVE_STATE_ANALYSIS_COMPLETE",role="OFFLINE_DIAGNOSTIC",
        expected_tasks=len(manifest["tasks"]),task_rows=len(coverage),snapshot_rows=len(snaprows),
        raw_target_rows=len(targets),unique_physical_evaluations=len(evaluations),
        original_Physical_LU=sum(r["physical_LU"]for r in evaluations),
        original_complex_RHS=sum(r["complex_RHS"]for r in evaluations),
        original_adjoint_RHS=sum(r["adjoint_RHS"]for r in evaluations),
        CPU_seconds=time.process_time()-START_CPU,wall_seconds=time.perf_counter()-START_WALL,
        authorization=str(a.compareauth.resolve()),diagnostic_seal=str(sealpath),
        primary_seals=locations,diagnostic_manifest=str(a.manifest.resolve()),
        platform=platform.platform(),model_source=str(Path(physics.__file__).resolve()),
        truth_NPZ_read=False,GT_labels_source=str(score_file)if scores else None,
        new_optimizer=False,new_projection=False,new_identity_verification=False,
        new_candidates_or_primary_selection=False,identity_scope="Use defect decomposition only; no residual identity error test",
        counterfactual_scope="Ordinary CIE/onehop fixed-current gradients, not G19 native gradients",
        current_transport_scope="Fixed postcurrent J/Q/A, beta_current, saved native coefficient gradient is PRE-current; no new update/SVD/Maxwell solve for this component",
        cache_scope="Exact np.array_equal material reuse within the same case/public operator/prior; LU not retained across unique materials",
        cost_scope="New same-finite-physics offline forward/adjoint diagnostics, not independent validation; all setup/failure/export costs paid. Nested timers never added twice",
        timer_scope="Whole analysis process from imports/authorization through tables/figures; final receipt write/interpreter exit require external whole-child CPU",
        issues=len(issues))
    dump(out/"RECEIPT.json",receipt)
    print(json.dumps(receipt,indent=2))
    return 0


if __name__=="__main__":raise SystemExit(main())


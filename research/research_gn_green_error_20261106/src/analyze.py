"""Frozen GN endpoint error geometry. CPU, one BLAS thread, no inverse reruns."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
    os.environ[key]="1"
from pathlib import Path
from collections import defaultdict
import csv,json,sys,time
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parents[1]
REPO=HERE.parents[1]
sys.path.insert(0,str(REPO))
from research.research_z_alg_r4_diversity.src.physics import Model

BANDS=("strong","weak","complement")
TAUS=(.03,.1,.3)
RANK_RTOL=1e-12
INT_TOL=1e-9
PROJ_TOL=1e-10
cost=[]

def energy(x):
    return float(np.vdot(x,x).real)

def timed_svd(a,label,full=False):
    t=time.perf_counter()
    u,d,vh=svd(a,full_matrices=full,check_finite=False,lapack_driver="gesdd")
    cost.append(dict(operation=label,rows=a.shape[0],columns=a.shape[1],
                     wall_seconds=time.perf_counter()-t,LU=0,source_RHS=0))
    return u,d,vh

def spectrum(a,label):
    # Wide receiver matrices need the full right basis, including structural nulls.
    u,d,vh=timed_svd(a,"fixed_spectrum_"+label,full=a.shape[0]<a.shape[1])
    ss=np.zeros(a.shape[1]);ss[:len(d)]=d
    return dict(matrix=a,s=ss,V=vh.conj().T,label=label)

def split(op,tau):
    rel=op["s"]/op["s"][0]
    ids=[np.flatnonzero(rel>=tau),
         np.flatnonzero((rel<tau)&(rel>RANK_RTOL)),
         np.flatnonzero(rel<=RANK_RTOL)]
    return dict(**op,ids=ids,bases=[op["V"][:,i] for i in ids],tau=tau)

def components(op,x):
    # Two groups plus the orthogonal remainder avoid costly full-null products.
    norms=[b.shape[1] for b in op["bases"]]
    omitted=int(np.argmax(norms))
    parts=[None,None,None]
    for k,b in enumerate(op["bases"]):
        if k!=omitted:
            parts[k]=b@(b.conj().T@x)
    parts[omitted]=x-sum(parts[k] for k in range(3) if k!=omitted)
    return parts

def joint_spaces(inn,ob):
    caps={};projected={};meta=[]
    n=inn["V"].shape[0]
    for a,va in enumerate(inn["bases"]):
        for b,vb in enumerate(ob["bases"]):
            ra,rb=va.shape[1],vb.shape[1]
            if not ra or not rb:
                caps[a,b]=np.empty((n,0),dtype=inn["V"].dtype)
                projected[a,b]=np.empty((n,0),dtype=inn["V"].dtype)
                meta.append(dict(in_band=BANDS[a],ob_band=BANDS[b],intersection_rank=0,projected_rank=0))
                continue
            cross=vb.conj().T@va
            u,d,vh=timed_svd(cross,f"projected_span_{inn['label']}_{inn['tau']}_{a}{b}")
            qz=vb@u[:,d>PROJ_TOL]
            other=np.concatenate([ob["ids"][k] for k in range(3) if k!=b])
            constraints=ob["V"][:,other].conj().T@va
            if constraints.shape[0]:
                _,sins,nh=timed_svd(constraints,f"intersection_{inn['label']}_{inn['tau']}_{a}{b}",full=True)
                rk=int(np.count_nonzero(sins>INT_TOL))
                qc=va@nh[rk:].conj().T
                rank_sensitivity={str(t):int(ra-np.count_nonzero(sins>t)) for t in (1e-10,1e-9,1e-8)}
            else:
                qc=va.copy();rank_sensitivity={str(t):ra for t in (1e-10,1e-9,1e-8)}
            caps[a,b]=qc;projected[a,b]=qz
            clipped=np.clip(d,0,1)
            comm=float(np.sqrt(np.sum(2*clipped**2*(1-clipped**2))))
            meta.append(dict(in_band=BANDS[a],ob_band=BANDS[b],intersection_rank=qc.shape[1],
                             projected_rank=qz.shape[1],intersection_rank_tolerance_sensitivity=rank_sensitivity,
                             projector_commutator_frobenius=comm,
                             principal_cosines=d.tolist()))
    return dict(inn=inn,ob=ob,caps=caps,projected=projected,meta=meta)

def share(num,den):
    return float(num/den) if den else None

def writecsv(path,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def aggregate(rows,keys,energy_fields):
    groups=defaultdict(list)
    for r in rows:groups[tuple(r[k] for k in keys)].append(r)
    out=[]
    for key,rs in groups.items():
        row=dict(zip(keys,key));row["cases"]=len(rs)
        for field in energy_fields:row[field]=sum(float(r[field]) for r in rs)
        row["pooled_fraction"]=share(row["component_energy"],row["total_energy"])
        row["mean_case_fraction"]=float(np.mean([r["fraction"] for r in rs if r["fraction"] is not None]))
        out.append(row)
    return out

def main():
    begin=time.perf_counter()
    out=HERE/"results";out.mkdir(exist_ok=True)
    old=REPO/"research/research_z_alg_r4_diversity"
    meta={r["task_id"]:r for r in json.loads((old/"RESULTS_SUMMARY.json").read_text())["rows"]}
    public=dict(np.load(old/"public_data/operators.npz",allow_pickle=False))
    G,S=public["G"],public["S"]
    model=Model(public)
    cache={}
    archive=np.load(old/"predictions/endpoints.npz",allow_pickle=False)
    select=[i for i in range(len(archive["m"]))
            if str(archive["backend"][i])=="cpu" and str(archive["method"][i]) in ("native_trf","reciprocal_trf")]
    select.sort(key=lambda i:str(archive["task_id"][i]))
    cases=[];mes=[];jes=[];truths=[];estimates=[];ids=[]
    t=time.perf_counter()
    for i in select:
        tid=str(archive["task_id"][i]);r=meta[tid]
        mhat=archive["m"][i].copy()
        with np.load(old/"truth_after_freeze"/(r["observation_id"]+".npz"),allow_pickle=False) as z:
            gt=z["m"].copy()
        key=gt.tobytes()
        if key not in cache:
            state=model.forward(gt)
            cache[key]={"j":state["j"].copy(),"y":state["y"].copy(),"Q":state["Q"].copy()}
        true=cache[key]
        state=model.forward(mhat)
        dm=mhat-gt;dj=state["j"]-true["j"]
        row={k:r[k] for k in ("task_id","observation_id","method","shape","amplitude","background","noise","joint_gate","data_gate","material_gate")}
        row.update(material_RMSE=float(np.linalg.norm(dm)/np.sqrt(model.n)),
                   material_NRMSE=float(np.linalg.norm(dm)/np.sqrt(model.n)/r["amplitude"]),
                   current_error_energy=energy(dj),material_error_energy=energy(dm),
                   internal_field_error_energy=energy(G@dj),true_data_error_energy=energy(S@dj),
                   relative_true_data_error=float(np.linalg.norm(S@dj)/np.linalg.norm(true["y"])),
                   stored_observed_data_residual=r["relative_data_residual"],
                   nrmse_old=r["NRMSE"])
        cases.append(row);mes.append(dm);jes.append(dj);truths.append(gt);estimates.append(mhat);ids.append(tid)
    physics_wall=time.perf_counter()-t
    cost.append(dict(operation="exact_currents_all_GN_and_unique_GT",wall_seconds=physics_wall,
                     LU=model.counts["LU"],source_RHS=model.counts["source_RHS"],
                     actual_complex_RHS=model.counts["actual_complex_RHS"],unique_GT=len(cache),
                     endpoint_count=len(cases)))
    np.savez_compressed(out/"error_arrays.npz",task_id=np.array(ids),delta_m=np.array(mes),
                        delta_j=np.array(jes),truth_m=np.array(truths),estimate_m=np.array(estimates))
    writecsv(out/"case_metrics.csv",cases)
    print(json.dumps(dict(phase="currents",cases=len(cases),unique_GT=len(cache),
                          seconds=physics_wall,LU=model.counts["LU"])),flush=True)
    with threadpool_limits(limits=1):
        inc=spectrum(G,"complex_G_in");obc=spectrum(S,"complex_G_ob")
        inr=spectrum(np.vstack((G.real,G.imag)),"real_G_in")
        obr=spectrum(np.vstack((S.real,S.imag)),"real_G_ob")
        frames={"complex":(inc,obc),"real":(inr,obr)}
        np.savez_compressed(out/"fixed_spectra.npz",in_s=inc["s"],ob_s=obc["s"],in_V=inc["V"],ob_V=obc["V"],
                            real_in_s=inr["s"],real_ob_s=obr["s"],real_in_V=inr["V"],real_ob_V=obr["V"])
        bundles={};space_meta=[];spectra_rows=[]
        for family,(gi,go) in frames.items():
            for tau in TAUS:
                ii,oo=split(gi,tau),split(go,tau)
                bundle=joint_spaces(ii,oo);bundles[family,tau]=bundle
                space_meta.append(dict(family=family,tau=tau,groups=bundle["meta"]))
                for op in (ii,oo):
                    for band,idx in zip(BANDS,op["ids"]):
                        spectra_rows.append(dict(family=family,tau=tau,operator=op["label"],band=band,
                            rank=len(idx),dimension=op["matrix"].shape[1],sigma_max=float(op["s"][0]),
                            min_nonzero_retained=float(op["s"][op["s"]>op["s"][0]*RANK_RTOL][-1]),
                            operator_frobenius_energy_fraction=share(float(np.sum(op["s"][idx]**2)),float(np.sum(op["s"]**2)))))
        print(json.dumps(dict(phase="fixed_subspaces",seconds=time.perf_counter()-begin)),flush=True)
        single=[];two=[];remainders=[];cumulative=[]
        for case,dm,dj in zip(cases,mes,jes):
            base={k:case[k] for k in ("task_id","method","shape","amplitude","noise","joint_gate")}
            for quantity,x,family in (("material",dm[:,None],"complex"),
                                      ("current",dj,"complex"),("material_real",dm[:,None],"real")):
                total=energy(x)
                for tau in TAUS:
                    z=bundles[family,tau];inn,ob=z["inn"],z["ob"]
                    common=dict(**base,quantity=quantity,tau=tau,total_energy=total)
                    for opname,op in (("G_in",inn),("G_ob",ob)):
                        c=op["V"].conj().T@x
                        mode=np.sum(np.abs(c)**2,axis=1)
                        output=float(np.sum(op["s"]**2*mode))
                        for band,idx in zip(BANDS,op["ids"]):
                            ee=float(np.sum(mode[idx]));oe=float(np.sum(op["s"][idx]**2*mode[idx]))
                            single.append(dict(**common,operator=opname,band=band,rank=len(idx),
                                component_energy=ee,fraction=share(ee,total),
                                component_RMS=float(np.sqrt(ee/x.size)),output_component_energy=oe,
                                output_total_energy=output,output_fraction=share(oe,output)))
                        if tau==.1:
                            cumulative.append(dict(task_id=case["task_id"],quantity=quantity,operator=opname,
                                mode_energy=mode))
                    pini=components(inn,x);pobo=components(ob,x)
                    outden={"in":energy(inn["matrix"]@x),"ob":energy(ob["matrix"]@x)}
                    for order in ("ob_then_in","in_then_ob"):
                        matrix=np.zeros((3,3));oein=0.;oeob=0.
                        for first in range(3):
                            parts=components(inn,pobo[first]) if order=="ob_then_in" else components(ob,pini[first])
                            for second,part in enumerate(parts):
                                a,b=(second,first) if order=="ob_then_in" else (first,second)
                                ee=energy(part);matrix[a,b]=ee
                                extra={}
                                if tau==.1:
                                    fi=energy(inn["matrix"]@part);fo=energy(ob["matrix"]@part)
                                    oein+=fi;oeob+=fo
                                    extra=dict(in_output_energy=fi,ob_output_energy=fo,
                                               in_output_fraction_of_actual=share(fi,outden["in"]),
                                               ob_output_fraction_of_actual=share(fo,outden["ob"]))
                                two.append(dict(**common,definition=order,in_band=BANDS[a],ob_band=BANDS[b],
                                                component_energy=ee,fraction=share(ee,total),**extra))
                        remainders.append(dict(**common,definition=order,input_allocated_energy=float(matrix.sum()),
                            input_unassigned_energy=total-float(matrix.sum()),
                            in_output_cross_energy=(outden["in"]-oein if tau==.1 else ""),
                            ob_output_cross_energy=(outden["ob"]-oeob if tau==.1 else ""),
                            in_output_total=outden["in"],ob_output_total=outden["ob"]))
                    for name,spaces in (("true_intersection",z["caps"]),("tsom_projected_span",z["projected"])):
                        eeall=0.
                        for (a,b),q in spaces.items():
                            ee=energy(q.conj().T@x);eeall+=ee
                            two.append(dict(**common,definition=name,in_band=BANDS[a],ob_band=BANDS[b],
                                            rank=q.shape[1],component_energy=ee,fraction=share(ee,total)))
                        remainders.append(dict(**common,definition=name,input_allocated_energy=eeall,
                                               input_unassigned_energy=(total-eeall if name=="true_intersection" else ""),
                                               note=("orthogonal remainder retained" if name=="true_intersection" else "spaces overlap; sum is not a partition")))
        writecsv(out/"single_operator.csv",single)
        writecsv(out/"twofold.csv",two)
        writecsv(out/"twofold_remainder.csv",remainders)
        writecsv(out/"operator_groups.csv",spectra_rows)
        np.savez_compressed(out/"cumulative_mode_energies.npz",task_id=np.array([r["task_id"] for r in cumulative]),
                            quantity=np.array([r["quantity"] for r in cumulative]),operator=np.array([r["operator"] for r in cumulative]),
                            mode_energy=np.array([r["mode_energy"] for r in cumulative]))
        (out/"subspace_metadata.json").write_text(json.dumps(space_meta,indent=2))
        singleagg=aggregate(single,["method","amplitude","noise","quantity","tau","operator","band"],
                             ["component_energy","total_energy","output_component_energy","output_total_energy"])
        for r in singleagg:r["pooled_output_fraction"]=share(r["output_component_energy"],r["output_total_energy"])
        twoagg=aggregate(two,["method","amplitude","noise","quantity","tau","definition","in_band","ob_band"],
                          ["component_energy","total_energy"])
        writecsv(out/"aggregate_single.csv",singleagg)
        writecsv(out/"aggregate_twofold.csv",twoagg)
        rg=defaultdict(list)
        for r in remainders:
            rg[tuple(r[k] for k in ("method","amplitude","noise","quantity","tau","definition"))].append(r)
        rag=[]
        for key,rs in rg.items():
            a=dict(zip(("method","amplitude","noise","quantity","tau","definition"),key))
            den=sum(r["total_energy"] for r in rs)
            a.update(total_energy=den,cases=len(rs),allocated_fraction=share(sum(r["input_allocated_energy"] for r in rs),den))
            if key[-1]!="tsom_projected_span":
                a["unassigned_fraction"]=share(sum(r["input_unassigned_energy"] for r in rs),den)
            if key[-1] in ("ob_then_in","in_then_ob") and key[-2]==.1:
                for op in ("in","ob"):
                    a[op+"_output_cross_fraction"]=share(sum(r[op+"_output_cross_energy"] for r in rs),
                                                            sum(r[op+"_output_total"] for r in rs))
            rag.append(a)
        writecsv(out/"aggregate_remainder.csv",rag)
    writecsv(out/"COST_LEDGER.csv",cost)
    summary=dict(status="COMPLETED_OFFLINE_DISCLOSED_GT_DIAGNOSTIC",client_date="2026-11-06",
                 cases=len(cases),unique_GT=len(cache),methods=["native_trf","reciprocal_trf"],
                 physical_LU=model.counts["LU"],source_RHS=model.counts["source_RHS"],
                 new_GN_optimizations=0,total_wall=time.perf_counter()-begin,
                 primary_tau=.1,sensitivity_taus=[.03,.3],rank_rtol=RANK_RTOL,
                 primary_single=[r for r in singleagg if r["tau"]==.1],
                 primary_twofold=[r for r in twoagg if r["tau"]==.1],
                 primary_remainders=[r for r in rag if r["tau"]==.1],
                 operator_groups=[r for r in spectra_rows if r["tau"]==.1],
                 scope="Fixed-G right-domain geometry; material coordinate alignment is not the nonlinear material Jacobian. Current differences give exact field/data differences.")
    (out/"SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps({k:v for k,v in summary.items() if k not in ("primary_single","primary_twofold","primary_remainders","operator_groups")}),flush=True)

if __name__=="__main__":
    with threadpool_limits(limits=1):main()


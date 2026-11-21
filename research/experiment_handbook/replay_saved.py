"""Read-only replay of exposed Maxwell experiments. No optimization or data generation."""
import argparse
import csv
import json
import os
import platform
import sys
import time
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS"):
    os.environ[key] = "1"
sys.dont_write_bytecode = True
import numpy as np


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def norm(a):
    return float(np.linalg.norm(a))


def reference(root, env, short, worktree):
    if env in os.environ:
        path = Path(os.environ[env])
        require(path.is_dir(), f"{env} must be a restored study directory")
        return path
    choices = [root/"research"/short, root/"research/delegated/max_struct_r1"/worktree/"research"/short]
    for path in choices:
        if path.is_dir():
            return path
    raise FileNotFoundError(f"Restore {short}, or set {env}")


def half_replay(study, out):
    cases = sorted((c for c in read(study/"data/cases.json")
                    if (c["cohort"], c["array"]) == ("confirmation", "PW180_10x20")), key=lambda c: c["case"])
    require(len(cases) == len({c["object_id"] for c in cases}) == 24, "Expected 24 distinct primary objects")
    methods = ("dual", "srap")
    views = ("raw", "light")
    complete = {}
    # Read all paid endpoint arrays before touching any postseal truth.
    for c in cases:
        for method in methods:
            path = study/"results/predictions"/f"{c['case']}__{method}.npz"
            with np.load(path, allow_pickle=False) as z:
                for view in views:
                    a = z[f"m_30_{view}"]
                    require(a.shape == (1024,) and np.all(np.isfinite(a)), f"Incomplete {path.name}")
                    complete[c["case"], method, view] = a
    rows = []
    for c in cases:
        with np.load(study/"data/truth"/(c["case"]+".npz"), allow_pickle=False) as z:
            truth = z["m"]
            require(truth.shape == (1024,) and str(z["object_id"]) == c["object_id"], "Mismatched exposed GT")
        inside = truth > 0
        true_support = truth >= .5 * truth.max()
        den = norm(truth)
        for view in views:
            for method in methods:
                a = complete[c["case"], method, view]
                err = a-truth
                support = a >= .5*truth.max()
                rows.append(dict(case=c["case"], object_id=c["object_id"], view=view, method=method, budget=30,
                                 rel_L2=norm(err)/den, inside_L2=norm(err[inside])/den,
                                 background_L2=norm(err[~inside])/den,
                                 half_IoU=np.count_nonzero(support & true_support)/np.count_nonzero(support | true_support)))
    # Open comparison tables only after independent material metrics are complete.
    hist_path = study/"results/confirmation_metrics.csv"
    with hist_path.open(newline="", encoding="utf-8") as f:
        history = {(r["case"], r["view"], r["method"]): r for r in csv.DictReader(f)
                   if r["array"] == "PW180_10x20" and r["budget"] == "30" and r["method"] in methods}
    require(len(history) == 96, "Historical export has a different denominator")
    summary = {"objects": len(cases), "completed_prediction_files": 48, "material_rows": len(rows)}
    discrepancies = {}
    for key in ("rel_L2", "inside_L2", "background_L2", "half_IoU"):
        discrepancies[key] = max(abs(r[key] - float(history[r["case"], r["view"], r["method"]][key])) for r in rows)
    summary["max_abs_error_vs_historical_rows"] = discrepancies
    require(max(discrepancies.values()) < 1e-12, "Saved curves and historical scores disagree")
    ref = read(study/"results/paired_summary.json")["confirmation"]
    pairs = []
    for view in views:
        values = {m: [r for r in rows if r["view"] == view and r["method"] == m] for m in methods}
        gains = []
        for first, second in zip(values["dual"], values["srap"]):
            require(first["case"] == second["case"], "Mismatched pairing")
            gain = (first["rel_L2"]-second["rel_L2"])/first["rel_L2"]
            gains.append(gain)
            pairs.append(dict(object_id=first["object_id"], view=view, rel_L2_fractional_gain=gain,
                              **{f"{m}_{k}": r[k] for m,r in zip(methods,(first,second))
                                 for k in ("rel_L2","inside_L2","background_L2","half_IoU")}))
        g = np.asarray(gains)
        rng = np.random.default_rng(20261121)
        ci = np.percentile(np.median(g[rng.integers(0,len(g),(10000,len(g)))], axis=1), [2.5,97.5])
        result = dict(wins=int(sum(g>1e-10)), ties=int(sum(abs(g)<=1e-10)), losses=int(sum(g< -1e-10)),
                      median_gain=float(np.median(g)), ci95=ci.tolist())
        for m,prefix in zip(methods, ("baseline", "srap")):
            result[prefix+"_median"] = float(np.median([r["rel_L2"] for r in values[m]]))
            result[m+"_better_air"] = sum(r["rel_L2"] < 1-1e-10 for r in values[m])
            for metric,name in (("inside_L2","inside"),("background_L2","background"),("half_IoU","iou")):
                result[f"{name}_{m}"] = float(np.median([r[metric] for r in values[m]]))
        historical = ref[f"PW180_10x20 30 {view}"]
        require(all(np.allclose(v, historical[k], rtol=0, atol=1e-12) for k,v in result.items()), "Paired summary discrepancy")
        summary[view] = result
    with (out/"half24_pairs.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(pairs[0]), lineterminator="\n"); w.writeheader(); w.writerows(pairs)
    return summary


def r11_replay(root, r12):
    p = r12/"supplement_paired"
    seal = read(p/"ALL44_PREDICTIONS_SEALED.json")
    require(seal["root_authorized_offline_scoring"] and seal["tasks"] == 44, "Historical cohort seal missing")
    cases = [c for c in read(p/"inputs/manifest.json")["cases"] if c["cohort"] == "A"]
    rows=[]
    for case in cases:
        key=case["case_id"]
        with np.load(p/case["input_file"], allow_pickle=False) as z:
            Y, E, G, S, m0 = [z[k] for k in ("Y","E0","G","S","m_init")]
        require((Y.shape,E.shape,G.shape,S.shape,m0.shape)==((48,12),(576,12),(576,576),(48,576),(576,)), "A contract")
        # The old manifest contains machine-specific provenance. Resolve by relative task id.
        r10dir = os.environ.get("R10_DIR")
        orig = [Path(r10dir)] if r10dir else [root/"research/research_max_struct_r10", root/"research/delegated/max_struct_r1/r10_worktree/research/research_max_struct_r10"]
        candidates = [d/"RESULTS/twod"/case["parent_task_id"]/"prediction.npz" for d in orig]
        parent_checked = False
        for path in candidates:
            if path.is_file():
                with np.load(path, allow_pickle=False) as z:
                    require(np.array_equal(z["m"],m0), "Original parent/start material differ")
                parent_checked=True
                break
        for method in ("g19","ap","dr","direct"):
            op=p/"RESULTS"/f"A__{key}__{method}"
            with np.load(op/"initial.npz", allow_pickle=False) as z:
                require(np.array_equal(z["m"],m0), "Arms do not share original initial m")
            with np.load(op/"primary.npz", allow_pickle=False) as z:
                m,j,q,yhat=[z[k] for k in ("m","j","Q","Y_prediction")]
            receipt=read(op/"primary.json")
            v=norm(yhat-Y)/norm(Y)
            require(abs(v-receipt["relative_residual"])<1e-9, "A receiver residual mismatch")
            field=norm(q-E-G@j)/norm(q); port=norm(S@j-yhat)/norm(yhat)
            require(field<1e-10 and port<1e-10, "Saved fields violate acquisition contract")
            rows.append(dict(object=key,parent=case["parent_method"],method=method,original_parent_array_checked=parent_checked,
                             m0_shape=list(m0.shape),Y_shape=list(Y.shape),data_rel=v,field_identity_error=field,
                             receive_identity_error=port,selected_phi=receipt["Phi"],paid_candidate_cpu=receipt["CPU_seconds"]))
    return dict(sealed_children=44, same_start_verified_arms=len(rows), records=rows, truth_reads=0)


def austria_replay(r8):
    base=r8/"RESULTS/paper_benchmarks";case="austria_eps4_400MHz";data=base/case
    receipt=read(data/"DATA_RECEIPT.json")
    with np.load(data/"PUBLIC_DATA.npz",allow_pickle=False) as z:
        Y=z["Y"]; freq=float(z["frequency"])
    require(Y.shape==(40,20) and freq==400e6 and receipt["forward_side"]==150, "Austria input dimensions")
    rows=[]
    for side, folder in ((32,data/"cie_dual"),(64,base/"original64"/case/"cie_dual")):
        require((folder/"PREDICTION_SEAL.json").is_file(),"Missing Austria seal")
        with np.load(folder/"PREDICTION.npz",allow_pickle=False) as z:
            m,yp=z["m"],z["Y_prediction"]
        r=read(folder/"RECEIPT.json");value=norm(yp-Y)/norm(Y)
        require(m.shape==(side**2,) and yp.shape==Y.shape and abs(value-r["relative_data_residual"])<1e-12,"Austria residual discrepancy")
        rows.append(dict(inverse_side=side,material_unknowns=m.size,data_rel=value,cpu_s=r["cpu_seconds"],wall_s=r["wall_seconds"],
                         min_chi=float(m.min()),max_chi=float(m.max())))
    return dict(case=case,forward_side=150,noise_fraction=receipt["noise_relative"],Y_shape=list(Y.shape),frequency_Hz=freq,method="cie_dual",rows=rows,truth_reads=0)


def fresnel_replay(root,study,r8):
    # Import only the documented pure data reader and empty-target calibration.
    sys.path.insert(0,str(r8.parents[1]))
    from research.research_max_struct_r8.src.fresnel_adapter import load_file,incident_fit
    raw_root=Path(os.environ.get("FRESNEL_2005_DIR",root/"research/delegated/a2_literature/data/2005_iop_21_6_S09"))
    measurement=Path(os.environ.get("FRESNEL_INPUTS",study/"data"))
    raw,Y,I=load_file(raw_root/"FoamTwinDielTM.exp",8.)
    require(raw.shape==(4338,7) and Y.shape==I.shape==(360,18),"Twin raw data contract")
    a,b=raw[:,0].astype(int)-1,raw[:,1].astype(int)-1
    require(np.array_equal(Y[b,a],((raw[:,3]-raw[:,5])-1j*(raw[:,4]-raw[:,6]))),"Raw field convention")
    require(np.array_equal(I[b,a],raw[:,5]-1j*raw[:,6]),"Empty-field phase")
    def coordinates(n):
        angles=np.deg2rad(np.arange(n)*360/n)
        return 1.67*np.stack([np.cos(angles),np.sin(angles)],1)
    rows=[];selectors=[]
    for layout in ("contiguous","interleaved"):
        case="fresnel_twin_8GHz__"+layout
        with np.load(measurement/"public"/(case+".npz"),allow_pickle=False) as z:
            meta=json.loads(str(z["case_json"])); ids=z["rx_ids"]; mask=z["valid_mask"]; yt=z["Y"]
            tx=np.array(meta["train_sources"])-1
            if not rows:
                e,cal,_=incident_fit(I,coordinates(18),coordinates(360),z["centers"],float(z["k"]))
            cal_error=norm(e[:,tx]-z["E"])/norm(z["E"])
            require(z["S"].shape==(360,3600) and ids.shape==mask.T.shape and cal_error<1e-10, "Twin incident/receiver contract")
        pairs=set()
        for j,s in enumerate(tx):
            receivers=ids[j,mask[:,j]]
            require(np.array_equal(Y[receivers,s],yt[mask[:,j],j]),"Source-specific packed channel mismatch")
            pairs.update((int(s),int(r)) for r in receivers)
        with np.load(measurement/"heldout"/(case+".npz"),allow_pickle=False) as z:
            held_tx=z["source_ids"];seen=z["source_ids_train_hemi"];y1=z["Y"];y2=z["Y_train_opposite_rx"]
            require(np.array_equal(Y[:,held_tx],y1,equal_nan=True),"Held source field disagrees with raw")
            require(np.array_equal(Y[:,seen][np.isfinite(y2)],y2[np.isfinite(y2)]),"Held receive field disagrees with raw")
            h1={(int(held_tx[j]),int(r)) for r,j in np.argwhere(np.isfinite(y1))}
            h2={(int(seen[j]),int(r)) for r,j in np.argwhere(np.isfinite(y2))}
        require(not (set(tx)&set(held_tx) or pairs & (h1|h2) or h1&h2), "Leakage in Twin masks")
        selectors.append((pairs,h1,h2))
        rows.append(dict(layout=layout,Y_packed_shape=list(yt.shape),source_dependent_mask_shape=list(mask.shape),training_pairs=len(pairs),
                         held_source_pairs=len(h1),held_receiver_pairs=len(h2),training_source_ids_zero_based=tx.tolist(),
                         held_source_ids_zero_based=held_tx.tolist(),shared_seen_source_ids_zero_based=seen.tolist(),
                         calibrated_incident_rel_diff=cal_error,empty_validation_median=float(np.median([r["validation_relative"] for r in cal]))))
    require(selectors[0][1:]==selectors[1][1:] and all(row["training_pairs"]==733 for row in rows),"Unmatched Fresnel test")
    return dict(raw_rows=raw.shape[0],frequency_GHz=8,original_source_count=18,receiver_slots=360,finite_pairs_per_source=np.isfinite(Y).sum(0).tolist(),variants=rows,material_truth=False)


def vector_replay(study):
    v=study/"vector"
    seal=read(v/"COHORT_SEAL.json")
    require(len(seal["arms"])==48 and "hemi_00_SRAP" in seal["arms"],"Vector cohort incomplete")
    with np.load(v/"data/hemi_00.npz",allow_pickle=False) as z:
        Y=z["Y"];delta=float(z["delta"])
    with np.load(v/"sealed_truth/00.npz",allow_pickle=False) as z: # Published 12-case cohort only.
        coarse=z["chi"];fine=z["chi24"]
    require(fine.size==24**3 and coarse.size==12**3 and np.array_equal(coarse,fine.reshape(12,2,12,2,12,2).mean((1,3,5)).ravel()),"Vector coarse scoring truth")
    acq=read(v/"acquisition.json")["hemi"]
    for role,n in (("transmit",8),("receive",16)):
        pol=np.asarray(acq[role]["polarizations"])
        require(pol.shape==(n,2,3) and np.allclose(pol@pol.transpose(0,2,1),np.eye(2)),"Vector polarization")
    with (v/"objects.csv").open(newline="",encoding="utf-8") as f:
        historical={r["method"]:r for r in csv.DictReader(f) if (r["object"],r["array"],r["budget"],r["view"])==("0","hemi","180","light")}
    scores={}
    for method in ("Dual","SRAP"):
        p=v/"predictions"/("hemi_00_"+method)
        r=read(p/"SEALED.json")
        with np.load(p/"materials.npz",allow_pickle=False) as z:m=z["light_180"]
        err=norm(m-coarse)/norm(coarse)
        require(abs(err-float(historical[method]["rel_L2"]))<1e-12 and Y.shape==(32,16),"Vector score mismatch")
        scores[method]=dict(relative_material_L2=err,physical_phi=r["selected"]["light_180"]["phi"],route_rank=r["route"].get("route_rank",0))
    return dict(case="hemi_00",source_excitations=16,receive_components=32,observed_complex_entries=Y.size,
                material_cells=coarse.size,current_components_per_cell=3,raw_field_data_shape=list(Y.shape),delta=delta,light180=scores)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("scene",choices=("all","half24","r11","austria","fresnel","vector"))
    ap.add_argument("--root",type=Path,default=Path.cwd(),help="Inv_SLAM or restored public repository root")
    ap.add_argument("--output",type=Path,default=Path(".artifacts/experiment_handbook_replay"))
    args=ap.parse_args();out=args.output;out.mkdir(exist_ok=True,parents=True)
    root=args.root.resolve()
    def get8():return reference(root,"R8_DIR","research_max_struct_r8","r8_worktree")
    def get12():return reference(root,"R12_DIR","research_max_struct_r12","r12_worktree")
    study=Path(os.environ.get("SRAP_DIR",root/"research/srap_cie_array_study_20261121"))
    tasks={"r11":lambda:r11_replay(root,get12()),"austria":lambda:austria_replay(get8()),"fresnel":lambda:fresnel_replay(root,study,get8()),"vector":lambda:vector_replay(study),"half24":lambda:half_replay(study,out)}
    for name,task in tasks.items():
        if args.scene not in ("all",name):continue
        start=time.perf_counter();cpu=time.process_time();data=task()
        record=dict(scene=name,wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,
                    python=platform.python_version(),numpy=np.__version__,new_inversions=0,new_Maxwell_solves=0,result=data)
        (out/(name+"_replay.json")).write_text(json.dumps(record,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")
        print(f"{name}: {record['wall_s']:.4f}s wall, {record['cpu_s']:.4f}s CPU",flush=True)
if __name__=="__main__":main()

#!/usr/bin/env python3
"""Research T E0: isolated exposed saved-array audit; no physical solve.

This executable deliberately imports no research solver or online proposal
module.  It never writes upstream data.  Truth is read only for this offline
audit, not exported as search initialization. Matrix applications below test
already-saved states; they are not new state, tangent, or adjoint solves.
"""
from __future__ import annotations
import csv
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

START, CPU = perf_counter(), process_time()
import numpy as np
from scipy import ndimage
from scipy.stats import chi2
import scipy.linalg
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
GAUSSIAN = ROOT.parents[1]
A30 = GAUSSIAN / "A30/a30-maxwell-material-manifold"
REL = A30 / "inputs/a29_release"
OLD_A = ROOT / "research/research_a_20261015/results/E2"
OUT = ROOT / "research/research_t_20261015/e0"
A_CASES = [f"p{p}_{b}_a1_snr20" for p in (250065, 250069, 250073, 250077) for b in ("weak", "ring")]
CASES = [f"p{p}_{b}_a1_snr20" for p in range(250065, 250081) for b in ("weak", "ring")]
TAU = float(np.sqrt(chi2.ppf(.95, 512)))


def forbidden(*args, **kwargs):
    raise RuntimeError("E0_FORBIDS_ANY_LINEAR_SOLVE_OR_FACTORIZATION")


# Defense against accidental future additions calling a physical solver.
for obj, names in ((np.linalg, ("solve", "inv", "pinv", "lstsq", "cholesky")),
                   (scipy.linalg, ("solve", "lu_solve", "lu_factor", "cho_solve", "inv", "lstsq"))):
    for name in names:
        setattr(obj, name, forbidden)


def read(path, keys=None):
    path = Path(path)
    if any("locked" in p.lower() for p in path.parts):
        raise ValueError("LOCKED_PATH_FORBIDDEN")
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in (z.files if keys is None else keys)}


def dump(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


def csv_write(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v) if isinstance(v, (list, dict)) else v for k, v in row.items()})


def topo(mask):
    mask = np.asarray(mask, bool).reshape(24, 24)
    _, components = ndimage.label(mask, ndimage.generate_binary_structure(2, 1))
    labels, count = ndimage.label(~mask, ndimage.generate_binary_structure(2, 2))
    edge = np.r_[labels[0], labels[-1], labels[:, 0], labels[:, -1]]
    exterior = set(map(int, edge)) - {0}
    return dict(components=int(components), holes=int(count)-len(exterior), pixels=int(mask.sum()))


def relative(delta, scale):
    return float(np.linalg.norm(delta) / max(float(np.linalg.norm(scale)), 1e-300))


def compatible(chi):
    return bool(np.isfinite(chi).all() and np.all(chi.real >= -.5-1e-10)
                and np.all(chi.real <= 4+1e-10) and np.all(chi.imag >= -1e-10)
                and np.all(chi.imag <= .5+1e-10))


def checks(material, current, field, data, bg):
    # c=chi*e, e=incident+G*c, y=S*c.  No state reconstruction.
    prediction = bg["S"] @ current
    field_check = bg["incident"] + bg["G"] @ current
    constitutive = current - material[:, None] * field
    rhs = material[:, None] * bg["incident"]
    equation = current - material[:, None] * (bg["G"] @ current) - rhs
    return dict(constitutive_relative=relative(constitutive, current),
                field_relative=relative(field-field_check, field),
                observed_prediction_relative=relative(data-prediction, data),
                backward_relative=relative(equation, rhs))


def count_masks(chi, truth, background):
    mag = np.abs(chi-background)
    return {str(frac): topo(mag > frac * np.max(np.abs(truth-background))) for frac in (.15, .25, .35)}


def phase_table(chi):
    # Rounded phase frequencies describe saved arrays, not a new material prior.
    vals, counts = np.unique(np.round(chi, 8), return_counts=True)
    order = np.argsort(-counts, kind="stable")[:12]
    return [{"chi_re": float(vals[i].real), "chi_im": float(vals[i].imag), "pixels": int(counts[i])} for i in order]


def make_figure(case, truth, n3, base, background, norms):
    fig, axes = plt.subplots(3, 5, figsize=(17, 10.0), constrained_layout=True)
    labels = ("OFFLINE truth", "Historical N3", "Research A base (one GN step)")
    ims = []
    for row, (name, chi, residual) in enumerate(zip(labels, (truth, n3, base), norms)):
        images = ((1+chi.real, .5, 5, "viridis", "Absolute Re(epsilon_r)=1+Re(chi)"),
                  (chi.imag, 0, .5, "magma", "Absolute Im(epsilon_r)=Im(chi)"),
                  ((np.abs(chi) <= 1e-10).astype(float), 0, 1, "gray_r", "Strict vacuum: |chi| <= 1e-10"),
                  (np.abs(chi-background), 0, 1.2, "viridis", "Historical magnitude |chi-background|"),
                  ((np.abs(chi-background) > .25*np.max(np.abs(truth-background))).astype(float),
                   0, 1, "gray_r", "OFFLINE increment mask (25% truth max)"))
        for col, (image, lo, hi, cmap, title) in enumerate(images):
            ax = axes[row, col]
            im = ax.imshow(image.reshape(24, 24).T, origin="lower", extent=(-.5,.5,-.5,.5),
                           vmin=lo, vmax=hi, cmap=cmap, interpolation="nearest")
            if row == 0:
                ax.set_title(title, fontsize=9)
            if col == 0:
                ax.set_ylabel(name + f"\nobserved whitened norm={residual:.4f}\ny [m]", fontsize=9)
            ax.set_xlabel("x [m]")
            ax.set_aspect("equal")
            if row == 2:
                fig.colorbar(im, ax=axes[:, col].tolist(), shrink=.55, location="bottom", pad=.03)
    fig.suptitle(f"{case}: absolute constitutive phase vs increment support\n"
                 f"All columns use fixed scales; air is not dark increment. Frozen tau={TAU:.6f}", fontsize=13)
    for suffix in ("png", "pdf"):
        fig.savefig(OUT / f"PHASE_AUDIT_{case}.{suffix}", dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    targets = read(REL / "results/a29/OFFLINE_TARGETS.npz")
    validation = read(REL / "inputs/a25/validation.npz")
    arrays = read(OLD_A / "MATERIAL_AND_PREDICTION_ARRAYS.npz")
    histories = json.loads((OLD_A / "SEARCH_HISTORY.json").read_text())
    protocol = json.loads((OLD_A / "PROTOCOL.json").read_text())
    original_geom = json.loads(str(validation["geometry_descriptor"]))
    t_ids, v_ids = list(targets["case_id"].astype(str)), list(validation["case_id"].astype(str))
    rows, candidate_rows, topologies, missing, phase_rows = [], [], [], [], []
    offline_arrays = {}
    case_map = {}
    for case in CASES:
        i, vi = t_ids.index(case), v_ids.index(case)
        truth, background = targets["truth"][i], targets["background"][i]
        pool_path = REL / "results/a29/pools" / (case + ".npz")
        cache_path = REL / "inputs/a28/results/a28/exact_states" / (case + ".npz")
        anatomy_path = REL / "inputs/a28/results/a28/anatomy" / (case + ".npz")
        pool, cache, anatomy = read(pool_path), read(cache_path), read(anatomy_path)
        bgname = case.split("_")[1]
        bg = read(REL / f"inputs/a25/background/{bgname}.npz",
                  ("G", "S", "incident", "yb", "chi_b", "grid", "volume", "geometry_descriptor"))
        y, clean, sigma = pool["observation"], validation["clean"][vi], float(pool["sigma"])
        n3 = pool["n3"]
        base = arrays[case + "__base"] if case in A_CASES else None
        white = np.sqrt(2) / sigma
        residual = lambda pred: float(np.linalg.norm((pred-y) * white))
        n3_norm = residual(cache["data"][0])
        base_norm = residual(arrays[case + "__base_prediction"]) if base is not None else None
        noise_norm = residual(clean)
        rng = np.random.default_rng(np.random.SeedSequence(validation["noise_seed"][vi].tolist()))
        replay_noise = sigma/np.sqrt(2) * (rng.standard_normal(y.shape)+1j*rng.standard_normal(y.shape))
        truth_names = anatomy["TRUTH_CENTERED_ORACLE_names"].astype(str).tolist()
        ti = truth_names.index("CENTER")
        t_checks = checks(truth, anatomy["TRUTH_CENTERED_ORACLE_current"][ti],
                          anatomy["TRUTH_CENTERED_ORACLE_field"][ti], anatomy["TRUTH_CENTERED_ORACLE_F"][ti], bg)
        n_checks = checks(n3, cache["current"][0], cache["field"][0], cache["data"][0], bg)
        geom_same = original_geom == json.loads(str(bg["geometry_descriptor"]))
        identity = dict(original_N3_matches_pool0=bool(np.array_equal(n3, pool["chi"][0])),
                        original_N3_matches_offline=bool(np.array_equal(n3, targets["N3"][i])),
                        original_N3_matches_cache=bool(np.array_equal(n3, cache["chi"][0])),
                        original_N3_cache_maxabs=float(np.max(np.abs(n3-cache["chi"][0]))),
                        original_N3_cache_numerically_matches=bool(np.allclose(n3,cache["chi"][0],rtol=0,atol=1e-14)),
                        A_N3_matches_original=bool(np.array_equal(n3, arrays[case+"__N3"])) if base is not None else None,
                        A_base_matches_original_N3=bool(np.array_equal(n3, base)) if base is not None else None,
                        A_base_minus_N3_maxabs=float(np.max(np.abs(base-n3))) if base is not None else None,
                        observation_A_matches_original=bool(np.array_equal(y, arrays[case+"__observation"])) if base is not None else None,
                        observation_matches_validation=bool(np.array_equal(y, validation["data"][vi])),
                        sigma_matches_validation=bool(sigma == float(validation["sigma"][vi])),
                        truth_matches_validation=bool(np.array_equal(truth, validation["chi"][vi])),
                        truth_matches_anatomy=bool(np.array_equal(truth, anatomy["TRUTH_CENTERED_ORACLE_chi"][ti])),
                        generator_clean_matches_anatomy_F=bool(np.array_equal(clean, anatomy["TRUTH_CENTERED_ORACLE_F"][ti])),
                        background_matches_original=bool(np.array_equal(background, bg["chi_b"])),
                        geometry_equal=geom_same)
        clean_diff = float(np.linalg.norm((clean-anatomy["TRUTH_CENTERED_ORACLE_F"][ti])*white))
        statuses = []
        if base_norm is not None and base_norm > TAU: statuses.append("INVALID_BASE")
        if base_norm is not None and base_norm > TAU: statuses.append("OPTIMIZATION_FAILURE_SUSPECTED_NOT_PROVEN")
        if case.startswith("p250065_"): statuses.append("TOPOLOGY_MISINTERPRETATION")
        row = dict(case_id=case, parent_id=int(targets["parent_id"][i]), family=str(targets["family"][i]),
                   background_id=bgname, sigma_complex_RMS=sigma, sigma_per_real=sigma/np.sqrt(2),
                   tau=TAU, tau_definition="sqrt(chi2.ppf(.95,512)); fixed by ResearchA, not an A30 original reconstruction gate",
                   N3_observed_residual=n3_norm, A_base_observed_residual=base_norm,
                   N3_box_legal=compatible(n3), A_base_box_legal=compatible(base) if base is not None else None, truth_box_legal=compatible(truth),
                   N3_status="VALID_NUMERICAL" if n3_norm <= TAU and compatible(n3) else "INVALID",
                   A_base_data_status="NOT_AUDITED" if base is None else ("VALID_NUMERICAL" if base_norm <= TAU and compatible(base) else "INVALID"),
                   Research_A_status="RUN_EXPOSED" if base is not None else "NOT_RUN",
                   A_base_state_audit_status="NOT_AUDITED", A_base_state_audit_reason="Full endpoint prediction saved; base fields/backward residual not retained per endpoint",
                   truth_noise_only_norm=noise_norm, truth_noise_chi2_upper_tail=float(chi2.sf(noise_norm**2,512)),
                   truth_noise_exceeds_tau=bool(noise_norm > TAU),
                   metadata_noise_replay_maxabs=float(np.max(np.abs(y-(clean+replay_noise)))),
                   original_sigma_formula=float(np.sqrt(np.mean(np.abs(clean-bg["yb"])**2))*.1),
                   generator_vs_inversion_clean_white_discrepancy=clean_diff,
                   truth_constitutive_relative=t_checks["constitutive_relative"], truth_backward_relative=t_checks["backward_relative"],
                   N3_constitutive_relative=n_checks["constitutive_relative"], N3_backward_relative=n_checks["backward_relative"],
                   source_geometry_equal=geom_same, identity=identity, root_cause_categories=statuses,
                   INVALID_BASE=bool(base_norm > TAU) if base is not None else None,
                   TOPOLOGY_MISINTERPRETATION="CONFIRMED_AIR_VS_INCREMENT_RISK" if case.startswith("p250065_") else "NOT_ESTABLISHED_PER_CASE",
                   OPTIMIZATION_FAILURE="SUSPECTED_NOT_PROVEN" if base_norm is not None and base_norm > TAU else "NOT_ESTABLISHED",
                   REPRESENTATION_FAILURE="OFFLINE_CAPACITY_GAP_NOT_ONLINE_CAUSE_PROOF",
                   POSSIBLE_AMBIGUITY="NOT_ESTABLISHED_BY_SINGLE_BASE",
                   actual_output_error_bound=None, actual_output_error_bound_status="NOT_AUDITED_NO_STABILITY_BOUND",
                   material_source=str(pool_path.relative_to(A30)), observation_source="inputs/a25/validation.npz and A29 pool exact identity",
                   N3_state_source=str(cache_path.relative_to(A30)),
                   original_geometry=dict(grid=24,cell_volume=1/576,side_m=1, k=2*np.pi,sources=8,receivers=32))
        row["MODEL_MISMATCH"] = "NO_CACHED_EVIDENCE; truth clean arrays and saved state equations agree"
        rows.append(row); case_map[case] = dict(row=row, truth_checks=t_checks, N3_checks=n_checks)
        materials=[("truth_OFFLINE",truth),("N3",n3),("background",background)]
        if base is not None: materials.append(("A_base",base))
        for name, image in materials:
            masks = count_masks(image, truth, background)
            phase_rows.append(dict(case_id=case, material=name, real_min=float(image.real.min()), real_max=float(image.real.max()),
                imag_min=float(image.imag.min()),imag_max=float(image.imag.max()),
                strict_vacuum_pixels=int(np.count_nonzero(np.abs(image)<=1e-10)),
                phase_frequencies_rounded_8decimals=phase_table(image)))
            for frac, result in masks.items():
                topologies.append(dict(case_id=case, material=name, threshold_fraction=frac, **result,
                                       definition="abs(chi-background)>fraction*max_abs_truth_increment; OFFLINE ONLY",foreground_connectivity=4,background_connectivity=8))
            offline_arrays[case+"__"+name+"__chi"]=image
        # Same eight observations, every retained A endpoint and old N3.
        percase = [h for h in histories if h["case_id"] == case]
        for h in [dict(arm="historical_N3", trajectory=-1,step=-1)] + percase:
            if h["arm"] == "historical_N3":
                image, pred, name, backward = n3, cache["data"][0], "historical_N3", n_checks["backward_relative"]
                source_ref = str(cache_path.relative_to(A30)); receipt_norm=None
            elif h["arm"] == "base":
                image,pred,name,backward=base,arrays[case+"__base_prediction"],"A_base",None
                source_ref=f"{OLD_A.relative_to(ROOT)}/MATERIAL_AND_PREDICTION_ARRAYS.npz:{case}__base"
                receipt_norm=float(h["residual"])
            else:
                name=h["array_prefix"];image,pred=arrays[name+"__chi"],arrays[name+"__prediction"]
                backward=h.get("field_backward_residual");source_ref=f"{OLD_A.relative_to(ROOT)}/MATERIAL_AND_PREDICTION_ARRAYS.npz:{name}"
                receipt_norm=float(h["residual"])
            rn=residual(pred);legal=compatible(image);valid=legal and rn<=TAU
            full_status = "INVALID" if not valid else ("VALID_NUMERICAL" if backward is not None and backward <=1e-10 else "NOT_AUDITED")
            candidate_rows.append(dict(case_id=case,candidate_id=name,arm=h["arm"],trajectory=h.get("trajectory"),step=h.get("step"),
                material_array_ref=source_ref, matching_prediction_saved=True, material_shared_all_8_sources=True,
                calibration_gain=1.,calibration_scope="fixed known calibration; original observations unchanged",
                material_box_legal=legal, observed_residual_norm=rn,tau=TAU,data_compatible=bool(rn<=TAU),
                validation_status=full_status,backward_residual=backward,
                constitutive_residual=n_checks["constitutive_relative"] if h["arm"]=="historical_N3" else None,
                constitutive_status="CHECKED_CACHED_ARRAYS" if h["arm"]=="historical_N3" else "NOT_AUDITED_FIELDS_NOT_SAVED",
                output_error_upper_bound=None,residual_receipt_absolute_difference=None if receipt_norm is None else abs(rn-receipt_norm),
                minimum_re=float(image.real.min()),maximum_re=float(image.real.max()),minimum_im=float(image.imag.min()),maximum_im=float(image.imag.max()),
                active_real_lower=int(np.count_nonzero(image.real<=-.5+1e-8)),active_real_upper=int(np.count_nonzero(image.real>=4-1e-8)),
                active_imag_lower=int(np.count_nonzero(image.imag<=1e-8)),active_imag_upper=int(np.count_nonzero(image.imag>=.5-1e-8)),
                strict_vacuum_pixels=int(np.count_nonzero(np.abs(image)<=1e-10)),
                increment_mask_pixels=int(np.count_nonzero(np.abs(image-background)>.18)),
                increment_changed_cells_vs_A_base=int(np.count_nonzero((np.abs(image-background)>.18)!=(np.abs(base-background)>.18))) if base is not None else None,
                historical_topology_025=count_masks(image,truth,background)["0.25"],
                source_recorded_new_states=h.get("fresh_states"),new_T_audit_states=0,new_T_audit_RHS=0))
        missing.append(dict(case_id=case,missing=[
            "Research A base/candidate fields/current unavailable, or A not run for this case",
            "rigorous output/stability error bound for all saved endpoints",
            "A30 Sparse64 full-forward endpoint predictions (offline capacity only)"],
            available=["A25 clean truth/current+observation+noise seed", "A28 exact N3 state+LU+prediction", "A28 exposed anatomy truth+field+current+F", "Research A all endpoint material+full predictions"]))
        if case == CASES[0]:
            make_figure(case,truth,n3,base,background,(noise_norm,n3_norm,base_norm))
    csv_write(OUT/"A30_CASE_AUDIT.csv",rows)
    csv_write(OUT/"CANDIDATE_AUDIT.csv",candidate_rows)
    csv_write(OUT/"HISTORICAL_MASK_TOPOLOGY.csv",topologies)
    csv_write(OUT/"ABSOLUTE_PHASE_AUDIT.csv",phase_rows)
    np.savez_compressed(OUT/"OFFLINE_AUDIT_MATERIAL_ARRAYS.npz",**offline_arrays)
    dump(OUT/"MISSING_CACHE_MANIFEST.json",dict(status="NO_CORE_CONTRACT_BLOCKER",cases=missing,
        missing_is_not_zero=True, new_solves_to_fill_gaps=0))
    original_pair_rows=list(csv.DictReader((A30/"results/a30/compensation_corrected/COMPENSATION_PAIRS.csv").open()))
    csv_write(OUT/"HISTORICAL_PAIR_DISTANCE_CONTEXT.csv",[r for r in original_pair_rows if r["case_id"] in CASES])
    capacity=json.loads((A30/"results/a30/capacity/CAPACITY_SUMMARY.json").read_text())
    capacity_rows=list(csv.DictReader((A30/"results/a30/capacity/CAPACITY_CASES.csv").open()))
    csv_write(OUT/"FROZEN_CAPACITY_ROWS.csv",[r for r in capacity_rows if r["case_id"] in CASES])
    capimages=read(A30/"results/a30/capacity/CAPACITY_IMAGES.npz",("case_id","truth","background","N3","sparse_64"))
    capacity_bycase=[]
    for case in CASES:
        j=list(capimages["case_id"].astype(str)).index(case)
        sparse, truth,bg,n3=[capimages[k][j] for k in ("sparse_64","truth","background","N3")]
        capacity_bycase.append(dict(case_id=case,sparse64_delta_nrmse=relative(sparse-truth,truth-bg),
            N3_delta_nrmse=relative(n3-truth,truth-bg),sparse64_box_legal=compatible(sparse),
            observed_data_residual=None,full_endpoint_available=False,
            fitting_objective="OFFLINE truth-target quadratic under frozen N3/TRAIN-selected Haar support; not observed-data inversion"))
    csv_write(OUT/"SPARSE64_CAPACITY_VS_RECOVERY.csv",capacity_bycase)
    n3_valid=sum(r["N3_status"]=="VALID_NUMERICAL" for r in rows)
    base_valid=sum(r["A_base_data_status"]=="VALID_NUMERICAL" for r in rows)
    numeric=max(v for c in case_map.values() for key in ("truth_checks","N3_checks") for v in c[key].values())
    failedbase=[r["case_id"] for r in rows if r["A_base_data_status"]=="INVALID"]
    summary=dict(status="E0_COMPLETE_ZERO_NEW_PHYSICS",new_full_state_attempts=0,new_complex_full_solver_RHS=0,
        new_factorizations=0,solver_modules_imported=[],locked_truth_reads=0,old_results_modified=False,
        original_revision_A30="717d3d26a290248bd2b9aeb55f3d03694b407b84",Research_A_revision="f1b674eb4893e39126701ed48f3d3dcef3098c3c",
        cases=32,Research_A_cases=8,Research_A_not_run_cases=24,original_N3_data_qualified=n3_valid,
        N3_cache_material_maxabs=max(r["identity"]["original_N3_cache_maxabs"] for r in rows),
        fixed8_original_N3_data_qualified=sum(r["N3_status"]=="VALID_NUMERICAL" for r in rows if r["case_id"] in A_CASES),
        A_base_data_qualified=base_valid,A_base_unqualified=8-base_valid,
        no_original_N3_is_Research_A_base=True,truth_noise_outliers=[r["case_id"] for r in rows if r["truth_noise_exceeds_tau"]],
        noise_metadata_replay_maxabs=max(r["metadata_noise_replay_maxabs"] for r in rows),
        max_cached_physical_identity_relative_error=numeric,generator_clean_prediction_discrepancy_max=0.,
        candidate_rows=len(candidate_rows),Research_A_rows=len(histories),
        Research_A_data_valid_endpoint_count=sum(r["data_compatible"] and r["material_box_legal"] for r in candidate_rows if r["arm"]!="historical_N3"),
        sparse64_offline_complex_gain=capacity["sparse_64"]["complex"]["point"],
        model_mismatch_status="NOT_OBSERVED_IN_CACHED_SAME_DATA_IDENTITIES",core_contract_recoverable=True,
        paid_experiment_blocker=False,qualified_handoff_caveat="New qualification must retain matching endpoint state/check receipt; saved A base fields unavailable",
        cached_matrix_vector_actions=dict(G=4*len(CASES)*8,S=2*len(CASES)*8),
        failed_base_cases=failedbase,cases_detail=case_map,
        truth_information_firewall="This executable is separate offline scoring only. Do not import or use OFFLINE_AUDIT_MATERIAL_ARRAYS or figures as online edit coordinates.",
        wall_seconds=perf_counter()-START,process_cpu_seconds=process_time()-CPU)
    dump(OUT/"E0_AUDIT.json",summary)
    write_reports(rows,summary,topologies)
    receipt=dict(status="COMPLETE",pid=os.getpid(),new_full_state_attempts=0,
        new_complex_full_solver_RHS=0,new_factorizations=0,GPU_seconds=0,old_ledgers_modified=False,
        implementation="No physics import; solve/factorization entry points disabled; saved matrix applications only",
        wall_seconds=perf_counter()-START,process_cpu_seconds=process_time()-CPU)
    receipt_path=OUT/"EXECUTION_RECEIPT.json"
    ledger_path=OUT/"E0_EXECUTION_RECEIPTS.jsonl"
    if receipt_path.exists() and not ledger_path.exists():
        with ledger_path.open("a") as f:
            f.write(json.dumps(json.loads(receipt_path.read_text()))+"\n")
    with ledger_path.open("a") as f:
        f.write(json.dumps(receipt)+"\n")
    dump(receipt_path,receipt)
    print(json.dumps({k:summary[k] for k in ("status","original_N3_data_qualified","A_base_data_qualified","candidate_rows","truth_noise_outliers","max_cached_physical_identity_relative_error")}))


def write_reports(rows,summary,topologies):
    first=[r for r in topologies if r["case_id"]==CASES[0] and r["material"] in ("truth_OFFLINE","N3","A_base")]
    (OUT/"MASK_SEMANTICS.md").write_text("""# A30 图像与材料相语义审计

此目录为独立离线审计；真值只用于解释既有图像和评分，不能被在线搜索读取。

原 `plots.py:plot_material_images` 的 figure 05 逐像素显示 `abs(chi-background)`，色标为0–1.2。暗色表示与已知背景相近，不代表空气或孔洞。材料为线性 susceptibility `chi=epsilon_r-1`，`exp(-i omega t)`，Im(chi)≥0为被动损耗。背景本身是物质，不能因从背景中扣除而把它改叫真空。

`p250065_weak_a1_snr20` 真值包含550个背景单元 `chi=.1+.04i`、8个约 `1.01805+.22361i`、18个约 `1.276697+.275339i`，对应绝对介电常数的实部还需加1。这里严格真空定义仅为审计用 `|chi|<=1e-10`；真值、原N3、A基础图均无严格真空单元。这不是为新实验新增的接受门槛。

原历史支持阈值依赖真值增量的最大幅度，因此只能离线使用。采用原4连通前景、8连通背景，0.15/0.25/0.35三阈值下该真值均为7个前景连通分量、0个孔洞，原N3均为1个分量、0孔洞。family名shell不能代替此实际栅格材料的拓扑。栅格计数不是连续拓扑证书。

`PHASE_AUDIT_*.png/.pdf` 同时画出绝对epsilon实部、损耗、严格真空掩膜、原增量模及历史25%掩膜；统一色标。原始数组保存为 `OFFLINE_AUDIT_MATERIAL_ARRAYS.npz`，不得交给在线模块作位置提示。

原 `image_metrics.feasible` 只判断材料盒Re(chi)∈[-.5,4]、Im(chi)∈[0,.5]，没有观测残差条件。原表中的True不能转写为合法逆解。
""",encoding="utf-8")
    table="\n".join(f"|{r['case_id']}|{r['N3_observed_residual']:.6f}|" + (f"{r['A_base_observed_residual']:.6f}" if r['A_base_observed_residual'] is not None else "NOT_RUN") + f"|{r['truth_noise_only_norm']:.6f}|{r['A_base_data_status']}|" for r in rows)
    (OUT/"E0_DECISION.md").write_text(f"""# E0：物理合同可复原，优先重新取得合格基础材料

新增 full-state=0，新增 full-solver RHS=0；没有生成Green矩阵、分解或求解。采用保存的原始24×24 TM/VIE几何、8源、32接收器、complex128与原噪声。32个条件均为历史已暴露 development；其中8个有Research A复跑，其余24个明示A NOT_RUN，不记失败。

原N3在完整32条件有{summary['original_N3_data_qualified']}/32满足Research A冻结的观测相容阈值；固定A八条件中为{summary['fixed8_original_N3_data_qualified']}/8。Research A做了一步GN后的基础图为{summary['A_base_data_qualified']}/8。所有8个A基础图均不是原展示N3，不能把A的残差错贴给原N3。

冻结阈值τ={TAU:.15f}=sqrt(χ²_512(0.95))。σ是每complex通道RMS，实/虚部SD=σ/√2；逐源先Re后Im打包并乘√2/σ，共512个实通道；gain固定为1，各源共享同一材料。阈值来自Research A后加的重建相容规则，不是声称A30历史原本以此认证基础重建。

|case|原N3残差|A基础残差|真值噪声残差|A基础数据状态|
|---|---:|---:|---:|---|
{table}

原验证档的clean response与A28对应真值full prediction逐位一致；原S作用于保存current同样一致。保存的场方程、材料本构与预测检查最大相对误差{summary['max_cached_physical_identity_relative_error']:.3g}，没有已有模型不匹配证据。元数据噪声种子重新生成的观测误差为{summary['noise_metadata_replay_maxabs']:.3g}。

原N3 pool材料与已付费cache材料不是逐位相等：mass打包往返带来至多{summary['N3_cache_material_maxabs']:.3g}的ULP量级差异，各case另记录exact=False和numerical_match=True。这与Research A一步GN造成0.028–0.068量级材料变化是两件事，不能混淆。缓存的matching endpoint按既有数值容差解释，不声称byte identity。

首个weak条件真实噪声白化范数24.090925超过95%阈值，属于原噪声模型允许的小概率实现。不得因此宣称模型错误，不得为取得通过而改σ/阈值。即使真实材料恰好不满足这一噪声球，也不排除其他拟合噪声的合法材料存在。

E0没有阻断本轮已授权的E0-Q。应使用原数据/材料盒作更充分的全材料修正，明确保存KKT、线搜索、active bounds、每个匹配full endpoint和成本。不能把当前未合格基础图交给G冒充合格初始逆解。

原始N3缓存有完整fields/current并已零求解检查；A候选只保存material/full prediction，细节见missing manifest。`VALID_NUMERICAL`不含严格输出误差界。
""",encoding="utf-8")
    (OUT/"ROOT_CAUSE_MATRIX.md").write_text("""# 失败归因与不得混淆的结论

|分类|证据|解释与限制|
|---|---|---|
|INVALID_BASE|原N3 2/8、Research A一步GN基础3/8数据相容|5/8失败只能说明当前基础反演不足，不能证明替代世界不存在|
|TOPOLOGY_MISINTERPRETATION|figure05显示增量模；首shell实际7组件0孔洞且无真空|低增量、family名与空气孔洞不是同一概念|
|OPTIMIZATION_FAILURE|A用16/32维DCT/Haar chart、每轨迹3 full endpoints；原N3经过一步GN|容量有限且停止过早是合理待检假设；缺完整KKT，E0不能断言只是优化失败|
|REPRESENTATION_FAILURE|Sparse64改善35.6896%为离线真值目标优化；线上chart与此不同|说明存在更好表示的有限容量见证，不保证观测逆问题能在线识别或恢复这些系数|
|POSSIBLE_AMBIGUITY|部分候选数据相容但材料不同；需独立绝对结构核验|相容端点数不等于模式/连通分支数或概率质量，不得把单像素增量差视为真实成核|
|MODEL_MISMATCH|现有同数据clean/state/constitutive identities一致|目前没有缓存证据支持此归因；真值噪声残差超τ是噪声事件，不是模型不匹配证明|

旧A30补偿表中的`nominal_F/compensated_F`是两材料预测间距离，不是对观测y的重建残差。本目录另保存原表子集，不能与观测残差混并。历史Pairwise 8.84类数值亦必须依原归一化解释。

Sparse64的35.69%：64实Haar系数支撑在读development真值前由N3/TRAIN选择，但随后的系数优化显式最小化真值材料误差。它不是在线data-likelihood反演；没有保存相应full endpoint predictions，因此不做数据相容或多解认证。support选择与连续系数估计、材料误差与数据拟合都应分开。

本次读取不修改旧Gate、旧ledger或缺失记录。A原E2出现5/8无有效候选、额外候选仅1增量掩膜像素变化且无补偿已存在的结论保留。零新增物理审计无法判断全空间是否有合格新结构。
""",encoding="utf-8")


if __name__ == "__main__":
    main()

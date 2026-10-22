"""Compare completed surrogate-force artifacts; zero new physics and no GT.

Only explicit eight chains' scalar traces and saved reference material means
are decompressed. Original errors/flags remain untouched. Source diagnostic
statuses are quoted, never upgraded; final interpretation belongs to parent.
"""
from pathlib import Path
import json
import sys
import zipfile

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from p3l.core import OUT, Book, write

SEQUENCE = OUT/"schedules/surrogate_force_hmc_fix1"
CAP = 32*1024*1024
KEYS = ("Phi", "accepted", "divergences", "large_energy_error", "force_calls",
        "proxy_Phi_error", "proxy_white_response_RMS", "requested_leapfrog_steps", "leapfrog_steps", "warmup")


def read(path):
    if path.stat().st_size > 16*1024*1024:
        raise ValueError("unexpectedly large JSON summary")
    return json.loads(path.read_text())


def bounded(path, keys):
    with zipfile.ZipFile(path) as archive:
        sizes = [archive.getinfo(key+".npy").file_size for key in keys]
        if max(sizes, default=0) > CAP or sum(sizes) > CAP:
            raise MemoryError("requested scalar/reference NPZ members exceed 32MiB")
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in keys}


def statistics(value):
    values = np.asarray(value, float).ravel()
    finite = values[np.isfinite(values)]
    result = dict(count=len(values), finite_count=len(finite), nonfinite_count=int((~np.isfinite(values)).sum()),
                  nan_count=int(np.isnan(values).sum()), positive_infinity_count=int(np.isposinf(values).sum()),
                  negative_infinity_count=int(np.isneginf(values).sum()),
                  status="ALL_FINITE" if len(finite) == len(values) else "NONFINITE_RAW_VALUES_RETAINED_IN_SOURCE")
    for name, function in (("mean", np.mean), ("median", np.median), ("minimum", np.min), ("maximum", np.max)):
        result[name] = float(function(finite)) if len(finite) else None
    result["q05_q95"] = np.quantile(finite, [.05, .95]).tolist() if len(finite) else None
    result["statistic_scope"] = "finite entries only; all omitted nonfinite counts exposed, no raw values changed"
    return result


def reference(arm):
    folder = OUT/f"C/references/DEV_030/surrogate_force_{arm}"
    source = read(folder/"REFERENCE_STATS.json")
    result = {key: source.get(key) for key in ("status", "reference_qualified", "cross_mechanism_status",
        "target_family", "warmup_removed", "chain_count", "production_draws_per_chain",
        "full_grid_mean_rms_mcse", "gates", "thresholds", "summaries")}
    result.update(source=str((folder/"REFERENCE_STATS.json").relative_to(OUT)),
                  qualification_action="SOURCE_STATUS_QUOTED_ONLY_NO_UPGRADE",
                  numerical_events_source=str((folder/"NUMERICAL_EVENT_AUDIT.json").relative_to(OUT)))
    mean = bounded(folder/"REFERENCE_ARRAYS.npz", ["material_mean"])["material_mean"]
    if mean.shape != (576,) or not np.all(np.isfinite(mean)):
        raise ValueError("saved reference material mean shape/finite contract failed")
    return result, mean


def arm_results(arm):
    rows, traces = [], []
    audit = read(OUT/f"C/jobs/audit_surrogate_force_{arm}_DEV030/CLOSED.json")
    if audit.get("status") != "COMPLETE":
        raise RuntimeError("reference audit must complete")
    for chain in range(4):
        folder = OUT/f"C/jobs/surrogate_force_{arm}_DEV030_c{chain}"
        closed, transition = read(folder/"CLOSED.json"), read(folder/"TRANSITION_SUMMARY.json")
        if closed.get("status") != "COMPLETE" or closed.get("error") is not None:
            raise RuntimeError("all four arm chains must complete")
        arrays = bounded(folder/"CHAIN.npz", KEYS)
        warmup = int(arrays["warmup"])
        if warmup != 1024 or any(arrays[k].shape != (5120,) for k in KEYS if k != "warmup"):
            raise ValueError("fixed 1024 warmup plus 4096 production scalar layout required")
        prod = slice(warmup, None)
        rows.append(dict(chain=chain, source=str((folder/"CHAIN.npz").relative_to(OUT)),
            raw_transition_summary=transition, closure_status=closed["status"],
            actual_all_cost_RHS=transition["cost"]["RHS"],
            actual_all_controlled_wall_seconds=closed["controlled_wall_seconds"],
            all_force_calls=transition["total_proxy_force_calls"],
            all_proxy_force_seconds=transition["proxy_force_seconds"],
            step_size=transition["step_size"], production_acceptance=float(arrays["accepted"][prod].mean()),
            production_force_calls=int(arrays["force_calls"][prod].sum()),
            warmup_divergences=int(arrays["divergences"][:warmup].sum()),
            production_divergences=int(arrays["divergences"][prod].sum()),
            warmup_large_energy_errors=int(arrays["large_energy_error"][:warmup].sum()),
            production_large_energy_errors=int(arrays["large_energy_error"][prod].sum()),
            production_Phi=statistics(arrays["Phi"][prod]),
            production_signed_proxy_Phi_error=statistics(arrays["proxy_Phi_error"][prod]),
            production_absolute_proxy_Phi_error=statistics(np.abs(arrays["proxy_Phi_error"][prod])),
            production_white_proxy_response_RMS=statistics(arrays["proxy_white_response_RMS"][prod]),
            production_requested_leapfrog_steps=statistics(arrays["requested_leapfrog_steps"][prod]),
            production_completed_leapfrog_steps=statistics(arrays["leapfrog_steps"][prod])))
        traces.append({key: arrays[key][prod] for key in ("Phi", "proxy_Phi_error", "proxy_white_response_RMS")})
    ref, mean = reference(arm)
    aggregate = dict(chain_count=4, warmup_per_chain=1024, production_per_chain=4096,
        actual_all_cost_RHS=sum(r["actual_all_cost_RHS"] for r in rows),
        chain_controlled_wall_seconds=sum(r["actual_all_controlled_wall_seconds"] for r in rows),
        reference_audit_controlled_wall_seconds=audit["controlled_wall_seconds"],
        all_controlled_wall_seconds=sum(r["actual_all_controlled_wall_seconds"] for r in rows)+audit["controlled_wall_seconds"],
        all_force_calls=sum(r["all_force_calls"] for r in rows),
        all_proxy_force_seconds=sum(r["all_proxy_force_seconds"] for r in rows),
        production_acceptance=float(np.mean([r["production_acceptance"] for r in rows])),
        production_divergences=sum(r["production_divergences"] for r in rows),
        warmup_divergences=sum(r["warmup_divergences"] for r in rows),
        production_large_energy_errors=sum(r["production_large_energy_errors"] for r in rows),
        warmup_large_energy_errors=sum(r["warmup_large_energy_errors"] for r in rows))
    for name, key, absolute in (("production_Phi", "Phi", False),
            ("production_absolute_proxy_Phi_error", "proxy_Phi_error", True),
            ("production_white_proxy_response_RMS", "proxy_white_response_RMS", False)):
        values = np.concatenate([r[key] for r in traces])
        aggregate[name] = statistics(np.abs(values) if absolute else values)
    return dict(chains=rows, aggregate=aggregate, reference=ref), traces, mean


def figures(arms, traces):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder = OUT/"C/surrogate_comparison_figures"
    folder.mkdir(exist_ok=True)
    outputs = []
    def finish(fig, name):
        path = folder/(name+".png")
        fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)
        outputs.append(str(path.relative_to(OUT)))
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)
    for ax, arm in zip(axes, ("linear", "material")):
        for chain, values in enumerate(traces[arm]):
            ax.plot(np.arange(1,4097), values["Phi"], linewidth=.55, label=f"chain {chain}")
        ax.set(title=arm+" force: actual four-chain production Phi", xlabel="Production iteration (warmup removed)", ylabel="Full Phi")
        ax.legend()
    fig.suptitle("Same frozen surrogate pilot budgets; full-target traces do not establish mixing")
    finish(fig, "four_chain_production_Phi")
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
    labels = ["linear", "material"]
    for ax, key, title in zip(axes[0], ("actual_all_cost_RHS", "all_controlled_wall_seconds", "all_force_calls"),
            ("Actual all-cost RHS", "Chain + reference-audit controlled wall (s)", "All proxy force calls")):
        ax.bar(labels, [arms[a]["aggregate"][key] for a in labels]); ax.set_title(title)
    axes[1,0].bar(labels, [arms[a]["reference"]["full_grid_mean_rms_mcse"] for a in labels])
    axes[1,0].set_title("Saved material-mean RMS MCSE (not GT error)")
    for ax, key, title in ((axes[1,1], "proxy_Phi_error", "Absolute proxy Phi error at paid proposal endpoints"),
                           (axes[1,2], "proxy_white_response_RMS", "White-response proxy RMS at paid endpoints")):
        data = []
        for arm in labels:
            values = np.concatenate([v[key] for v in traces[arm]])
            values = np.abs(values) if key == "proxy_Phi_error" else values
            data.append(values[np.isfinite(values)])
        ax.boxplot(data, tick_labels=labels, showfliers=False)
        ax.set_yscale("symlog", linthresh=1e-6); ax.set_title(title)
    fig.suptitle("Matched warmup/production policy; all failures and paid setup included; no posterior-accuracy extrapolation")
    finish(fig, "matched_policy_cost_error_MCSE")
    return outputs


def markdown(report):
    lines = ["# Surrogate-force HMC: saved two-arm comparison", "",
        "This analysis reads completed artifacts only. It does not qualify a reference, issue a GO/NO-GO or novelty judgment, or estimate posterior accuracy against truth.", "",
        "The original parser failure remains preserved as ERROR with zero RHS. The valid fix1 sequence uses the same paid job IDs and the minimal production-policy fix recorded by the parent.", "",
        "| Arm | All RHS | Controlled wall, chain + audit s | Force calls | Production acceptance | Step-size range | Material mean RMS MCSE | Saved reference status |",
        "|---|---:|---:|---:|---:|---:|---:|---|"]
    for arm in ("linear", "material"):
        item = report["arms"][arm]; a = item["aggregate"]; r = item["reference"]
        steps = [row["step_size"] for row in item["chains"]]
        lines.append(f'| {arm} | {a["actual_all_cost_RHS"]} | {a["all_controlled_wall_seconds"]:.6g} | {a["all_force_calls"]} | {a["production_acceptance"]:.6g} | {min(steps):.6g}..{max(steps):.6g} | {r["full_grid_mean_rms_mcse"]:.6g} | {r["status"]} |')
    lines += ["", f'Between-arm saved mean-material RMS difference: {report["between_arm_saved_mean_material_RMS_difference"]:.9g}. This is a disagreement of two empirical references, not error against a true posterior or GT.', "",
        "Every chain has 1024 warmup + 4096 production iterations. Costs include each chain's actual charged setup/warmup/production and its reference-audit wall. Frozen production step sizes differ by arm/chain. Proxy discrepancies were measured after MH at already-paid proposal endpoints, including rejected proposals; they are not current-state Phi errors.", "",
        "JSON retains each original transition summary, numerical event counts, latent/property/material summaries and unassessed flags. Nonfinite proxy values remain in original NPZs and are counted explicitly; finite-only descriptive statistics never replace those raw failures.", "",
        "Ordinary Full HMC is external context with unequal warmup/production lengths and potentially different setup. Its saved diagnostics are linked separately; no equal-budget or speedup claim is made.", "",
        "## Figures", ""]
    lines.extend(f'![{Path(path).stem}]({path.removeprefix("C/")})' for path in report["figures"])
    return "\n".join(lines)+"\n"


def main():
    sequence = read(SEQUENCE/"CLOSED.json")
    if sequence.get("status") != "COMPLETE":
        raise RuntimeError("valid fix1 surrogate sequence must close COMPLETE before analysis")
    if (OUT/"C/SURROGATE_COMPARISON.json").exists():
        raise FileExistsError("comparison already exists; no implicit overwrite or replay")
    book = Book("C", "analyze_surrogate_force_fix1_comparison", max_wall=600, analysis_only=True)
    status, error = "ERROR", None
    try:
        arms, traces, means = {}, {}, {}
        for arm in ("linear", "material"):
            arms[arm], traces[arm], means[arm] = arm_results(arm)
        external_path = OUT/"A/references/DEV_030/hmc_pilot/REFERENCE_STATS.json"
        external = dict(status="NOT_RUN")
        if external_path.exists():
            saved = read(external_path)
            external = {k: saved.get(k) for k in ("status", "reference_qualified", "warmup_removed", "production_draws_per_chain", "chain_count", "summaries", "full_grid_mean_rms_mcse")}
            external.update(source=str(external_path.relative_to(OUT)), comparison_scope="EXTERNAL_UNEQUAL_WARMUP_LENGTH_SETUP_NOT_MATCHED_COST")
        old_path = OUT/"schedules/surrogate_force_hmc/CLOSED.json"
        report = dict(status="COMPLETE_ARTIFACT_COMPARISON_PARENT_REVIEW_REQUIRED", parent="DEV_030",
            valid_sequence_source=str((SEQUENCE/"CLOSED.json").relative_to(OUT)), valid_sequence_closure=sequence,
            original_failed_sequence=read(old_path) if old_path.exists() else dict(status="NOT_RUN"),
            arms=arms, external_Full_HMC_context=external,
            between_arm_saved_mean_material_RMS_difference=float(np.sqrt(np.mean(np.abs(means["linear"]-means["material"])**2))),
            new_physics_RHS=0, GT_read=False, Target_constructed=False,
            scientific_judgment="PARENT_ONLY", reference_qualification_changed=False, novelty_claim=False,
            semantics="Conventional frozen surrogate-force leapfrog with Full endpoint Hamiltonian MH; diagnostics only, no posterior-accuracy extrapolation")
        report["figures"] = figures(arms, traces)
        write(OUT/"C/SURROGATE_COMPARISON.json", report)
        (OUT/"C/SURROGATE_COMPARISON.md").write_text(markdown(report))
        status = "COMPLETE"
        print(json.dumps(dict(status=report["status"], new_physics_RHS=0,
                             arms={arm: value["aggregate"] for arm, value in arms.items()})))
    except BaseException as exc:
        error = repr(exc)
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()

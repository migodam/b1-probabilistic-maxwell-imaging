"""Offline DEV image risks of a declared profile-property GN Gaussian.

No physics/backend/Target/Book imports or solver calls. This script is intended
to be wrapped by the parent's stage-B resource owner. It must run in a dedicated
offline process after the paid worker exits. It does not qualify a posterior,
teacher, probability, credible interval, or held-out/NN result.
"""
from pathlib import Path
import argparse
import json
import sys
import zipfile

import numpy as np
from scipy.linalg import cho_solve, solve_triangular

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from p1_material.prior import MaterialPrior, saturation
from p3l.estimators import complex_geometric_median, material_diagnostics
from p3l.moments import gaussian_material_moments

DEFAULT_OUT = ROOT / "research/p3l_local_posterior_imaging"
LABEL = "PROFILE_PROPERTY_PLUGIN_GN_GAUSSIAN_REFERENCE"
ARRAY_CAP = 128 * 1024 * 1024


def write(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def save(path, **arrays):
    for value in arrays.values():
        if np.asarray(value).nbytes > ARRAY_CAP:
            raise MemoryError("individual array exceeds frozen 128 MiB limit")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
    temporary.replace(path)


def read_arrays(path, names):
    # Check zip member length before decompressing; do not load unrelated LU,
    # currents, truth latent states, predictions, or other cache entries.
    with zipfile.ZipFile(path) as archive:
        for name in names:
            if archive.getinfo(name + ".npy").file_size > ARRAY_CAP:
                raise MemoryError("requested array exceeds frozen 128 MiB limit")
    with np.load(path, allow_pickle=False) as archive:
        result = {name: archive[name] for name in names}
    if any(v.nbytes > ARRAY_CAP or not np.all(np.isfinite(v)) for v in result.values()):
        raise ValueError("nonfinite or oversized array")
    return result


def frozen_configuration():
    return dict(
        revision="P3L_B_GAUSSIAN_IMAGE_ESTIMATORS_001", reference=LABEL,
        parents=[dict(parent=f"DEV_{i:03d}", role="train" if i < 24 else "validation")
                 for i in range(32)],
        parent_selection="ALL_32_DEVELOPMENT_NO_SELECTION_NO_DROPPING_FAILURES",
        anchors="A/anchors/{parent}.npz; q, field, receiver_T required",
        covariance="Cz=(I64+Jz.T@Jz)^-1; stable Cholesky; fixed anchor properties",
        field_mean="TARGET.field_mean + MaterialPrior.A @ anchor_q[:64]",
        field_variance="diag(A @ Cz @ A.T)",
        whitening="source-major Re32 then Im32; multiply sqrt(2)/sigma",
        gaussian_center="paid optimized anchor; no new mode optimization or linear recentering",
        estimates=["gaussian_material_mean", "optimized_anchor", "original_prior_material_mean",
                   "decoder_at_latent_mean", "gaussian_material_geometric_median_mc"],
        geometric_median=dict(enabled=True, samples=2048, seed=2026102281,
                              parent_rng="SeedSequence([seed,parent_index])",
                              tolerance=1e-10, max_iterations=2000,
                              interpretation="GAUSSIAN_REFERENCE_MC_NOT_FULL_POSTERIOR"),
        risks=dict(masks=["full_grid", "ROI"], metrics=["complex_MSE", "complex_MAE", "complex_RMSE"],
                   comparisons=["optimized_anchor", "original_prior_material_mean"],
                   difference="estimator risk minus baseline risk; negative favors estimator",
                   aggregation=["all_development", "train", "validation"],
                   bootstrap=dict(replicates=2000, seed=2026102295, interval=[.025, .975],
                                  method="paired parent resampling within reported role; mean differences",
                                  missing="retain all parents; full-cohort result null if any missing; no success-only primary")),
        phases="freeze -> all material arrays and generation manifest -> offline DEV chi reads -> evaluation",
        truth_access="only inputs/offline_only/DEV_000..031.npz chi in evaluation phase",
        new_Maxwell_RHS=0, new_LU=0, new_forward_states=0,
        array_limit_bytes=ARRAY_CAP, heldout_access=False, NN_training=False,
        online_proposals=False, GT_selection=False, qualified_teacher=False,
        qualification="No full/local posterior, probabilities, credible certification or Gaussian adequacy claim")


def register(out):
    config = frozen_configuration()
    path = out / "configs/B_GAUSSIAN_IMAGE_ESTIMATORS.json"
    if path.exists():
        if json.loads(path.read_text()) != config:
            raise RuntimeError("frozen Gaussian estimator config differs; explicit revision required")
    else:
        # Exclusive creation prevents concurrent re-registration or overwrite.
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as stream:
            stream.write(json.dumps(config, indent=2, allow_nan=False) + "\n")
    return config


def stationarity(row):
    if row is None:
        return dict(status="ANCHOR_PANEL_ROW_MISSING", qualified_stationary_mode=False)
    kkt = row.get("KKT")
    return dict(status="PAID_ANCHOR_KKT_REPORTED_NOT_REQUALIFIED", KKT=kkt,
                noise_ball_qualified=row.get("noise_ball_qualified"),
                method=row.get("method"), objective=row.get("objective"),
                qualified_stationary_mode=False,
                interpretation="Anchor metadata only; no stationary Laplace or posterior certificate")


def construct(out, prior, cfg, row, config, prior_mean):
    parent = row["parent"]
    paid = read_arrays(out / f"A/anchors/{parent}.npz", ("q", "field", "receiver_T"))
    q, field, transfer = (paid[k] for k in ("q", "field", "receiver_T"))
    if q.shape != (66,) or not prior.legal(q):
        raise ValueError("paid anchor q must be legal native Gaussian64 x uniform2")
    if field.shape != (576, 8) or transfer.shape != (32, 576):
        raise ValueError("exact paid field(576,8) and receiver_T(32,576) required; no fallback")
    anchor, aux = prior.decode(q)
    Dg = prior.derivative(q)
    Cchi = np.einsum("qp,ps->qsp", transfer, field).reshape(256, 576)
    response = Cchi @ Dg[:, :64]
    # response row ordering is receiver-major/source-minor. Recover receivers,
    # sources, coordinates before source-major Re32/Im32 whitening.
    response = response.reshape(32, 8, 64).transpose(1, 0, 2)
    Jz = np.concatenate((response.real, response.imag), axis=1).reshape(512, 64)
    Jz *= np.sqrt(2.) / float(cfg["noise"]["sigma"])
    H = np.eye(64) + Jz.T @ Jz
    L = np.linalg.cholesky(H)
    covariance = cho_solve((L, True), np.eye(64))
    # A L^-T has covariance A H^-1 A.T; sum squares avoids diagonal cancellation.
    marginal_root = solve_triangular(L, prior.A.T, lower=True).T
    field_mean = prior.A @ q[:64] + cfg["field_mean"]
    field_sd = np.sqrt(np.sum(marginal_root * marginal_root, axis=1))
    moments = gaussian_material_moments(field_mean, field_sd, aux["inclusion_increment"],
                                        host=prior.host, mask=prior.mask,
                                        halfwidth=prior.temperature)
    arrays = dict(q_anchor=q, gaussian_field_mean=field_mean, gaussian_field_sd=field_sd,
                  gaussian_latent_covariance=covariance,
                  gaussian_material_mean=moments.mean, gaussian_material_variance=moments.variance,
                  gaussian_material_second_abs_moment=moments.second_abs_moment,
                  optimized_anchor=anchor, original_prior_material_mean=prior_mean,
                  decoder_at_latent_mean=prior.decode(q)[0])
    median_info = dict(status="NOT_RUN")
    mc = config["geometric_median"]
    if mc["enabled"]:
        rng = np.random.default_rng(np.random.SeedSequence([mc["seed"], int(parent[-3:])]))
        deviations = solve_triangular(L.T, rng.normal(size=(64, mc["samples"])), lower=False).T
        field_draws = (q[None, :64] + deviations) @ prior.A.T + cfg["field_mean"]
        materials = prior.host + aux["inclusion_increment"] * prior.mask[None, :] * saturation(field_draws / prior.temperature)
        median = complex_geometric_median(materials, tolerance=mc["tolerance"], max_iterations=mc["max_iterations"])
        diagnostics = material_diagnostics(materials)
        arrays.update(gaussian_material_geometric_median_mc=median.value,
                      median_converged=median.converged, median_iterations=median.iterations,
                      median_normalized_subgradient_residual=median.normalized_subgradient_residual)
        median_info = dict(status="GAUSSIAN_REFERENCE_MC_NUMERICAL_ACTION",
                           samples=mc["samples"], seed=[mc["seed"], int(parent[-3:])],
                           converged_pixels=int(np.count_nonzero(median.converged)),
                           maximum_subgradient_residual=float(np.max(median.normalized_subgradient_residual)),
                           material_spread_rms=diagnostics["material_spread_rms"],
                           interpretation="Gaussian-reference MC only; no full-posterior teacher or convergence claim")
    dest = out / f"B/gaussian_images/{parent}.npz"
    save(dest, **arrays)
    return dict(status="MATERIAL_ARRAYS_SAVED_BEFORE_TRUTH", arrays=str(dest.relative_to(out)),
                median=median_info, field_covariance_scope="full correlated latent covariance; pixel marginal moments",
                decoder_at_latent_mean_equals_anchor=bool(np.array_equal(arrays["decoder_at_latent_mean"], anchor)))


def loss(estimate, truth, mask):
    difference = np.abs(estimate[mask] - truth[mask])
    mse = float(np.mean(difference ** 2))
    return dict(complex_MSE=mse, complex_MAE=float(np.mean(difference)), complex_RMSE=float(np.sqrt(mse)))


def aggregate(rows, config):
    result = {}
    bootstrap = config["risks"]["bootstrap"]
    for role in config["risks"]["aggregation"]:
        cohort = [r for r in rows if role == "all_development" or r["role"] == role]
        rng = np.random.default_rng(bootstrap["seed"])
        # Same paired parent indices across every estimator/metric in a role.
        indices = rng.integers(0, len(cohort), size=(bootstrap["replicates"], len(cohort)))
        group = dict(parent_count=len(cohort), parents=[r["parent"] for r in cohort], comparisons={})
        for estimator in config["estimates"]:
            for baseline in config["risks"]["comparisons"]:
                if estimator == baseline:
                    continue
                for mask in config["risks"]["masks"]:
                    for metric in config["risks"]["metrics"]:
                        key = f"{estimator}_minus_{baseline}/{mask}/{metric}"
                        values = [r.get("differences", {}).get(estimator, {}).get(baseline, {}).get(mask, {}).get(metric) for r in cohort]
                        missing = [r["parent"] for r, value in zip(cohort, values) if value is None]
                        item = dict(parent_count=len(cohort), missing_parents=missing,
                                    difference_mean=None, paired_bootstrap_95_interval=None,
                                    status="INCOMPLETE_COHORT_NO_SUCCESS_ONLY_PRIMARY" if missing else "COMPLETE_DESCRIPTIVE_DEV_COHORT")
                        if not missing:
                            values = np.asarray(values)
                            replicates = values[indices].mean(axis=1)
                            item.update(difference_mean=float(values.mean()),
                                        paired_bootstrap_95_interval=np.quantile(replicates, bootstrap["interval"]).tolist())
                        group["comparisons"][key] = item
        result[role] = group
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    out = args.out.resolve()
    config = register(out)  # Frozen before any anchor generation or truth read.
    parents = json.loads((out / "configs/PARENTS.json").read_text())
    declared = [dict(parent=r["parent"], role=r["role"]) for r in parents if r["group"] == "development"]
    if declared != config["parents"]:
        raise RuntimeError("all32 development identities/order/roles differ from frozen config")
    cfg = json.loads((out / "configs/TARGET.json").read_text())
    geometry = json.loads((ROOT / "research/p1_material_target_20261018/inputs/geometry.json").read_text())
    prior = MaterialPrior(cfg, geometry["points"])
    if cfg["field_mean"] != 0.:
        raise ValueError("the exact original prior mean .5 mixture requires centered symmetric field")
    mean_increment = complex(np.mean(prior.real_interval), np.mean(prior.imag_interval))
    prior_mean = np.full(576, prior.host, complex) + .5 * mean_increment * prior.mask
    panel = json.loads((out / "A/ANCHOR_PANEL.json").read_text())
    panel_rows = {r["parent"]: r for r in panel["parents"]}
    rows = []
    # Generation phase: every estimator array is persisted for every successful
    # parent before the first offline truth path is opened. Failures stay rows.
    for row in config["parents"]:
        result = dict(**row, reference=LABEL, stationarity=stationarity(panel_rows.get(row["parent"])))
        try:
            result.update(construct(out, prior, cfg, row, config, prior_mean))
        except Exception as error:
            result.update(status="ESTIMATOR_CONSTRUCTION_FAILED_RETAINED", error=repr(error))
        rows.append(result)
    generation = dict(phase="ALL_ESTIMATORS_PERSISTED_OR_FAILURE_RECORDED_BEFORE_TRUTH", truth_read=False,
                      all32_retained=len(rows) == 32, rows=rows)
    write(out / "B/GAUSSIAN_IMAGE_GENERATION.json", generation)
    # Separate evaluation phase. Read only the DEV chi array, never its latent
    # truth or any train/validation/locked pool. No estimates are changed here.
    for result in rows:
        if result["status"] != "MATERIAL_ARRAYS_SAVED_BEFORE_TRUTH":
            continue
        parent = result["parent"]
        try:
            truth = read_arrays(out / f"inputs/offline_only/{parent}.npz", ("chi",))["chi"]
            if truth.shape != (576,):
                raise ValueError("DEV truth chi must have shape (576,)")
            estimates = read_arrays(out / result["arrays"], config["estimates"])
            risks = {name: {mask_name: loss(value, truth, mask) for mask_name, mask in
                     (("full_grid", np.ones(576, bool)), ("ROI", prior.mask))} for name, value in estimates.items()}
            differences = {name: {baseline: {mask: {metric: risks[name][mask][metric] - risks[baseline][mask][metric]
                            for metric in config["risks"]["metrics"]} for mask in config["risks"]["masks"]}
                            for baseline in config["risks"]["comparisons"] if baseline != name} for name in estimates}
            result.update(evaluation_status="OFFLINE_DEV_RISKS_ONLY", risks=risks, differences=differences)
        except Exception as error:
            result.update(evaluation_status="OFFLINE_EVALUATION_FAILED_RETAINED", evaluation_error=repr(error))
    report = dict(reference=LABEL, estimator_config="configs/B_GAUSSIAN_IMAGE_ESTIMATORS.json",
                  all32_retained=True, no_GT_selection=True, new_Maxwell_RHS=0,
                  teacher_qualified=False, full_sampler_qualified=False, rows=rows,
                  aggregation=aggregate(rows, config),
                  interpretation="Profile property plug-in GN Gaussian image reference only; descriptive DEV risks, no full/local posterior, credible, probabilities, held-out or NN claim")
    write(out / "B/GAUSSIAN_IMAGE_RISK.json", report)
    print(json.dumps(dict(status="OFFLINE_GAUSSIAN_IMAGE_EVALUATION_RECORDED", parents=len(rows),
                          construction_failures=sum("error" in r for r in rows),
                          evaluation_failures=sum(r.get("evaluation_status") != "OFFLINE_DEV_RISKS_ONLY" for r in rows),
                          new_Maxwell_RHS=0)))


if __name__ == "__main__":
    main()

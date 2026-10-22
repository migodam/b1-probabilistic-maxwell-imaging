"""Paid transformed-coordinate GN Gaussian construction; no truth access.

The P1 target is unchanged. The probit-coordinate mode is a coordinate MAP,
not the native Gaussian64 x Uniform2 MAP. No Gaussian, GH refinement or IS ESS
alone qualifies the full posterior or a teacher. Offline scoring is separate.
"""
from pathlib import Path
import argparse
import gc
import json
import sys

import numpy as np
from scipy.linalg import cho_solve, solve_triangular
from scipy.optimize import least_squares
from scipy.special import logsumexp, ndtr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from p3l.core import (OUT, Book, Target, BudgetStop, write, save, plain,
                      to_probit, from_probit, probit_derivative)
from p3l.moments import gaussian_material_moments


class StateCap(RuntimeError):
    pass


def freeze(importance):
    config = dict(
        revision="P3L_B_JOINT_PROBIT_GAUSSIAN_001",
        reason="Fixed-property Gaussian mean worsened DEV MSE by about 49.7%; test joint property compensation and exact boundary transport without truth selection",
        parents=[f"DEV_{i:03d}" for i in range(32)],
        roles={f"DEV_{i:03d}": "train" if i < 24 else "validation" for i in range(32)},
        reference="FULL_CORRELATED_JOINT_PROBIT_GN_GAUSSIAN_REFERENCE_NOT_FULL_POSTERIOR",
        target_changed=False, prior_changed=False, noise_changed=False,
        target="P1 likelihood and prior transported exactly to z64 x eta2; u=normalCDF(eta)",
        energy="Phi(Tx)+.5*||x||^2; normalized log density=-energy-33log(2pi)",
        initialization="to_probit(paid native A anchor, initialization=True); inward property move initialization only",
        optimizer=dict(method="TRF", residual="[white(F-y), x66]",
                       jacobian="[Jq*probit_derivative(x);I66]", fresh_state_cap=64,
                       ftol=1e-9, xtol=1e-9, gtol=1e-7, x_scale="jac",
                       retain="minimum transformed energy evaluated Point; failed/capped parents retained"),
        mode_interpretation="Probit coordinate MAP differs from native q MAP; neither is a material Bayes mean",
        covariance="H^-1, H=I66+Jx.T Jx; GN reference, not exact Hessian",
        conditioning="B=-Hzz^-1 Hzeta; Cz|eta=Hzz^-1; Seta=Hetaeta-Hetaz Hzz^-1 Hzeta; Ceta=Seta^-1",
        integration=dict(method="standard-normal tensor Gauss-Hermite over eta2", orders=[5, 9, 17],
                         inner="gaussian_material_moments with conditional field mean/SD and fixed u=normalCDF(eta)",
                         outputs="Echi and E|chi|²; variance=E|chi|²-|Echi|²; successive-order refinement diagnostics only"),
        controls=["native_paid_anchor", "decoder_at_probit_mode", "fixed_property_mode_gaussian_mean",
                  "joint_property_correlated_gaussian_mean"],
        importance=dict(enabled=bool(importance), parents=["DEV_000", "DEV_028", "DEV_030", "DEV_026"],
                        draws_per_parent=64, seed=2026102282,
                        draw="x~N(probit_mode,H^-1)", physics="64 full forward states, 8 RHS each; no adjoints/Jacobians",
                        weights="exp(-Phi-.5*x²-33log2pi)/normalized_N(mode,H^-1)",
                        interpretation="full-target IS diagnostic; ESS alone never qualification"),
        serial=True, transfer_cache_limit=4, maximum_individual_array_bytes=128*1024*1024,
        book=dict(stage="B", job="joint_probit_gaussian", max_wall=7200),
        GT_read=False, network=False, qualified_teacher=False,
        scoring="separate offline process after all controls saved; all32 failure parents retained")
    path = OUT / "configs/B_JOINT_GAUSSIAN_REVISION.json"
    if path.exists():
        if json.loads(path.read_text()) != config:
            raise RuntimeError("joint Gaussian config already frozen differently; explicit new revision required")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x") as stream:
            stream.write(json.dumps(config, indent=2, allow_nan=False) + "\n")
    return config


def optimize(target, y, parent, folder, native):
    """Retain only current and best full Points; cap every fresh attempt."""
    settings = dict(ftol=1e-9, xtol=1e-9, gtol=1e-7)
    x0 = to_probit(native.q, initialization=True)
    best = last = None
    best_x = last_x = None
    best_energy = np.inf
    fresh = 0
    history = []
    trace_path = folder / "MODE_HISTORY.jsonl"
    if trace_path.exists():
        raise FileExistsError("mode history already exists; no implicit paid replay")
    folder.mkdir(parents=True, exist_ok=True)
    trace = trace_path.open("x")

    def log(event, **fields):
        row = plain(dict(event=event, parent=parent, fresh_state_attempts=fresh, **fields))
        history.append(row)
        trace.write(json.dumps(row, allow_nan=False) + "\n")
        trace.flush()

    def get(x):
        nonlocal best, last, best_x, last_x, best_energy, fresh
        x = np.array(x, float, copy=True)
        if x.shape != (66,) or not np.all(np.isfinite(x)):
            raise ValueError("finite probit coordinates required")
        if last_x is not None and np.array_equal(x, last_x):
            return last
        if best_x is not None and np.array_equal(x, best_x):
            last, last_x = best, best_x.copy()
            return best
        q = from_probit(x)
        if np.array_equal(q, native.q):
            point = native
        else:
            if fresh >= 64:
                raise StateCap("64 fresh mode state attempts exhausted")
            fresh += 1
            log("FRESH_FORWARD_ATTEMPT", x=x, q=q)
            point = target.evaluate(q, parent + "_joint_probit_mode")
        phi = target.phi(point, y)
        energy = phi + .5 * float(x @ x)
        if not np.isfinite(energy):
            raise FloatingPointError("nonfinite transformed energy")
        native_objective = phi + .5 * float(q[:64] @ q[:64])
        log("EVALUATED", endpoint=point.key, x=x, q=q, Phi=phi,
            transformed_energy=energy, native_objective=native_objective)
        last, last_x = point, x.copy()
        if energy < best_energy:
            best, best_x, best_energy = point, x.copy(), energy
            write(folder / "BEST_MODE.json", dict(x=best_x, q=best.q, endpoint=best.key,
                  transformed_energy=best_energy, native_objective=native_objective,
                  coordinate_mode="PROBIT_COORDINATE_MAP_NOT_NATIVE_MAP"))
        return point

    def fun(x):
        p = get(x)
        return np.r_[target.white(p.state.data - y), x]

    def jac(x):
        p = get(x)
        Jx = target.jacobian(p) * probit_derivative(x)[None, :]
        return np.vstack((Jx, np.eye(66)))

    optimizer = dict(status="NOT_STARTED", success=False)
    try:
        result = least_squares(fun, x0, jac=jac, method="trf", tr_solver="exact",
                               x_scale="jac", max_nfev=64, **settings)
        optimizer = dict(status="SCIPY_FINISHED", success=bool(result.success),
                         scipy_status=int(result.status), message=str(result.message),
                         nfev=int(result.nfev), njev=int(result.njev),
                         scipy_optimality=float(result.optimality), terminal_cost=float(result.cost))
    except BudgetStop:
        raise
    except StateCap as error:
        optimizer = dict(status="STRICT_FRESH_STATE_CAP", success=False, error=repr(error))
    except Exception as error:
        optimizer = dict(status="OPTIMIZER_FAILURE_RETAINED", success=False, error=repr(error))
    finally:
        trace.close()
        write(folder / "MODE_OPTIMIZATION.json", dict(optimizer=optimizer,
              fresh_state_attempts=fresh, history=history,
              initialization_x=x0, initialization_q=from_probit(x0),
              boundary_initialization_only=True))
    if best is None:
        raise RuntimeError("no finite transformed mode Point available; preserve history")
    # This may pay a missing receiver Jacobian, but never an extra forward.
    Jq = target.jacobian(best)
    Jx = Jq * probit_derivative(best_x)[None, :]
    residual = np.r_[target.white(best.state.data - y), best_x]
    Jaug = np.vstack((Jx, np.eye(66)))
    gradient = Jaug.T @ residual
    H = Jaug.T @ Jaug
    scale = np.sqrt(np.diag(H))
    stationarity = dict(gradient=gradient, gradient_inf=float(np.max(np.abs(gradient))),
                        scaled_gradient_inf=float(np.max(np.abs(gradient) / scale)),
                        optimizer=optimizer, fresh_state_attempts=fresh,
                        coordinate_mode="PROBIT_COORDINATE_MAP_NOT_NATIVE_MAP",
                        qualification="numerical stationarity diagnostic only; no exact Hessian/posterior qualification")
    target.persist_cache(best, folder / "RAW_MODE_CACHE.npz")
    save(folder / "JOINT_MODE_GEOMETRY.npz", x_mode=best_x, q_mode=best.q,
         Jq=Jq, Jx=Jx, H_GN=H, transformed_gradient=gradient)
    write(folder / "STATIONARITY.json", stationarity)
    return best, best_x, H, stationarity


def integrate(prior, mode, x, H, folder):
    Hzz, Hzeta = H[:64, :64], H[:64, 64:]
    Lz = np.linalg.cholesky(Hzz)
    solve_cross = cho_solve((Lz, True), Hzeta)
    B = -solve_cross
    Schur = H[64:, 64:] - Hzeta.T @ solve_cross
    Schur = .5 * (Schur + Schur.T)
    Leta = np.linalg.cholesky(Schur)
    Ceta = cho_solve((Leta, True), np.eye(2))
    Cz = cho_solve((Lz, True), np.eye(64))
    eta_root = solve_triangular(Leta.T, np.eye(2), lower=False)
    field_root = solve_triangular(Lz, prior.A.T, lower=True).T
    field_sd = np.sqrt(np.sum(field_root ** 2, axis=1))
    _, aux = prior.decode(mode.q)
    fixed = gaussian_material_moments(prior.config["field_mean"] + prior.A @ x[:64], field_sd,
              aux["inclusion_increment"], host=prior.host, mask=prior.mask, halfwidth=prior.temperature)
    save(folder / "CONTROLS.npz", decoder_at_probit_mode=mode.state.chi,
         fixed_property_mode_gaussian_mean=fixed.mean,
         fixed_property_mode_gaussian_variance=fixed.variance,
         conditional_field_sd=field_sd, conditional_latent_covariance=Cz,
         eta_covariance=Ceta, eta_precision_Schur=Schur, conditional_regression=B)
    previous = None
    refinements = []
    for order in (5, 9, 17):
        nodes, weights = np.polynomial.hermite.hermgauss(order)
        nodes = np.sqrt(2.) * nodes
        weights = weights / np.sqrt(np.pi)
        mean = np.zeros(576, complex)
        second = np.zeros(576)
        for i, a in enumerate(nodes):
            for j, b in enumerate(nodes):
                delta_eta = eta_root @ np.array([a, b])
                u = ndtr(x[64:] + delta_eta)
                increment = complex(prior.real_interval[0] + np.diff(prior.real_interval)[0] * u[0],
                                    prior.imag_interval[0] + np.diff(prior.imag_interval)[0] * u[1])
                field_mean = prior.config["field_mean"] + prior.A @ (x[:64] + B @ delta_eta)
                moments = gaussian_material_moments(field_mean, field_sd, increment,
                          host=prior.host, mask=prior.mask, halfwidth=prior.temperature)
                weight = weights[i] * weights[j]
                mean += weight * moments.mean
                second += weight * moments.second_abs_moment
        raw_variance = second - np.abs(mean) ** 2
        save(folder / f"JOINT_GH_ORDER_{order}.npz", chi_mean=mean,
             chi_variance=np.maximum(0., raw_variance), chi_second_abs_moment=second,
             raw_complex_modulus_variance=raw_variance)
        entry = dict(order=order, minimum_raw_variance=float(raw_variance.min()),
                     status="DETERMINISTIC_GAUSSIAN_REFERENCE_QUADRATURE_NOT_POSTERIOR_CERTIFICATE")
        if previous is not None:
            old_order, old_mean, old_second = previous
            entry.update(previous_order=old_order,
                         mean_refinement_complex_RMS=float(np.sqrt(np.mean(np.abs(mean - old_mean) ** 2))),
                         mean_refinement_complex_max=float(np.max(np.abs(mean - old_mean))),
                         second_abs_moment_refinement_max=float(np.max(np.abs(second - old_second))))
        refinements.append(entry)
        previous = (order, mean.copy(), second.copy())
        write(folder / "GH_REFINEMENT.json", dict(orders=refinements,
              interpretation="successive quadrature differences only, no Gaussian adequacy/full-posterior proof"))
    return refinements


def importance(target, y, parent, mode, x, H, folder):
    rng = np.random.default_rng(np.random.SeedSequence([2026102282, int(parent[-3:])]))
    L = np.linalg.cholesky(H)
    draws = x + solve_triangular(L.T, rng.normal(size=(66, 64)), lower=False).T
    qs, predictions, materials, logweights, phis = [], [], [], [], []
    logdet_half = float(np.log(np.diag(L)).sum())
    # Every completed draw is preserved atomically even if a later attempt fails.
    for index, draw in enumerate(draws):
        q = from_probit(draw)
        p = target.evaluate(q, parent + "_joint_gaussian_full_target_IS")
        phi = target.phi(p, y)
        delta = draw - x
        logproposal = logdet_half - 33*np.log(2*np.pi) - .5*float(delta @ H @ delta)
        logtarget = -phi - .5*float(draw @ draw) - 33*np.log(2*np.pi)
        logweight = logtarget - logproposal
        save(folder / "IS_DRAWS" / f"draw_{index:03d}.npz", x=draw, q=q,
             prediction=p.state.data, chi=p.state.chi, Phi=np.array(phi),
             log_proposal=np.array(logproposal), log_target=np.array(logtarget),
             logweight=np.array(logweight))
        qs.append(q.copy()); predictions.append(p.state.data.copy())
        materials.append(p.state.chi.copy()); logweights.append(logweight); phis.append(phi)
        del p
    logweights = np.asarray(logweights)
    weights = np.exp(logweights - logsumexp(logweights))
    materials = np.asarray(materials)
    mean = weights @ materials
    save(folder / "FULL_TARGET_IS.npz", x=draws, q=np.asarray(qs), weights=weights,
         logweights=logweights, Phi=np.asarray(phis), prediction=np.asarray(predictions),
         chi=materials, weighted_material_mean=mean,
         weighted_material_variance=weights @ np.abs(materials - mean) ** 2)
    report = dict(draws=64, ESS=float(1 / (weights @ weights)),
                  max_normalized_weight=float(weights.max()),
                  log_mean_importance_weight=float(logsumexp(logweights) - np.log(64)),
                  forward_RHS=512, receiver_or_adjoint_RHS=0,
                  status="FULL_TARGET_IMPORTANCE_DIAGNOSTIC_NOT_QUALIFIED_POSTERIOR_OR_TEACHER")
    write(folder / "FULL_TARGET_IS.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--importance", action="store_true", help="also pay 64 full-target IS draws on each fixed four-parent pilot")
    args = parser.parse_args()
    config = freeze(args.importance)  # Before any fresh physical query.
    rows = json.loads((OUT / "configs/PARENTS.json").read_text())
    dev = [r for r in rows if r["group"] == "development"]
    if [r["parent"] for r in dev] != config["parents"] or any(r["role"] != config["roles"][r["parent"]] for r in dev):
        raise RuntimeError("frozen all32 DEV identities/roles differ")
    if args.importance:
        deep = json.loads((OUT / "A/DEEP_REFERENCE_PARENTS.json").read_text())["parents"]
        if deep[:4] != config["importance"]["parents"]:
            raise RuntimeError("y-only first-four panel differs from frozen IS panel")
    book = Book("B", "joint_probit_gaussian", max_wall=7200)
    status, error = "COMPLETE", None
    results = [dict(parent=r["parent"], role=r["role"], status="NOT_RUN") for r in dev]
    summary_path = OUT / "B/JOINT_GAUSSIAN_REFERENCE.json"
    try:
        target = Target(book, keep_fields=True)
        for result in results:
            parent = result["parent"]
            folder = OUT / f"B/joint_gaussian/{parent}"
            native = mode = None
            try:
                y = target.observation(parent)
                native = target.load_cache(OUT / f"A/anchors/{parent}.npz", owner="paid y-only A native anchor; probit Gaussian revision")
                save(folder / "NATIVE_ANCHOR_CONTROL.npz", q=native.q, chi=native.state.chi,
                     native_objective=np.array(target.objective(native, y)))
                mode, x, H, stationarity = optimize(target, y, parent, folder, native)
                result.update(status="PROBIT_MODE_AND_GN_GEOMETRY_SAVED", stationarity=plain(stationarity))
                result["quadrature"] = integrate(target.prior, mode, x, H, folder)
                result["status"] = "JOINT_GAUSSIAN_REFERENCE_SAVED_NOT_POSTERIOR_QUALIFIED"
                if args.importance and parent in config["importance"]["parents"]:
                    result["importance"] = importance(target, y, parent, mode, x, H, folder)
            except BudgetStop:
                result.update(status="BUDGET_STOP_PARTIAL_ARTIFACTS_RETAINED")
                raise
            except Exception as exc:
                result.update(status="PARENT_FAILURE_RETAINED", error=repr(exc))
                write(folder / "FAILURE.json", result)
            finally:
                native = mode = None
                if target.model is not None:
                    target.model._transfers.clear()
                gc.collect()
                write(summary_path, dict(config="configs/B_JOINT_GAUSSIAN_REVISION.json", rows=results,
                      all32_retained=True, GT_read=False, target_changed=False,
                      status="DECLARED_GAUSSIAN_REFERENCES_NOT_FULL_POSTERIOR_OR_TEACHER"))
    except BaseException as exc:
        status, error = "ERROR", repr(exc)
        raise
    finally:
        book.close(status, error)


if __name__ == "__main__":
    main()

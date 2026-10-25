#!/usr/bin/env python3
"""Reproducible linear-algebra checks for Research W-SLAM.

This script runs no PDE and does not establish Maxwell identifiability.
Run: python deliverables/code/observability_toy.py
"""
from pathlib import Path
import json
import numpy as np


def blockdiag(blocks):
    out = np.zeros((sum(b.shape[0] for b in blocks), sum(b.shape[1] for b in blocks)))
    i = j = 0
    for b in blocks:
        out[i:i + b.shape[0], j:j + b.shape[1]] = b
        i += b.shape[0]
        j += b.shape[1]
    return out


def sym(a):
    return (a + a.T) / 2


def invsqrt(a):
    vals, vecs = np.linalg.eigh(sym(a))
    assert vals.min() > 0
    return (vecs * vals ** -0.5) @ vecs.T


def schur_information(a, b, lam=None):
    if lam is None:
        lam = np.zeros((b.shape[1], b.shape[1]))
    h = sym(b.T @ b + lam)
    return sym(a.T @ a - a.T @ b @ np.linalg.pinv(h, rcond=1e-12) @ b.T @ a)


def dynamic_information_full_rank(a, b, d, q):
    c = np.linalg.pinv(b, rcond=1e-12) @ a
    r = a - b @ c
    g = b.T @ b
    return sym(r.T @ r + c.T @ d.T @ np.linalg.solve(q + d @ np.linalg.solve(g, d.T), d @ c))


def dynamic_information_singular(a, b, d, q):
    _, sv, vt = np.linalg.svd(b, full_matrices=True)
    rank = int((sv > max(b.shape) * np.finfo(float).eps * sv.max()).sum())
    v = vt[:rank].T
    n = vt[rank:].T
    c = np.linalg.pinv(b, rcond=1e-12) @ a
    r = a - b @ c
    qi = invsqrt(q)
    dn = qi @ d @ n
    p = sym(np.eye(d.shape[0]) - dn @ np.linalg.pinv(dn, rcond=1e-12))
    dt = p @ qi @ d @ v
    ct = p @ qi @ d @ c
    gr = v.T @ b.T @ b @ v
    return sym(r.T @ r + ct.T @ np.linalg.solve(np.eye(d.shape[0]) + dt @ np.linalg.solve(gr, dt.T), ct)), rank


def main():
    rng = np.random.default_rng(20261025)
    results = {"scope": "finite-dimensional linear-Gaussian algebra only; no Maxwell/PDE run", "seed": 20261025}
    # Three times, two target coordinates per time, five observations per time.
    b = blockdiag([rng.normal(size=(5, 2)) for _ in range(3)])
    a = rng.normal(size=(15, 4))
    d = np.zeros((4, 6))
    d[:2, :2] = -np.eye(2)
    d[:2, 2:4] = np.eye(2)
    d[2:, 2:4] = -np.eye(2)
    d[2:, 4:] = np.eye(2)
    z = rng.normal(size=(4, 4))
    q = z @ z.T + 0.3 * np.eye(4)
    direct = schur_information(a, b, d.T @ np.linalg.solve(q, d))
    derived = dynamic_information_full_rank(a, b, d, q)
    err = float(np.linalg.norm(direct - derived, 2))
    results["dynamic_full_rank"] = {"operator_norm_error": err, "minimum_schur_eigenvalue": float(np.linalg.eigvalsh(direct).min())}
    assert err < 1e-10

    # Singular per-time target sensing; first eliminate target-null motion freedom.
    b = blockdiag([rng.normal(size=(3, 2)) @ rng.normal(size=(2, 4)) for _ in range(3)])
    a = rng.normal(size=(9, 4))
    d = rng.normal(size=(8, 12))
    z = rng.normal(size=(8, 8))
    q = z @ z.T + np.eye(8)
    direct = schur_information(a, b, d.T @ np.linalg.solve(q, d))
    derived, rank = dynamic_information_singular(a, b, d, q)
    err = float(np.linalg.norm(direct - derived, 2))
    results["dynamic_singular"] = {"target_sensing_rank": rank, "target_dimension": 12, "operator_norm_error": err}
    assert err < 1e-9

    # Explicit false-recovery counterexample if G^{-1} is replaced by G^dagger.
    b = np.array([[1., 0., 0., 0.], [0., 0., 1., 0.]])
    a = np.array([[1.], [-1.]])
    d = np.array([[-1., -1., 1., 0.], [0., -1., 0., 1.]])
    q = np.eye(2)
    c = np.linalg.pinv(b) @ a
    naive = c.T @ d.T @ np.linalg.solve(q + d @ np.linalg.pinv(b.T @ b) @ d.T, d @ c)
    exact = schur_information(a, b, d.T @ d)
    corrected, _ = dynamic_information_singular(a, b, d, q)
    compensation = np.array([[-1.], [2.], [1.], [2.]])
    results["singular_pseudoinverse_counterexample"] = {
        "wrong_information_from_naive_G_pinv": float(naive.item()),
        "joint_schur_information": float(exact.item()),
        "corrected_formula_information": float(corrected.item()),
        "compensating_target_state_for_wall_delta_1": compensation[:, 0].tolist(),
        "data_residual_norm": float(np.linalg.norm(a + b @ compensation)),
        "motion_residual_norm": float(np.linalg.norm(d @ compensation))}
    assert abs(exact.item()) < 1e-12 and abs(corrected.item()) < 1e-12
    assert abs(naive.item() - 4 / 3) < 1e-12

    # Two scalar observations with a stochastic motion model x2=f*x1+eta.
    c1, c2, f, g1, g2, process_variance = 1., -1., 1., 4., 9., 0.3
    b = np.diag(np.sqrt([g1, g2]))
    a = b @ np.array([[c1], [c2]])
    d = np.array([[-f, 1.]])
    exact = schur_information(a, b, d.T @ d / process_variance).item()
    closed = (c2 - f * c1) ** 2 / (process_variance + f ** 2 / g1 + 1 / g2)
    results["two_time_scalar"] = {"joint_schur": float(exact), "closed_form": float(closed), "absolute_error": float(abs(exact - closed))}
    assert abs(exact - closed) < 1e-10

    # A new measurement after conditional nuisance integration.
    nw, nu = 5, 3
    z = rng.normal(size=(nw + nu, nw + nu))
    h = z.T @ z + np.eye(nw + nu)
    j = rng.normal(size=(4, nw + nu))
    a, b = j[:, :nw], j[:, nw:]
    huu, huw = h[nw:, nw:], h[nw:, :nw]
    old = h[:nw, :nw] - h[:nw, nw:] @ np.linalg.solve(huu, huw)
    hn = h + j.T @ j
    new = hn[:nw, :nw] - hn[:nw, nw:] @ np.linalg.solve(hn[nw:, nw:], hn[nw:, :nw])
    c = a - b @ np.linalg.solve(huu, huw)
    noise = np.eye(4) + b @ np.linalg.solve(huu, b.T)
    update = c.T @ np.linalg.solve(noise, c)
    err = float(np.linalg.norm(new - old - update, 2))
    results["conditional_innovation"] = {"operator_norm_error": err, "minimum_information_increment_eigenvalue": float(np.linalg.eigvalsh(sym(new - old)).min())}
    assert err < 1e-10

    # Exact bridge from map marginal precision to a linear task covariance.
    t = rng.normal(size=(2, nw + nu))
    tw, tu = t[:, :nw], t[:, nw:]
    effective_task = tw - tu @ np.linalg.solve(huu, huw)
    task_full = t @ np.linalg.solve(h, t.T)
    task_split = tu @ np.linalg.solve(huu, tu.T) + effective_task @ np.linalg.solve(old, effective_task.T)
    err = float(np.linalg.norm(task_full - task_split, 2))
    results["task_covariance_bridge"] = {"operator_norm_error": err}
    assert err < 1e-10

    # Shared nuisance can recover a direction invisible to either experiment.
    single1 = schur_information(np.array([[1.]]), np.array([[1.]])).item()
    single2 = schur_information(np.array([[1.]]), np.array([[-1.]])).item()
    joint = schur_information(np.array([[1.], [1.]]), np.array([[1.], [-1.]])).item()
    results["shared_compensation_recovery"] = {"single1": single1, "single2": single2, "joint": joint, "status": "old/standard result, not a new contribution"}
    assert abs(single1) < 1e-12 and abs(single2) < 1e-12 and abs(joint - 2) < 1e-12

    # Absolute information is nondecreasing, normalized retention need not be.
    kold = schur_information(np.array([[1.]]), np.array([[0.]])).item()
    knew = schur_information(np.array([[1.], [1.]]), np.array([[0.], [1.]])).item()
    results["retention_nonmonotonicity"] = {"old_absolute_information": kold, "new_absolute_information": knew, "old_retention": kold, "new_retention": knew / 2, "status": "old/standard result, not a new contribution"}
    assert abs(kold - 1) < 1e-12 and abs(knew - 1) < 1e-12
    results["all_checks_passed"] = True
    output = Path(__file__).resolve().parents[1] / "results" / "observability_toy.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

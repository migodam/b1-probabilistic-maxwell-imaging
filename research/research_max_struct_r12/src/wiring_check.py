"""Small scalar integration and shared-vector-material wiring checks only.

No study, optimization, derivative check or truth-scoring experiment is run.
"""
from pathlib import Path
import importlib.util
import json
import tempfile

import numpy as np

from field_current import (DenseFieldCurrentProjectors, PhysicalAudit,
                           UnreachableDataBall, run_projection)


def main():
    source = Path(__file__).parents[1] / "sources/pro/MAX_STRUCT_RESEARCH/feasibility/consensus_reflection_toy.py"
    spec = importlib.util.spec_from_file_location("pro_toy_unmodified", source)
    toy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(toy)
    _, G, E, S = toy.setup(n=3, ns=2, nr=5)
    # Public arbitrary data; no truth argument is exposed to the core.
    Y = np.full((5, 2), .03 + .02j)
    physics = toy.Physics(G, E, S, Y)

    def audit(m):
        before = physics.lu_count, physics.rhs_count
        data, J, Q = physics.forward(m)
        prior = .001 * float(m @ m)
        return PhysicalAudit(data + prior, J, Q,
                             metrics={"data_objective": data, "prior": prior},
                             costs={"physical_lu": physics.lu_count - before[0],
                                    "physical_rhs": physics.rhs_count - before[1]})

    report = {"scope": "small wiring only", "scalar": {}}
    for mode in ("ap", "dr"):
        p = DenseFieldCurrentProjectors(G, E, S, Y, delta=.01)
        snapshots = []
        run = run_projection(mode, p, np.full(9, .1), audit,
                             steps=6, audit_every=2,
                             checkpoint_callback=lambda row: snapshots.append(row.iteration))
        assert snapshots == [0, 2, 4, 6]
        assert run.objective == min(r.objective for r in run.history)
        assert run.physical_costs == {"physical_lu": 4., "physical_rhs": 8.}
        assert run.costs["H_iteration_rhs"] == 12
        assert p.costs["H_cholesky"] == 1 and p.costs["H_setup_rhs"] == 5
        for row in run.history:
            assert row.objective == row.metrics["data_objective"] + row.metrics["prior"]
            assert row.states["X"].shape == E.shape
        with tempfile.TemporaryDirectory() as directory:
            run.save(directory)
            arrays = np.load(Path(directory) / "states.npz", allow_pickle=False)
            assert np.array_equal(arrays["material"], run.material)
            saved = json.loads((Path(directory) / "history.json").read_text())
            assert len(saved["history"]) == 4
        report["scalar"][mode] = {"iterations": run.iterations,
                                   "audits": len(run.history),
                                   "physical_costs": run.physical_costs,
                                   "H_setup_rhs": p.costs["H_setup_rhs"]}
    # Components deliberately request distinct slopes; result must remain ONE
    # scalar per cell, selected using all components and both illuminations.
    ev = np.ones((6, 2), complex)
    vector = DenseFieldCurrentProjectors(np.zeros((6, 6)), ev,
                                        np.eye(6), np.zeros((6, 2)), 100.,
                                        components=3, scale=1., lower=-2., upper=3.)
    X = np.array([[.2, .4], [1., .3], [2., .8],
                  [-.6, -.1], [-1., -1.3], [-.2, -.9]], complex)
    Q = ev.copy()
    cx, cq, m = vector.material(X, Q)
    expanded = np.repeat(m, 3)
    assert m.shape == (2,)
    assert np.linalg.norm(cx - expanded[:, None] * cq) < 1e-13
    # Independent tiny eigensystem readout for these two mixed-component cells.
    errors = []
    for cell in range(2):
        xi, qi = X[cell * 3:(cell + 1) * 3], Q[cell * 3:(cell + 1) * 3]
        aa, bb = float(np.vdot(xi, xi).real), float(np.vdot(qi, qi).real)
        cc = float(np.vdot(xi, qi).real)
        _, vectors = np.linalg.eigh([[aa, -cc], [-cc, bb]])
        candidates = [-2., 3.]
        if abs(vectors[0, 0]) > 1e-20:
            root = vectors[1, 0] / vectors[0, 0]
            if -2 <= root <= 3:
                candidates.append(root)
        reference = min(np.linalg.norm((xi - a * qi) / np.hypot(1, a)) ** 2
                        for a in candidates)
        actual = (np.linalg.norm(cx[cell * 3:(cell + 1) * 3] - xi) ** 2
                  + np.linalg.norm(cq[cell * 3:(cell + 1) * 3] - qi) ** 2)
        errors.append(abs(actual - reference))
    assert max(errors) < 1e-12
    report["vector_shared_material"] = {"cells": 2, "components": 3, "sources": 2,
                                         "material": m.tolist(),
                                         "max_objective_difference": max(errors)}
    # Minimal receiver-contract edge wiring: different masks and whitening,
    # one global ball, plus unreachable and pseudoinverse-limit branches.
    p = DenseFieldCurrentProjectors(np.zeros((2, 2)), np.ones((2, 2)),
                                    np.eye(2), np.ones((2, 2)), .3,
                                    receiver_masks=[np.array([True, False]), np.array([0, 1])],
                                    whiteners=[np.array([2.]), np.diag([1., .5])])
    ax, _ = p.affine_ball(np.zeros((2, 2)), np.ones((2, 2)))
    residual = np.sqrt(sum(np.linalg.norm(port.T @ ax[:, j] - y) ** 2
                           for j, (port, y) in enumerate(zip(p.ports, p.targets))))
    assert abs(residual - .3) < 1e-12
    zeros = np.zeros((1, 1))
    try:
        DenseFieldCurrentProjectors(zeros, np.ones((1, 1)), zeros,
                                    np.ones((1, 1)), .5)
    except UnreachableDataBall:
        pass
    else:
        raise AssertionError("unreachable data ball accepted")
    p = DenseFieldCurrentProjectors(zeros, np.ones((1, 1)),
                                    np.array([[1.], [0.]]), np.array([[2.], [1.]]),
                                    1., scale=1.)
    ax, _ = p.affine_ball(zeros, np.ones((1, 1)))
    assert p.last_A_info["pseudoinverse_limit"] and np.allclose(ax, 2.)
    report["receiver_contract"] = {"masked_whitened_global_residual": float(residual),
                                   "unreachable_raised": True,
                                   "pseudoinverse_limit": True}
    # A zero remaining budget retains the initial physical audit/incumbent.
    budget_run = run_projection("dr", p, np.zeros(1),
                                lambda m: PhysicalAudit(1., zeros.copy(), np.ones((1, 1))),
                                steps=10, cpu_budget_seconds=0.)
    assert budget_run.iterations == 0 and budget_run.stop_reason == "cpu_budget"
    report["zero_budget"] = {"iterations": budget_run.iterations,
                              "stop_reason": budget_run.stop_reason,
                              "initial_incumbent_retained": True}
    (Path(__file__).parent / "wiring_check_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

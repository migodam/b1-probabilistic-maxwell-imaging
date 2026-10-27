"""Zero-PDE driver helper checks.  Never call g0_operators.main()."""
from pathlib import Path
import json

import numpy as np
from scipy.sparse import csc_matrix

import g0_operators as driver


def main():
    out = Path(__file__).resolve().parent
    rng = np.random.default_rng(271027)
    raw = rng.normal(size=(24, 4)) + 1j * rng.normal(size=(24, 4))
    snapshots = np.column_stack([raw[:, 0], 2 * raw[:, 0], raw[:, 1],
                                 raw[:, 2], raw[:, 3], np.zeros(24)])
    q, qr = driver.two_pass_joint_qr(snapshots, [f"synthetic_{j}" for j in range(6)])
    assert q.shape == (24, 4), q.shape
    assert np.linalg.norm(q.conj().T @ q - np.eye(4)) < 1e-12
    assert [row["accepted"] for row in qr["columns"]] == [True, False, True, True, True, False]

    class SyntheticBook:
        def __init__(self):
            self.rows = []

        def dense(self, *args, **kwargs):
            self.rows.append(("dense", args, kwargs))

        def event(self, *args, **kwargs):
            self.rows.append(("event", args, kwargs))

    class SyntheticCache:
        pass

    class SyntheticAssembly:
        pass

    cache = SyntheticCache()
    cache.U = np.arange(7)
    cache.C = np.eye(7, dtype=complex) * .2
    cache.uw = rng.normal(size=(7, 2)) + 1j * rng.normal(size=(7, 2))
    assembly = SyntheticAssembly()
    assembly.Vtarget = csc_matrix(np.diag(np.linspace(.05, .3, 7) + .01j))
    assembly.bX = rng.normal(size=(7, 2)) + 1j * rng.normal(size=(7, 2))
    book = SyntheticBook()
    diagnostic, arrays = driver.inverse_core_diagnostic(cache, assembly, book)
    expected = 1 / np.min(np.abs(1 - .2 * assembly.Vtarget.diagonal()))
    assert abs(diagnostic["inverse_one_norm_estimate"] - expected) < 1e-12
    assert diagnostic["actual_equation_backward"]["backward"] < 1e-14
    assert len([row for row in book.rows if row[0] == "dense"]) == 1
    assert np.linalg.norm(arrays["core_probe_residual"]) < 1e-12

    _, gates, cases, _ = driver.frozen_plan()
    memory = driver.owned_memory_plan(28419, 2257, 1_500_000_000)
    try:
        driver.owned_memory_plan(60000, 2257, 1_500_000_000)
    except MemoryError:
        larger_rejected = True
    else:
        raise AssertionError("larger memory plan must be refused")

    # Exact failures must retain raw output.  A ROM with a correct reduced
    # solve must preserve its inaccurate full-equation result as accuracy.
    case = {"id": "synthetic_gate"}
    case_dir = out / case["id"]
    case_dir.mkdir(exist_ok=True)
    matrix = csc_matrix(np.eye(3, dtype=complex))
    rhs = np.array([[1.], [2.], [3.]], complex)
    reference = rhs.copy(), rhs[:1].T.copy()
    bad = rhs + .01, reference[1] + .01, {"backward": 0.}
    row = driver.save_method(out, case, "woodbury", bad, reference,
                             matrix, rhs, {}, gates)
    assert row["strict_algebra_status"] == "FAIL"
    assert (case_dir / "woodbury.npz").exists()
    basis = np.array([[1.], [0.], [0.]], complex)
    rom = np.array([[1.], [0.], [0.]], complex), reference[1] + .1, {"backward": .2}
    rom_row = driver.save_method(out, case, "goal_rom", rom, reference,
                                 matrix, rhs, {}, gates, basis)
    assert rom_row["strict_failures"] == []
    assert rom_row["full_equation_backward"]["backward"] > 1e-10
    assert rom_row["reduced_equation_backward"]["backward"] == 0.

    report = dict(
        status="PASS_ZERO_PDE_HELPER_CHECKS",
        new_FEM_mesh_assembly_sparse_LU_or_physical_rhs=0,
        driver_main_executed=False,
        physical_results_directory_created=driver.OUT.exists(),
        joint_QR_retained_rank=qr["retained_columns"],
        joint_QR_gram_error=qr["gram_identity_frobenius_error"],
        synthetic_core_inverse_norm_relative_error=abs(diagnostic["inverse_one_norm_estimate"] / expected - 1),
        synthetic_core_extra_dense_LU_count=1,
        strict_failure_arrays_preserved=True,
        ROM_full_residual_is_accuracy_only=True,
        frozen_plan_first_seven_cases_validated=len(cases),
        named_owned_array_estimate_for_representative_free_count=memory,
        larger_memory_plan_rejected=larger_rejected,
        scope="synthetic algebra and static manifest checks only; no Maxwell qualification",
        earlier_attempt=dict(interpreter="system python3", status="DEPENDENCY_UNAVAILABLE",
                             error="ModuleNotFoundError: No module named numpy", physical_work=0))
    assert report["physical_results_directory_created"] is False
    (out / "ZERO_PDE_CHECKS.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

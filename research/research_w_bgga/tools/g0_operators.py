"""Seven-condition, same-discrete G0 operator comparison (explicit physical work).

Run once through ``tools/run_job.py G0 <unique-job-name> python tools/g0_operators.py``.
The wrapper is the sole owner of aggregate process-CPU accounting.  This script
records measured kernel CPU/wall time without charging aggregate CPU twice.

The full, resolvent, paid-cache Schur, and preconditioned Krylov comparisons are
ordinary exact algebra.  The fixed goal ROM is an accuracy experiment.  Neither
this driver nor successful same-space comparisons complete the G0 stress panel.
"""
from __future__ import annotations

from pathlib import Path
import gc
import json
import os
import time
import warnings

import numpy as np
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from scipy.sparse import save_npz
from scipy.sparse.linalg import LinearOperator, onenormest

from bgga.accounting import BudgetBook, _default
from bgga.fixed_mesh import TargetAssembler, public_fixed_mesh, public_mesh_key
from bgga.operators import FixedSpaceCache, goal_galerkin, relative_backward
from g0_qualify import background_system
from wslam.fem_backend import FEM, line_green
from wslam.world import World


ROOT = Path(__file__).resolve().parents[1]
PHASE = "G0"
OUT = ROOT / "results" / "g0_operators"
CASE_IDS = ("weak_center", "strong_center", "left_small", "right_small",
            "left_large", "right_large", "low_loss_high_contrast")
TRAINING_IDS = ("weak_center", "left_small")
HELD_OUT_IDS = tuple(case for case in CASE_IDS if case not in TRAINING_IDS)
EXACT_METHODS = ("full", "woodbury", "schur_from_paid_ports", "krylov")
EXPECTED_PORTS = 2257
DEFLATION_TOL = 1e-12
INVERSE_ESTIMATOR = dict(t=2, itmax=5)
INVERSE_ESTIMATOR_SEED = 271027
THREAD_VARIABLES = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS")


class StrictQualificationFailure(RuntimeError):
    """An artifact has been saved, and the strict algebraic run must stop."""


def write_json(path, value):
    """Update only this run's status files; never replace another run directory."""
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, default=_default, indent=2,
                                    allow_nan=False) + "\n")
    temporary.replace(path)


def measured(call):
    cpu, wall = time.process_time(), time.perf_counter()
    # SciPy reports a singular dense LU through a warning; it is a hard failure.
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", LinAlgWarning)
            value = call()
    except Exception as error:
        error.bgga_kernel_timing = dict(measured_kernel_cpu_seconds=time.process_time() - cpu,
                                        measured_kernel_wall_seconds=time.perf_counter() - wall)
        raise
    return value, dict(measured_kernel_cpu_seconds=time.process_time() - cpu,
                       measured_kernel_wall_seconds=time.perf_counter() - wall)


def finite_arrays(*arrays):
    return all(np.isfinite(array).all() for array in arrays)


def relative_difference(value, reference):
    absolute = float(np.linalg.norm(value - reference))
    denominator = float(np.linalg.norm(reference))
    relative = absolute / denominator if denominator else (0. if not absolute else None)
    return dict(absolute=absolute, reference_norm=denominator, relative=relative)


def world_from_case(case):
    return World(np.asarray(case["support"], bool), case["epsinf"], case["ell"],
                 case["pose"], radius=case["radius"], contrast=case["contrast"],
                 target_loss=case["target_loss"])


def frozen_plan():
    contract = json.loads((ROOT / "CONTRACT.json").read_text())
    gates = json.loads((ROOT / "GATES.json").read_text())[PHASE]
    cases = json.loads((ROOT / "G0_CASES.json").read_text())["qualification"][:7]
    if tuple(case["id"] for case in cases) != CASE_IDS:
        raise ValueError("The first seven predeclared conditions changed; no adaptive selection")
    for case in cases:
        if (np.asarray(case["support"]).shape != (2, 2)
                or not np.array_equal(case["support"], np.ones((2, 2), int))
                or case["epsinf"] != 3.75 or case["ell"] != .125):
            raise ValueError("All seven conditions must share the frozen full-2x2 background")
    if (gates["relative_same_discrete_identity"] != 1e-9
            or gates["linear_backward_residual"] != 1e-10):
        raise ValueError("G0 strict thresholds changed; this driver does not select new thresholds")
    fs, physics = contract["fixed_space"], contract["physics"]
    if (fs["h"] != .17 or fs["wall_subdivision"] != 2
            or fs["target_step"] != .125 or fs["circle_polygon_n"] != 384
            or physics["threads"] != 1 or physics["relative_frequencies"] != [1.0]
            or physics["noise_std_each_real_imag"] != .01):
        raise ValueError("This driver implements only the frozen public P3 discretization")
    qualification_path = ROOT / "results" / "g0_qualification" / "RESULTS.json"
    qualified = {row["case"]: row for row in json.loads(qualification_path.read_text())}
    for case_id in CASE_IDS:
        row = qualified.get(case_id)
        if (row is None or row.get("numerical_qualification") is not True
                or row.get("target_ports") != EXPECTED_PORTS
                or row["algebraic_backward"]["backward"] > gates["linear_backward_residual"]):
            raise ValueError(f"Prior fixed-space prerequisite is absent/failed: {case_id}")
    return contract, gates, cases, dict(
        schema="BGGA_G0_ORDINARY_OPERATOR_COMPARISON_V1",
        case_ids=CASE_IDS, training_ids=TRAINING_IDS, held_out_ids=HELD_OUT_IDS,
        background=dict(support=[[1, 1], [1, 1]], epsinf=3.75, ell=.125),
        background_sparse_lu=1, background_primal_rhs=2,
        public_ports=EXPECTED_PORTS, public_port_setup_calls=1,
        public_port_rhs=EXPECTED_PORTS, paid_background_adjoint_rhs=16,
        full_baseline_target_lu=7, full_baseline_target_rhs=14,
        known_non_krylov_physical_rhs=2 + EXPECTED_PORTS + 16 + 14,
        krylov_rhs="Every actual background-preconditioner application is separately reserved",
        schur="schur_from_paid_ports; no independent exterior LU or extra 2259 RHS",
        goal_basis=dict(snapshot_order=["background primal (2)", "weak_center full (2)",
                                        "left_small full (2)", "background receiver adjoint (16)"],
                        candidate_columns=22, orthogonalization="two-pass modified Gram-Schmidt QR",
                        deflation_relative_tolerance=DEFLATION_TOL,
                        rank_rule="algebraic deflation only; no rank sweep or held-out-field selection"),
        core_stability=dict(case="low_loss_high_contrast", core="I - C V",
                            dimension=EXPECTED_PORTS, additional_dense_lu=1,
                            inverse_one_norm_estimator=INVERSE_ESTIMATOR,
                            estimator_seed=INVERSE_ESTIMATOR_SEED,
                            interpretation="numerical estimate only; not an inf-sup lower bound"),
        prior_qualification_source=str(qualification_path.relative_to(ROOT)),
        ordinary_exact_methods=EXACT_METHODS, rom_method="goal_rom",
        strict_identity_threshold=gates["relative_same_discrete_identity"],
        strict_backward_threshold=gates["linear_backward_residual"],
        rom_full_equation_residual="accuracy result; no exactness threshold or automatic fallback",
        held_out_policy="Retain held-out accuracy even when inaccurate; never modify Q/rank",
        perturbation_jacobian_and_information_stress_status="NOT_RUN_BY_THIS_DRIVER",
        whole_G0_status="PARTIAL; independent G0 requirements remain pending",
        novelty="ordinary operator baselines; no BGGA novelty inference",
        aggregate_cpu_authority="tools/run_job.py completion receipt, once")


def two_pass_joint_qr(snapshots, labels, tolerance=DEFLATION_TOL):
    """Fixed-order two-pass QR with only predeclared algebraic deflation.

    Each input is normalized before QR, so unit differences between primal and
    adjoint fields do not select a scientific rank.  The fixed 1e-12 test only
    removes a column already represented to that relative numerical tolerance.
    """
    if snapshots.ndim != 2 or snapshots.shape[1] != len(labels):
        raise ValueError("one frozen provenance label per joint snapshot is required")
    if not finite_arrays(snapshots):
        raise ValueError("nonfinite training/background/adjoint snapshot")
    q = np.empty(snapshots.shape, dtype=np.complex128)
    accepted, rows = 0, []
    for column, label in enumerate(labels):
        source = snapshots[:, column]
        original_norm = float(np.linalg.norm(source))
        candidate = source / original_norm if original_norm else np.zeros_like(source)
        projections = []
        for _ in range(2):
            coeff = q[:, :accepted].conj().T @ candidate
            projections.append(float(np.linalg.norm(coeff)))
            candidate = candidate - q[:, :accepted] @ coeff
        norm = float(np.linalg.norm(candidate))
        keep = bool(original_norm and norm > tolerance)
        rows.append(dict(input_column=column, source=label, original_norm=original_norm,
                         normalized_remainder_norm=norm, accepted=keep,
                         projection_norm_by_pass=projections))
        if keep:
            q[:, accepted] = candidate / norm
            accepted += 1
    q = q[:, :accepted].copy()
    if not accepted:
        raise StrictQualificationFailure("All fixed goal snapshots were algebraically deflated")
    orthogonality = float(np.linalg.norm(q.conj().T @ q - np.eye(accepted)))
    # This is a numerical construction check, not a selected model rank.
    if not np.isfinite(orthogonality) or orthogonality > 1e-10:
        raise StrictQualificationFailure("Two-pass goal QR orthogonality check failed")
    return q, dict(method="fixed-order two-pass modified Gram-Schmidt QR",
                   normalized_inputs=True, tolerance=tolerance, input_columns=len(labels),
                   retained_columns=accepted, algebraic_deflation_only=True,
                   gram_identity_frobenius_error=orthogonality, columns=rows,
                   held_out_fields_used=False, rank_sweep=False)


def factor_storage_bytes(factor):
    return int(sum(array.nbytes for matrix in (factor.L, factor.U)
                   for array in (matrix.data, matrix.indices, matrix.indptr))
               + factor.perm_r.nbytes + factor.perm_c.nbytes)


def owned_memory_plan(n, m, max_bytes):
    # Keep build_ports' original guard and cover Schur's additional live dense
    # copies.  This bounds named owned arrays, not process RSS or opaque BLAS/LU
    # workspace; the wrapper separately reports peak child RSS.
    port_lower = 16 * (n * m + 3 * m * m)
    comparison_estimate = 16 * (n * m + 5 * m * m + n * 60 + 16 * m)
    plan = dict(port_cache_original_estimate_bytes=int(port_lower),
                named_comparison_array_estimate_bytes=int(comparison_estimate),
                max_bytes=int(max_bytes),
                model="complex128 n*m + five m*m dense buffers + 60 n-vectors + receiver ports",
                includes="shared Green images, restricted C, Schur/core copies, snapshots and fields",
                excludes="FEM/Basis storage, sparse LU storage and opaque library workspace",
                process_RSS_guarantee=False)
    if max(port_lower, comparison_estimate) > max_bytes:
        raise MemoryError(f"Owned-array guard refuses allocation: {plan}")
    return plan


def inverse_core_diagnostic(cache, assembly, book):
    """One predeclared 2257-core estimate, with a single additional dense LU.

    Dense solves are algebra on the already-paid port cache.  They have no new
    PDE RHS charge; every solve is still logged with its dense dimensions/RHS.
    """
    m = len(cache.U)
    v = assembly.Vtarget[cache.U][:, cache.U].toarray()
    core = np.eye(m, dtype=complex) - cache.C @ v
    del v
    book.dense(PHASE, "one predeclared low-loss/high-contrast core diagnostic LU", m)
    factor = lu_factor(core)
    counts = dict(solve_calls=0, right_hand_sides=0, adjoint_solve_calls=0,
                  adjoint_right_hand_sides=0)

    def apply(rhs, adjoint=False):
        columns = 1 if rhs.ndim == 1 else rhs.shape[1]
        counts["adjoint_solve_calls" if adjoint else "solve_calls"] += 1
        counts["adjoint_right_hand_sides" if adjoint else "right_hand_sides"] += columns
        book.event(PHASE, "dense inverse-norm-estimator triangular solve",
                   dense_dimension=m, dense_rhs=int(columns), transpose="H" if adjoint else "N",
                   new_physical_rhs=0)
        return lu_solve(factor, rhs, trans=2 if adjoint else 0)

    inverse = LinearOperator((m, m), dtype=np.complex128,
                             matvec=apply, matmat=apply,
                             rmatvec=lambda rhs: apply(rhs, True),
                             rmatmat=lambda rhs: apply(rhs, True))
    # onenormest may resample sign vectors.  Fix only this diagnostic's random
    # stream and restore the caller's stream; no case/basis is selected with it.
    random_state = np.random.get_state()
    try:
        np.random.seed(INVERSE_ESTIMATOR_SEED)
        inverse_norm = float(onenormest(inverse, **INVERSE_ESTIMATOR))
    finally:
        np.random.set_state(random_state)
    # Probe uses the actual resolvent RHS of this predeclared case.
    rhs = cache.C @ (assembly.bX + assembly.Vtarget @ cache.uw)[cache.U]
    book.event(PHASE, "dense core diagnostic actual-equation solve",
               dense_dimension=m, dense_rhs=int(rhs.shape[1]), new_physical_rhs=0)
    solution = lu_solve(factor, rhs)
    defect = relative_backward(core, solution, rhs)
    if not finite_arrays(core, solution) or not np.isfinite(inverse_norm):
        raise StrictQualificationFailure("Nonfinite predeclared core stability diagnostic")
    core_norm = float(np.max(np.sum(np.abs(core), axis=0)))
    result = dict(case="low_loss_high_contrast", description="low-loss/high-contrast attempt",
                  core="I-CV", dimension=m, dense_factorizations=1,
                  core_one_norm=core_norm, inverse_one_norm_estimate=inverse_norm,
                  condition_one_norm_estimate_product=core_norm * inverse_norm,
                  inverse_estimator=INVERSE_ESTIMATOR, estimator_seed=INVERSE_ESTIMATOR_SEED,
                  estimator_work=counts,
                  actual_equation_backward=defect,
                  interpretation="finite numerical inverse-norm estimate; not an inf-sup lower bound",
                  singular_values="NOT_COMPUTED", resonance_claim="NONE",
                  jitter=False, pseudoinverse=False)
    return result, dict(core_probe_rhs=rhs, core_probe_solution=solution,
                        core_probe_residual=core @ solution - rhs)


def save_method(out, case, method, result, reference, matrix, rhs, timing, gates, q=None):
    """Persist fields and diagnostics before deciding a strict stop."""
    u, data, returned = result
    full_u, full_data = reference
    verify_cpu, verify_wall = time.process_time(), time.perf_counter()
    residual = matrix @ u - rhs
    finite = finite_arrays(u, data, residual, full_u, full_data)
    # Nonfinite values cannot be JSON numbers.  Preserve their raw arrays first.
    defect = relative_backward(matrix, u, rhs) if finite else None
    field_error = relative_difference(u, full_u) if finite else None
    data_error = relative_difference(data, full_data) if finite else None
    exact = method in EXACT_METHODS
    failures = []
    if not finite:
        failures.append("NONFINITE_OUTPUT")
    elif exact:
        for name, error in (("field", field_error), ("response", data_error)):
            if error["relative"] is None or error["relative"] > gates["relative_same_discrete_identity"]:
                failures.append(f"{name.upper()}_SAME_DISCRETE_IDENTITY")
        if defect["backward"] > gates["linear_backward_residual"]:
            failures.append("FULL_EQUATION_BACKWARD")
        if not np.isfinite(returned["backward"]) or returned["backward"] > gates["linear_backward_residual"]:
            failures.append("METHOD_REPORTED_BACKWARD")
        if "core_backward" in returned and (
                not np.isfinite(returned["core_backward"]["backward"])
                or returned["core_backward"]["backward"] > gates["linear_backward_residual"]):
            failures.append("METHOD_CORE_BACKWARD")
    projected = None
    if not exact and finite:
        core, reduced_rhs = q.conj().T @ (matrix @ q), q.conj().T @ rhs
        coefficients = q.conj().T @ u
        projected = relative_backward(core, coefficients, reduced_rhs)
        if projected["backward"] > gates["linear_backward_residual"]:
            failures.append("ROM_REDUCED_SOLVE_BACKWARD")
    verification = dict(verification_cpu_seconds=time.process_time() - verify_cpu,
                        verification_wall_seconds=time.perf_counter() - verify_wall)
    destination = out / case["id"] / f"{method}.npz"
    io_cpu, io_wall = time.process_time(), time.perf_counter()
    np.savez(destination, response=data, scattered_field=u, equation_residual=residual,
             field_difference_to_full=u - full_u, response_difference_to_full=data - full_data)
    row = dict(case=case["id"], method=method, response_shape=list(data.shape),
               field_shape=list(u.shape), ordinary_exact_same_discrete=exact,
               strict_algebra_status="FAIL" if failures else ("PASS" if exact else "REDUCED_SOLVE_PASS"),
               strict_failures=failures, full_equation_backward=defect,
               method_reported=returned if finite else "NONFINITE; raw arrays saved",
               field_difference_to_full=field_error, response_difference_to_full=data_error,
               full_response_difference_whitened_l2=float(np.linalg.norm(data - full_data) / .01) if finite else None,
               reduced_equation_backward=projected,
               rom_partition=("TRAINING" if case["id"] in TRAINING_IDS else "HELD_OUT") if not exact else None,
               rom_accuracy_status="MEASURED; no accuracy gate, basis change or fallback" if not exact else None,
               solve_timing=timing, **verification,
               artifact_cpu_seconds=time.process_time() - io_cpu,
               artifact_wall_seconds=time.perf_counter() - io_wall,
               arrays=str(destination.relative_to(out)))
    write_json(destination.with_suffix(".json"), row)
    return row


def main():
    # No PDE work or output directory is started without the fixed grant and
    # single-thread wrapper.  Re-running never overwrites a previous result.
    if os.environ.get("BGGA_STAGE") != PHASE or any(os.environ.get(key) != "1" for key in THREAD_VARIABLES):
        raise RuntimeError("Run through tools/run_job.py G0 so one-thread execution and wrapper CPU accounting apply")
    contract, gates, cases, plan = frozen_plan()
    book = BudgetBook()
    authorization = book.authorization(PHASE)
    OUT.mkdir(parents=True, exist_ok=False)
    rows, setups = [], []
    context = dict(step="initialization", case=None, method=None)
    status = dict(operator_comparison_status="RUNNING", whole_G0_status=plan["whole_G0_status"],
                  stress_panel_status="NOT_RUN_BY_THIS_DRIVER", exact_methods=EXACT_METHODS,
                  goal_rom_held_out_ids=HELD_OUT_IDS, aggregate_cpu_authority=plan["aggregate_cpu_authority"])
    write_json(OUT / "PLAN.json", plan)
    write_json(OUT / "STATUS.json", status)

    def persist():
        write_json(OUT / "RESULTS.json", rows)
        write_json(OUT / "SETUP.json", setups)
        write_json(OUT / "STATUS.json", status)

    def setup(label, call):
        context.update(step=label, case=None, method=None)
        value, timing = measured(call)
        setups.append(dict(step=label, **timing))
        persist()
        return value

    def record(case, method, result, reference, matrix, rhs, timing, q=None):
        row = save_method(OUT, case, method, result, reference, matrix, rhs, timing, gates, q)
        rows.append(row)
        persist()
        print(json.dumps(row, default=_default, allow_nan=False), flush=True)
        if row["strict_failures"]:
            raise StrictQualificationFailure(f"{case['id']}/{method}: {row['strict_failures']}")
        return row

    try:
        tx = np.array([[-1.5, 0.], [1.5, 0.]])
        rx = np.column_stack((np.linspace(-2, 2, 16), np.zeros(16)))
        p = setup("one public fixed P3 FEM assembly", lambda: FEM(
            1., tx, rx, book, h=.17, order=3, mesh=public_fixed_mesh(), phase=PHASE,
            residual_tolerance=1e-10, max_dofs=contract["physics"]["max_dofs"],
            max_factor_bytes=contract["physics"]["max_factor_bytes"]))
        assembler = setup("public target assembler", lambda: TargetAssembler(p.basis, p._free_dofs))
        ports = assembler.public_roi_free_indices
        if len(ports) != EXPECTED_PORTS:
            raise StrictQualificationFailure(f"Public port layout changed: {len(ports)} != {EXPECTED_PORTS}")
        max_bytes = contract["physics"]["max_factor_bytes"]
        memory = owned_memory_plan(len(p._free_dofs), len(ports), max_bytes)
        write_json(OUT / "MEMORY_PLAN.json", memory)
        R = p._receiver_probe[:, p._free_dofs]
        direct = line_green(rx, tx, p.k).T
        Lw, bw = setup("shared full-2x2 eps=3.75 ell=.125 background assembly",
                       lambda: background_system(p, world_from_case(cases[0])))
        cache = setup("one background sparse LU and two primal RHS",
                      lambda: FixedSpaceCache(Lw, bw, R, direct, ports, book, PHASE))
        actual_factor_bytes = factor_storage_bytes(cache.factor)
        write_json(OUT / "BACKGROUND.json", dict(
            full_dofs=p.ndofs, free_dofs=len(p._free_dofs), public_ports=len(ports),
            sparse_factor_storage_bytes=actual_factor_bytes, max_factor_bytes=max_bytes,
            background_backward=relative_backward(Lw, cache.uw, bw),
            public_mesh_key=public_mesh_key(), paid_primal_rhs=2))
        np.savez(OUT / "background.npz", rhs=bw, scattered_field=cache.uw,
                 response=cache.wall_data, direct=direct, free_dofs=p._free_dofs,
                 public_ports=ports, equation_residual=Lw @ cache.uw - bw)
        save_npz(OUT / "background_operator.npz", Lw)
        save_npz(OUT / "receiver_operator.npz", R.tocsr())
        if actual_factor_bytes > max_bytes:
            raise MemoryError("Actual background SuperLU storage exceeds frozen 1.5GB guard")
        if (not finite_arrays(cache.uw, cache.wall_data)
                or relative_backward(Lw, cache.uw, bw)["backward"] > gates["linear_backward_residual"]):
            raise StrictQualificationFailure("Background primal backward residual failed")
        port_info = setup("one complete paid public Green-port setup",
                          lambda: cache.build_ports(max_bytes=max_bytes))
        write_json(OUT / "PUBLIC_PORTS.json", dict(**port_info, setup_calls=1,
                   shared_by=["woodbury", "schur_from_paid_ports"],
                   independent_exterior_schur_status="NOT_RUN; would require extra 2259 RHS"))

        case_inputs = {}

        def get_inputs(case):
            if case["id"] not in case_inputs:
                context.update(step="target assembly", case=case["id"], method=None)
                assembly, timing = measured(lambda: assembler.assemble(
                    case["pose"], case["radius"], p.k,
                    case["contrast"] + 1j * case["target_loss"],
                    lambda points: line_green(points, tx, p.k), source_count=2))
                destination = OUT / case["id"]
                destination.mkdir(exist_ok=False)
                save_npz(destination / "target_update.npz", assembly.Vtarget)
                np.savez(destination / "target_input.npz", bx=assembly.bX,
                         touched_free_indices=assembly.touched_free_indices,
                         quadrature_points=assembly.quadrature.points,
                         quadrature_weights=assembly.quadrature.weights)
                write_json(destination / "CASE.json", dict(case=case, assembly=assembly.metadata,
                           assembly_timing=timing, same_public_ports=len(ports),
                           rom_partition="TRAINING" if case["id"] in TRAINING_IDS else "HELD_OUT"))
                case_inputs[case["id"]] = assembly
            return case_inputs[case["id"]]

        # Exactly one full solve per condition.  The two frozen training full
        # solves occur first; their later baseline rows reuse saved paid solves.
        training_fields = []
        for case_id in TRAINING_IDS:
            case = next(case for case in cases if case["id"] == case_id)
            assembly = get_inputs(case)
            context.update(step="training full baseline", case=case_id, method="full")
            result, timing = measured(lambda: cache.full(assembly.Vtarget, assembly.bX))
            record(case, "full", result, result[:2], Lw - assembly.Vtarget, bw + assembly.bX, timing)
            training_fields.append(result[0])
        context.update(step="sixteen paid background receiver adjoint RHS", case=None, method=None)

        def solve_adjoints():
            book.reserve(PHASE, "goal ROM paid background receiver adjoints", rhs=16,
                         scope="all fixed receivers; no held-out target field")
            rhs = R.conj().T.toarray()
            fields = cache.factor.solve(rhs, trans="H")
            return fields, relative_backward(Lw.conj().T, fields, rhs), rhs

        (adjoints, adjoint_defect, adjoint_rhs) = setup("sixteen paid background receiver adjoint RHS", solve_adjoints)
        np.savez(OUT / "background_adjoints.npz", fields=adjoints, rhs=adjoint_rhs,
                 equation_residual=Lw.conj().T @ adjoints - adjoint_rhs)
        write_json(OUT / "BACKGROUND_ADJOINTS.json", dict(rhs=16, transpose="H",
                   backward=adjoint_defect, factorization="same paid background LU"))
        if (not finite_arrays(adjoints) or adjoint_defect["backward"] > gates["linear_backward_residual"]):
            raise StrictQualificationFailure("Background adjoint backward residual failed")
        labels = ([f"background_primal_source_{j}" for j in range(2)]
                  + [f"{case_id}_full_source_{j}" for case_id in TRAINING_IDS for j in range(2)]
                  + [f"background_receiver_adjoint_{j}" for j in range(16)])
        snapshots = np.column_stack([cache.uw, *training_fields, adjoints])
        q, qr = setup("fixed joint two-pass goal QR", lambda: two_pass_joint_qr(snapshots, labels))
        np.savez(OUT / "goal_basis.npz", Q=q)
        write_json(OUT / "GOAL_BASIS.json", dict(**qr, training_ids=TRAINING_IDS,
                   held_out_ids=HELD_OUT_IDS, frozen_before_held_out_full_solves=True,
                   training_full_solves_are_baseline_solves=True,
                   background_adjoint_setup_rhs=16, automatic_fallback=False))
        del snapshots, training_fields, adjoints, adjoint_rhs
        gc.collect()

        for case in cases:
            assembly = get_inputs(case)
            matrix, rhs = (Lw - assembly.Vtarget).tocsc(), bw + assembly.bX
            context.update(step="ordinary operator comparison", case=case["id"], method="full")
            if case["id"] in TRAINING_IDS:
                with np.load(OUT / case["id"] / "full.npz") as saved:
                    reference = saved["scattered_field"].copy(), saved["response"].copy()
            else:
                result, timing = measured(lambda: cache.full(assembly.Vtarget, assembly.bX))
                reference = result[:2]
                record(case, "full", result, reference, matrix, rhs, timing)
            for method in EXACT_METHODS[1:]:
                context["method"] = method
                result, timing = measured(lambda: getattr(cache, method)(assembly.Vtarget, assembly.bX))
                record(case, method, result, reference, matrix, rhs, timing)
                del result
            context["method"] = "goal_rom"
            result, timing = measured(lambda: goal_galerkin(matrix, rhs, R, direct, q, book, PHASE))
            record(case, "goal_rom", result, reference, matrix, rhs, timing, q)
            del result, reference, matrix, rhs
            if case["id"] == "low_loss_high_contrast":
                context.update(step="one low-loss/high-contrast core stability estimate", method="core_diagnostic")
                (diagnostic, arrays), timing = measured(lambda: inverse_core_diagnostic(cache, assembly, book))
                write_json(OUT / "CORE_STABILITY.json", dict(**diagnostic, **timing))
                np.savez(OUT / "core_stability.npz", **arrays)
                if diagnostic["actual_equation_backward"]["backward"] > gates["linear_backward_residual"]:
                    raise StrictQualificationFailure("Predeclared core diagnostic backward residual failed")
            gc.collect()
        status.update(operator_comparison_status="COMPLETED", condition_count=7,
                      method_result_count=len(rows), full_baseline_solves=7,
                      rom_basis_columns=q.shape[1], held_out_accuracy_preserved=True,
                      stress_panel_status="NOT_RUN_BY_THIS_DRIVER",
                      remaining="finite perturbation/Jacobian, information/posterior stress and final scientific review")
        write_json(OUT / "LEDGER_SNAPSHOT.json", dict(
            phase=PHASE, authorization=authorization, current=book.snapshot(PHASE),
            aggregate_cpu_note="Current wrapper completion CPU is not yet charged; use its final receipt"))
        persist()
        return 0
    except Exception as error:
        status.update(operator_comparison_status="STOPPED_FAILURE", last_context=context.copy(),
                      error_type=type(error).__name__, error=str(error),
                      failed_kernel_timing=getattr(error, "bgga_kernel_timing", None),
                      automatic_retry=False, fallback=False,
                      whole_G0_status="PARTIAL/FAILED_SUBCHECK; no whole-G0 pass")
        write_json(OUT / "STOP.json", dict(**status, saved_method_results=len(rows)))
        persist()
        print(json.dumps(status, default=_default, allow_nan=False), flush=True)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())

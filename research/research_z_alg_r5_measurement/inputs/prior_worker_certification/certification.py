"""Directed-interval LOCAL Krawczyk certification for scalar-TM real roots.

The certified statement is deliberately narrow.  At a supplied binary64 point
``m`` it proves, at most, that there is one real material root in the reported
small box.  The constants ``G``, ``E``, ``S`` and the target ``y`` are treated
as exact binary64 numbers.  No continuum Green-function or PDE error is
certified.
"""
from __future__ import annotations

from contextlib import contextmanager
import math
import time
import traceback
from typing import Any, Sequence

import mpmath as mp


STATUS_CERTIFIED = "CERTIFIED_LOCAL_REAL_ROOT"
STATUS_UNKNOWN = "UNKNOWN"
PRECISION_BITS = (64, 128, 256)
RADIUS_MULTIPLIERS = (1.0e-7, 1.0e-6, 1.0e-5)

SCOPE = (
    "LOCAL_REAL_ROOT_OF_THE_ALGEBRAIC_SCALAR_TM_MAP_ONLY; constants G,E,S,y "
    "are exact binary64 inputs; certification is valid only in the reported "
    "small real box X"
)

EXCLUSIONS = (
    "NO_GLOBAL_ROOT_COVERAGE_OR_EXHAUSTIVENESS",
    "NO_COVERAGE_OUTSIDE_REPORTED_X",
    "NO_COMPLEX_MATERIAL_ROOTS",
    "NO_HANKEL_INTEGRAL_OR_GEOMETRY_ACCURACY_CERTIFICATION",
    "NO_PDE_OR_DISCRETIZATION_ERROR_CERTIFICATION",
    "NO_PHYSICAL_MATERIAL_INTERPRETATION_BEYOND_BINARY64_SCALAR_TM",
    "NO_CERTIFICATION_OF_A_DIFFERENT_TARGET_OR_PERTURBED_INPUT",
    "NO_CLAIM_THAT_OTHER_ROOTS_OUTSIDE_X_DO_NOT_EXIST",
)


class IntervalFailure(RuntimeError):
    """An interval enclosure could not be made rigorous."""


def _nested(value: Any) -> Any:
    if hasattr(value, "tolist"):
        return value.tolist()
    return value


def _real_binary64(value: Any) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("All real inputs must be finite binary64 numbers")
    return value


def _complex_binary64(value: Any) -> complex:
    z = complex(value)
    real, imag = float(z.real), float(z.imag)
    if not math.isfinite(real) or not math.isfinite(imag):
        raise ValueError("All complex inputs must have finite binary64 parts")
    return complex(real, imag)


def _matrix(value: Any, rows: int, cols: int, name: str) -> list[list[Any]]:
    data = _nested(value)
    try:
        result = [list(row) for row in data]
    except TypeError as exc:
        raise ValueError(f"{name} must be a {rows}x{cols} matrix") from exc
    if len(result) != rows or any(len(row) != cols for row in result):
        raise ValueError(f"{name} must be a {rows}x{cols} matrix")
    return result


def _vector(value: Any, length: int, name: str) -> list[Any]:
    data = _nested(value)
    try:
        result = list(data)
    except TypeError as exc:
        raise ValueError(f"{name} must have length {length}") from exc
    if len(result) != length:
        raise ValueError(f"{name} must have length {length}")
    return result


def _prepare_inputs(
    G: Any, E: Any, S: Any, target: Any, m: Any
) -> tuple[
    list[list[complex]],
    list[list[complex]],
    list[list[complex]],
    list[complex],
    list[float],
    int,
    int,
]:
    m_values = [_real_binary64(value) for value in _vector(m, len(_nested(m)), "m")]
    n = len(m_values)
    if n not in (4, 6):
        raise ValueError("This API supports n=2s with n in {4,6}")
    s = n // 2

    target_rows = _matrix(target, 1, s, "target")
    y = [_complex_binary64(value) for value in target_rows[0]]
    G_rows = _matrix(G, n, n, "G")
    E_rows = _matrix(E, n, s, "E")
    S_rows = _matrix(S, 1, n, "S")

    G_value = [[_complex_binary64(value) for value in row] for row in G_rows]
    E_value = [[_complex_binary64(value) for value in row] for row in E_rows]
    S_value = [[_complex_binary64(value) for value in row] for row in S_rows]
    return G_value, E_value, S_value, y, m_values, n, s


def _b64_solve_matrix(
    A: Sequence[Sequence[complex]], B: Sequence[Sequence[complex]]
) -> list[list[complex]]:
    """Solve A X=B in binary64 using row-pivoted Gauss-Jordan operations."""
    n = len(A)
    k = len(B[0]) if B else 0
    if n == 0 or len(B) != n or any(len(row) != k for row in B):
        raise ValueError("Incompatible solve shapes")
    aug = [list(A[i]) + list(B[i]) for i in range(n)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(aug[row][col]))
        if aug[pivot][col] == 0:
            raise ValueError("Binary64 point matrix is singular")
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        divisor = aug[col][col]
        for j in range(col, n + k):
            aug[col][j] /= divisor
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if factor:
                for j in range(col, n + k):
                    aug[row][j] -= factor * aug[col][j]
    return [row[n:] for row in aug]


def _binary64_preconditioner(
    G: Sequence[Sequence[complex]],
    E: Sequence[Sequence[complex]],
    S: Sequence[Sequence[complex]],
    m: Sequence[float],
) -> tuple[list[list[float]], dict[str, Any]]:
    """Construct a binary64 inverse-real-Jacobian left preconditioner C."""
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    n, s = len(m), len(E[0])
    A = [
        [complex(1 if i == j else 0) - m[i] * G[i][j] for j in range(n)]
        for i in range(n)
    ]
    b = [[m[i] * E[i][source] for source in range(s)] for i in range(n)]
    at = [[A[j][i] for j in range(n)] for i in range(n)]
    current_columns = _b64_solve_matrix(A, b)
    current = [
        [current_columns[i][source] for i in range(n)] for source in range(s)
    ]
    transfer_column = [
        row[0]
        for row in _b64_solve_matrix(
            at, [[S[0][i]] for i in range(n)]
        )
    ]
    transfer = transfer_column
    E_total = [
        [
            E[i][source]
            + sum(G[i][j] * current[source][j] for j in range(n))
            for source in range(s)
        ]
        for i in range(n)
    ]
    J_complex = [
        [transfer[i] * E_total[i][source] for i in range(n)]
        for source in range(s)
    ]
    J_real = []
    for source in range(s):
        J_real.append([value.real for value in J_complex[source]])
        J_real.append([value.imag for value in J_complex[source]])

    # Invert the binary64 2s x n real Jacobian.  Any nonsingular C is a valid
    # Krawczyk preconditioner; being the point inverse affects width, not proof.
    identity_rhs = [
        [1.0 if i == j else 0.0 for j in range(n)] for i in range(n)
    ]
    C = _b64_solve_matrix(J_real, identity_rhs)
    metadata = {
        "arithmetic": "BINARY64_POINT_PRECONDITIONER_ONLY",
        "complex_LU_equivalents": 2,
        "real_inverse_LU_equivalents": 1,
        "complex_RHS": s + 1,
        "real_RHS": n,
        "cpu_seconds": time.process_time() - start_cpu,
        "wall_seconds": time.perf_counter() - start_wall,
    }
    return C, metadata


def _iv(value: Any) -> mp.iv.mpf:
    """Exact binary64-to-interval conversion at the active interval precision."""
    return mp.iv.mpf(float(value))


def _iv_complex(value: complex) -> mp.iv.mpc:
    return mp.iv.mpc(float(value.real), float(value.imag))


def _contains_zero(value: mp.iv.mpf | mp.iv.mpc) -> bool:
    if isinstance(value, mp.iv.mpc):
        return _contains_zero(value.real) and _contains_zero(value.imag)
    return bool(value.a.a <= 0 <= value.b.b)


def _abs_upper(value: mp.iv.mpf) -> mp.mpf:
    return max(abs(value.a.a), abs(value.a.b), abs(value.b.a), abs(value.b.b))


def _abs_lower(value: mp.iv.mpf | mp.iv.mpc) -> mp.mpf:
    if isinstance(value, mp.iv.mpc):
        return max(_abs_lower(value.real), _abs_lower(value.imag))
    if value.a.a <= 0 <= value.b.b:
        return mp.mpf("0")
    return min(abs(value.a.a), abs(value.b.b))


def _interval_inverse(
    A: Sequence[Sequence[mp.iv.mpc]],
) -> tuple[list[list[mp.iv.mpc]], list[dict[str, Any]]]:
    """Interval Gauss-Jordan inverse with every divisor kept strictly off 0."""
    n = len(A)
    zero, one = mp.iv.mpc(0), mp.iv.mpc(1)
    aug = [
        list(A[i]) + [one if i == j else zero for j in range(n)]
        for i in range(n)
    ]
    pivots: list[dict[str, Any]] = []
    for col in range(n):
        candidates = [row for row in range(col, n) if not _contains_zero(aug[row][col])]
        if not candidates:
            raise IntervalFailure(
                f"NO_NONZERO_INTERVAL_PIVOT_AT_STEP_{col + 1}"
            )
        pivot = max(candidates, key=lambda row: _abs_lower(aug[row][col]))
        if pivot != col:
            aug[col], aug[pivot] = aug[pivot], aug[col]
        divisor = aug[col][col]
        lower = _abs_lower(divisor)
        if lower <= 0:
            raise IntervalFailure(f"PIVOT_{col + 1}_CONTAINS_ZERO")
        for j in range(2 * n):
            aug[col][j] = aug[col][j] / divisor
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            for j in range(2 * n):
                aug[row][j] = aug[row][j] - factor * aug[col][j]
        pivots.append(
            {
                "step": col + 1,
                "source_row": pivot + 1,
                "lower_modulus_bound": lower,
            }
        )
    return [row[n:] for row in aug], pivots


def _interval_forward_jacobian(
    X: Sequence[mp.iv.mpf],
    G: Sequence[Sequence[complex]],
    E: Sequence[Sequence[complex]],
    S: Sequence[Sequence[complex]],
    target: Sequence[complex],
) -> tuple[
    list[mp.iv.mpf],
    list[list[mp.iv.mpf]],
    list[list[mp.iv.mpc]],
    list[dict[str, Any]],
]:
    """Enclose F(X) and the real material Jacobian over X."""
    n, s = len(X), len(E[0])
    A: list[list[mp.iv.mpc]] = []
    for i in range(n):
        row = []
        for j in range(n):
            value = _iv_complex(G[i][j]) * X[i]
            row.append(mp.iv.mpc(1 if i == j else 0) - value)
        A.append(row)
    invA, pivot_log = _interval_inverse(A)

    current: list[list[mp.iv.mpc]] = []
    for source in range(s):
        rhs = [X[j] * _iv_complex(E[j][source]) for j in range(n)]
        col = []
        for i in range(n):
            col.append(sum((invA[i][j] * rhs[j] for j in range(n)), mp.iv.mpc(0)))
        current.append(col)

    observed: list[mp.iv.mpc] = []
    for source in range(s):
        observed.append(
            sum(
                (_iv_complex(S[0][i]) * current[source][i] for i in range(n)),
                mp.iv.mpc(0),
            )
        )
    F_flat: list[mp.iv.mpf] = []
    for source in range(s):
        F_flat.extend(
            [observed[source].real - _iv(target[source].real),
             observed[source].imag - _iv(target[source].imag)]
        )

    E_total: list[list[mp.iv.mpc]] = []
    for i in range(n):
        values = []
        for source in range(s):
            scattered = sum(
                (_iv_complex(G[i][j]) * current[source][j] for j in range(n)),
                mp.iv.mpc(0),
            )
            values.append(_iv_complex(E[i][source]) + scattered)
        E_total.append(values)

    transfer: list[mp.iv.mpc] = []
    for i in range(n):
        transfer.append(
            sum(
                (_iv_complex(S[0][row]) * invA[row][i] for row in range(n)),
                mp.iv.mpc(0),
            )
        )
    J: list[list[mp.iv.mpf]] = []
    for source in range(s):
        J.append([(transfer[i] * E_total[i][source]).real for i in range(n)])
        J.append([(transfer[i] * E_total[i][source]).imag for i in range(n)])
    return F_flat, J, E_total, pivot_log


def _krawczyk(
    X: Sequence[mp.iv.mpf],
    x: Sequence[mp.iv.mpf],
    F_at_x: Sequence[mp.iv.mpf],
    J_X: Sequence[Sequence[mp.iv.mpf]],
    C: Sequence[Sequence[float]],
) -> tuple[list[mp.iv.mpf], list[mp.iv.mpf], list[mp.mpf]]:
    n, rows = len(X), len(J_X)
    CF = [
        sum((_iv(C[i][j]) * F_at_x[j] for j in range(n)), mp.iv.mpf(0))
        for i in range(rows)
    ]
    M = [
        [
            mp.iv.mpf(1 if i == j else 0)
            - sum((_iv(C[i][k]) * J_X[k][j] for k in range(rows)), mp.iv.mpf(0))
            for j in range(n)
        ]
        for i in range(rows)
    ]
    K = []
    row_sums = []
    for i in range(rows):
        displacement = sum(
            (M[i][j] * (X[j] - x[j]) for j in range(n)), mp.iv.mpf(0)
        )
        K.append(x[i] - CF[i] + displacement)
        row_sums.append(sum((_abs_upper(M[i][j]) for j in range(n)), mp.mpf(0)))
    return K, M, row_sums


def _report_interval(value: mp.iv.mpf) -> dict[str, float]:
    """Widen once outward to binary64 only for JSON reporting."""
    lower_endpoint, upper_endpoint = value.a.a, value.b.b
    if not (mp.isfinite(lower_endpoint) and mp.isfinite(upper_endpoint)):
        raise IntervalFailure("Nonfinite interval endpoint")
    lower_float = float(lower_endpoint)
    upper_float = float(upper_endpoint)
    if not (math.isfinite(lower_float) and math.isfinite(upper_float)):
        raise IntervalFailure("Interval endpoint is outside binary64 reporting range")
    lower = math.nextafter(lower_float, -math.inf)
    upper = math.nextafter(upper_float, math.inf)
    return {"lo": lower, "hi": upper}


def _populate_box(
    m: Sequence[float], radius: float, legal_box: tuple[float, float]
) -> list[mp.iv.mpf]:
    lower, upper = legal_box
    boxes = []
    for value in m:
        center = mp.iv.mpf(value)
        halfwidth = mp.iv.mpf(radius)
        legal_lower = mp.iv.mpf(lower)
        legal_upper = mp.iv.mpf(upper)
        lower_value = max((center - halfwidth).a.a, legal_lower.a.a)
        upper_value = min((center + halfwidth).b.b, legal_upper.b.b)
        boxes.append(mp.iv.mpf([lower_value, upper_value]))
    return boxes


def _strictly_inside(
    K: Sequence[mp.iv.mpf], X: Sequence[mp.iv.mpf]
) -> tuple[bool, list[dict[str, Any]]]:
    rows = []
    strict = True
    for i, (k, x) in enumerate(zip(K, X)):
        lower_margin = k.a - x.a
        upper_margin = x.b - k.b
        row_strict = bool(k.a.a > x.a.b and x.b.a > k.b.b)
        strict = strict and row_strict
        rows.append(
            {
                "coordinate": i + 1,
                "strict": row_strict,
                "lower_margin": lower_margin,
                "upper_margin": upper_margin,
            }
        )
    return strict, rows


@contextmanager
def _interval_precision(bits: int):
    old = mp.iv.prec
    mp.iv.prec = bits
    try:
        yield
    finally:
        mp.iv.prec = old


def certify(
    G: Any,
    E: Any,
    S: Any,
    target: Any,
    m: Any,
    box: tuple[float, float] = (0.5, 30.0),
    timeout: float = 30.0,
) -> dict[str, Any]:
    """Certify at most one real scalar-TM root in a small box around ``m``.

    G, E, S and target are interpreted as binary64 constants.  The returned
    status is either ``CERTIFIED_LOCAL_REAL_ROOT`` or ``UNKNOWN``.
    """
    wall_start, cpu_start = time.perf_counter(), time.process_time()
    deadline = wall_start + float(timeout)
    result: dict[str, Any] = {
        "status": STATUS_UNKNOWN,
        "scope": SCOPE,
        "exclusions": list(EXCLUSIONS),
        "precision_schedule_bits": list(PRECISION_BITS),
        "radius_multipliers": list(RADIUS_MULTIPLIERS),
        "timeout_seconds": float(timeout),
        "attempts": [],
        "counter": {
            "interval_matrix_inversions": 0,
            "interval_point_evaluations": 0,
            "interval_box_evaluations": 0,
            "preconditioner_complex_LU_equivalents": 0,
            "preconditioner_real_LU_equivalents": 0,
            "preconditioner_complex_RHS": 0,
            "preconditioner_real_RHS": 0,
        },
    }
    try:
        G_value, E_value, S_value, y, m_value, n, s = _prepare_inputs(
            G, E, S, target, m
        )
        lower, upper = _real_binary64(box[0]), _real_binary64(box[1])
        if not lower < upper:
            raise ValueError("box must have lower < upper")
        result["shapes"] = {
            "G": [n, n],
            "E": [n, s],
            "S": [1, n],
            "target": [1, s],
            "m": [n],
            "real_equation_rows": [2 * s],
            "real_row_order": "SOURCE-MAJOR; REAL THEN IMAGINARY",
        }
        result["input_arithmetic"] = (
            "BINARY64_CONSTANTS_CONVERTED_EXACTLY_TO_DIRECTED_INTERVALS"
        )
        if any(not (lower <= value <= upper) for value in m_value):
            result["failure_reasons"] = ["CENTER_OUTSIDE_GLOBAL_BOX"]
            return result

        C, preconditioner_cost = _binary64_preconditioner(
            G_value, E_value, S_value, m_value
        )
        result["preconditioner"] = {
            "method": "BINARY64_INVERSE_REAL_POINT_JACOBIAN",
            "C": C,
            **preconditioner_cost,
        }
        result["counter"].update(
            {
                "preconditioner_complex_LU_equivalents": 2,
                "preconditioner_real_LU_equivalents": 1,
                "preconditioner_complex_RHS": s + 1,
                "preconditioner_real_RHS": n,
            }
        )

        norm_inf = max(1.0, max(abs(value) for value in m_value))
        for bits in PRECISION_BITS:
            if time.perf_counter() >= deadline:
                result["failure_reasons"] = ["TIMEOUT_BEFORE_PRECISION"]
                break
            with _interval_precision(bits):
                x = [_iv(value) for value in m_value]
                point_start_wall, point_start_cpu = (
                    time.perf_counter(),
                    time.process_time(),
                )
                try:
                    F_x, _, _, point_pivots = _interval_forward_jacobian(
                        x, G_value, E_value, S_value, y
                    )
                    inversion_ok = True
                except (IntervalFailure, ZeroDivisionError) as exc:
                    inversion_ok = False
                    point_failure = str(exc)
                    F_x, point_pivots = [], []
                result["counter"]["interval_matrix_inversions"] += 1
                result["counter"]["interval_point_evaluations"] += 1
                point_cpu = time.process_time() - point_start_cpu
                point_wall = time.perf_counter() - point_start_wall

                if not inversion_ok:
                    for multiplier in RADIUS_MULTIPLIERS:
                        radius = float(multiplier * norm_inf)
                        X = _populate_box(m_value, radius, (lower, upper))
                        result["attempts"].append(
                            {
                                "precision_bits": bits,
                                "radius": radius,
                                "status": STATUS_UNKNOWN,
                                "x": [_report_interval(value) for value in X],
                                "reasons": [
                                    "POINT_INTERVAL_NONSINGULARITY_FAILED",
                                    point_failure,
                                ],
                                "point_cpu_seconds": point_cpu,
                                "point_wall_seconds": point_wall,
                            }
                        )
                    continue

                for multiplier in RADIUS_MULTIPLIERS:
                    if time.perf_counter() >= deadline:
                        result["failure_reasons"] = ["TIMEOUT_BEFORE_RADIUS_TRIAL"]
                        break
                    radius = float(multiplier * norm_inf)
                    attempt_start_wall, attempt_start_cpu = (
                        time.perf_counter(),
                        time.process_time(),
                    )
                    reasons: list[str] = []
                    X = _populate_box(m_value, radius, (lower, upper))
                    if any(value.a >= value.b for value in X):
                        attempt = {
                            "precision_bits": bits,
                            "radius": radius,
                            "status": STATUS_UNKNOWN,
                            "reasons": ["RADIUS_TRIAL_DEGENERATE_OR_OUTSIDE_BOX"],
                            "x": [_report_interval(value) for value in X],
                        }
                        result["attempts"].append(attempt)
                        continue
                    try:
                        _, J_X, _, box_pivots = _interval_forward_jacobian(
                            X, G_value, E_value, S_value, y
                        )
                        result["counter"]["interval_matrix_inversions"] += 1
                        result["counter"]["interval_box_evaluations"] += 1
                    except (IntervalFailure, ZeroDivisionError) as exc:
                        result["counter"]["interval_matrix_inversions"] += 1
                        result["attempts"].append(
                            {
                                "precision_bits": bits,
                                "radius": radius,
                                "status": STATUS_UNKNOWN,
                                "reasons": ["BOX_INTERVAL_NONSINGULARITY_FAILED", str(exc)],
                                "x": [_report_interval(value) for value in X],
                                "cpu_seconds": time.process_time() - attempt_start_cpu,
                                "wall_seconds": time.perf_counter() - attempt_start_wall,
                            }
                        )
                        continue

                    K, M, row_sums = _krawczyk(X, x, F_x, J_X, C)
                    strict, coordinate_margins = _strictly_inside(K, X)
                    contraction = all(value < 1 for value in row_sums)
                    if not strict:
                        reasons.append("STRICT_INCLUSION_FAILED")
                    if not contraction:
                        reasons.append("ROW_SUM_CONTRACTION_NOT_STRICT")
                    if time.perf_counter() > deadline:
                        reasons.append("TIMEOUT_AFTER_ATTEMPT")

                    certified = not reasons
                    attempt = {
                        "precision_bits": bits,
                        "radius": radius,
                        "status": STATUS_CERTIFIED if certified else STATUS_UNKNOWN,
                        "reasons": reasons,
                        "x": [_report_interval(value) for value in X],
                        "point_F": [_report_interval(value) for value in F_x],
                        "point_nonsingularity_pivots": [
                            {
                                **entry,
                                "lower_modulus_bound": float(
                                    entry["lower_modulus_bound"]
                                ),
                            }
                            for entry in point_pivots
                        ],
                        "box_nonsingularity_pivots": [
                            {
                                **entry,
                                "lower_modulus_bound": float(
                                    entry["lower_modulus_bound"]
                                ),
                            }
                            for entry in box_pivots
                        ],
                        "K": [_report_interval(value) for value in K],
                        "strict_inclusion": strict,
                        "inclusion_coordinate_margins": [
                            {
                                **entry,
                                "lower_margin": _report_interval(
                                    entry["lower_margin"]
                                ),
                                "upper_margin": _report_interval(
                                    entry["upper_margin"]
                                ),
                            }
                            for entry in coordinate_margins
                        ],
                        "contraction_row_sums": [
                            _report_interval(mp.iv.mpf(value))
                            for value in row_sums
                        ],
                        "contraction_row_sum_upper": [float(value) for value in row_sums],
                        "contraction_strict": contraction,
                        "J_interval": [
                            [_report_interval(value) for value in row] for row in J_X
                        ],
                        "M_interval": [
                            [_report_interval(value) for value in row] for row in M
                        ],
                        "cpu_seconds": time.process_time() - attempt_start_cpu,
                        "wall_seconds": time.perf_counter() - attempt_start_wall,
                    }
                    result["attempts"].append(attempt)
                    if certified:
                        result.update(
                            {
                                "status": STATUS_CERTIFIED,
                                "precision_bits": bits,
                                "radius": radius,
                                "intervals": {
                                    "X": attempt["x"],
                                    "F_at_x": attempt["point_F"],
                                    "J_over_X": attempt["J_interval"],
                                    "K": attempt["K"],
                                    "contraction_matrix_M": attempt["M_interval"],
                                },
                                "uniqueness": (
                                    "STRICT_KRAWCZYK_INCLUSION_AND_STRICT_REAL_ROW_SUM_CONTRACTION"
                                ),
                            }
                        )
                        break
                if result["status"] == STATUS_CERTIFIED:
                    break
        if result["status"] != STATUS_CERTIFIED:
            result.setdefault("failure_reasons", ["NO_RADIUS_AND_PRECISION_ATTEMPT_CERTIFIED"])
        result["point_evaluation_timing_last"] = {
            "cpu_seconds": locals().get("point_cpu", 0.0),
            "wall_seconds": locals().get("point_wall", 0.0),
        }
        return result
    except Exception as exc:
        result["status"] = STATUS_UNKNOWN
        result["failure_reasons"] = [f"INPUT_OR_INTERVAL_ERROR: {type(exc).__name__}: {exc}"]
        result["exception_traceback"] = traceback.format_exc()
        return result
    finally:
        result["cpu_seconds"] = time.process_time() - cpu_start
        result["wall_seconds"] = time.perf_counter() - wall_start
        result["interval_reporting"] = (
            "operations and all proof checks use mpmath.iv; endpoints are widened "
            "outward to binary64 only when serializing this report"
        )


__all__ = [
    "certify",
    "STATUS_CERTIFIED",
    "STATUS_UNKNOWN",
    "PRECISION_BITS",
    "RADIUS_MULTIPLIERS",
    "SCOPE",
    "EXCLUSIONS",
]

"""Conditional, separately accounted A30 24x24 TM falsification experiment.

Importing this module performs no physics.  ``run`` requires the parent's E1
scientific entry decision.  Original A30 data, decisions and ledgers are read
only.  These are exposed development cases, not independent validation.
"""
from __future__ import annotations

from collections import Counter
from contextlib import contextmanager
from pathlib import Path
import argparse
import csv
import importlib
import json
import os
import subprocess
import sys
import time

import numpy as np
from scipy import ndimage
from scipy.optimize import LinearConstraint, minimize
from scipy.stats import chi2


SOURCE = Path(os.environ.get('A30_ROOT',str(Path(__file__).resolve().parents[4] / "A30/a30-maxwell-material-manifold")))
STAGE_REL = Path("research/research_a_20261015")
PARENTS = (250065, 250069, 250073, 250077)
ARMS = ("finite_history", "multistart", "weak_null", "same_edit_no_compensation")
EXPECTED_A30_COMMIT = "717d3d26a290248bd2b9aeb55f3d03694b407b84"


def _plain(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        return {str(k): _plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_plain(v) for v in x]
    return x


def _write(path, obj):
    Path(path).write_text(json.dumps(_plain(obj), indent=2, allow_nan=False) + "\n")


class E2BudgetStop(RuntimeError):
    pass


class ResearchABook:
    """Append-only, nonresettable new-stage accounting required by A30 adapter."""

    def __init__(self, output, state_limit=384, rhs_limit=9216, wall_limit=3600.):
        self.root = Path(output).resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.ledger = self.root / "COST_LEDGER.jsonl"
        self.counts = Counter()
        self.state_limit = min(384, int(state_limit))
        self.rhs_limit = min(9216, int(rhs_limit))
        self.wall_limit = min(3600., float(wall_limit))
        self.started = time.perf_counter()
        self.cpu_started = time.process_time()
        self.prior_wall = sum(r.get("execution_wall_seconds", 0.) for r in self.rows()
                              if r.get("event") == "EXECUTION_END")
        pending = [r for r in self.rows() if r.get("event") == "EXECUTION_BEGIN"]
        ended = [r for r in self.rows() if r.get("event") == "EXECUTION_END"]
        if len(pending) != len(ended):
            raise E2BudgetStop("UNCLOSED_PREVIOUS_EXECUTION_REQUIRES_ACCOUNTING_REPAIR")
        self.append({"event": "EXECUTION_BEGIN", "pid": os.getpid(),
                     "prior_execution_wall_seconds": self.prior_wall})

    def rows(self):
        return [json.loads(s) for s in self.ledger.read_text().splitlines()] if self.ledger.exists() else []

    def append(self, row):
        row = dict(row, timestamp_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
        with self.ledger.open("a") as f:
            f.write(json.dumps(_plain(row), allow_nan=False) + "\n")
            f.flush()
            os.fsync(f.fileno())

    def check(self):
        if self.prior_wall + time.perf_counter() - self.started >= self.wall_limit:
            raise E2BudgetStop("E2_CUMULATIVE_EXECUTION_WALL_LIMIT")

    def _reserve(self, event, amount, limit, **metadata):
        self.check()
        if int(amount) != amount or amount < 1:
            raise ValueError("POSITIVE_INTEGER_RESERVATION_REQUIRED")
        lock = self.root / "accounting.lock"
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        try:
            os.write(fd, str(os.getpid()).encode())
            os.close(fd)
            used = sum(r.get("amount", 0) for r in self.rows() if r.get("event") == event)
            if used + amount > limit:
                raise E2BudgetStop(event + "_HARD_LIMIT")
            self.append(dict(event=event, amount=int(amount), cumulative=used + amount,
                             limit=limit, **metadata))
        finally:
            lock.unlink(missing_ok=True)

    def reserve_state(self, material, *, resolution, purpose):
        self._reserve("STATE_RESERVED", 1, self.state_limit, resolution=resolution,
                      purpose=purpose, material_shape=list(np.shape(material)))

    def reserve_rhs(self, count, role):
        self._reserve("RHS_RESERVED", int(count), self.rhs_limit, role=role)

    @contextmanager
    def span(self, label, **counts):
        self.check()
        self.counts.update(counts)
        start, cpu = time.perf_counter(), time.process_time()
        status = "FAILED"
        try:
            yield
            status = "COMPLETE"
        finally:
            self.append(dict(event="ACTION", label=label, status=status, counts=counts,
                             wall_seconds=time.perf_counter() - start,
                             process_cpu_seconds=time.process_time() - cpu))

    def finish(self, status):
        self.append(dict(event="EXECUTION_END", status=status,
                         execution_wall_seconds=time.perf_counter() - self.started,
                         process_cpu_seconds=time.process_time() - self.cpu_started,
                         aggregate_counts=dict(self.counts)))

    def receipt(self):
        rows = self.rows()
        return dict(new_full_state_attempts=sum(r.get("amount", 0) for r in rows if r.get("event") == "STATE_RESERVED"),
                    new_complex_LU_RHS=sum(r.get("amount", 0) for r in rows if r.get("event") == "RHS_RESERVED"),
                    controlled_execution_wall_seconds=self.prior_wall + time.perf_counter() - self.started,
                    state_limit=self.state_limit, rhs_limit=self.rhs_limit,
                    wall_limit_seconds=self.wall_limit, GPU_seconds=0,
                    old_A30_ledger_modified=False, old_P2_credit_used=0)


def _modules(a30_root):
    source = Path(a30_root).resolve() / "src"
    existing = sys.modules.get("a30")
    if existing is not None and not Path(existing.__file__).resolve().is_relative_to(source):
        raise ImportError("FOREIGN_A30_MODULE_COLLISION")
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
    physics = importlib.import_module("a30.physics")
    representations = importlib.import_module("a30.representations")
    for module in (physics, representations):
        if not Path(module.__file__).resolve().is_relative_to(source):
            raise ImportError("UNFROZEN_A30_IMPORT")
    return physics, representations


def pack(data):
    """A30: each source's Re receivers then Im receivers (not B1 packing)."""
    y = np.asarray(data)
    return np.concatenate((y.real.T, y.imag.T), axis=1).reshape(-1)


def material_distance(a, b):
    return float(np.sqrt(np.mean(np.abs(np.asarray(a) - b) ** 2)) / 1.2)


def _legal_step(image, update):
    """Contract one material-space ray; never clip pixels or fill protected holes."""
    x, d = np.r_[image.real, image.imag], np.r_[update.real, update.imag]
    lo, hi = np.r_[np.full(576, -.5), np.zeros(576)], np.r_[np.full(576, 4.), np.full(576, .5)]
    # Remove only optimizer roundoff normal to a binding box face. A finite
    # outward displacement still contracts tozero and is recorded as failure.
    tiny_outward = (((x <= lo+1e-12) & (d < 0)) | ((x >= hi-1e-12) & (d > 0))) & (np.abs(d) <= 1e-10)
    d[tiny_outward] = 0.
    alpha = 1.
    positive, negative = d > 0, d < 0
    if positive.any():
        alpha = min(alpha, float(np.min((hi[positive] - x[positive]) / d[positive])))
    if negative.any():
        alpha = min(alpha, float(np.min((lo[negative] - x[negative]) / d[negative])))
    contracted = max(0., alpha) * (1 - 1e-12) * d
    return contracted[:576] + 1j * contracted[576:]


def _bases(representations):
    yy, xx = np.meshgrid(np.arange(24) + .5, np.arange(24) + .5, indexing="ij")
    modes = sorted(((a, b) for a in range(5) for b in range(5)), key=lambda x: (sum(x), x))[0:8]
    low = np.array([(np.cos(np.pi * a * yy / 24) * np.cos(np.pi * b * xx / 24)).ravel()
                    for a, b in modes]).T
    low /= np.maximum(np.linalg.norm(low, axis=0), 1e-300)
    transform, _ = representations.haar24_transform()
    # The fixed multiscale order is representation-only, independent of truth.
    haar = np.asarray(transform).T[:, :16]
    if haar.shape != (576, 16):
        raise ValueError("A30_HAAR_BASIS_LAYOUT")
    return (np.c_[low, 1j * low], np.c_[haar, 1j * haar])


def _jacobian(model, state, basis, white):
    transfer = model.receiver_transfer(state)
    raw = np.einsum("qn,nk,ns->qsk", transfer, basis, state.field)
    return np.column_stack([pack(raw[:, :, k]) for k in range(basis.shape[1])]) * white


def _gn_proposal(model, state, basis, observation, white, trust=.35, damping=1e-3):
    jac = _jacobian(model, state, basis, white)
    residual = pack(state.data - observation) * white
    scale = max(float(np.linalg.norm(jac, 2)) ** 2, 1e-12)
    hessian = (jac.T @ jac) / scale + damping * np.eye(jac.shape[1])
    gradient = (jac.T @ residual) / scale
    material_basis = np.r_[basis.real, basis.imag]
    material = np.r_[state.chi.real, state.chi.imag]
    low = np.r_[np.full(576, -.5), np.zeros(576)] - material
    high = np.r_[np.full(576, 4.), np.full(576, .5)] - material
    active = np.linalg.norm(material_basis, axis=1) > 1e-14
    constraint = LinearConstraint(material_basis[active], np.maximum(low[active], -trust),
                                  np.minimum(high[active], trust))
    fit = minimize(lambda z: float(.5*z@hessian@z + gradient@z),
                   np.zeros(jac.shape[1]), jac=lambda z: hessian@z+gradient,
                   constraints=constraint, method="SLSQP", options=dict(maxiter=60, ftol=1e-10))
    update = basis @ fit.x
    model.book.append(dict(event="MATERIAL_QP", converged=bool(fit.success), status=str(fit.message),
                           iterations=int(fit.nit), real_parameters=basis.shape[1], trust=trust,
                           predicted_objective=float(fit.fun)))
    return state.chi + _legal_step(state.chi, update)


def _project_material(image):
    """Declared pixel-box projection for global starts; not a current null claim."""
    return np.clip(image.real, -.5, 4.) + 1j * np.clip(image.imag, 0., .5)


def _backend_audit(model, anchor, book):
    """One new3-cell real-Maxwell state; finite secant/Schur/dual checks.

    Cost:1freshstate and51LU RHS (8forward +8secant +8tangent +8remainder
    +3anchor columns +8smallcore +8dual). Existing N3 receiver transfer reused.
    """
    import scipy.linalg as la
    rng = np.random.default_rng(261015299)
    chosen = np.array([122, 276, 425])
    h = np.zeros(576, dtype=np.complex128)
    # Strictly inward legal finite changes, independent oftruth/support.
    h[chosen] = .10 * np.where(anchor.chi[chosen].real < 1.75, 1., -1.) + \
        .02j * np.where(anchor.chi[chosen].imag < .25, 1., -1.)
    second = model.full_state(anchor.chi + h, purpose="ResearchA_E2_backend_finite3cell_audit")
    rhs = h[:, None] * anchor.field
    secant = model.solve_anchor(rhs, second, role="finite_secant_audit")
    tangent = model.solve_anchor(rhs, anchor, role="material_tangent_audit")
    remainder = model.solve_anchor(h[:, None] * (model.G @ tangent), second, role="finite_remainder_audit")
    columns = np.eye(576, dtype=np.complex128)[:, chosen]
    inverse_columns = model.solve_anchor(columns, anchor, role="local_Schur_anchor_columns")
    core = np.eye(3) - h[chosen, None] * (model.G[chosen] @ inverse_columns)
    with book.span("local_exact_Schur_core_factor", small_core_LU_factorizations=1, small_core_dimension=3):
        factor = la.lu_factor(core)
    book.reserve_rhs(8, "local_Schur_small_core_forward")
    with book.span("local_exact_Schur_core_solve", small_core_RHS=8):
        coefficients = la.lu_solve(factor, h[chosen, None] * anchor.field[chosen])
    schur_secant = inverse_columns @ coefficients
    # A tiny second-order remainder loses digits when subtracted from the
    # full output. Test the dual identity at a resolvable arbitrary primal
    # residual; do not relax its tolerance or alter a physical observation.
    perturbation=(rng.standard_normal(second.current.shape)+1j*rng.standard_normal(second.current.shape))
    perturbation*=.01*np.linalg.norm(second.current)/np.linalg.norm(perturbation)
    approximate_current = anchor.current + tangent + perturbation
    approximate_data = model.S @ approximate_current
    primal_residual = second.chi[:, None] * model.incident - second.matrix @ approximate_current
    load = rng.standard_normal((32, 8)) + 1j * rng.standard_normal((32, 8))
    dual = model.solve_anchor(model.S.conj().T @ load, second, adjoint=True, role="dual_correction_audit")
    lhs = np.vdot(second.data - approximate_data, load)
    rhs_dual = np.vdot(primal_residual, dual)
    relative = lambda a,b: float(np.linalg.norm(a-b) / max(1e-15, np.linalg.norm(b)))
    direction = (rng.standard_normal(576) + 1j * rng.standard_normal(576)) * .001
    dot = model.cached_material_dot_identity(direction, load, anchor)
    errors = dict(finite_current_secant=relative(second.current-anchor.current, secant),
        finite_data_secant=relative(second.data-anchor.data, model.S @ secant),
        exact_nonlinear_remainder=relative(second.data-anchor.data-model.S@tangent, model.S@remainder),
        exact_local_Schur=relative(schur_secant, secant),
        dual_correction=relative(lhs, rhs_dual), material_mass_adjoint=dot["relative_error"],
        material_realizability=float(np.linalg.norm(second.current-second.chi[:,None]*second.field) /
                                      max(1e-15,np.linalg.norm(second.current))),
        backward_residual=float(second.residual))
    result = dict(status="PASS" if max(errors.values()) < 1e-8 else "FAIL", errors=errors,
        physics="Actual saved24x24TM geometry; one newlycharged3cell changed state",
        source_law="chi linear; c=sqrt(v)j, field=sqrt(v)E, v=1/576",
        affected_cells=chosen.tolist(), newlycharged_states=1, newlycharged_RHS=51,
        synthetic_fixture=False, gain="fixed1", whitening="sqrt(2)/sigma; no renoising")
    if result["status"] != "PASS":
        raise RuntimeError("E2_ACTUAL_BACKEND_IDENTITY_FAILURE:"+json.dumps(result))
    return result


def _hole_edit(parent, background, rng, index):
    magnitude = np.abs(parent - background).reshape(24, 24)
    allowed = np.argwhere(magnitude > .15 * 1.2)
    if len(allowed) == 0:
        allowed = np.array([[12, 12]])
    # Rotate quantile-ranked eligible sites, then jitter; no oracle mask/family.
    order = np.argsort(magnitude[allowed[:, 0], allowed[:, 1]], kind="stable")
    position = allowed[order[int((.2 + .2 * (index % 4)) * (len(order) - 1))]]
    position = np.clip(position + rng.integers(-2, 3, 2), 2, 21)
    radius = (1.5, 2.2, 2.2, 1.5)[index % 4]
    yy, xx = np.mgrid[:24, :24]
    mask = ((yy - position[0]) ** 2 + (xx - position[1]) ** 2 <= radius ** 2).ravel()
    edited = np.array(parent, copy=True)
    edited[mask] = background[mask]
    return edited, mask, dict(center_pixels=position.tolist(), radius_pixels=radius,
                             protected_cells=int(mask.sum()), edit="object-removal-to-known-background")


def _material_metrics(image, truth, background):
    error = image - truth
    delta = truth - background
    roi = np.abs(delta) > .25 * float(np.max(np.abs(delta)))
    predicted = np.abs(image-background) > .25 * float(np.max(np.abs(delta)))
    actual2, predicted2 = roi.reshape(24,24), predicted.reshape(24,24)
    structure4 = ndimage.generate_binary_structure(2, 1)
    edge_actual = actual2 & ~ndimage.binary_erosion(actual2, structure=structure4, border_value=0)
    edge_predicted = predicted2 & ~ndimage.binary_erosion(predicted2, structure=structure4, border_value=0)
    cavity = ndimage.binary_fill_holes(actual2) & ~actual2
    union = np.count_nonzero(roi | predicted)
    def holes(mask):
        filled = ndimage.binary_fill_holes(mask)
        labels, count = ndimage.label(filled & ~mask, structure=structure4)
        areas = [int(np.count_nonzero(labels==i)) for i in range(1,count+1)]
        return sorted(areas, reverse=True)
    boundary_distance = None
    if edge_actual.any() and edge_predicted.any():
        distance_actual, distance_predicted = ndimage.distance_transform_edt(~edge_actual), ndimage.distance_transform_edt(~edge_predicted)
        boundary_distance = float((distance_actual[edge_predicted].mean()+distance_predicted[edge_actual].mean())/48)
    return dict(full_complex_MAE=float(np.mean(np.abs(error))), full_Re_MAE=float(np.mean(np.abs(error.real))),
                full_Im_MAE=float(np.mean(np.abs(error.imag))),
                complex_permittivity_NMAE=float(np.mean(np.abs(error)) / max(np.mean(np.abs(delta)), 1e-300)),
                object_ROI_complex_MAE=float(np.mean(np.abs(error[roi]))) if roi.any() else None,
                cavity_ROI_complex_MAE=float(np.mean(np.abs(error.reshape(24,24)[cavity]))) if cavity.any() else None,
                support_IoU=float(np.count_nonzero(roi & predicted)/union) if union else 1.,
                symmetric_boundary_distance_domain_units=boundary_distance,
                hole_areas_pixels=holes(predicted2), truth_hole_areas_pixels=holes(actual2),
                material_distance_truth=material_distance(image, truth))


def _select(states, base, k):
    if not states:
        return []
    first = min(range(len(states)), key=lambda i: material_distance(states[i]["chi"], base))
    chosen = [first]
    while len(chosen) < min(k, len(states)):
        candidates = [i for i in range(len(states)) if i not in chosen]
        chosen.append(max(candidates, key=lambda i: min(material_distance(states[i]["chi"], states[j]["chi"])
                                                       for j in chosen)))
    return chosen


def run(root, *, e1_gate_path, a30_root=SOURCE):
    """Run once after parent records ``E2_allowed: true`` in its E1 decision.

    Maximum planned cost: 329 new states and 9139 new complex LU RHS, including
    eight cached-LU audits, coarse endpoints, all rejected corrections and four
    matched uncompensated edits per case.  Caps stay384/9216/3600s.
    """
    root, a30_root = Path(root).resolve(), Path(a30_root).resolve()
    gate_path = Path(e1_gate_path).resolve()
    gate = json.loads(gate_path.read_text())
    if gate.get("E2_allowed") is not True:
        raise E2BudgetStop("PARENT_E1_SCIENTIFIC_ENTRY_NOT_GRANTED")
    output = root / STAGE_REL / "results/E2"
    output.mkdir(parents=True, exist_ok=True)
    if (output / "RESULTS.json").exists():
        raise E2BudgetStop("COMPLETED_E2_MUST_NOT_BE_RERUN")
    book = ResearchABook(output)
    status, rows, arrays, failures = "FAILED", [], {}, []
    start = time.perf_counter()
    try:
        physics, representations = _modules(a30_root)
        basis_pair = _bases(representations)
        release = a30_root / "inputs/a29_release"
        cases = [f"p{p}_{bg}_a1_snr20" for p in PARENTS for bg in ("weak", "ring")]
        metadata = {x["case_id"]: x for x in csv.DictReader((release / "results/a29/COVERAGE_CASES.csv").open())}
        _write(output / "PROTOCOL.json", dict(protocol="Research-A-E2-A30-TM-v1", source_commit=EXPECTED_A30_COMMIT,
            parent_entry_decision=str(gate_path), cases=cases, arms=ARMS,
            measurement="Original saved complex observations; no renoising or rescaling",
            white="sqrt(2)/case_sigma", gain="fixed1; original known calibration",
            compatibility_radius=float(np.sqrt(chi2.ppf(.95, 512))),
            material_dedup_rms_normalized_1p2=.03, coarse="frozen N3 plus1bounded material-GN step",
            corrected_arm_states_max=12, trajectories_per_arm=4, full_states_per_trajectory=3,
            rounds="two trajectories in16real DCT coordinates then two in32real Haar coordinates",
            uncompensated="four exactly matched finite edits;4paidstates, not gratuitously repeated12times",
            planned_fresh_states=329, planned_complex_LU_RHS=9139, limits=dict(states=384,RHS=9216,executionwall=3600),
            posterior_probability=None, exposure="HISTORICALLY_EXPOSED_DEVELOPMENT",
            E3="NOT_AUTHORIZED", GPU=0, truth="separate final offline scorer only"))
        adapters = {}
        threshold = float(np.sqrt(chi2.ppf(.95, 512)))
        for case_i, case in enumerate(cases):
            book.check()
            bgid = metadata[case]["background_id"]
            if bgid not in adapters:
                whitelist = [release / "inputs/a28/results/a28/exact_states" / (c + ".npz")
                             for c in cases if metadata[c]["background_id"] == bgid]
                adapters[bgid] = physics.FrozenMaxwellAdapter.from_background_npz(
                    release / "inputs/a25/background" / (bgid + ".npz"), book,
                    allowed_cache_paths=whitelist, purpose_prefix="Research_A_E2")
            model = adapters[bgid]
            with np.load(release / "results/a29/pools" / (case + ".npz"), allow_pickle=False) as z:
                observation, sigma, background = z["observation"].copy(), float(z["sigma"]), z["background"].copy()
                n3 = z["n3"].copy()
            white = np.sqrt(2.) / sigma
            anchor = model.cached_anchor(n3, release / "inputs/a28/results/a28/exact_states" / (case + ".npz"))
            anchor_jac = [_jacobian(model, anchor, basis, white) for basis in basis_pair]
            if case_i == 0:
                _write(output / "ACTUAL_BACKEND_AUDIT.json", _backend_audit(model, anchor, book))
            proposal = _gn_proposal(model, anchor, basis_pair[1], observation, white)
            coarse_trial = model.full_state(proposal, purpose=case + ":coarse_GN")
            norm = lambda state: float(np.linalg.norm(pack(state.data - observation) * white))
            coarse = coarse_trial if norm(coarse_trial) < norm(anchor) else anchor
            arrays[case + "__observation"], arrays[case + "__background"] = observation, background
            arrays[case + "__N3"], arrays[case + "__base"] = n3, np.array(coarse.chi)
            arrays[case + "__base_prediction"] = np.array(coarse.data)
            rows.append(dict(case_id=case, arm="base", trajectory=-1, step=0, residual=norm(coarse),
                             compatible=norm(coarse) <= threshold, fresh_states=1))
            finite_edits = []
            for arm_i, arm in enumerate(ARMS[:3]):
                rng = np.random.default_rng(261015201 + 1000 * case_i + arm_i)
                archive = []
                if norm(coarse) <= threshold:
                    archive.append(dict(chi=np.array(coarse.chi), data=np.array(coarse.data), residual=norm(coarse), name="base"))
                for trajectory in range(4):
                    round_id = trajectory // 2
                    basis = np.array(basis_pair[round_id], copy=True)
                    protected = np.zeros(576, dtype=bool)
                    if arm == "finite_history":
                        parent = coarse.chi if not archive else archive[trajectory % len(archive)]["chi"]
                        proposal, protected, edit = _hole_edit(parent, background, rng, trajectory)
                        basis[protected] = 0
                        finite_edits.append((proposal.copy(), protected.copy(), edit))
                    elif arm == "multistart":
                        # Independent global coefficient restart; no oracle chart/support.
                        coefficients = rng.standard_normal(basis.shape[1])
                        update = basis @ coefficients
                        update *= .8 / max(float(np.max(np.abs(update))), 1e-300)
                        proposal = _project_material(coarse.chi + update)
                        edit = dict(edit="independent_global_multiscale_restart")
                    else:
                        # Exact SVD of the paid N3 material tangent; finite trial stillfull.
                        _, singular, vh = np.linalg.svd(anchor_jac[round_id], full_matrices=False)
                        count = max(2, basis.shape[1] // 3)
                        coefficients = vh[-count:].T @ rng.standard_normal(count)
                        update = basis @ coefficients
                        update *= .8 / max(float(np.max(np.abs(update))), 1e-300)
                        proposal = _project_material(coarse.chi + update)
                        edit = dict(edit="N3-anchored_material_weak_direction_with_box_projection",
                                    weakest_singular=float(singular[-1]),
                                    box_projection_defect=float(np.linalg.norm(proposal-coarse.chi-update)))
                    current = None
                    trust = .35
                    for step in range(3):
                        book.check()
                        trial = model.full_state(proposal, purpose=f"{case}:{arm}:trajectory{trajectory}:step{step}")
                        trial_norm = norm(trial)
                        improved = current is None or trial_norm < norm(current)
                        if improved:
                            current = trial
                        else:
                            trust *= .3
                        key = f"{case}__{arm}__t{trajectory}__s{step}"
                        arrays[key + "__chi"], arrays[key + "__prediction"] = np.array(trial.chi), np.array(trial.data)
                        compatible = trial_norm <= threshold
                        distinct = compatible and all(material_distance(trial.chi, item["chi"]) >= .03 for item in archive)
                        if distinct:
                            archive.append(dict(chi=np.array(trial.chi), data=np.array(trial.data), residual=trial_norm, name=key))
                        rows.append(dict(case_id=case, arm=arm, round=round_id, real_parameters=basis.shape[1],
                            trajectory=trajectory, step=step, residual=trial_norm, compatible=compatible,
                            distinct_discovery=distinct, archive_size=len(archive), nonlinear_step_accepted=improved,
                            field_backward_residual=float(trial.residual), protected_cells=int(protected.sum()),
                            protected_defect=float(np.max(np.abs(trial.chi[protected] - background[protected]))) if protected.any() else 0.,
                            edit=edit, array_prefix=key, fresh_states=1))
                        if step < 2:
                            proposal = _gn_proposal(model, current, basis, observation, white, trust=trust)
                    np.savez_compressed(output / "MATERIAL_AND_PREDICTION_ARRAYS.npz", **arrays)
                    _write(output / "SEARCH_HISTORY.json", rows)
            # Same proposals, no compensation; costsarepaid instead ofhidden replays.
            # The failed first audit spent91RHS. Drop one redundant matched
            # no-comp endpoint in each of the firsttwo cases, preserving the
            # hard9216 cumulative cap including that retained failure.
            unmatched_limit=3 if case_i<2 else 4
            for trajectory, (proposal, protected, edit) in enumerate(finite_edits[:unmatched_limit]):
                trial = model.full_state(proposal, purpose=f"{case}:no_compensation:t{trajectory}")
                key = f"{case}__same_edit_no_compensation__t{trajectory}__s0"
                arrays[key + "__chi"], arrays[key + "__prediction"] = np.array(trial.chi), np.array(trial.data)
                rows.append(dict(case_id=case, arm=ARMS[3], trajectory=trajectory, step=0,
                    residual=norm(trial), compatible=norm(trial) <= threshold, array_prefix=key,
                    protected_cells=int(protected.sum()), edit=edit, fresh_states=1))
            np.savez_compressed(output / "MATERIAL_AND_PREDICTION_ARRAYS.npz", **arrays)
            _write(output / "SEARCH_HISTORY.json", rows)
            print(json.dumps(dict(case=case, status="FULL_PHYSICS_SEARCH_COMPLETE", **book.receipt())), flush=True)
        # Search and proposal state is frozen before the clean offline process.
        # This online process never opens OFFLINE_TARGETS or receives its truth.
        np.savez_compressed(output / "MATERIAL_AND_PREDICTION_ARRAYS.npz", **arrays)
        _write(output / "SEARCH_HISTORY.json", rows)
        environment = dict(os.environ)
        source_directory = str(Path(__file__).resolve().parents[1])
        environment["PYTHONPATH"] = source_directory + os.pathsep + environment.get("PYTHONPATH", "")
        subprocess.run([sys.executable, "-m", "research_a.e2_offline", "--output", str(output),
                        "--source", str(a30_root)], env=environment, check=True)
        status = "COMPLETE"
    except BaseException as exc:
        failures.append(dict(error=type(exc).__name__, message=str(exc), rows_complete=len(rows),
                             execution_wall_seconds=time.perf_counter() - start))
        _write(output / "FAILURES.json", failures)
        np.savez_compressed(output / "MATERIAL_AND_PREDICTION_ARRAYS.npz", **arrays)
        _write(output / "SEARCH_HISTORY.json", rows)
        raise
    finally:
        receipt = book.receipt()
        book.finish(status)
        _write(output / "RESULTS.json", dict(status=status, resources=receipt, histories=len(rows),
            failures=failures, scientific_verdict="PARENT_REVIEW_REQUIRED",
            probabilities="CANDIDATE_ONLY_NULL", original_A30_NO_GO="UNCHANGED",
            repeated_campaign="REFUSE_RUN_IF_RESULTS_EXISTS"))
    return dict(status=status, resources=receipt, output=str(output))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--a30-root", type=Path, default=SOURCE)
    parser.add_argument("--e1-gate", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, a30_root=args.a30_root, e1_gate_path=args.e1_gate)))

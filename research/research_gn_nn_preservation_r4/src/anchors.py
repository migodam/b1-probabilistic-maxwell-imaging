"""Public-observation-only mature full-data TRF anchors and fixed features.

R4 Model, regularizer and material prolongation are reused without modification.
No loader here accepts a truth path, family, amplitude or generator seed.
"""
from __future__ import annotations

from contextlib import contextmanager
import json
import platform,socket,sys
from pathlib import Path
import time
import numpy as np
from scipy.optimize import least_squares
from threadpoolctl import threadpool_limits, threadpool_info
from research.research_z_alg_r4_diversity.src.physics import Model, Regularizer, material_basis, resample_material

DEFAULTS = dict(wall_seconds=45., max_lu=500, max_rhs=30000,
                lower=.02, upper=20., initial_material=.1,
                stages=[8, 12, 24], stage_fractions=[.20, .25, .55],
                ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=2000, discrepancy_stop=False)


def read_npz(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def dump(path, value):
    def safe(obj):
        if isinstance(obj, dict): return {str(k): safe(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)): return [safe(v) for v in obj]
        if isinstance(obj, np.ndarray): return safe(obj.tolist())
        if isinstance(obj, np.generic): return safe(obj.item())
        if isinstance(obj, float) and not np.isfinite(obj): return None
        return obj
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(safe(value), indent=2, allow_nan=False)+"\n", encoding="utf-8")


def prepare_reference(operators_path, reference_path):
    """Shared public m=.1 Jacobian, exact source-major interleaved packing."""
    reference_path = Path(reference_path)
    if reference_path.exists():
        return {"status": "REUSED_CACHE", "path": str(reference_path)}
    started = time.perf_counter()
    model = Model(read_npz(operators_path))
    with threadpool_limits(limits=1):
        state = model.forward(np.full(model.n, .1), label="public_reference", jac=True)
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(reference_path, J_ref=state["J"], Y_ref=state["y"], m_ref=state["m"])
    receipt = dict(status="PREPARED_PUBLIC_REFERENCE", wall_seconds=time.perf_counter()-started,
                   counts=dict(model.counts), reference_material=.1,
                   packing="source-major receiver-major real-imag interleaved",
                   shared_setup=True, truth_access="NONE")
    dump(reference_path.with_suffix(".json"), receipt)
    return receipt


class BudgetStop(RuntimeError):
    pass


class PhaseStop(RuntimeError):
    pass


class DiscrepancyStop(PhaseStop):
    pass


class CountedModel(Model):
    @contextmanager
    def operation(self, kind, rhs=0, **cost):
        if getattr(self, "budget_owner", None) is not None:
            self.budget_owner.check(int(rhs), int(cost.get("LU", 0)))
        with super().operation(kind, rhs=rhs, **cost):
            yield


class AnchorRun:
    def __init__(self, operators, observation, weakreg=0., config=None):
        self.cfg = {**DEFAULTS, **(config or {})}
        self.started = time.perf_counter()
        self.model = CountedModel(operators)
        self.model.budget_owner = self
        self.Y = np.asarray(observation["Y"], complex)
        self.noise = float(observation["noise"])
        self.discrepancy_threshold = .0011 if self.noise > 0 else 1e-6
        self.scale2 = max(float(np.vdot(self.Y, self.Y).real), 1e-300)
        if self.Y.shape != (48, 12) or self.model.n != 576:
            raise ValueError("Full R4 geometry required")
        if weakreg not in (0., 1e-5):
            raise ValueError("Only base or predesignated weakreg=1e-5 anchors")
        self.weakreg = weakreg
        self.reg = Regularizer((24, 24), smooth=weakreg, tv=weakreg, epsilon=.01)
        self.phase_wall = self.phase_lu = self.phase_rhs = np.inf
        self.events, self.trace = [], []
        self.best = None
        self.status, self.stop_reason = "COMPLETED_NUMERICAL", None
        self.endpoint = None
        self.all_unknowns_released = False

    def elapsed(self):
        return time.perf_counter()-self.started

    def check(self, rhs=0, lu=0):
        for name, current, cap in (("WALL_SECONDS", self.elapsed(), self.cfg["wall_seconds"]),
                                   ("MAX_LU", self.model.counts["LU"]+lu, self.cfg["max_lu"]),
                                   ("MAX_COMPLEX_RHS", self.model.counts["actual_complex_RHS"]+rhs, self.cfg["max_rhs"])):
            if (current >= cap if name == "WALL_SECONDS" else current > cap):
                raise BudgetStop(name)
        if self.elapsed() >= self.phase_wall: raise PhaseStop("PHASE_SECONDS")
        if self.model.counts["LU"]+lu > self.phase_lu: raise PhaseStop("PHASE_LU")
        if self.model.counts["actual_complex_RHS"]+rhs > self.phase_rhs: raise PhaseStop("PHASE_COMPLEX_RHS")

    def snapshot(self, state, coefficients, stage):
        relative = float(np.linalg.norm(state["y"]-self.Y)/np.sqrt(self.scale2))
        objective = .5*relative**2+self.reg.evaluate(state["m"])[0]
        row = dict(stage=stage, relative_residual=relative, objective=objective,
                   state_relative_residual=state["state_relative_residual"], seconds=self.elapsed())
        saved = dict(**row, m=state["m"].copy(), j=state["j"].copy(), Q=state["Q"].copy(), y=state["y"].copy(), c=np.array(coefficients).copy())
        self.trace.append(row)
        if self.best is None or objective < self.best["objective"]:
            self.best = saved
        return saved

    def stage(self, coefficients, stage, fraction):
        self.phase_wall = self.cfg["wall_seconds"]*fraction
        self.phase_lu = int(round(self.cfg["max_lu"]*fraction))
        self.phase_rhs = int(round(self.cfg["max_rhs"]*fraction))
        B = material_basis((24, 24), stage)
        cached_x, cached_state, cached_saved = None, None, None
        accepted = None
        first = len(self.trace)
        started, before = time.perf_counter(), dict(self.model.counts)

        def state_at(x):
            nonlocal cached_x, cached_state, cached_saved
            self.check()
            if cached_x is None or not np.array_equal(cached_x, x):
                cached_state = self.model.forward(B@x, label=f"anchor_stage{stage}")
                cached_x = np.array(x).copy()
                cached_saved = self.snapshot(cached_state, x, stage)
            return cached_state

        def residual(x):
            state = state_at(x)
            data = self.model.pack(state["y"]-self.Y)/np.sqrt(self.scale2)
            if not self.weakreg:
                return data
            d = self.reg.D@state["m"]
            q = np.sqrt(d*d+self.reg.epsilon**2)
            return np.concatenate((data, np.sqrt(2*self.reg.smooth/self.reg.edges)*d,
                                   np.sqrt(2*self.reg.tv/self.reg.edges)*d/np.sqrt(q+self.reg.epsilon)))

        def jacobian(x):
            nonlocal accepted
            state = state_at(x)
            # SciPy requests a new Jacobian after accepting a TRF step.
            # Saving here prevents rejected residual trials becoming endpoints.
            accepted = cached_saved
            self.endpoint = accepted
            if stage == 24:
                self.all_unknowns_released = True
            J = self.model.jacobian(state)@B/np.sqrt(self.scale2)
            if stage == 24 and self.cfg["discrepancy_stop"] and accepted["relative_residual"] <= self.discrepancy_threshold:
                raise DiscrepancyStop("PUBLIC_DISCREPANCY_THRESHOLD")
            if not self.weakreg:
                return J
            d = self.reg.D@state["m"]
            q = np.sqrt(d*d+self.reg.epsilon**2)
            DB = self.reg.D@B
            td = np.sqrt(2*self.reg.tv/self.reg.edges)*(1/np.sqrt(q+self.reg.epsilon)-d*d/(2*q*(q+self.reg.epsilon)**1.5))
            return np.vstack((J, np.sqrt(2*self.reg.smooth/self.reg.edges)*DB, td[:, None]*DB))

        status, reason, fit_info = "FINISHED", None, {}
        try:
            fit = least_squares(residual, np.clip(coefficients, self.cfg["lower"], self.cfg["upper"]),
                                jac=jacobian, method="trf", bounds=(self.cfg["lower"], self.cfg["upper"]),
                                x_scale="jac", tr_solver="exact", ftol=self.cfg["ftol"],
                                xtol=self.cfg["xtol"], gtol=self.cfg["gtol"], max_nfev=self.cfg["max_nfev"])
            # fit.x is the optimizer endpoint, and is already a computed physical state.
            if cached_x is not None and np.array_equal(fit.x, cached_x):
                accepted = cached_saved
            elif accepted is None or not np.array_equal(fit.x, accepted["c"]):
                state_at(fit.x)
                accepted = cached_saved
            self.endpoint = accepted
            fit_info = dict(optimizer_status=int(fit.status), optimizer_success=bool(fit.success),
                            message=str(fit.message), nfev=int(fit.nfev), njev=int(fit.njev), optimality=float(fit.optimality))
        except DiscrepancyStop as exc:
            status, reason = "DISCREPANCY_STOP", str(exc)
            self.status, self.stop_reason = status, reason
        except PhaseStop as exc:
            status, reason = "PHASE_BUDGET_STOP", str(exc)
        except BudgetStop as exc:
            status, reason = "GLOBAL_BUDGET_STOP", str(exc)
            self.status, self.stop_reason = "BUDGET_STOP", reason
        except np.linalg.LinAlgError as exc:
            # Preserve the accepted physical endpoint for this observed SciPy
            # gesdd failure. Do not switch LAPACK drivers, optimizers or limits.
            if "svd did not converge" not in str(exc).lower():
                raise
            status, reason = "NUMERICAL_STOP", "TRF_SVD_NONCONVERGENCE"
            self.status, self.stop_reason = status, reason
            fit_info = dict(optimizer_success=False, message=str(exc),
                            numerical_error_type=type(exc).__name__,
                            numerical_error_operation="SciPy dense exact TRF SVD",
                            fallback_or_restart_used=False)
        self.events.append(dict(stage=stage, parameter_dimension=stage*stage, status=status,
                                stop_reason=reason, seconds=time.perf_counter()-started,
                                start_counts=before, end_counts=dict(self.model.counts),
                                evaluated_points=len(self.trace)-first, accepted_endpoint_available=accepted is not None,
                                **fit_info))
        return np.array(coefficients) if accepted is None else accepted["c"].copy()

    def run(self):
        c, old_shape = np.full(64, self.cfg["initial_material"]), (8, 8)
        for stage, fraction in zip(self.cfg["stages"], np.cumsum(self.cfg["stage_fractions"])):
            if stage != old_shape[0]:
                c = resample_material(c, old_shape, (stage, stage))
            c = self.stage(c, stage, fraction)
            old_shape = (stage, stage)
            if self.stop_reason is not None:
                break
        if self.endpoint is None:
            raise RuntimeError("No accepted valid anchor endpoint")
        return self.endpoint


def generate_anchor(operators_path, public_path, output_path, reference_path, weakreg=0., config=None):
    """Only public files enter this API. Return scalar receipt; write NPZ states."""
    cold = time.perf_counter()
    operators, observation, reference = read_npz(operators_path), read_npz(public_path), read_npz(reference_path)
    if set(observation) != {"Y", "noise"}:
        raise ValueError("Observation payload must contain only Y and noise")
    io_seconds = time.perf_counter()-cold
    with threadpool_limits(limits=1):
        run = AnchorRun(operators, observation, weakreg, config)
        endpoint = run.run()
        solve_seconds = run.elapsed()
        counts_before = dict(run.model.counts)
        features_started = time.perf_counter()
        gb = reference["J_ref"].T@run.model.pack(run.Y-endpoint["y"])/run.scale2
        feature_seconds = time.perf_counter()-features_started
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, m0=endpoint["m"], j0=endpoint["j"], Q0=endpoint["Q"],
                        Y=run.Y, gb=gb, noise=np.array(run.noise),
                        best_valid_m=run.best["m"], best_valid_j=run.best["j"], best_valid_Q=run.best["Q"])
    receipt = dict(schema="gn-nn-r4-anchor-v1", execution_origin=dict(hostname=socket.gethostname(),platform=platform.platform(),machine=platform.machine(),python=sys.version.split()[0],backend="NumPy/SciPy CPU BLAS1",hardware_label="Apple Silicon M4 CPU" if sys.platform=="darwin" else "Windows i9-14900HX CPU"), status=run.status, stop_reason=run.stop_reason,
                   public_path=str(public_path), output_path=str(output_path), truth_access="NONE_PUBLIC_Y_AND_OPERATORS_ONLY",
                   weakreg=weakreg, config=run.cfg, discrepancy_threshold=run.discrepancy_threshold, solve_wall_seconds=solve_seconds,
                   feature_wall_seconds=feature_seconds, feature_counts={"LU": 0, "actual_complex_RHS": 0, "fixed_J_transpose_matvec": 1},
                   anchor_plus_feature_wall_seconds=solve_seconds+feature_seconds,
                   input_io_seconds=io_seconds, process_task_wall_seconds=time.perf_counter()-cold,
                   counts=counts_before, final_stage=endpoint["stage"], all_unknowns_released=run.all_unknowns_released,
                   endpoint_semantics="actual last-stage accepted TRF point; rejected trials excluded",
                   final_relative_residual=endpoint["relative_residual"], best_valid_relative_residual=run.best["relative_residual"],
                   state_relative_residual=endpoint["state_relative_residual"], events=run.events, trace=run.trace,
                   budget_overshoot_seconds=max(0., solve_seconds+feature_seconds-run.cfg["wall_seconds"]),
                   lu_budget_overshoot=max(0, counts_before.get("LU", 0)-run.cfg["max_lu"]),
                   rhs_budget_overshoot=max(0, counts_before.get("actual_complex_RHS", 0)-run.cfg["max_rhs"]),
                   shared_reference_setup_charged_separately=True, feature_formula="J_ref.T pack(Y-F(m0))/||Y||^2", threadpools=threadpool_info())
    dump(output_path.with_suffix(".json"), receipt)
    return receipt

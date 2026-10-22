from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from .core import (OUT,Book,Target,BudgetStop,write,save,plain,to_probit,from_probit,probit_derivative)
class StateCap(RuntimeError): pass

def probit_optimize(target, y, parent, folder, native, cap=64):
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
            if fresh >= cap:
                raise StateCap("fresh mode state attempts exhausted")
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
                               x_scale="jac", max_nfev=cap, **settings)
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

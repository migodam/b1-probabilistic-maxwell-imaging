"""Independent replay of the frozen HC-CIE disk-cell operators and objectives."""
from dataclasses import dataclass
import json
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import least_squares, root
from scipy.special import hankel1, j1


RANK_RTOL = 1e-13
BOX_UPPER = 30.0
MAX_EVALS = 180
SUCCESS_RTOL = 1e-5


def pack(z):
    return np.r_[z.real.ravel(), z.imag.ravel()]


def real_columns(columns):
    return np.concatenate([columns.real, columns.imag], axis=0)


def unpack(z, shape):
    size = np.prod(shape)
    return (z[:size] + 1j*z[size:]).reshape(shape)


def complex_array(record, key):
    return np.asarray(record[key]["real"]) + 1j*np.asarray(record[key]["imag"])


def complex_json(z):
    return {"real": z.real.tolist(), "imag": z.imag.tolist()}


def wedge(j, e):
    return j[..., 0]*e[..., 1] - j[..., 1]*e[..., 0]


def disk_operators(points, k, radius, receiver_angles, receiver_radius):
    distance = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    np.fill_diagonal(distance, 1.)
    factor = .5j*np.pi*k*radius*j1(k*radius)
    green = factor*hankel1(0, k*distance)
    np.fill_diagonal(green, .5j*np.pi*k*radius*hankel1(1, k*radius) - 1)
    receiver = receiver_radius*np.column_stack(
        [np.cos(receiver_angles), np.sin(receiver_angles)])
    distance = np.linalg.norm(receiver[:, None] - points[None, :], axis=2)
    return green, factor*hankel1(0, k*distance)


@dataclass
class Model:
    green: np.ndarray
    source: np.ndarray
    receiver: np.ndarray
    data: np.ndarray
    truth: np.ndarray
    points: np.ndarray
    radius: float

    def field(self, j):
        return self.source + self.green @ j

    def forward(self, chi):
        return np.linalg.solve(np.eye(len(chi)) - chi[:, None]*self.green,
                               chi[:, None]*self.source)

    def box_readout(self, j):
        e = self.field(j)
        denominator = np.sum(abs(e)**2, axis=1)
        numerator = np.sum(e.conj()*j, axis=1)
        raw = np.divide(numerator, denominator, out=np.zeros_like(numerator),
                        where=denominator > 0)
        return raw, np.clip(raw.real, 0., BOX_UPPER)

    def audit(self, j, truth_current):
        audit_start, audit_cpu_start = perf_counter(), process_time()
        raw, chi = self.box_readout(j)
        checked_current = self.forward(chi)
        record = dict(
            current_relative_error=float(np.linalg.norm(j-truth_current)
                                         / np.linalg.norm(truth_current)),
            raw_material=complex_json(raw),
            raw_material_relative_error=float(np.linalg.norm(raw-self.truth)
                                              / np.linalg.norm(self.truth)),
            box_material=chi.tolist(),
            box_material_relative_error=float(np.linalg.norm(chi-self.truth)
                                              / np.linalg.norm(self.truth)),
            auxiliary_data_relative=float(np.linalg.norm(self.receiver@j-self.data)
                                          / np.linalg.norm(self.data)),
            physical_data_relative=float(np.linalg.norm(self.receiver@checked_current-self.data)
                                         / np.linalg.norm(self.data)),
            raw_minors_norm=float(np.linalg.norm(wedge(j, self.field(j)))),
            real_box_closure_norm=float(np.linalg.norm(j-chi[:, None]*self.field(j))),
            min_field_row_norm=float(np.min(np.linalg.norm(self.field(j), axis=1))),
        )
        return dict(**record, audit_lu=1, audit_rhs=int(self.source.shape[1]),
                    audit_elapsed_seconds=perf_counter()-audit_start,
                    audit_cpu_seconds=process_time()-audit_cpu_start)


def frozen_model(directory, case):
    if case == "soft3":
        name = "tm_disjoint_elimination_counterexample.json"
    else:
        name = "tm_hardfiber_counterexample.json"
    source = json.loads((directory/name).read_text())
    points = np.asarray(source["positions"])
    radius = np.sqrt(source["area"]/np.pi)
    if case == "soft3":
        g, s = disk_operators(points, source["k"], radius, source["receiver_angles"], 2.)
        b = complex_array(source, "B")
        chi = np.asarray(source["chi"])
        truth_j = np.linalg.solve(np.eye(chi.size)-chi[:, None]*g, chi[:, None]*b)
        model = Model(g, b, s, s@truth_j, chi, points, radius)
    else:
        model = Model(*(complex_array(source, key) for key in ("G", "B", "S", "Y")),
                      np.asarray(source["chi"]), points, radius)
    return model, source


class CurrentProblem:
    def __init__(self, model, basis, base, include_data=False, weight=1.):
        self.model = model
        self.basis = basis
        self.base = base
        self.shape = (basis.shape[1], model.source.shape[1])
        self.include_data = include_data
        self.weight = weight
        self.complex_directions = np.array(
            [basis@e.reshape(self.shape) for e in np.eye(np.prod(self.shape))])
        self.directions = np.concatenate([self.complex_directions,
                                          1j*self.complex_directions])
        self.green_directions = model.green@self.directions
        self.port_directions = model.receiver@self.directions
        self.second_minors = self.weight*(
            wedge(self.directions[:, None], self.green_directions[None, :])
            + wedge(self.directions[None, :], self.green_directions[:, None]))

    def current(self, z):
        return self.base + self.basis@unpack(z, self.shape)

    def coords(self, j):
        return pack(self.basis.conj().T@(j-self.base))

    def residual_and_jac(self, z, form="minors", projector=None):
        model = self.model
        j = self.current(z)
        e = model.field(j)
        dirs, dfield = self.directions, self.green_directions
        if form == "minors":
            residual = wedge(j, e)
            derivatives = wedge(dirs, e) + wedge(j, dfield)
        else:
            raw, chi = model.box_readout(j)
            denominator = np.sum(abs(e)**2, axis=1)
            derivatives_chi = np.zeros_like(denominator + dirs[..., 0].real)
            interior = (raw.real > 0) & (raw.real < BOX_UPPER) & (denominator > 0)
            dnum = np.sum((dfield.conj()*j + e.conj()*dirs).real, axis=2)
            dden = 2*np.sum((e.conj()*dfield).real, axis=2)
            derivatives_chi[:, interior] = ((dnum - chi*dden)[:, interior]
                                             / denominator[interior])
            residual = j - chi[:, None]*e
            derivatives = dirs - chi[:, None]*dfield - derivatives_chi[..., None]*e
        if projector is not None:
            residual = projector@residual
            derivatives = derivatives@projector.T
        residual, derivatives = self.weight*residual.ravel(), self.weight*derivatives.reshape(len(dirs), -1)
        if self.include_data:
            residual = np.r_[(model.receiver@j-model.data).ravel(), residual]
            derivatives = np.hstack([self.port_directions.reshape(len(dirs), -1), derivatives])
        return pack(residual), real_columns(derivatives.T)

    def gradient(self, z, form="minors"):
        f, a = self.residual_and_jac(z, form)
        return a.T@f

    def analytic_hessian(self, z):
        j = self.current(z)
        r = self.weight*wedge(j, self.model.field(j))
        _, a = self.residual_and_jac(z)
        return a.T@a + (self.second_minors@r.conj()).real

    def hessian_box(self, z):
        j = self.current(z)
        e = self.model.field(j)
        raw, chi = self.model.box_readout(j)
        f, a = self.residual_and_jac(z, "real_box")
        closure = j-chi[:, None]*e
        dirs, dfield = self.directions, self.green_directions
        denominator = np.sum(abs(e)**2, axis=1)
        interior = (raw.real > 0) & (raw.real < BOX_UPPER) & (denominator > 0)
        dnum = np.sum((dfield.conj()*j + e.conj()*dirs).real, axis=2)
        dden = 2*np.sum((e.conj()*dfield).real, axis=2)
        dchi = np.zeros_like(dnum)
        dchi[:, interior] = ((dnum-chi*dden)[:, interior]/denominator[interior])
        # The interior normal equation eliminates the term containing D²chi.
        contracted = (dchi @ np.sum(dfield.conj()*closure, axis=2).real.T)
        return a.T@a - self.weight**2*(contracted+contracted.T)

    def refine_bad(self, z):
        sol = root(self.gradient, z, jac=self.analytic_hessian, method="hybr",
                   options={"xtol": 1e-11})
        return sol.x, dict(root_flag=bool(sol.success), evaluations=sol.nfev,
                           message=str(sol.message))

    def solve(self, start, form, projector=None):
        def value(z):
            return self.residual_and_jac(z, form, projector)[0]

        def jacobian(z):
            return self.residual_and_jac(z, form, projector)[1]

        t0, c0 = perf_counter(), process_time()
        sol = least_squares(value, start, jac=jacobian, method="lm", max_nfev=MAX_EVALS,
                            ftol=1e-12, xtol=1e-12, gtol=1e-12)
        elapsed, cpu = perf_counter()-t0, process_time()-c0
        return sol.x, dict(status=int(sol.status), nfev=int(sol.nfev), njev=int(sol.njev),
                           residual_norm=float(np.linalg.norm(sol.fun)),
                           gradient_norm=float(np.linalg.norm(sol.grad)),
                           elapsed_seconds=elapsed, cpu_seconds=cpu, message=sol.message)


class MaterialProblem:
    def __init__(self, model):
        self.model = model
        self.cache_x = None
        self.factorizations = self.rhs_solved = 0

    def current(self, x, direction=None):
        from scipy.linalg import lu_factor, lu_solve

        if not np.array_equal(x, self.cache_x):
            self.cache_x = x.copy()
            self.factor = lu_factor(np.eye(len(x)) - x[:, None]*self.model.green)
            self.j = lu_solve(self.factor, x[:, None]*self.model.source)
            self.factorizations += 1
            self.rhs_solved += self.model.source.shape[1]
            self.e = self.model.field(self.j)
        if direction is None:
            return self.j
        self.rhs_solved += len(direction)*self.model.source.shape[1]
        directions = direction*self.e[None, :, :]
        rhs = directions.transpose(1, 0, 2).reshape(len(x), -1)
        return lu_solve(self.factor, rhs).reshape(len(x), len(direction), -1).transpose(1, 0, 2)

    def residual(self, x):
        return pack(self.model.receiver@self.current(x) - self.model.data)

    def jac(self, x):
        basis = np.eye(len(x))[:, :, None]
        d = self.model.receiver@self.current(x, basis)
        return real_columns(d.reshape(len(x), -1).T)

    def objective_hessian(self, x, direction):
        from scipy.linalg import lu_solve

        f = self.residual(x)
        da = self.current(x, direction[None, :, None])[0]
        db = lu_solve(self.factor, 2*direction[:, None]*(self.model.green@da))
        d1 = pack(self.model.receiver@da)
        return float(d1@d1 + f@pack(self.model.receiver@db))

    def full_hessian(self, x):
        basis = np.eye(len(x))
        diagonal = np.array([self.objective_hessian(x, h) for h in basis])
        hessian = np.diag(diagonal)
        for row in range(len(x)):
            for col in range(row):
                mixed = self.objective_hessian(x, basis[row]+basis[col])
                hessian[row, col] = hessian[col, row] = (mixed-diagonal[row]-diagonal[col])/2
        return hessian

    def refine_box_stationary(self, candidate, truth_current):
        x = candidate.copy()
        lower = x <= 1e-8
        upper = x >= BOX_UPPER-1e-8
        x[lower] = 0.
        x[upper] = BOX_UPPER
        free = ~(lower | upper)

        def assembled(z):
            y = x.copy()
            y[free] = z
            return y

        def gradient(y):
            return self.jac(y).T@self.residual(y)

        if free.any():
            stationary = root(lambda z: gradient(assembled(z))[free], x[free],
                              jac=lambda z: self.full_hessian(assembled(z))[np.ix_(free, free)],
                              options={"xtol": 1e-11})
            x[free] = stationary.x
            root_flag = bool(stationary.success)
        else:
            root_flag = True
        g = gradient(x)
        h = self.full_hessian(x)
        return dict(
            material=x.tolist(), root_flag=root_flag,
            distance_from_solver_endpoint=float(np.linalg.norm(x-candidate)),
            free_gradient_norm=float(np.linalg.norm(g[free])), gradient=g.tolist(),
            interior_active_set_clear=bool(np.all((x[free]>0)&(x[free]<BOX_UPPER))),
            strict_lower_kkt=bool(np.all(g[lower]>0)),
            strict_upper_kkt=bool(np.all(g[upper]<0)),
            lower_active_cells=np.flatnonzero(lower).tolist(),
            upper_active_cells=np.flatnonzero(upper).tolist(),
            free_hessian_eigenvalues=np.linalg.eigvalsh(h[np.ix_(free, free)]).tolist(),
            **self.model.audit(self.model.forward(x), truth_current),
        )

    def finite_response(self, left, right):
        from scipy.linalg import lu_solve

        jl = self.current(left).copy()
        el = self.model.field(jl)
        jr = self.current(right)
        difference = lu_solve(self.factor, (right-left)[:, None]*el)
        return float(np.linalg.norm(self.model.receiver@(jr-jl-difference)))

    def solve(self, start, method):
        t0, c0 = perf_counter(), process_time()
        self.cache_x = None
        self.factorizations = self.rhs_solved = 0
        if method == "trf_gn":
            sol = least_squares(self.residual, start, jac=self.jac, bounds=(0., BOX_UPPER),
                                method="trf", max_nfev=MAX_EVALS,
                                ftol=1e-12, xtol=1e-12, gtol=1e-12)
            x, f, a = sol.x, sol.fun, sol.jac
            info = dict(status=int(sol.status), nfev=int(sol.nfev), njev=int(sol.njev),
                        message=sol.message)
        else:
            x, f, a, info = self.projected_lm(start)
        elapsed, cpu = perf_counter()-t0, process_time()-c0
        g = a.T@f
        return x, dict(**info, residual_norm=float(np.linalg.norm(f)),
                       gradient_norm=float(np.linalg.norm(g)),
                       projected_gradient_norm=float(np.linalg.norm(x-np.clip(x-g, 0., BOX_UPPER))),
                       lu=self.factorizations, rhs=self.rhs_solved,
                       elapsed_seconds=elapsed, cpu_seconds=cpu)

    def projected_lm(self, start):
        x = start.copy()
        f, a = self.residual(x), self.jac(x)
        nfev = njev = 1
        damping = max(1e-12, 1e-3*np.linalg.norm(a, 2)**2)
        for _ in range(MAX_EVALS-1):
            g = a.T@f
            if np.linalg.norm(x-np.clip(x-g, 0., BOX_UPPER), np.inf) <= 1e-12:
                return x, f, a, dict(status=1, nfev=nfev, njev=njev,
                                     message="projected gradient tolerance")
            free = ~(((x <= 0) & (g >= 0)) | ((x >= BOX_UPPER) & (g <= 0)))
            sub = a[:, free]
            step = np.zeros_like(x)
            step[free] = np.linalg.solve(sub.T@sub + damping*np.eye(np.sum(free)), -g[free])
            trial = np.clip(x+step, 0., BOX_UPPER)
            if np.linalg.norm(trial-x) <= 1e-13*(1+np.linalg.norm(x)):
                return x, f, a, dict(status=2, nfev=nfev, njev=njev,
                                     message="projected step tolerance")
            f_trial = self.residual(trial)
            nfev += 1
            if np.linalg.norm(f_trial) < np.linalg.norm(f):
                x, f = trial, f_trial
                a = self.jac(x)
                njev += 1
                damping = max(1e-15, damping/3)
            else:
                damping *= 10
        return x, f, a, dict(status=0, nfev=nfev, njev=njev,
                              message="evaluation limit")

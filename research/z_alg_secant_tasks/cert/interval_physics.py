"""Directed-rounding disk-pulse operators and finite-dimensional Krawczyk proofs."""

from time import perf_counter, process_time

from flint import __version__ as flint_version
from flint import acb, acb_mat, arb, arb_mat, ctx
import numpy as np

from array_secant import RADIUS, SOURCE_ANGLES, RECEIVER_ANGLES, known_host


PRECISION = 192
BOX_HALF_WIDTHS = (1e-11, 3e-11, 1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 3e-8, 1e-7, 1e-6)
BASE_PORTS = ((0, 0), (1, 0), (0, 1))  # (receiver, source)
ANGLE_PORTS = tuple((receiver, source + 2) for receiver, source in BASE_PORTS)
ctx.prec = PRECISION
ctx.threads = 1


def interval_environment():
    return {"python_flint": flint_version, "precision_bits": ctx.prec, "flint_threads": ctx.threads,
            "nominal_half_widths": BOX_HALF_WIDTHS}


def geom_decimal(x):
    return arb(str(float(x)))


def point(x):
    return arb(float(x))


def bounds(x):
    return {
        "lo": float(np.nextafter(float(x.lower()), -np.inf)),
        "hi": float(np.nextafter(float(x.upper()), np.inf)),
    }


def complex_bounds(z):
    return {"real": bounds(z.real), "imag": bounds(z.imag)}


def identity(n, matrix=acb_mat):
    return matrix([[int(i == j) for j in range(n)] for i in range(n)])


def column(vector):
    return arb_mat([[x] for x in vector])


def diag(vector):
    return acb_mat([[v if i == j else 0 for j in range(len(vector))]
                    for i, v in enumerate(vector)])


def ball_norm(vector):
    lo, hi = arb(0), arb(0)
    for x in vector:
        lower, upper = x.abs_lower(), x.abs_upper()
        lo += lower * lower
        hi += upper * upper
    return lo.sqrt().union(hi.sqrt())


def dist(first, second):
    return ball_norm(x - y for x, y in zip(first, second, strict=True))


def hankel(order, x):
    x = acb(x)
    return x.bessel_j(order) + acb(0, 1) * x.bessel_y(order)


def polar(radius, angles):
    polar_points = []
    for angle in angles:
        theta = angle * arb.pi() / 180
        polar_points.append((arb(radius) * theta.cos(), arb(radius) * theta.sin()))
    return polar_points


class ReducedFrequency:
    """An Arb enclosure of the exact decimal and angle specified pulse physics."""

    def __init__(self, geometry, n, frequency_ratio, with_angle):
        self.n = n
        centers = [[geom_decimal(v) for v in xy] for xy in geometry]
        sources = polar("0.95", SOURCE_ANGLES if with_angle else SOURCE_ANGLES[:2])
        receivers = polar("1.07", RECEIVER_ANGLES[:2])
        nc, ns = len(centers), len(sources)
        radius, pi = geom_decimal(RADIUS), arb.pi()
        k, imag = 2 * pi * arb(frequency_ratio), acb(0, 1)
        ka = k * radius
        factor = imag * pi * ka * acb(ka).bessel_j(1) / 2
        diagonal = imag * pi * ka * hankel(1, ka) / 2 - 1
        G = acb_mat(nc, nc)
        for i, x in enumerate(centers):
            G[i, i] = diagonal
            for j, y in enumerate(centers[:i]):
                d = dist(x, y)
                if not d > 2 * radius:
                    raise ArithmeticError("Disk disjointness was not proved.")
                G[i, j] = G[j, i] = factor * hankel(0, k * d)
        B = acb_mat([[imag * hankel(0, k * dist(x, s)) / 4 for s in sources] for x in centers])
        S = acb_mat([[factor * hankel(0, k * dist(r, x)) for x in centers] for r in receivers])
        if not all(dist(a, x) > radius for a in sources + receivers for x in centers):
            raise ArithmeticError("A transducer was not proved exterior to every disk.")
        fixed = [arb(0)] * n + [geom_decimal(x) for x in known_host()[1]]
        A0 = identity(nc) - diag(fixed) * G
        X = A0.solve(acb_mat([
            [fixed[i] * B[i, source] for source in range(ns)] + [int(i == j) for j in range(n)]
            for i in range(nc)
        ]))
        host = acb_mat([[X[i, source] for source in range(ns)] for i in range(nc)])
        embed = acb_mat([[X[i, ns + j] for j in range(n)] for i in range(nc)])
        Hfull, Bfull = G * embed, B + G * host
        self.H = acb_mat([[Hfull[i, j] for j in range(n)] for i in range(n)])
        self.B = acb_mat([[Bfull[i, source] for source in range(ns)] for i in range(n)])
        self.R = S * embed
        self.host = S * host
        self.solve_count = 0

    def response(self, material):
        self.solve_count += 1
        D = diag(material)
        system = identity(self.n) - D * self.H
        # A dyadic preconditioner makes the ball matrix nearly diagonal.
        # Inverting its product certifies the invertibility of both factors.
        preconditioner = system.mid().inv().mid()
        inverse = (preconditioner * system).inv() * preconditioner
        W = inverse * D * self.B
        E = self.B + self.H * W
        P = self.R * inverse
        response = self.host + self.R * W
        return response, P, E

    def channels(self, material, ports):
        response, P, E = self.response(material)
        return [response[r, s] for r, s in ports], [
            [P[r, j] * E[j, s] for j in range(self.n)] for r, s in ports
        ]


class IntervalModel:
    def __init__(self, geometry, n, with_frequency=False):
        self.f0 = ReducedFrequency(geometry, n, "1", with_angle=with_frequency)
        self.f1 = ReducedFrequency(geometry, n, "1.2", with_angle=False) if with_frequency else None

    @property
    def solve_count(self):
        return self.f0.solve_count + (self.f1.solve_count if self.f1 is not None else 0)

    def response(self, material, arm):
        ports_f0 = BASE_PORTS + ANGLE_PORTS if arm == "angle" else BASE_PORTS
        observations, jacobians = self.f0.channels(material, ports_f0)
        if arm == "frequency":
            more, more_jacobian = self.f1.channels(material, BASE_PORTS)
            power_scale = arb("1.2").sqrt()
            observations.extend(power_scale * x for x in more)
            jacobians.extend([power_scale * x for x in row] for row in more_jacobian)
        exposure = 1 / arb(2).sqrt() if arm != "base" else arb(1)
        J = arb_mat([[exposure * x.real for x in row] for row in jacobians]
                    + [[exposure * x.imag for x in row] for row in jacobians])
        return [exposure * x for x in observations], J

    def host(self, arm):
        ports = BASE_PORTS + ANGLE_PORTS if arm == "angle" else BASE_PORTS
        observations = [self.f0.host[r, s] for r, s in ports]
        if arm == "frequency":
            observations.extend(arb("1.2").sqrt() * self.f1.host[r, s] for r, s in BASE_PORTS)
        return [x / arb(2).sqrt() for x in observations] if arm != "base" else observations


def make_real(values):
    return column([x.real for x in values] + [x.imag for x in values])


def certify_chord(model, center, t, v, U, root, midpoint_jacobian, roi,
                  arm="base", single=False, task_target=None):
    """Enclose a corrected pair, or one endpoint sharing the exact reference."""
    start_wall, start_cpu = perf_counter(), process_time()
    solve_before = model.solve_count
    dimension = U.shape[1]
    n = len(center)
    nuisance = arb_mat(U.tolist())
    z = [point(x) for x in center]
    slope = [point(x) for x in v]
    center_u = [point(x) for x in root]
    cutoff = [arb("0.5"), arb(30)]
    roi_indices = np.flatnonzero(roi)
    if any(sum(int(U[i, j]) for i in roi_indices) != 0 for j in range(dimension)):
        raise ArithmeticError("The claimed chart does not preserve the declared task exactly.")
    task_slope = sum(slope[i] for i in roi_indices) / len(roi_indices)
    task_shift = point(t) * task_slope if task_target is None else arb(task_target)
    step = point(t) if task_target is None else task_shift / task_slope
    base = [step * x for x in slope]
    attempted = []

    def endpoints(u):
        change = [base[i] + sum(int(U[i, j]) * u[j] for j in range(dimension))
                  for i in range(n)]
        return [m + h for m, h in zip(z, change, strict=True)], [
            m - h for m, h in zip(z, change, strict=True)
        ]

    def evaluate(u):
        plus, minus = endpoints(u)
        yp, Jp = model.response(plus, arm)
        if single:
            return make_real([a - b for a, b in zip(yp, y0, strict=True)]), Jp * nuisance
        ym, Jm = model.response(minus, arm)
        return make_real([(p - m) / 2 for p, m in zip(yp, ym, strict=True)]), ((Jp + Jm) * nuisance) * arb("0.5")

    try:
        y0, _ = model.response(z, arm)
        if arm == "base":
            baseline_ref = y0
        else:
            # Shared noise RMS is keyed to the unweighted baseline at the same physical center.
            baseline_ref, _ = model.response(z, "base")
        baseline_signal = ball_norm(y - host for y, host in zip(baseline_ref, model.host("base"), strict=True))
        noise = arb("0.1") * baseline_signal * (arb(len(y0)) / len(baseline_ref)).sqrt()
        fmid, _ = evaluate(center_u)
        inverse = np.linalg.inv(midpoint_jacobian)
        preconditioner = arb_mat([[point(x) for x in row] for row in inverse])
        determinant = preconditioner.det()
        if determinant.contains(0):
            raise ArithmeticError("Preconditioner nonsingularity was not proved.")
        residual = preconditioner * fmid
        I = identity(dimension, arb_mat)
        for half_width in BOX_HALF_WIDTHS:
            widths = [arb(0, half_width)] * dimension
            box = [u + r for u, r in zip(center_u, widths, strict=True)]
            plus_box, minus_box = endpoints(box)
            checked = plus_box if single else plus_box + minus_box
            if not all(cutoff[0] < x < cutoff[1] for x in checked):
                attempted.append({"half_width": half_width, "reason": "physical_box_not_interior"})
                continue
            try:
                _, Jbox = evaluate(box)
            except (ValueError, ZeroDivisionError) as error:
                attempted.append({"half_width": half_width, "reason": str(error)})
                continue
            error_map = I - preconditioner * Jbox
            contraction = max(sum(error_map[i, j].abs_upper() for j in range(dimension)).upper()
                              for i in range(dimension))
            image = error_map * column(widths) - residual
            image_rows = [image[i, 0] for i in range(dimension)]
            if not contraction < arb(1) or not all(w.contains_interior(img)
                                                    for w, img in zip(widths, image_rows, strict=True)):
                attempted.append({"half_width": half_width, "reason": "strict_Krawczyk_inclusion_failed"})
                continue
            # Any fixed point in the tested box lies in its Krawczyk image.
            root_enclosure = [u + img for u, img in zip(center_u, image_rows, strict=True)]
            plus, minus = endpoints(root_enclosure)
            common, _ = model.response(plus, arm)
            # In the one-sided problem, its exact solution shares the declared reference.
            if single:
                common = y0
            delta_y = arb(0) if single else ball_norm(a - b for a, b in zip(common, y0, strict=True))
            clearances = [gap for x in (plus if single else plus + minus)
                          for gap in (x - cutoff[0], cutoff[1] - x)]
            clearance = min(x.lower() for x in clearances).union(min(x.upper() for x in clearances))
            return {
                "status": "certified_exact_root_in_chart", "arm": arm,
                "half_width": half_width, "contraction_infinity_norm": bounds(contraction),
                "Krawczyk_image_offset_inf_norm_upper": bounds(max(x.abs_upper() for x in image_rows)),
                "preconditioner_determinant": bounds(determinant),
                "chart_root_enclosure": [bounds(x) for x in root_enclosure],
                "material_plus_enclosure": [bounds(x) for x in plus],
                "material_minus_enclosure": None if single else [bounds(x) for x in minus],
                "material_box_min_clearance": bounds(clearance),
                "task_shift_plus_from_reference": bounds(task_shift),
                "task_gap_plus_minus": None if single else bounds(2 * task_shift),
                "common_complex_readout": [complex_bounds(x) for x in common],
                "reference_complex_readout": [complex_bounds(x) for x in y0],
                "readout_reference_distance": bounds(delta_y),
                "noise_10pct_radius": bounds(noise),
                "readout_reference_distance_over_10pct": bounds(delta_y / noise),
                "reference_observation_admits_endpoint_or_pair": bool(delta_y < noise),
                "failed_box_attempts": attempted,
                "interval_state_solves": model.solve_count - solve_before,
                "wall_seconds": perf_counter() - start_wall, "cpu_seconds": process_time() - start_cpu,
            }
    except (ArithmeticError, ValueError, ZeroDivisionError, np.linalg.LinAlgError) as error:
        attempted.append({"reason": str(error)})
    return {
        "status": "not_certified", "arm": arm, "failed_box_attempts": attempted,
        "interval_state_solves": model.solve_count - solve_before,
        "wall_seconds": perf_counter() - start_wall, "cpu_seconds": process_time() - start_cpu,
    }

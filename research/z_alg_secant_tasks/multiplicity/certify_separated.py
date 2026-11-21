"""Enclose roots and distinct cross ports for the exact translated Hankel model."""

import argparse
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter, process_time

from separated import HERE, FIXED, write_json
import numpy as np
from flint import acb, acb_mat, arb, arb_mat, ctx, __version__ as flint_version
from threadpoolctl import threadpool_limits


WIDTHS = (1e-11, 3e-11, 1e-10, 3e-10, 1e-9, 3e-9, 1e-8, 1e-7, 1e-6)
ctx.prec = 192
ctx.threads = 1


def binary(x):
    return arb(float(x))


def ball_range(x):
    return [float(np.nextafter(float(x.lower()), -np.inf)),
            float(np.nextafter(float(x.upper()), np.inf))]


def complex_range(z):
    return {"real": ball_range(z.real), "imag": ball_range(z.imag)}


def identity(n, mat=acb_mat):
    return mat([[int(i == j) for j in range(n)] for i in range(n)])


def diagonal(values):
    return acb_mat([[x if i == j else 0 for j in range(len(values))] for i, x in enumerate(values)])


def column(values):
    return arb_mat([[x] for x in values])


def distance(x, y):
    a, b = x[0]-y[0], x[1]-y[1]
    return (a*a+b*b).sqrt()


def ball_norm(values):
    lower, upper = arb(0), arb(0)
    for value in values:
        lo, hi = value.abs_lower(), value.abs_upper()
        lower += lo * lo
        upper += hi * hi
    return lower.sqrt().union(upper.sqrt())


def hankel(order, value):
    argument = acb(value)
    return argument.bessel_j(order) + acb(0, 1) * argument.bessel_y(order)


class ExactModel:
    """Known-host reduction built from ball-certified all-to-all coefficients."""

    def __init__(self, source, count, R):
        exact_xy = {
            name: [(arb(x) + component*R, arb(y)) for component in range(count) for x, y in local]
            for name, local in source["geometry_decimal"].items()
        }
        disks, sources, receivers = (exact_xy[name] for name in ("centers", "sources", "receivers"))
        n, ns = len(disks), len(sources)
        fixed_local = source["known_material_contrast_decimal"]
        fixed = [arb(x) for _ in range(count) for x in fixed_local]
        active_indices = [j + i*len(fixed_local) for i in range(count) for j in source["active_indices"]]
        self.n = len(active_indices)
        self.nr, self.ns = len(receivers), ns
        shots_per_component = len(source["geometry_decimal"]["sources"])
        self.local = tuple((i, i*shots_per_component+b) for i in range(self.nr)
                           for b in range(shots_per_component))
        self.cross = tuple((i, s) for i in range(self.nr) for s in range(ns) if (i, s) not in self.local)
        self.target = [acb(binary(re), binary(im)) for re, im in source["local_observation_binary64"]] * count

        radius, k, pi = arb(source["disk_radius_decimal"]), 2*arb.pi(), arb.pi()
        i = acb(0, 1)
        ka = k*radius
        cell = i*pi*ka*acb(ka).bessel_j(1)/2
        diag = i*pi*ka*hankel(1, ka)/2-1
        G = acb_mat(n, n)
        for row, x in enumerate(disks):
            G[row, row] = diag
            for col, y in enumerate(disks[:row]):
                d = distance(x, y)
                if not d > 2*radius:
                    raise ArithmeticError("The exact disks are not proved disjoint.")
                G[col, row] = G[row, col] = cell*hankel(0, k*d)
        if not all(distance(x, q) > radius for x in disks for q in sources+receivers):
            raise ArithmeticError("A transducer is not proved outside the disks.")
        B = acb_mat([[i*hankel(0, k*distance(x, s))/4 for s in sources] for x in disks])
        S = acb_mat([[cell*hankel(0, k*distance(r, x)) for x in disks] for r in receivers])
        A = identity(n)-diagonal(fixed)*G
        solved = A.solve(acb_mat([[fixed[j]*B[j, source] for source in range(ns)]
                                   + [int(j == a) for a in active_indices] for j in range(n)]))
        incident_host_current = acb_mat([[solved[j, s] for s in range(ns)] for j in range(n)])
        embedding = acb_mat([[solved[j, ns+p] for p in range(self.n)] for j in range(n)])
        GE, Ehost = G*embedding, B+G*incident_host_current
        self.G = acb_mat([[GE[j, p] for p in range(self.n)] for j in active_indices])
        self.B = acb_mat([[Ehost[j, s] for s in range(ns)] for j in active_indices])
        self.S = S*embedding
        self.host_readout = S*incident_host_current
        self.state_solve_count = 0

    def state(self, active):
        self.state_solve_count += 1
        D = diagonal(active)
        A = identity(self.n)-D*self.G
        left = A.mid().inv().mid()
        Ainv = (left*A).inv()*left
        j = Ainv*D*self.B
        observed, internal = self.S*Ainv, self.B+self.G*j
        y = self.host_readout+self.S*j
        J = [[observed[r, i]*internal[i, s] for i in range(self.n)] for r, s in self.local]
        return y, arb_mat([[x.real for x in row] for row in J]+[[x.imag for x in row] for row in J])

    def local_residual(self, y):
        f = [y[r, s]-target for (r, s), target in zip(self.local, self.target, strict=True)]
        return column([x.real for x in f]+[x.imag for x in f])


def certify_point(proof, seed, estimate, lo, hi):
    start, cpu = perf_counter(), process_time()
    before = proof.state_solve_count
    center = [binary(x) for x in estimate]
    y, J = proof.state(center)
    approx_inverse = J.inv().mid()
    det = approx_inverse.det()
    if det.contains(0):
        raise ArithmeticError("The preconditioner was not proved nonsingular.")
    dy = approx_inverse*proof.local_residual(y)
    bad = []
    for w in WIDTHS:
        error = arb(0, w)
        X = [x+error for x in center]
        if not all(lo < point < hi for point in X):
            bad.append([w, "not_inside_legal_material_box"])
            continue
        try:
            _, Jbox = proof.state(X)
        except (ValueError, ZeroDivisionError) as issue:
            bad.append([w, str(issue)])
            continue
        M = identity(proof.n, arb_mat)-approx_inverse*Jbox
        contraction = max(sum(M[i, j].abs_upper() for j in range(proof.n)).upper() for i in range(proof.n))
        K = M*column([error]*proof.n)-dy
        image = [K[i, 0] for i in range(proof.n)]
        if not contraction < 1 or not all(error.contains_interior(interval) for interval in image):
            bad.append([w, "strict_Krawczyk_inclusion_failed"])
            continue
        root = [x+interval for x, interval in zip(center, image, strict=True)]
        predicted, _ = proof.state(root)
        legal_guards = [gap for x in root for gap in (x-lo, hi-x)]
        min_clearance = min(x.lower() for x in legal_guards).union(min(x.upper() for x in legal_guards))
        return {
            "branch": seed, "status": "certified_simple_real_root",
            "Krawczyk_radius": w, "contraction_inf": ball_range(contraction),
            "material_intervals": [ball_range(x) for x in root],
            "min_legal_clearance": ball_range(min_clearance),
            "all_complex_port_intervals": [[complex_range(predicted[r, s]) for s in range(proof.ns)]
                                           for r in range(proof.nr)],
            "previous_tested_boxes": bad, "state_solve_count": proof.state_solve_count-before,
            "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
        }, root, predicted
    return {
        "branch": seed, "status": "not_certified", "previous_tested_boxes": bad,
        "state_solve_count": proof.state_solve_count-before,
        "wall_seconds": perf_counter()-start, "cpu_seconds": process_time()-cpu,
    }, None, None


def cross_proofs(model, roots):
    pairs = []
    for (a, za, ya), (b, zb, yb) in combinations(roots, 2):
        material_gap = ball_norm(x-z for x, z in zip(za, zb, strict=True))
        cross = [ya[i, j]-yb[i, j] for i, j in model.cross]
        distance = ball_norm(cross)
        pairs.append({"first": a, "second": b, "material_l2": ball_range(material_gap),
                      "cross_ports_l2": ball_range(distance),
                      "cross_separation_by_complex_channel": [ball_range(abs(c)) for c in cross],
                      "cross_ports_strictly_separate": bool(distance > 0)})
    one_channel = []
    if pairs:
        for i, (rx, source) in enumerate(model.cross):
            low = min(row["cross_separation_by_complex_channel"][i][0] for row in pairs)
            if low > 0:
                one_channel.append({"receiver_index": rx, "source_index": source,
                                    "min_pairwise_separation_lower_bound": low})
    return pairs, one_channel


def run(parameters, numeric, output):
    started, cpu = perf_counter(), process_time()
    attempts = []
    for label, numeric_dir in numeric.items():
        s, c = perf_counter(), process_time()
        metadata = json.loads((numeric_dir/"numeric_ladder.json").read_text())
        selected = next(row for row in metadata["rows"] if row["id"] == label)
        proof = ExactModel(parameters, selected["components"], selected["spacing_reference_wavelengths"])
        preparation = {"wall_seconds": perf_counter()-s, "cpu_seconds": process_time()-c,
                       "known_host_full_ball_solve": 1,
                       "known_host_full_ball_rhs": proof.ns+proof.n}
        assert len(selected["unique_branch_representatives"]) >= 2
        check, roots = [], []
        with np.load(numeric_dir/"materials_and_ports.npz", allow_pickle=False) as archive:
            for branch in selected["unique_branch_representatives"]:
                initial = archive[f"{label}__{branch}__material"]
                interval, box, fields = certify_point(proof, branch, initial, *parameters["legal_material_range"])
                check.append(interval)
                if box is not None:
                    roots.append((branch, box, fields))
                print(label, branch, interval["status"], interval.get("Krawczyk_radius"), flush=True)
        pairs, one_port = cross_proofs(proof, roots)
        attempts.append({
            "id": label, "components": selected["components"], "spacing": selected["spacing_reference_wavelengths"],
            "root_certificates": check, "certified_root_count": len(roots),
            "verified_all_certified_materials_distinct": all(pair["material_l2"][0] > 0 for pair in pairs),
            "cross_ports_receiver_and_source_zero_based": proof.cross,
            "pairwise_separation_certificates": pairs,
            "cross_ports_individually_separating_all_these_roots": one_port,
            "preparation": preparation, "state_enclosure_count": proof.state_solve_count,
        })
    write_json(output, {
        "proof_scope": "exact decimal local geometry + integer x shifts, k=2pi, true outgoing Hankel via Arb; saved exact binary64 target",
        "exact_target_complex": parameters["local_observation_binary64"],
        "root_box_proof": "strict Krawczyk inclusion and infinity norm contraction in the legal material box",
        "precision_bits": ctx.prec, "flint_version": flint_version, "flint_threads": ctx.threads,
        "wall_seconds": perf_counter()-started, "cpu_seconds": process_time()-cpu,
        "results": attempts,
    })


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--inputs", type=Path, default=FIXED)
    cli.add_argument("--output", type=Path, default=HERE/"results/krawczyk.json")
    args = cli.parse_args()
    inputs = json.loads(args.inputs.read_text())
    with threadpool_limits(limits=1, user_api="blas"):
        run(inputs, {"L2_R40000": HERE/"results/intermediate", "L3_R200000": HERE/"results"}, args.output)

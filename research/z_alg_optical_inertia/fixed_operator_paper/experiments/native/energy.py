"""Certify dictionary rank and a source-only energy threshold on native boxes."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh
from flint import arb, acb, ctx
from threadpoolctl import threadpool_limits

from source_index import (
    HERE, STUDY, CONTROL_ID, load, ProofGrid, diag, adj, ball, diagonal_times,
    large_sign, positive_dimension, block, gersh,
)


def clocks():
    return time.perf_counter(), time.process_time()


def elapsed(start):
    w, c = clocks()
    return dict(wall_seconds=w-start[0], cpu_seconds=c-start[1])


def endpoints(a):
    return [float(np.nextafter(float(a.lower()), -np.inf)),
            float(np.nextafter(float(a.upper()), np.inf))]


def physical_source_graph(ell, upper, states):
    """Exhibit a real exact constitutive coefficient admitting every source.

    If one constant is zero the current is zero. Otherwise use an interval
    of legal constants: det(I-mG) is a degree <= n polynomial, equal to one
    at zero, so finitely many of the constants can be singular.
    """
    lo, hi = float(np.max(ell)), float(np.min(upper))
    if lo <= 0 <= hi:
        return dict(proved=True, method="vacuum current graph J=0", constant_interval=[0., 0.])
    if lo < hi:
        return dict(proved=True, method="legal interval of uniform material; nonzero determinant polynomial",
                    constant_interval=[lo, hi], determinant_degree_upper=states,
                    determinant_at_zero=1)
    return dict(proved=False, method="constant interval test inconclusive", constant_interval=[lo, hi])


def source_form(g, b, ell, upper, w):
    pairs = [(arb(float(l)), arb(float(u))) for l, u in zip(ell, upper, strict=True)]
    mid = [(l+u)/2 for l, u in pairs]
    prod = [l*u for l, u in pairs]
    mg = diagonal_times(mid, g)
    positive = diagonal_times(prod, g)
    physical = acb(w)*(diag([1]*len(ell))-mg-adj(mg)+adj(g)*positive)
    nv = acb(w)*(diagonal_times(mid, b)-adj(g)*diagonal_times(prod, b))
    cb = acb(w)*adj(b)*diagonal_times(prod, b)
    return cb-adj(nv)*physical.solve(nv), physical


def arb_frame_signs(matrix):
    """Count signs using a high-precision frame and verify its restriction."""
    values, vectors = matrix.eig(right=True, algorithm='approx')
    counts = []
    for polarity in (1, -1):
        indices = [i for i, value in enumerate(values) if polarity*value.real > 0]
        indices.sort(key=lambda i: polarity*values[i].real, reverse=True)
        if not indices:
            counts.append(0)
            continue
        frame = block(vectors, range(matrix.nrows()), indices)
        projected = polarity*(adj(frame)*matrix*frame)
        while indices and not gersh(projected, upper=False, n=len(indices)) > 0:
            indices.pop()
        counts.append(len(indices))
    return counts


class Pencil:
    def __init__(self, source, power):
        self.source = source
        self.power = power
        self.n = source.nrows()
        self.tested = {}
        self.refined = set()

    def signs(self, threshold, refine=False):
        if threshold in self.tested and (not refine or threshold in self.refined):
            return self.tested[threshold]
        matrix = self.source if threshold == 0 else self.source-acb(threshold)*self.power
        measured = self.tested.get(threshold)
        if measured is None:
            measured = positive_dimension(matrix)
        if refine and sum(measured) < self.n:
            refined = arb_frame_signs(matrix)
            measured = [max(x, y) for x, y in zip(measured, refined, strict=True)]
        if sum(measured) > self.n:
            raise ArithmeticError('Verified positive and negative counts overlap.')
        self.tested[threshold] = measured
        if refine:
            self.refined.add(threshold)
        return measured

    def eigenvalue_intervals(self, nroots):
        # The generalized midpoint calculation supplies proposals only.
        s, g = (np.array(x.tolist(), dtype=complex) for x in (self.source, self.power))
        approx64 = eigh((s+s.conj().T)/2, (g+g.conj().T)/2, eigvals_only=True)[::-1]
        roots = self.power.solve(self.source).eig(algorithm='approx')
        approximate = sorted((float(x.real) for x in roots), reverse=True)
        brackets = []
        for i in range(nroots):
            estimate = float(approximate[i])
            lower, upper = None, None
            if estimate > 0.:
                refine = abs(estimate) < abs(approximate[0])*1e-9
                for frac in (1e-8, 1e-6, 1e-4, 1e-2, .1, .5):
                    lo, hi = estimate*(1-frac), estimate*(1+frac)
                    if lower is None and lo > 0 and self.signs(lo, refine=refine)[0] >= i+1:
                        lower = lo
                    if upper is None and self.signs(hi, refine=refine)[1] >= self.n-i:
                        upper = hi
                    if lower is not None and upper is not None:
                        break
            if upper is None:
                # A positive shifted sign bound remains valid for a mode
                # whose sign at zero cannot be resolved.
                proposed = max(abs(estimate)*2, abs(approximate[0])*2.**-120, 1e-45)
                for _ in range(50):
                    if self.signs(proposed, refine=True)[1] >= self.n-i:
                        upper = proposed
                        break
                    proposed *= 8
            brackets.append(dict(index=i+1, lower=lower, upper=upper, proposal=estimate,
                                 binary64_proposal=float(approx64[i])))
        return brackets


def threshold_bound(q, p_low, p_high, brackets, radius):
    """Infimum energy for q strictly positive modes after subtracting r^2 I."""
    if p_high < q:
        return dict(status="IMPOSSIBLE_WITH_FINITE_ENERGY", interval=None,
                    reason="strict upper bound on dictionary-positive directions is below q")
    if not q:
        return dict(status="CERTIFIED_FINITE_THRESHOLD", interval=[0., 0.])
    floor = None
    if all(b["upper"] is not None for b in brackets):
        floor = endpoints(arb(radius)**2*sum((1/arb(b["upper"]) for b in brackets), arb(0)))[0]
    if p_low >= q and all(b["lower"] is not None for b in brackets):
        ceiling = endpoints(arb(radius)**2*sum((1/arb(b["lower"]) for b in brackets), arb(0)))[1]
        return dict(status="CERTIFIED_FINITE_THRESHOLD", interval=[floor, ceiling])
    return dict(status="ENERGY_NECESSITY_ONLY" if floor is not None else "UNRESOLVED",
                interval=[floor, None] if floor is not None else None,
                reason="a strict positive lower endpoint is missing for at least one required mode")


def store(result, path):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def main(path):
    start = clocks()
    spec = json.loads((HERE/'contract.json').read_text())
    disk = load(STUDY/'ip_guided_v4/boxes/results', [CONTROL_ID])[CONTROL_ID]
    with np.load(STUDY/'ip_final_v5/strong_control/strong_ring_arrays.npz', allow_pickle=False) as rings:
        ids, profile, ring_radius = rings['ids'], rings['c_solid'], float(rings['radius'][0])
    published = json.loads((HERE/'source_index.json').read_text())['rows']
    _, (g, _, b) = ProofGrid(disk.geometry).restricted(ids)
    # Fluence over all 9216 cells, in units of beam 0 fluence.
    bd = ball(disk.geometry.B)
    gpower = adj(bd)*bd
    gpower = gpower/acb(gpower[0, 0].real)
    pd = positive_dimension(gpower)
    if pd[0] != b.ncols():
        raise ArithmeticError("Incident fluence metric is unresolved.")
    out = dict(scope="source-only pencil; old-error upper contract r^2 I; arbitrary dictionary coefficients",
               operator_contract="EXACT_STORED_FLOAT of padded spectrum and source/receiver arrays",
               precision_bits=ctx.prec,
               arb_threads=ctx.threads, blas_threads=2, dictionary_size=b.ncols(),
               metric="incident fluence summed over the full 2 m square, in units of original beam 0",
               metric_verified_positive=pd[0], common_setup=elapsed(start), rows=[],
               total=elapsed(start), finished=False)
    store(out, path)
    for name, width, mu in spec['cases']:
        t = clocks()
        if name.startswith('solid'):
            ell, up = np.maximum(0, profile-width), profile+width
        else:
            ell, up = np.zeros(len(ids)), np.full(len(ids), float(width))
        row = dict(box=name, half_width_or_cap=width, weight=mu)
        try:
            source, physical = source_form(g, b, ell, up, mu)
            saved_counts = {x['physical_negative'] for x in published if x['box'] == name}
            if len(saved_counts) != 1:
                q = large_sign(physical, np.array(physical.tolist(), dtype=complex))
                provenance = "direct full physical inertia recheck"
            else:
                q = saved_counts.pop()
                provenance = "strict full physical inertia, same stored box and operator, source_index.json"
            del physical
            graph = physical_source_graph(ell, up, len(ids))
            pencil = Pencil(source, gpower)
            pos64, neg64 = pencil.signs(0.)
            pos, neg = pencil.signs(0., refine=True)
            hi = min(q, pencil.n-neg) if graph['proved'] else pencil.n-neg
            roots = pencil.eigenvalue_intervals(q) if hi >= q else []
            shifted_lower = max((k for x, (k, _) in pencil.tested.items() if x > 0), default=0)
            lo = max(pos, shifted_lower)
            if hi < lo:
                raise ArithmeticError("Source graph and verified spectral counts are inconsistent.")
            energy = {obj: dict(radius=radius, **threshold_bound(q, lo, hi, roots, radius))
                      for obj, radius in (('disk', disk.radius), ('hollow ring', ring_radius))}
            row.update(negative_material_modes=q, physical_mode_provenance=provenance,
                       physical_source_graph=graph, positive_dictionary_modes_interval=[lo, hi],
                       positive_dictionary_modes=lo if lo == hi else None,
                       dictionary_deficit_interval=[q-hi, q-lo], dictionary_deficit=q-lo if lo == hi else None,
                       binary64_frame_positive_at_least=pos64, binary64_frame_negative_at_least=neg64,
                       direct_unshifted_positive_at_least=pos, direct_unshifted_negative_at_least=neg,
                       unshifted_signs_unresolved=pencil.n-pos-neg,
                       positive_generalized_roots=roots, energy=energy,
                       count_probes=[dict(threshold=x, positive_at_least=p, negative_at_least=n,
                                          arb_precision_frame_tested=x in pencil.refined)
                                     for x, (p, n) in sorted(pencil.tested.items())],
                       energy_threshold=energy['disk']['interval'],
                       status="CERTIFIED_RANK" if lo == hi else "RANK_INTERVAL")
        except (ArithmeticError, ZeroDivisionError, ValueError, np.linalg.LinAlgError) as error:
            row.update(status="UNRESOLVED", error_type=type(error).__name__, reason=str(error))
        row['seconds'] = elapsed(t)
        out['rows'].append(row)
        out['total'] = elapsed(start)
        store(out, path)
        print(name, row.get('positive_dictionary_modes_interval'), row.get('energy'), row['seconds'], flush=True)
    out['finished'] = True
    out['total'] = elapsed(start)
    store(out, path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, default=192)
    parser.add_argument('--output', type=Path, default=HERE/'energy.json')
    args = parser.parse_args()
    ctx.prec, ctx.threads = args.precision, 3
    with threadpool_limits(limits=2):
        main(args.output)

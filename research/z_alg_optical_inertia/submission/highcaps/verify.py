#!/usr/bin/env python3
"""Prove at least ten negative directions in a local multiplier box."""
import itertools
import json
import sys
from time import perf_counter, process_time

import numpy as np
from flint import acb, acb_mat, arb, ctx
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from screen import APERTURE, FINE, HERE, N, SIDE, STUDY, TARGET_BOXES, blocks, hermitian

sys.path.insert(0, str(STUDY))
from core.certificate import ball, gershgorin, safe_float, star

PRECISION = 160


def arb_fft2_apply(input_vector, spectrum):
    """Apply the saved DFT multiplier with exact-input complex ball arithmetic."""
    p = len(spectrum)
    zero = acb(0)
    horizontal = [
        acb.dft([input_vector[y*SIDE+x] for x in range(SIDE)] + [zero]*(p-SIDE))
        for y in range(SIDE)
    ]
    column_ifft = []
    for x in range(p):
        column = acb.dft([horizontal[y][x] for y in range(SIDE)] + [zero]*(p-SIDE))
        column_ifft.append(acb.dft([column[y]*spectrum[y][x] for y in range(p)], inverse=True))
    return [value for y in range(SIDE)
            for value in acb.dft([column_ifft[x][y] for x in range(p)], inverse=True)[:SIDE]]


def block_rows(array, ids):
    return acb_mat(len(ids), array.ncols(), [array[int(i), j] for i in ids for j in range(array.ncols())])


def exact_forms(proposal, gallery, spec):
    started = perf_counter()
    s, v = ball(gallery), ball(proposal)
    sh = star(s)
    v = v - sh * (s*sh).solve(s*v)
    projection_residual = sum((abs(e).upper()**2 for e in (s*v).entries()), arb(0)).upper().sqrt()
    dtype = [[acb(float(z.real), float(z.imag)) for z in row] for row in spec]
    columns = [arb_fft2_apply([v[i, j] for i in range(N)], dtype) for j in range(v.ncols())]
    gv = acb_mat(N, v.ncols(), [col[i] for i in range(N) for col in columns])
    diagonal_forms = []
    sector_forms = []
    optical_forms = []
    for ids in blocks():
        vv = block_rows(v, ids)
        c = star(vv) * block_rows(gv, ids)
        diagonal_forms.append(star(vv) * vv)
        sector_forms.append((c + star(c))/2)
        optical_forms.append((c - star(c))/acb(0, 2))
    gram = sum(diagonal_forms, acb_mat(v.ncols(), v.ncols()))
    gram_lo = gershgorin(gram, upper=False)
    if not gram_lo > 0:
        raise ArithmeticError('Projected current vectors failed the rank check.')
    return (diagonal_forms, sector_forms, optical_forms), {
        'gram_floor_lo': safe_float(gram_lo, upper=False),
        'receive_projection_residual_hi': safe_float(projection_residual, upper=True),
        'forms_wall_seconds': perf_counter()-started,
    }


def strict_ceiling(form):
    midpoint = np.array([[complex(item) for item in row] for row in form.tolist()])
    vectors = ball(eigh(hermitian(midpoint), driver='evd', check_finite=False)[1])
    gram = star(vectors)*vectors
    if not gershgorin(gram, upper=False) > 0:
        raise ArithmeticError('Eigenvector congruence rank unresolved.')
    transformed_upper = gershgorin(star(vectors) * form * vectors, upper=True)
    denom = gershgorin(gram, upper=bool(transformed_upper < 0))
    return (transformed_upper / denom).upper()


def certify_vertices(forms, case):
    started = perf_counter()
    w_radius, t_radius, t_lo, t_hi = case
    signs = tuple(itertools.product((-1, 1), repeat=4))
    theta = [float(t_lo)] if t_hi == t_lo else [float(t_lo), float(t_hi)]
    best_worst = None
    worst_vertex = None
    checked = 0
    a, h, o = forms
    w_radius = arb(str(w_radius))
    t_radius = arb(str(t_radius))
    for u, t in itertools.product((1., 1.5), theta):
        sector = [aa - acb(u)*hh for aa, hh in zip(a, h, strict=True)]
        for ws, ts in itertools.product(signs, signs):
            form = acb_mat(10, 10)
            for i in range(4):
                form += acb(1 + ws[i]*w_radius)*sector[i] + acb(arb(t) + ts[i]*t_radius)*o[i]
            upper = strict_ceiling(form)
            if best_worst is None or upper > best_worst:
                best_worst = upper
                worst_vertex = {'u': u, 't': t, 'w_signs': ws, 'optical_signs': ts}
            checked += 1
    return dict(box=case, verified=bool(best_worst < 0),
                max_projected_eigenvalue_upper=safe_float(best_worst, upper=True),
                worst_corner=worst_vertex, corners_checked=checked, vertex_wall_seconds=perf_counter()-started)


def main():
    started = perf_counter(), process_time()
    ctx.prec, ctx.threads = PRECISION, 1
    with np.load(FINE, allow_pickle=False) as data:
        spec, receive = data['spectrum'], data[APERTURE+'__S']
    candidates = json.loads((HERE/'screen.json').read_text())['candidates']
    receipt = dict(aperture=APERTURE, state_cells=N, source_count=10, cap_range=[1., 1.5],
                   exact_operator='saved binary64 192^2 FFT multiplier and 20 native receiver rows',
                   precision_bits=PRECISION, method='exact kernel projection and 10D projected full-wave matrices',
                   boxes=[])
    for cfg in TARGET_BOXES:
        proposal = candidates[' '.join(map(str, cfg))]
        center = proposal['center']
        with np.load(HERE/'subspace_proposals.npz', allow_pickle=False) as data:
            v = data[f't_{center:g}']
        forms, form_cost = exact_forms(v, receive, spec)
        result = certify_vertices(forms, cfg)
        result.update(proposal_center=center, screen_maximum=proposal['score'], **form_cost)
        receipt['boxes'].append(result)
        receipt.update(wall_seconds=perf_counter()-started[0], cpu_seconds=process_time()-started[1])
        (HERE/'certificates.json').write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
        print(result, flush=True)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()

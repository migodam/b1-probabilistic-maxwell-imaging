"""Bracket the fixed disk hypothesis's physical cutoff and ordinary Neumann cutoff."""
import argparse
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx
from scipy.linalg import cholesky, eigh, solve
from scipy.sparse.linalg import svds
from threadpoolctl import threadpool_limits

from run import CONTROL, CONTROL_ID, HERE, clocks, coordinates, herm, load, read, star, timing, write
from verify import ProofGrid, ball, diagonal_times, frob, frob2, hi, lo, tr


def fraction_enclosure(v):
    x = float(v)
    midpoint = Fraction.from_float(x)
    return [float(np.nextafter(x, -np.inf)) if midpoint > v else x,
            float(np.nextafter(x, np.inf)) if midpoint < v else x]


def freeze_family(obs, result_dir, strict):
    start = clocks()
    record = strict['parameters']
    support_name = obs.id + '__rank_tail_disk__active'
    with np.load(result_dir / 'followup_disk_witnesses.npz', allow_pickle=False) as f:
        ids = f[support_name]
        l, u = (f[obs.id + '__rank_tail_disk__' + name] for name in ('lower', 'upper'))
    xy = coordinates(obs.geometry.side)
    radial_sq = np.sum((xy[ids] - record['center'])**2, axis=1)
    c = record['peak'] * (1. - .25 * radial_sq / record['radius_fraction']**2)
    with np.load(CONTROL / (CONTROL_ID + '.npz'), allow_pickle=False) as f:
        m = f['material'].ravel()
    if np.setdiff1d(np.flatnonzero(m), ids).size:
        raise ArithmeticError('The archived low disk has nonzero cells outside the hypothesis support.')
    if not strict['same_support_low_physical_witness']['agrees_with_own_data_ball']:
        raise ArithmeticError('The material must already have passed the exact original noise ball.')
    c_exact = [Fraction.from_float(float(x)) for x in c]
    m_exact = [Fraction.from_float(float(x)) for x in m[ids]]
    l_exact = [Fraction.from_float(float(x)) for x in l]
    u_exact = [Fraction.from_float(float(x)) for x in u]
    deviations = [abs(m_i-c_i) for m_i, c_i in zip(m_exact, c_exact, strict=True)]
    k = int(np.argmax(deviations))
    width_witness = deviations[k]
    if any(v < 0 for v in l_exact) or not all(b <= x <= d for x, b, d in zip(c_exact, l_exact, u_exact, strict=True)):
        raise ArithmeticError('The saved high box must include its nominal center.')
    # Binary rounded endpoints in the previous certificate differ slightly from
    # c +/- 2.5. Use the largest ideal family member provably inside that box.
    inclusion_width = min(
        min(d-x for d, x in zip(u_exact, c_exact, strict=True)),
        min(x-b for x, b in zip(c_exact, l_exact, strict=True) if b > 0),
    )
    return ids, c, inclusion_width, dict(
        hypothesis_center='binary64 per-cell profile from freeze, treated as exact',
        ideal_family='fixed vacuum outside active cells; max(0,c_i-w) <= m_i <= c_i+w on support',
        support_cells=len(ids), width_witness_entry=fraction_enclosure(width_witness),
        width_witness_entry_exact=f'{width_witness.numerator}/{width_witness.denominator}',
        exclusion_width_safe_from_saved_endpoints=fraction_enclosure(inclusion_width),
        exclusion_width_safe_exact=f'{inclusion_width.numerator}/{inclusion_width.denominator}',
        limiting_linear_pixel=int(ids[k]), limiting_row_column=[int(x) for x in np.unravel_index(int(ids[k]), (96, 96))],
        limiting_nominal_contrast=float(c[k]), limiting_physical_contrast=float(m[ids[k]]),
        witnessed_peak=float(np.max(m)),
        physical_first_fit_width_bracket=[fraction_enclosure(inclusion_width)[0],
                                          fraction_enclosure(width_witness)[1]],
        witness_residual_over_radius_upper=strict['same_support_low_physical_witness']['residual_over_radius_upper'],
        exclusion_residual_over_radius_lower=strict['distance_over_radius_lower'],
        **timing(start))


def verify_cutoffs(obs, ids, c, inclusion_width, width_receipt):
    start = clocks()
    grid = ProofGrid(obs.geometry)
    (G_np, _, _), (G, _, _) = grid.restricted(ids)
    n = len(ids)
    I_np = np.eye(n)
    I = ball(I_np)
    m = [arb(float(x)) for x in c]
    A = I - diagonal_times(m, G)
    A_np = I_np - c[:, None] * G_np
    T_np = solve(A_np.T, G_np.T, check_finite=False).T
    _, sig, vh = svds(T_np, k=1, which='LM', tol=1.e-11, rng=117)
    v = ball(solve(A_np, vh[0].conj(), check_finite=False)[:, None])
    norm_lower = frob(G * v).lower() / frob(A * v).upper()
    # A threshold close to the float proposal is verified by positive
    # definiteness of q^2 A^*A - G^*G, with a residual enclosure.
    begun = clocks()
    q0 = float(sig[0])
    upper_trial = q0 * (1 + 1.e-7)
    metric_np = herm(upper_trial**2 * (star(A_np) @ A_np) - star(G_np) @ G_np)
    margin_numeric = float(eigh(metric_np, subset_by_index=(0, 0), check_finite=False)[0][0])
    shift = .45 * margin_numeric
    F = ball(cholesky(metric_np - shift * I_np, lower=True, check_finite=False))
    metric = acb(arb(upper_trial) * arb(upper_trial)) * (tr(A) * A) - tr(G) * G
    defect = frob(metric - acb(shift) * I - F*tr(F)).upper()
    positive_margin = (arb(shift) - defect).lower()
    if not positive_margin > 0:
        raise ArithmeticError('Positive projection criterion for norm upper is unresolved.')
    width = arb(inclusion_width.numerator) / inclusion_width.denominator
    cutoff_upper, cutoff_lower = 1/norm_lower, 1/arb(upper_trial)
    norm_check = dict(unit_width_T_norm_lower=lo(norm_lower), unit_width_T_norm_upper=upper_trial,
                      sufficient_2norm_half_width_bracket=[lo(cutoff_lower), hi(cutoff_upper)],
                      metric_positive_gap_lower=lo(positive_margin), metric_cholesky_defect_upper=hi(defect),
                      excluded_width_over_norm_cutoff_lower=lo(width * norm_lower),
                      excluded_width_over_norm_cutoff_upper=hi(width * arb(upper_trial)),
                      witness_width_over_norm_cutoff_lower=lo(
                          arb(width_receipt['width_witness_entry'][0]) * norm_lower),
                      witness_width_over_norm_cutoff_upper=hi(
                          arb(width_receipt['width_witness_entry'][1]) * arb(upper_trial)),
                      full_Arb_gram_products=3, **timing(begun))
    print('2-norm cutoff', norm_check['sufficient_2norm_half_width_bracket'],
          'exclusion factor', norm_check['excluded_width_over_norm_cutoff_lower'], flush=True)
    # At this width choose every material cell at its upper endpoint. The
    # perturbation then has spectrum of w G A^-1. A certified solve and
    # |tr(T^4)| <= n rho(T)^4 prove actual divergence on this candidate.
    begun = clocks()
    T = acb(width) * A.transpose().solve(G.transpose()).transpose()
    T2 = T * T
    T4 = T2 * T2
    trace = sum((T4[i, i] for i in range(n)), acb(0))
    base = (abs(trace).lower()/n).lower()
    rho_lower = (base.log()/4).exp().lower()
    if not rho_lower > 1:
        raise ArithmeticError('Fourth spectral moment did not establish divergence.')
    spectral = dict(powers=4, exact_T_source_columns=n, full_Arb_products=2,
                    ideal_exclusion_width_operator_rho_lower=lo(rho_lower),
                    ideal_exclusion_width_operator_numeric_rho=float(np.max(abs(np.linalg.eigvals(T_np)))) *
                    fraction_enclosure(inclusion_width)[0],
                    abs_trace_divided_by_n_lower=lo(base),
                    unit_width_T_rho_lower=lo(rho_lower/width), **timing(begun))
    return dict(norm_threshold=norm_check, spectral_check=spectral, **timing(start))


def fixed_optical_obstruction(obs, ids, c, inclusion_width):
    """Recheck a dark witness after adding one ordinary global optical scalar."""
    start = clocks()
    prior = ctx.prec
    ctx.prec = max(256, prior)
    try:
        grid = ProofGrid(obs.geometry)
        (G_np, S_np, _), (G, S, _) = grid.restricted(ids)
        n = len(ids)
        A_np = np.eye(n) - c[:, None] * G_np
        A = ball(np.eye(n)) - diagonal_times([arb(float(x)) for x in c], G)
        width = arb(inclusion_width.numerator) / inclusion_width.denominator
        _, _, vh = np.linalg.svd(S_np, full_matrices=True)
        Z = star(vh[S_np.shape[0]:])
        AZ, GZ = A_np @ Z, float(inclusion_width) * (G_np @ Z)
        restricted = herm(star(AZ) @ AZ - star(GZ) @ GZ)
        _, trial = eigh(restricted, subset_by_index=(0, 0), check_finite=False)
        proposal = ball(Z @ trial)
        St = tr(S)
        v = proposal - St * (S * St).solve(S * proposal)
        gv = G * v
        sector = (frob2(A * v) - width * width * frob2(gv)).upper()
        optical = (tr(v) * gv)[0, 0].imag
        optical_abs = abs(optical).upper()
        if not sector < 0 or not optical_abs > 0:
            raise ArithmeticError('The fixed optical comparison was not resolved.')
        cutoff = (-sector / optical_abs).lower()
        return dict(status='negative_receiver_null_direction',
                    tested_state='S^*S + mu C + tau Im(G) with W=I',
                    compatible_normalization='mu>0 and |tau/mu| strictly below the reported safe cutoff',
                    dark_vector_norm_lower=lo(frob(v)), sector_quadratic_upper=hi(sector),
                    optical_quadratic_interval=[lo(optical), hi(optical)],
                    optical_quadratic_absolute_upper=hi(optical_abs),
                    guaranteed_absolute_optical_over_sector_multiplier_cutoff=lo(cutoff),
                    **timing(start))
    finally:
        ctx.prec = prior


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--threads', type=int, default=3)
    parser.add_argument('--optical-only', action='store_true',
                        help='append the scalar optical comparison to an existing disk-radius result')
    args = parser.parse_args()
    ctx.prec, ctx.threads = 128, args.threads
    start = clocks()
    obs = load(args.out, [CONTROL_ID])[CONTROL_ID]
    strict = read(args.out / 'followup_disk.json')['cases'][CONTROL_ID]['strict_checks']['rank_tail_disk']
    ids, c, width, width_receipt = freeze_family(obs, args.out, strict)
    result = dict(observation=CONTROL_ID, candidate=strict['parameters'], width_bracket=width_receipt,
                  arb_precision=ctx.prec, arb_threads=ctx.threads)
    if args.optical_only:
        result = read(args.out / 'disk_radius.json')
    else:
        print('physical cutoff bracket', width_receipt['physical_first_fit_width_bracket'], flush=True)
        write(args.out / 'disk_radius.json', dict(result, status='exact_width_bracket_verified'))
        result.update(verify_cutoffs(obs, ids, c, width, width_receipt))
        result['overall_cost'] = timing(start)
    result['fixed_global_optical_comparison'] = fixed_optical_obstruction(obs, ids, c, width)
    result['status'] = 'complete'
    result['scope'] = ('The first-fit infimum is bracketed. The best attainable quadratic lower endpoint and '
                       'the best unrestricted weight or optical formulation were not computed. '
                       'The all-positive-scalar obstruction applies to the recorded identity sector alone.')
    write(args.out / 'disk_radius.json', result)
    print('spectral lower', result['spectral_check']['ideal_exclusion_width_operator_rho_lower'], flush=True)


if __name__ == '__main__':
    with threadpool_limits(1):
        main()

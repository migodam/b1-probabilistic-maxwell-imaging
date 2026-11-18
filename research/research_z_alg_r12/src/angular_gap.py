"""Graf transfer and equal-power coherent illumination of archived feasible material pairs."""
from argparse import ArgumentParser
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

for name in ('OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from scipy.linalg import eigh, lu_factor, lu_solve
from scipy.spatial.distance import cdist
from scipy.special import hankel1, j0, jv

from weak_r4 import DEFAULT_ARCHIVE, HERE, PARENT_IDS, R4, read_npz

RINGS = (1.60, .72)
SIZES = (16, 24)
FEED_CAPS = (2., 4.)
DENSE_ANGLES = 256
ORDER = 64


def cost(start, cpu, **counts):
    return dict(wall_seconds=perf_counter()-start, cpu_seconds=process_time()-cpu, **counts)


def norm(value):
    return float(np.linalg.norm(value))


def relative(value, reference):
    return norm(value-reference)/max(norm(reference), 1e-300)


def modal_rms(array):
    return norm(array)/np.sqrt(array.shape[1])


def angle_grid(count):
    return 2*np.pi*np.arange(count)/count


def angles_to_xy(radius, angles):
    return radius*np.column_stack((np.cos(angles), np.sin(angles)))


def incident(centers, k, positions):
    return .25j*hankel1(0, k*cdist(centers, positions))


def receiver_map(ops, m):
    """L(m) = S [I - diag(m)G]⁻¹ diag(m), computed by receiver-transpose solves."""
    factor = lu_factor(np.eye(len(m))-m[:, None]*ops['G'], check_finite=False)
    return lu_solve(factor, ops['S'].T, trans=1, check_finite=False).T*m


def modes(centers, k):
    order = np.arange(-ORDER, ORDER+1)
    radius = norm_rows(centers)
    phi = np.arctan2(centers[:, 1], centers[:, 0])
    inside = (.25j*jv(order[None], k*radius[:, None]))*np.exp(1j*phi[:, None]*order)
    return order, inside, {R:hankel1(order, k*R) for R in RINGS}


def norm_rows(points):
    return np.linalg.norm(points, axis=1)


def synthesize(coefficients, order, angles):
    return coefficients @ np.exp(-1j*order[:, None]*angles[None])


def fourier(data, order):
    return np.fft.ifft(data, axis=1)[:, order % data.shape[1]]


def mode_groups(order, coeff, delta, count=8):
    groups = []
    denom = norm(coeff)**2
    for n in range(np.max(np.abs(order)) + 1):
        band = coeff[:, np.abs(order)==n]
        groups.append(dict(abs_order=n, fraction=norm(band)**2/denom,
                           rms_gap_over_noise=norm(band)/delta))
    return sorted(groups, key=lambda row:-row['fraction'])[:count]


def layouts(centers, k):
    output = {}
    for R in RINGS:
        for n in SIZES:
            angles = angle_grid(n)
            xy = angles_to_xy(R, angles)
            E0 = incident(centers, k, xy)
            gram = j0(k*cdist(xy, xy))
            eigen = np.linalg.eigvalsh(gram)
            output[(R,n)] = dict(angles=angles, xy=xy, E0=E0, gram=gram, eigen=eigen)
    return output


def hermitian(value):
    return .5*(value+value.conj().T)


def feed_menu(delta, gram, rng):
    """Finite regularized power-eigenvector menu plus single and random feeds."""
    n = gram.shape[0]
    directions = [(f'single:{i}', d) for i, d in enumerate(np.eye(n))]
    objective = hermitian(delta.conj().T @ delta)
    for ridge in np.r_[0., np.geomspace(1e-9, 1e2, 56)]:
        _, vectors = eigh(objective, gram + ridge*np.eye(n), subset_by_index=(n-3, n-1), check_finite=False)
        directions.extend((f'ridge:{ridge:.8g}:eig:{j+1}', w) for j, w in enumerate(vectors[:, ::-1].T))
    random = rng.standard_normal((48,n)) + 1j*rng.standard_normal((48,n))
    directions.extend((f'random:{i}', w) for i, w in enumerate(random))
    feeds = []
    for label, vector in directions:
        energy = float(np.vdot(vector, gram@vector).real)
        if energy <= 0:
            continue
        w = vector/np.sqrt(energy)
        feeds.append((norm(delta@w), norm(w), label, w))
    return feeds


def modal_wave(coeff, order, angles, feed):
    return coeff * (np.exp(-1j*order[:, None]*angles[None]) @ feed)[None]


def describe_wave(coeff, order, layout, feed, gap, delta_noise):
    modal = modal_wave(coeff, order, layout['angles'], feed)
    high = np.sum(modal[:, np.abs(order) >= 6], axis=1)
    low = np.sum(modal[:, np.abs(order) <= 5], axis=1)
    incident_norm = norm(layout['E0']@feed)
    power = float(np.vdot(feed, layout['gram']@feed).real)
    return dict(gap_over_noise=norm(gap)/delta_noise, incident_l2=incident_norm,
                feed_l2=norm(feed), emitted_power=power,
                norm_gap_per_incident_field=norm(gap)/incident_norm,
                modal_relative_closure=relative(low+high, gap), low_coherent_gap_over_noise=norm(low)/delta_noise,
                high_coherent_gap_over_noise=norm(high)/delta_noise,
                low_high_interference_over_noise2=2*float(np.vdot(low,high).real)/delta_noise**2,
                strongest_isolated_orders=mode_groups(order, modal, delta_noise))


def compare_array(case, coefficients, order, maps, noise, design):
    rng = np.random.default_rng([2026111906, int(case.split('p')[-1])])
    report, fields, predictions, feeds = {}, [], [], {}
    for n in SIZES:
        for R in RINGS:
            name = f'{R:.2f}_{n}'
            grid = design[(R,n)]
            y = maps@grid['E0']
            scores = norm_rows(y.T)
            i = int(np.argmax(scores))
            single = np.eye(n, dtype=complex)[i]
            menu = feed_menu(y, grid['gram'], rng)
            entry = dict(radius=R, emitters=n, gram_min=float(grid['eigen'][0]), gram_max=float(grid['eigen'][-1]),
                         candidates_in_menu=len(menu), designs={})
            arms = [('single', f'best independent angle {i}', single)]
            for cap in FEED_CAPS:
                goal, feed_norm, method, w = max((t for t in menu if t[1] <= cap+1e-12), key=lambda t:t[0])
                arms.append((f'coherent_l2_{cap:g}', method, w))
            # Waveform choice is based on the fixed endpoint-pair gap only.
            for arm, tag, w in arms:
                label = name+'_'+arm
                response = y@w
                fields.append(grid['E0']@w)
                predictions.append(response)
                feeds[label] = w
                entry['designs'][arm] = dict(selection=tag, **describe_wave(coefficients[R], order, grid, w, response, noise))
            report[name] = entry
        for arm in ('single', 'coherent_l2_2', 'coherent_l2_4'):
            far = report[f'1.60_{n}']['designs'][arm]
            near = report[f'0.72_{n}']['designs'][arm]
            factor = far['incident_l2']/near['incident_l2']
            near['near_far_gap_ratio_equal_emitted_power'] = near['gap_over_noise']/far['gap_over_noise']
            near['near_far_gap_ratio_equal_sampled_incident_l2'] = factor*near['gap_over_noise']/far['gap_over_noise']
            near['near_rescaled_emitted_power_for_far_incident_l2'] = factor**2
            near['near_rescaled_feed_l2_for_far_incident_l2'] = factor*near['feed_l2']
    return report, np.column_stack(fields), np.column_stack(predictions), feeds


def angular_test(order, coeff, direct, old, theta, noise, H):
    high = np.abs(order)>=6
    low = ~high
    old_n = old.shape[1]
    old_coef = fourier(old, order[low])
    alias = old_coef-coeff[RINGS[0]][:, low]
    # The Nyquist bin aliases ±6, so transport uses the 11 unique low orders.
    old_bin_order = np.arange(-old_n//2, old_n//2)
    old_all_coeff = fourier(old, old_bin_order)
    old_exact_alias_mask = ((order[None]-order[low, None]) % old_n==0) & high[None]
    old_exact_alias = coeff[RINGS[0]] @ old_exact_alias_mask.T
    details = dict(old_unique_orders=order[low].tolist(), nyquist_orders=[-6, 6],
                   old_nyquist_row_gap_over_noise=norm(old_all_coeff[:, 0])/noise,
                   old_alias_relative_to_exact_low=norm(alias)/norm(coeff[RINGS[0]][:, low]),
                   old_alias_prediction_relative_error=relative(old_exact_alias, alias),
                   radii={})
    for R in RINGS:
        result = direct[R]
        true_total = modal_rms(result)
        oracle_low = synthesize(coeff[R][:,low], order[low], theta)
        transported = old_coef*H[R][low][None]/H[RINGS[0]][low][None]
        transported_full = synthesize(transported, order[low], theta)
        full_coeff = coeff[R]
        fresh_high = synthesize(full_coeff[:,high], order[high], theta)
        error = result-transported_full
        alias_contribution = oracle_low-transported_full
        f_direct = fourier(result, order)
        hankel_pred = coeff[RINGS[0]]*H[R][None]/H[RINGS[0]][None]
        h_ratio = {str(n):float(abs(H[R][order==n][0]/H[RINGS[0]][order==n][0])) for n in [0, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48]}
        rec = dict(dense_data_rms_gap_over_noise=true_total/noise,
                   best_single_angle_gap_over_noise=float(max(norm_rows(result.T))/noise),
                   best_single_angle_deg=float(theta[np.argmax(norm_rows(result.T))]*180/np.pi),
                   exact_low_fraction=norm(full_coeff[:,low])**2/norm(full_coeff)**2,
                   high_power_fraction=norm(full_coeff[:,high])**2/norm(full_coeff)**2,
                   high_abs_ge_8_power_fraction=norm(full_coeff[:,abs(order)>=8])**2/norm(full_coeff)**2,
                   high_abs_ge_12_power_fraction=norm(full_coeff[:,abs(order)>=12])**2/norm(full_coeff)**2,
                   strongest_angles=mode_groups(order, full_coeff, noise),
                   graf_coefficients_against_dense_fft=relative(hankel_pred, f_direct),
                   graphic_truncation_relative={str(n):relative(synthesize(hankel_pred[:,abs(order)<=n], order[abs(order)<=n], theta),result) for n in (16,32,48,64)},
                   low_transport_error_rms_over_noise=modal_rms(error)/noise,
                   low_alias_error_rms_over_noise=modal_rms(alias_contribution)/noise,
                   high_missing_error_rms_over_noise=modal_rms(fresh_high)/noise,
                   low_alias_relative_to_oracle_low=relative(transported,full_coeff[:,low]),
                   low_alias_error_power_over_total=modal_rms(alias_contribution)**2/true_total**2,
                   low_total_error_power_over_total=modal_rms(error)**2/true_total**2,
                   error_parseval_relative_gap=abs(modal_rms(error)**2-modal_rms(alias_contribution)**2-modal_rms(fresh_high)**2)/true_total**2,
                   exact_parseval_relative_gap=abs(true_total**2-norm(full_coeff)**2)/true_total**2,
                   predicted_hankel_magnitude_ratios=h_ratio)
        details['radii'][f'{R:.2f}'] = rec
    # Train no material model: this stresses propagating a finite-sample far FFT
    # past roundoff rather than transferring exact physics coefficients.
    far_fft = fourier(direct[RINGS[0]], order)
    near = direct[RINGS[1]]
    details['dense_far_to_near_fft_transport'] = {
        str(n):relative(synthesize((far_fft*H[RINGS[1]][None]/H[RINGS[0]][None])[:,abs(order)<=n],order[abs(order)<=n],theta),near)
        for n in (5,12,16,24,32,48,64)}
    return details


def run(args):
    sys.path.insert(0, str(args.archive.resolve()))
    from research.research_z_alg_r4_diversity.src.physics import Model

    operators = args.archive/R4/'data/operators.npz'
    ops = read_npz(operators)
    centers, k = ops['centers'], float(ops['k'])
    old_angles = np.arctan2(ops['sources'][:,1],ops['sources'][:,0])
    if not np.allclose(np.exp(1j*old_angles), np.exp(1j*angle_grid(len(old_angles)))):
        raise ValueError('Old DFT requires the saved equal-angle, zero-offset array order')
    if not np.allclose(norm_rows(ops['sources']), RINGS[0]):
        raise ValueError('The historical array does not use the expected radius')
    if max(norm_rows(centers))+float(ops['radius']) >= RINGS[1]:
        raise ValueError('The Graf domain must fit inside the nearer line-source ring')
    tick, cpu0 = perf_counter(), process_time()
    theta = angle_grid(DENSE_ANGLES)
    order, inside, H = modes(centers, k)
    dense_inputs = {R:incident(centers,k,angles_to_xy(R,theta)) for R in RINGS}
    design = layouts(centers, k)
    results = dict(config=dict(cases=args.cases, old_samples=ops['E0'].shape[1], receivers=ops['S'].shape[0],
                               material_cells=ops['G'].shape[0], angular_validation_samples=DENSE_ANGLES,
                               max_angular_order=ORDER, radii=RINGS, array_sizes=SIZES, feed_l2_caps=FEED_CAPS,
                               fixed_radiated_power=1., region_radius=float(max(norm_rows(centers))+float(ops['radius'])),
                               observation_noise='old archived 0.001 relative complex RMS × ||Y|| divided by sqrt(12), fixed for each new full receiver shot',
                               power='J0(k |x_i-x_j|) Gram for 2D equal-k radiating monopoles; feeder norm cap and one coherent 48-receiver exposure per constraint',
                               angular_convention='theta sampled at 2πj/N; coefficients multiply exp(-inθ); 11 unique old DFT orders −5 to 5; old Nyquist mixes ±6',
                               truth_usage='None. Archived old noisy data used for noise scale only. Actual prediction uses exact independent solves on the archived two materials.'),
                   geometry_setup=cost(tick,cpu0), cases={})
    arrays = dict(angular_orders=order, dense_angles=theta)
    with np.load(args.materials, allow_pickle=False) as archived:
        for case in args.cases:
            plus = archived[case+'_plus']
            minus = archived[case+'_minus']
            obs = archived[case+'_old_Y']
            noise = .001*norm(obs)/np.sqrt(obs.shape[1])
            tick, cpu0 = perf_counter(), process_time()
            response_map = receiver_map(ops, plus)-receiver_map(ops, minus)
            map_cost = cost(tick,cpu0,full_LU=2,receiver_transpose_RHS=2*ops['S'].shape[0])
            tick, cpu0 = perf_counter(), process_time()
            raw_modes = response_map@inside
            coeff = {R:raw_modes*H[R][None] for R in RINGS}
            arrays[case+'_modal_coefficients'] = np.stack([coeff[R] for R in RINGS])
            layout, selected_fields, predicted, feeds = compare_array(case, coeff, order, response_map, noise, design)
            proposal_cost = cost(tick,cpu0,full_LU=0,complex_RHS=0)
            tick, cpu0 = perf_counter(), process_time()
            E0 = np.column_stack([ops['E0'], *(dense_inputs[R] for R in RINGS), selected_fields])
            model = Model(dict(ops,E0=E0))
            y = model.forward(plus)['y'] - model.forward(minus)['y']
            check_cost = cost(tick,cpu0,**model.counts)
            old_end = obs.shape[1]
            angle_end = old_end+len(RINGS)*DENSE_ANGLES
            old = y[:,:old_end]
            actual = y[:,angle_end:]
            direct = {R:y[:,old_end+i*DENSE_ANGLES:old_end+(i+1)*DENSE_ANGLES] for i,R in enumerate(RINGS)}
            tick, cpu0 = perf_counter(), process_time()
            angular = angular_test(order, coeff, direct, old, theta, noise, H)
            i = 0
            for name, row in layout.items():
                for arm, d in row['designs'].items():
                    d['actual_gap_over_noise'] = norm(actual[:,i])/noise
                    d['full_VIE_validation_relative_gap'] = relative(predicted[:,i],actual[:,i])
                    arrays[case+'_'+name+'_'+arm] = feeds[name+'_'+arm]
                    i += 1
            investigation_cost = cost(tick,cpu0)
            results['cases'][case] = dict(noise_per_view=noise, old_endpoint_frobenius_gap_over_noise=norm(old)/(noise*np.sqrt(obs.shape[1])),
                                          old_response_map_relative_validation=relative(response_map@ops['E0'],old),
                                          operator_cost=map_cost, design_cost=proposal_cost, validation_cost=check_cost,
                                          angular_postprocess_cost=investigation_cost, angular=angular, arrays=layout)
            print(case, 'modal-high near/far', *(round(angular['radii'][R]['high_power_fraction'],4) for R in ('0.72','1.60')),
                  'single near/far', *(round(layout[R+'_24']['designs']['single']['actual_gap_over_noise'],3) for R in ('0.72','1.60')),
                  'coherent2 near/far', *(round(layout[R+'_24']['designs']['coherent_l2_2']['actual_gap_over_noise'],3) for R in ('0.72','1.60')),
                  'coherent4 near/far', *(round(layout[R+'_24']['designs']['coherent_l2_4']['actual_gap_over_noise'],3) for R in ('0.72','1.60')), flush=True)
    args.output.with_suffix('.json').write_text(json.dumps(results, indent=2, allow_nan=False)+'\n')
    np.savez_compressed(args.output.with_suffix('.npz'), **arrays)


if __name__=='__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, default=Path(os.environ.get('ZALG_WEAK_ARCHIVE', DEFAULT_ARCHIVE)))
    parser.add_argument('--materials', type=Path, default=HERE/'results/weak_r4.npz')
    parser.add_argument('--output', type=Path, default=HERE/'results/angular_gap')
    parser.add_argument('--cases', nargs='+', default=list(PARENT_IDS))
    run(parser.parse_args())

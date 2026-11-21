"""Open the historical material references only after the online pair seal."""

import argparse
import csv
import json
from pathlib import Path

from run import OUT, ball_status, dump_json
from archives import AUS, R8, RING, SRAP_CASES, WORKSPACE, observations, rel, truth_path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def select_panel_pair(document):
    rows = [row for row in document['pairs'] if row['direction'] == 'projected_task']
    if not rows:
        return None
    return max(rows, key=lambda row: (row['archive_ball']['both_inside'],
                                      row['same_budget_center_ball']['both_inside'],
                                      row['task_gap']))


def materials(m, truth, support, q):
    diff = m - truth
    return {
        'material_rel_l2': float(np.linalg.norm(diff) / np.linalg.norm(truth)),
        'target_rel_l2': float(np.linalg.norm(diff[support]) / np.linalg.norm(truth[support])),
        'background_rmse': float(np.sqrt(np.mean(diff[~support] ** 2))),
        'target_mean_chi': float(np.mean(m[support])),
        'background_mean_chi': float(np.mean(m[~support])),
        'regional_mean_chi': float(q @ m),
        'regional_absolute_task_error': float(abs(q @ diff)),
        'peak_chi': float(np.max(m)),
    }


def original_austria_radius(Y):
    """Recover only the original noise norm, independently of any material label."""
    contract = R8 / 'delegated/paper_baselines/PAPER_CONTRACTS.json'
    table = load(contract)['NIE2016']['cases']
    index = next(i for i, case in enumerate(table) if case['id'] == AUS)
    rng = np.random.default_rng(8600 + index)
    unit = rng.normal(size=Y.shape) + 1j * rng.normal(size=Y.shape)
    unit /= np.linalg.norm(unit)
    fraction = .1
    # Y = Yclean + δu and δ = .1 ||Yclean|| give a quadratic in δ.
    a = 1 - fraction**2
    b = fraction**2 * np.vdot(unit, Y).real
    c = -fraction**2 * np.linalg.norm(Y)**2
    return float((-b + np.sqrt(b**2 - a*c)) / a)


def make_figure(path, gallery, include_reference=True):
    ncols = 5 if include_reference else 4
    fig, axes = plt.subplots(len(gallery), ncols, figsize=(2.6*ncols, 2.45*len(gallery) + 1), layout='constrained')
    fig.suptitle('Full-wave material pairs on archived Maxwell acquisitions', fontsize=19, x=.07, ha='left')
    maximum = max(1e-9, max(float(np.max(np.abs(material[2] - material[1])))
                           for _, _, _, material, _ in gallery))
    for ii, ((document, row, truth, material, side), axis) in enumerate(zip(gallery, axes, strict=True)):
        center, minus, plus = material
        imgs = [center, minus, plus]
        cols = ['Historical center', 'Minus', 'Plus']
        if include_reference:
            imgs.append(truth)
            cols.append('Coarse material reference')
        imgs.append(plus - minus)
        cols.append('Contrast difference')
        for jj, (ax, title, image) in enumerate(zip(axis, cols, imgs, strict=True)):
            if jj == ncols - 1:
                ax.imshow(image.reshape(side, side), origin='lower', cmap='RdBu_r', vmin=-maximum, vmax=maximum)
            else:
                ax.imshow(image.reshape(side, side), origin='lower', cmap='viridis', vmin=0, vmax=20)
            if ii == 0:
                ax.set_title(title, fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
            if jj == 0:
                status = 'archive ball' if row['archive_ball']['both_inside'] else (
                    'center ball' if row['same_budget_center_ball']['both_inside'] else 'outside both balls'
                )
                if document['case'] == AUS and status == 'archive ball':
                    status = '12% archive ball'
                ax.set_ylabel(f'{document["case"]}\nΔ task {row["task_gap"]:.3f} · {status}', fontsize=8)
    scale = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(0, 20), cmap='viridis')
    fig.colorbar(scale, ax=axes[:, :-1].ravel().tolist(), shrink=.45, location='bottom', label='Real contrast (common scale)')
    change_scale = matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(-maximum, maximum), cmap='RdBu_r')
    fig.colorbar(change_scale, ax=axes[:, -1].tolist(), shrink=.7, location='right', label='Contrast difference (common scale)')
    fig.savefig(path, dpi=160)
    plt.close(fig)


def csv_file(filename, records):
    keys = list(dict.fromkeys(key for row in records for key in row))
    with filename.open('w', encoding='utf-8', newline='') as file:
        writer = csv.DictWriter(file, keys)
        writer.writeheader()
        writer.writerows(records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    out = args.output
    seal = load(out / 'ENDPOINT_SEAL.json')
    if not seal['truth_read'] is False or seal['cases'] != [RING, AUS, *SRAP_CASES]:
        raise ValueError('The online endpoint seal is incomplete')
    documents = {case: load(out / f'{case}.json') for case in seal['cases']}
    panel_pairs = {case: select_panel_pair(document) for case, document in documents.items()}
    baseline_rows, pair_rows, report, gallery = [], [], {}, []
    for obs in observations():
        case, doc = obs.id, documents[obs.id]
        with np.load(truth_path(case), allow_pickle=False) as labels:
            truth = labels['m'].ravel()
            true_y = labels['Y_clean'] if 'Y_clean' in labels else None
        host = .1 if case == RING else 0.
        support = truth > host + .05
        if not np.any(support) or np.all(support):
            raise ValueError(f'Empty target/background support classification: {case}')
        field = obs.model.forward(truth, label='offline_coarse_reference')['y']
        y_norm = np.linalg.norm(obs.Y)
        r8radius = original_austria_radius(obs.Y) if case == AUS else None
        with np.load(out / f'{case}.npz', allow_pickle=False) as stored:
            q = stored['task_vector']
            note = {
                'truth_reference': rel(truth_path(case)),
                'truth_regional_chi': float(q @ truth),
                'truth_coarse_relative_residual_against_archive': float(np.linalg.norm(field - obs.Y) / y_norm),
                'truth_reference_fits_archive_ball': bool(np.linalg.norm(field - obs.Y) <= obs.radius),
                'target_support_threshold': host + .05,
                'offline_reference_full_wave_replay_count': dict(obs.model.counts),
            }
            if true_y is not None:
                note['truth_replay_difference_from_archived_clean_field'] = float(np.linalg.norm(field - true_y))
            if r8radius is not None:
                note['original_r8_noise_absolute_radius_reconstructed_offline'] = r8radius
                note['original_r8_noise_relative_to_observation'] = float(r8radius / y_norm)
            for arm, account in doc['historical_baselines'].items():
                metric = materials(stored[f'history__{arm}__m'], truth, support, q)
                entry = {'observation': case, 'historical_method': arm, 'archive_ball': account['archive_ball'],
                         'relative_residual': account['relative_residual'], **metric}
                if r8radius is not None:
                    entry['original_r8_noise_ball'] = ball_status(account['distance_to_archive_interval'], r8radius)
                baseline_rows.append(entry)
            for pair in doc['pairs']:
                prefix = pair['array_prefix']
                common = {
                    'observation': case, 'center': pair['center'], 'direction': pair['direction'],
                    'amplitude': pair['amplitude'], 'prefix': prefix, 'task_gap': pair['task_gap'],
                    'material_rms_between_ends': pair['material_rms_between_ends'],
                    'linearized_field_separation': pair['linearized_field_separation'],
                    'field_gap_low': pair['field_separation_interval'][0],
                    'field_gap_high': pair['field_separation_interval'][1],
                    'even_drift_low': pair['even_midpoint_drift_interval'][0],
                    'even_drift_high': pair['even_midpoint_drift_interval'][1],
                    'archive_ball_both': pair['archive_ball']['both_inside'],
                    'common_center_ball_both': pair['same_budget_center_ball']['both_inside'],
                    'midpoint_overlap': pair['midpoint_same_budget_overlap_proven'],
                }
                for end in ('minus', 'plus'):
                    mat = stored[f'{prefix}__m_{end}']
                    score = materials(mat, truth, support, q)
                    distance = pair[f'{end}_distance_to_archive_interval']
                    entry = {**common, 'end': end, 'archive_ball': pair['archive_ball'][end],
                             'archive_distance_low': distance[0], 'archive_distance_high': distance[1],
                             'archive_radius': obs.radius, **score}
                    if r8radius is not None:
                        entry['original_r8_noise_ball'] = ball_status(distance, r8radius)
                    pair_rows.append(entry)
            chosen = panel_pairs[case]
            if chosen is not None:
                prefix = chosen['array_prefix']
                gallery.append((doc, chosen, truth, [stored[f'history__{chosen["center"]}__m'],
                                stored[f'{prefix}__m_minus'], stored[f'{prefix}__m_plus']], obs.side))
                note['displayed_pair_selected_before_truth'] = prefix
        report[case] = note
    csv_file(out / 'BASELINES.csv', baseline_rows)
    csv_file(out / 'PAIRS.csv', pair_rows)
    dump_json(out / 'OFFLINE_SCORE.json', report)
    make_figure(out / 'HIGH_DIMENSIONAL_ENDPOINTS.png', gallery, include_reference=False)
    private_figures = WORKSPACE / '.artifacts/z_alg_secant_inputs'
    private_figures.mkdir(parents=True, exist_ok=True)
    make_figure(private_figures / 'HIGH_DIMENSIONAL_PAIRS.png', gallery)


if __name__ == '__main__':
    main()

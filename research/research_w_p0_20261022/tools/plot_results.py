#!/usr/bin/env python3
"""Render measured caches only. Does not import or execute the physics backend.

Usage: python tools/plot_results.py [--root ROOT]
Corrected source-average caches are used only when all configured seeds have
nx=128 and 192 records and files. Original raw results are never changed.
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

plt.rcParams.update({'font.size': 10, 'axes.titlesize': 11, 'figure.dpi': 120,
                     'savefig.dpi': 180, 'svg.fonttype': 'none'})


def read_json(path):
    return json.loads(path.read_text())


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main(root):
    cfg = read_json(root / 'w_config.json')
    out = root / 'figures'
    raw = out / 'raw'
    raw.mkdir(parents=True, exist_ok=True)
    entries, missing = [], []

    def save(fig, name, rows, source, note):
        fig.savefig(out / (name + '.png'), bbox_inches='tight')
        fig.savefig(out / (name + '.svg'), bbox_inches='tight')
        plt.close(fig)
        write_csv(raw / (name + '.csv'), rows)
        entries.append(dict(figure=name, source=source, note=note))

    original = root / 'results/qualification'
    corrected = root / 'results/qualification_source_average'
    def complete(folder):
        path = folder / 'WALL_GRIDS.json'
        if not path.exists():
            return False
        rows = read_json(path)
        return all(any(r['seed'] == seed and r['nx'] == nx for r in rows)
                   and (folder / f'WALL_seed{seed}_nx{nx}.npz').exists()
                   for seed in cfg['truth_seeds'] for nx in (128, 192))
    wall_dir = corrected if complete(corrected) else original
    wall_source = str(wall_dir.relative_to(root))
    wall_label = 'source-cell averaged' if wall_dir == corrected else 'original collocation'
    if wall_dir != corrected:
        missing.append('Corrected wall grid set incomplete: original caches shown explicitly; rerun after completion.')

    path = root / 'results/p0a/results.json'
    if path.exists():
        sweep = read_json(path)['bandwidth_sweep']
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        for key, label in [('wall_distance', 'Calibration wall distance'),
                           ('joint_shifted_target_distance', 'Joint distance (shifted target)')]:
            ax.plot([r['fractional_bandwidth'] for r in sweep], [r[key] for r in sweep], 'o-', label=label)
        ax.axhline(1, color='gray', ls=':', label='Distance 1')
        ax.set(xlabel='Fractional bandwidth', ylabel='Whitened mean-data distance',
               title='1D analytic regression: bandwidth dependence')
        ax.legend(); ax.grid(alpha=.2)
        save(fig, '01_1d_bandwidth', sweep, str(path.relative_to(root)),
             'Independent planar analytic regression; not a 2D witness. Non-nested frequency grids.')
    else:
        missing.append('01_1d_bandwidth NOT_DRAWN: missing analytic JSON.')

    # Export actual analytic complex arrays as long-form CSV, one NPZ at a time.
    analytic_npz = root / 'results/p0a/arrays.npz'
    if analytic_npz.exists():
        with np.load(analytic_npz, allow_pickle=False) as cache:
            records = []
            for key in cache.files:
                a = cache[key]
                for i, value in enumerate(a.ravel()):
                    records.append(dict(array=key, index=i, real=float(np.real(value)), imag=float(np.imag(value))))
            write_csv(raw / '01_1d_cached_arrays.csv', records)

    extent = [cfg['domain'][0][0], cfg['domain'][0][1], cfg['domain'][1][0], cfg['domain'][1][1]]
    geom_rows, field_rows = [], []
    fig_g, axes_g = plt.subplots(len(cfg['truth_seeds']), 3, figsize=(12, 7), squeeze=False)
    fig_f, axes_f = plt.subplots(len(cfg['truth_seeds']), 4, figsize=(14, 7), squeeze=False)
    loaded = 0
    for si, seed in enumerate(cfg['truth_seeds']):
        path = wall_dir / f'WALL_seed{seed}_nx192.npz'
        if not path.exists():
            missing.append(f'Geometry/field seed {seed} NOT_DRAWN: missing nx192 cache.')
            for ax in list(axes_g[si]) + list(axes_f[si]):
                ax.text(.5, .5, 'NOT_DRAWN: cache missing', ha='center', transform=ax.transAxes)
            continue
        with np.load(path, allow_pickle=False) as cache:
            eps = cache['epsilon_cal']
            ny, nx = eps.shape
            xs = extent[0] + (np.arange(nx) + .5) * (extent[1] - extent[0]) / nx
            ys = extent[2] + (np.arange(ny) + .5) * (extent[3] - extent[2]) / ny
            air = (np.abs(eps - 1) < 1e-12).astype(float)
            for ci, (data, label, cmap, norm) in enumerate([
                    (eps.real, 'Re(epsilon), calibration', 'viridis', Normalize(1, 6)),
                    (eps.imag, 'Im(epsilon), f/f0=1', 'magma', Normalize(0, .04)),
                    (air, 'Air cells (exterior + cavity)', 'gray_r', Normalize(0, 1))]):
                ax = axes_g[si, ci]
                im = ax.imshow(data, origin='lower', extent=extent, cmap=cmap, norm=norm, aspect='equal')
                ax.set(title=f'Seed {seed}\n{label}', xlabel='x / lambda0', ylabel='y / lambda0')
                fig_g.colorbar(im, ax=ax, shrink=.7)
            fields = {label: cache[f'field_{label}'][:, :, 0].copy() for label in ('cal', 'task')}
            maxmag = max(np.abs(a).max() for a in fields.values())
            for li, label in enumerate(('cal', 'task')):
                field = fields[label]
                for mi, (values, title, cmap, norm) in enumerate([
                        (np.abs(field), 'Magnitude', 'inferno', Normalize(0, maxmag)),
                        (np.angle(field), 'Phase [rad]', 'twilight', Normalize(-np.pi, np.pi))]):
                    ax = axes_f[si, li * 2 + mi]
                    im = ax.imshow(values, origin='lower', extent=extent, cmap=cmap, norm=norm, aspect='equal')
                    ax.contour(xs, ys, eps.real, levels=[1.05], colors='white', linewidths=.4)
                    ax.set(title=f'{seed} {label}: {title}', xlabel='x / lambda0', ylabel='y / lambda0')
                    fig_f.colorbar(im, ax=ax, shrink=.7)
            for iy, y in enumerate(ys):
                for ix, x in enumerate(xs):
                    geom_rows.append(dict(seed=seed, nx=nx, x=x, y=y, epsilon_real=eps[iy, ix].real,
                                          epsilon_imag=eps[iy, ix].imag, air_cell=int(air[iy, ix])))
                    field_rows.append(dict(seed=seed, nx=nx, x=x, y=y,
                        cal_real=fields['cal'][iy, ix].real, cal_imag=fields['cal'][iy, ix].imag,
                        cal_magnitude=abs(fields['cal'][iy, ix]), cal_phase=np.angle(fields['cal'][iy, ix]),
                        task_real=fields['task'][iy, ix].real, task_imag=fields['task'][iy, ix].imag,
                        task_magnitude=abs(fields['task'][iy, ix]), task_phase=np.angle(fields['task'][iy, ix])))
            loaded += 1
    fig_g.suptitle(f'Cached irregular dielectric geometries, nx=192 ({wall_label})')
    fig_f.suptitle(f'Actual total field: first source, f/f0=1, nx=192 ({wall_label})')
    fig_g.tight_layout(); fig_f.tight_layout()
    if loaded:
        save(fig_g, '02_geometry', geom_rows, wall_source, 'Calibration material; exact-air cells include cavity and exterior; coefficient boundary cells are area averaged.')
        save(fig_f, '03_first_source_field', field_rows, wall_source, 'Cached total field, source index 0. Calibration and task use the same source; per-seed shared magnitude scale.')
    else:
        plt.close(fig_g); plt.close(fig_f)

    mesh_path = wall_dir / 'WALL_GRIDS.json'
    wall_rows = read_json(mesh_path) if mesh_path.exists() else []
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    convergence_rows = []
    for folder, version in [(original, 'Original'), (corrected, 'Source averaged')]:
        cylinder_path = folder / 'CYLINDER.json'
        cylinder = read_json(cylinder_path) if cylinder_path.exists() else []
        if cylinder:
            axes[0].semilogy([r['nx'] for r in cylinder], [r['whitened_error'] for r in cylinder],
                             'o-' if folder == original else 's--', label=version + ' vs Mie')
            for r in cylinder:
                convergence_rows.append(dict(test='independent_Mie', source=str(cylinder_path.relative_to(root)), **r))
        else:
            missing.append(version + ' independent Mie convergence NOT_DRAWN: missing JSON.')
        history = folder / 'WALL_GRIDS.json'
        history = read_json(history) if history.exists() else []
        for si, seed in enumerate(cfg['truth_seeds']):
            sr = [r for r in history if r['seed'] == seed and r.get('mesh')]
            for ki, key in enumerate(sorted({k for r in sr for k in r['mesh']})):
                rr = [r for r in sr if key in r['mesh']]
                label = f"S{si+1} {key.replace('_mesh_whitened', '').replace('_', ' ')}"
                color = f'C{si*4+ki}'
                axes[1].semilogy([r['nx'] for r in rr], [r['mesh'][key] for r in rr],
                    color=color, marker='o' if folder == original else 's',
                    ls='-' if folder == original else '--', label=label if folder == original else None)
                for r in rr:
                    convergence_rows.append(dict(test=key, seed=seed, nx=r['nx'], previous_nx=r['previous_nx'],
                        whitened_error=r['mesh'][key], source=str(folder.relative_to(root))))
    for ax in axes:
        ax.axhline(cfg['numerics']['mesh_whitened_tolerance'], color='red', ls='--', label='Tolerance 0.1')
        ax.set(xlabel='Horizontal grid size nx', ylabel='Full-data whitened norm')
        ax.grid(alpha=.2); ax.legend(fontsize=7)
    axes[0].set_title('Independent cylinder / Mie reference')
    axes[1].set_title('Wall grids: original ladder + corrected finest pair')
    axes[1].text(.02, .10, 'Circles / solid: original; squares: source averaged\nS1/S2: configured seed order',
                 transform=axes[1].transAxes, fontsize=7,
                 bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
    fig.tight_layout()
    save(fig, '04_convergence', convergence_rows, 'results/qualification; results/qualification_source_average',
         'Full vector norms, not per-channel RMS. Adjacent mesh consistency is not a continuum certificate.')

    # Hold only the 2 small receiver-data sets; never hold two NPZ archives open.
    response = {}
    for seed in cfg['truth_seeds']:
        path = wall_dir / f'WALL_seed{seed}_nx192.npz'
        if path.exists():
            with np.load(path, allow_pickle=False) as cache:
                response[seed] = {k: cache[k].copy() for k in cache.files if k.startswith('data_')}
    if len(response) == 2:
        seed1, seed2 = cfg['truth_seeds']
        distance_rows = []
        sigma = cfg['noise']['each_real_quadrature_std']
        for band, ff in cfg['frequencies'].items():
            def data(seed, label):
                return np.stack([response[seed][f'data_{f:.8f}_{label}'] for f in ff])
            c1, c2 = data(seed1, 'cal'), data(seed2, 'cal')
            t1, t2 = data(seed1, 'task'), data(seed2, 'task')
            distance_rows.append(dict(band=band, fractional_bandwidth=max(ff)-min(ff), number_frequencies=len(ff),
                calibration_distance=np.linalg.norm(c1-c2)/sigma,
                fixed_target_response_distance=np.linalg.norm((t1-c1)-(t2-c2))/sigma,
                fixed_target_total_data_distance=np.linalg.norm(t1-t2)/sigma,
                target_centroid_difference=0, searched_witness=False, GO=False, nx=192, source=wall_source))
        distance_rows.sort(key=lambda r:r['fractional_bandwidth'])
        fig, ax = plt.subplots(figsize=(8, 4.8))
        for key, label in [('calibration_distance', 'Calibration'), ('fixed_target_response_distance', 'Same-target response (task minus cal)'),
                           ('fixed_target_total_data_distance', 'Same-target total data')]:
            ax.plot([r['fractional_bandwidth'] for r in distance_rows], [r[key] for r in distance_rows], 'o-', label=label)
        ax.set(xlabel='Fractional bandwidth (cached five frequencies)', ylabel='Whitened mean-data distance',
               title=f'Two-seed reference-pair diagnostic ({wall_label})')
        ax.text(.02, .78, 'Fixed same target: centroid difference = 0\nNo candidate search; no witness; no GO', transform=ax.transAxes, va='top', fontsize=9)
        ax.grid(alpha=.2); ax.legend(loc='center right', fontsize=8)
        save(fig, '05_reference_pair_bandwidth', distance_rows, wall_source,
             'Only two cached bandwidths; connecting lines guide the eye. Mean response distances, not sampled noisy data. No searched witness.')
    else:
        missing.append('05_reference_pair_bandwidth NOT_DRAWN: both nx192 receiver caches required.')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    cost_rows = []
    for seed in cfg['truth_seeds']:
        rows = [r for r in wall_rows if r['seed'] == seed]
        axes[0].plot([r['nx'] for r in rows], [r['state_wall_seconds'] for r in rows], 'o-', label=str(seed))
        for r in rows:
            worst = max(r.get('mesh', {}).values(), default=None)
            cost_rows.append(dict(source=wall_source, seed=seed, nx=r['nx'], measured_state_seconds=r['state_wall_seconds'],
                                  worst_adjacent_mesh_norm=worst, scope='state solves only; excludes setup and end-to-end inference'))
        rows = [r for r in rows if r.get('mesh')]
        axes[1].loglog([max(r['mesh'].values()) for r in rows], [r['state_wall_seconds'] for r in rows], 'o-', label=str(seed))
    bench_path = original / 'BENCHMARK.json'
    if bench_path.exists():
        bench = read_json(bench_path)
        cost_rows.append(dict(source=str(bench_path.relative_to(root)), nx=bench['nx'], measured_state_seconds=bench['wall_seconds'],
            setup_seconds=bench['setup_seconds'], scope='separate disk benchmark, one frequency; state includes Jacobian'))
        axes[0].text(.03, .96, f"Separate nx={bench['nx']} disk benchmark:\nsetup {bench['setup_seconds']:.4g} s\nstate + Jacobian {bench['wall_seconds']:.4g} s", transform=axes[0].transAxes, va='top', fontsize=8)
    axes[0].set(xlabel='Horizontal grid size nx', ylabel='Measured state wall seconds', title='Wall states: partial cost, all cached frequencies')
    axes[1].axvline(.1, color='red', ls='--', label='Tolerance 0.1')
    axes[1].set(xlabel='Worst adjacent full-data mesh norm', ylabel='Measured state wall seconds', title='Partial cost versus measured mesh consistency')
    for ax in axes:
        ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.suptitle('Measured partial timings: excludes wall setup, search, and end-to-end inference')
    fig.tight_layout()
    save(fig, '06_partial_cost', cost_rows, wall_source, 'No end-to-end inference cost or speedup measured. Wall setup timings are absent; separate benchmark is explicitly labeled.')

    four_path = root / 'results/four_blocks/RESULTS.json'
    if four_path.exists():
        rows = read_json(four_path)
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.7))
        # All eight measured cases are shown; x encodes frequency and labels encode contrast/seed.
        metrics = [('full_field_equivalence_relative', 'Discrete identity: relative field error', True),
                   ('born_induced_relative_error', 'Born induced-response relative error', False),
                   ('core_sigma_min', 'Core minimum singular value', False)]
        for ax, (key, title, log) in zip(axes, metrics):
            for si, seed in enumerate(cfg['truth_seeds']):
                for pi, probe in enumerate(sorted({r['probe_index'] for r in rows})):
                    rr = sorted([r for r in rows if r['seed'] == seed and r['probe_index'] == probe], key=lambda r:r['frequency'])
                    if rr:
                        ax.plot([r['frequency'] for r in rr], [r[key] for r in rr],
                                color=f'C{si}', marker='o' if pi == 0 else 's', ls='-' if pi == 0 else '--',
                                label=f'S{si+1}, target {probe}')
            if log:
                ax.set_yscale('log')
            ax.set(title=title, xlabel='Frequency f/f0')
            ax.grid(alpha=.2)
            ax.legend(fontsize=7)
        fig.suptitle('Four-block discrete diagnostic: task-restricted 21-cell ports')
        fig.text(.5, .01, 'Not a qualified wall witness or reduced-order model; near-resonance behavior is not established.',
                 ha='center', fontsize=9)
        fig.tight_layout(rect=(0, .045, 1, .94))
        save(fig, '07_fourblock_diagnostic', rows, str(four_path.relative_to(root)),
             'Task-restricted 21-cell ports; discrete identity only on an unqualified continuum grid. Not a wall witness/ROM; no near-resonance claim.')
    else:
        missing.append('07_fourblock_diagnostic NOT_DRAWN: missing RESULTS.json.')

    missing.extend(['Candidate search/witness: NOT_RUN; no candidate figure.', 'Posterior/NN/3D/held-out probe validation/inference speedup: NOT_RUN; no fabricated figures.'])
    manifest = dict(wall_source=wall_source, corrected_complete=wall_dir == corrected, figures=entries, missing_or_not_run=missing)
    (out / 'FIGURE_STATUS.json').write_text(json.dumps(manifest, indent=2) + '\n')
    mapping = root / 'research/delegated/w_figures/FIGURE_MAP.md'
    mapping.parent.mkdir(parents=True, exist_ok=True)
    text = ['# Research W P0 figure map', '', f'Wall cache used: `{wall_source}` ({wall_label}).', '',
            'All plots read saved results only. PNG and SVG share a stem; underlying CSVs are in `figures/raw/`.', '']
    for e in entries:
        text += [f"- `figures/{e['figure']}.png` / `.svg`: {e['note']} Source: `{e['source']}`."]
    text += ['', '## Missing / not run', ''] + ['- ' + m for m in missing]
    mapping.write_text('\n'.join(text) + '\n')
    print(json.dumps(dict(figures=len(entries), source=wall_source, status=str(out / 'FIGURE_STATUS.json'))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    main(parser.parse_args().root.resolve())

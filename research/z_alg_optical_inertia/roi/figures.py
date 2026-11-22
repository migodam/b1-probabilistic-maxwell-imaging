#!/usr/bin/env python3
"""Chart all archived ROI outcomes and the original-physics controls."""
from scan import CANONICAL, LEGACY, OUTPUT, STUDY, read

import argparse
import csv
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


def scans(path):
    result = {}
    for U in (20., 4.):
        scan = read(path / f'scan_U{U:g}.json')
        result[U] = {case: (item, grp) for grp in scan['geometries'] for case, item in grp['case_results'].items()}
    return result


def audits(path):
    result = {}
    for U in (20., 4.):
        source = path / f'verified_U{U:g}.json'
        doc = read(source) if source.exists() else {}
        result[U] = ({case:g for g in doc.get('shared_geometry_audits', ()) for case in g['cases']},
                     {(c['case'], c['region']):c for c in doc.get('selected_region_certificates', ())
                      if any(c['verification'].get('strictly_excludes', {}).values())})
    return result


def tabulation(results, verification, endpoints):
    rows = []
    for U, scan in results.items():
        audit, certificates = verification[U]
        for case, (record, group) in scan.items():
            sources = 20 if case == CANONICAL[0] else 10
            minimum_base = min(d['nu_lower_all_regions'] for d in group['base_inertias'])
            for i, region in enumerate(record['regions']):
                cert = certificates.get((case, i))
                proof = audit.get(case, {})
                verified_nu = proof.get('verified_negative_directions_each_region_lower')
                endpoint = endpoints[case]['regions'][i]
                for u in record['tested_caps']:
                    candidates = region['trials'].get(str(u), {})
                    numerical_nu = min((r['nu'] for r in candidates.values()), default=minimum_base)
                    if cert is not None and u <= cert['cap']:
                        verdict = 'excluded by interval verification'
                    elif verified_nu is not None and verified_nu >= sources:
                        verdict = 'source inertia exhausted by interval verification'
                    else:
                        verdict = 'no numeric exclusion'
                    rows.append(dict(U=U, case=case, row=region['row'], column=region['col'],
                                     cap=u, source_count=sources, numeric_evaluated_min_negative_count=numerical_nu,
                                     verified_negative_count_lower=verified_nu, verdict=verdict,
                                     numerical_frobenius_residual_lower=max((r['lower_bound_numeric']
                                        for r in candidates.values()), default=0.),
                                     certified_excluded_cap=cert['cap'] if cert is not None else None,
                                     reference_peak=endpoint['reference_peak'], reference_mean=endpoint['reference_mean'],
                                     Dual_CIE_peak=endpoint['endpoints']['dual_raw']['peak'],
                                     GN_peak=endpoint['endpoints']['physical_gn']['peak'],
                                     projection_zeroed_residual_U20=endpoint['literal_absence_projection']['residual']))
    return rows


def mechanism_table(ax, values, kind):
    if kind == 'exclusion':
        display = np.ma.array(values, mask=~np.isfinite(values))
        cmap = matplotlib.colormaps['viridis'].copy()
        cmap.set_bad('#e8edf0')
        ax.imshow(display, origin='lower', cmap=cmap, vmin=0, vmax=2)
    elif kind == 'inertia':
        ax.imshow(values, origin='lower', cmap='Blues', vmin=1, vmax=7.)
    elif kind == 'projection':
        ax.imshow(values, origin='lower', cmap='RdYlGn_r', norm=TwoSlopeNorm(vmin=0, vcenter=1, vmax=6))
    for row in range(4):
        for col in range(4):
            value = values[row, col]
            text = '–' if not np.isfinite(value) else f'{value:g}' if kind=='exclusion' else f'{value:.1f}'
            color = '#1c2830' if (kind=='exclusion' or value<5.) else 'white'
            ax.text(col, row, text, ha='center', va='center', fontsize=8.2, color=color)
    ax.set_xticks(np.arange(.5, 3., 1.), minor=True)
    ax.set_yticks(np.arange(.5, 3., 1.), minor=True)
    ax.tick_params(which='both', bottom=False, left=False)
    ax.grid(which='minor', color='white', lw=.8)


def plot_grid(path, results, verification, endpoints, arrays):
    fig, axes = plt.subplots(7, 6, figsize=(14.7, 16.2), layout='constrained')
    fig.suptitle('Local contrast exclusions on archived scattering data', fontsize=16, fontweight='bold')
    headings = ('Coarse reference', 'Historical Dual-CIE', 'Full-wave GN', 'Verified excluded cap\nU = 20',
                'Smallest evaluated ν / sources\nconditional U = 4', 'Zeroed-region fit /\nsample noise radius')
    for i, case in enumerate(CANONICAL+LEGACY):
        data = endpoints[case]
        a = axes[i]
        im = None
        for j, name in enumerate(('truth_coarse_m','dual_raw__m','physical_gn__m')):
            im = a[j].imshow(arrays[f'{case}__{name}'].reshape(32, 32), origin='lower', cmap='magma', vmin=0., vmax=4.)
            for edge in (7.5, 15.5, 23.5):
                a[j].axhline(edge, color='white', lw=.45, alpha=.65)
                a[j].axvline(edge, color='white', lw=.45, alpha=.65)
        for j, name in ((1,'dual_raw'),(2,'physical_gn')):
            resid = data['methods'][name]['port_relative_residual']
            a[j].set_xlabel(f'Port residual {resid:.1%}', fontsize=8)
        a[0].set_ylabel(data['label'].replace(' · ', '\n'), fontsize=9, rotation=0,
                        labelpad=9, va='center', ha='right')
        certs = verification[20.][1]
        verified_peak = np.array([certs.get((case,j),{}).get('cap',np.nan) for j in range(16)]).reshape(4,4)
        mechanism_table(a[3], verified_peak, 'exclusion')
        regions = results[4.][case][0]['regions']
        minimum = [min(t['nu'] for trials in r['trials'].values() for t in trials.values()) for r in regions]
        mechanism_table(a[4], np.array(minimum).reshape(4,4)/(20 if i==0 else 10), 'inertia')
        radius = min(v for name,v in data['radii'].items() if name!='conditional_noise_plus_coarse_proxy')
        scaled = [r['literal_absence_projection']['residual']/radius for r in data['regions']]
        mechanism_table(a[5], np.array(scaled).reshape(4,4), 'projection')
        for j, ax in enumerate(a):
            ax.set(xticks=[], yticks=[])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if i==0:
                ax.set_title(headings[j], fontsize=10)
    cbar = fig.colorbar(im, ax=axes[:, :3], orientation='horizontal', aspect=35, shrink=.82, pad=.017, extend='max')
    cbar.set_label('Electric contrast χ (values above 4 saturate)', fontsize=9)
    fig.supxlabel('“–” indicates no exclusion. The certificate requires ν < source count. '
                  'Zeroed-region fit is a numeric diagnostic, with feasibility only at values ≤ 1.', fontsize=10)
    fig.savefig(path / 'roi_archived_regions.png', dpi=150)
    fig.savefig(path / 'roi_archived_regions.pdf')
    plt.close(fig)


def plot_rank(path, results, verification):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.5), layout='constrained')
    fig.suptitle('Negative directions left outside each 8×8 test region', fontsize=15, fontweight='bold')
    cases = (CANONICAL[0], CANONICAL[1], CANONICAL[2], LEGACY[0])
    labels = ('Austria 20×40', 'SRAP 10×20 half', 'SRAP 10×20 full', 'Legacy 10×20 half')
    for p, U in enumerate((20.,4.)):
        ax=axes[p]
        for idx, (case, label) in enumerate(zip(cases,labels,strict=True)):
            _, group = results[U][case]
            ratios = sorted({r['ratio'] for r in group['base_inertias']})
            numeric = [min(r['nu_base'] for r in group['base_inertias'] if r['ratio']==eta) for eta in ratios]
            source_count = 20 if idx==0 else 10
            proof = verification[U][0].get(case,{})
            lower = proof.get('verified_negative_directions_each_region_lower')
            c=f'C{idx}'
            ax.plot(np.arange(len(ratios)), np.array(numeric)/source_count, color=c, alpha=.68,
                    marker=('o','s','^','D')[idx], label=label)
            if lower is not None:
                ax.axhline(lower/source_count, color=c, ls=':', lw=1.35)
        ax.axhline(1., color='black', ls='--', lw=.7, label='Source count')
        ax.set_xticks(np.arange(len(ratios)), [f'{x:g}' for x in ratios])
        ax.set(title=f'U = {U:g}' + (' (conditional)' if U==4. else ''),
               ylabel='Negative eigenvalues / source count', xlabel='Background optical weight λ/μ')
        ax.grid(alpha=.18)
        ax.legend(loc='best', fontsize=8, frameon=False)
    fig.supxlabel('Solid lines show the base numerical minimum over μ. Dotted lines are rigorously verified for every region.', fontsize=9)
    fig.savefig(path / 'roi_negative_directions.png', dpi=170)
    fig.savefig(path / 'roi_negative_directions.pdf')
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=type(OUTPUT),default=OUTPUT)
    parser.add_argument('--figures',type=type(OUTPUT),default=STUDY/'figures')
    args=parser.parse_args()
    args.figures.mkdir(parents=True,exist_ok=True)
    result, verification = scans(args.data), audits(args.data)
    endpoints = read(args.data/'endpoint_metrics.json')
    rows = tabulation(result,verification,endpoints)
    with (args.data/'regions.csv').open('w',newline='',encoding='utf8') as file:
        csv_writer=csv.DictWriter(file,fieldnames=list(rows[0]));csv_writer.writeheader();csv_writer.writerows(rows)
    with np.load(args.data/'endpoint_arrays.npz',allow_pickle=False) as arrays:
        plot_grid(args.figures,result,verification,endpoints,arrays)
    plot_rank(args.figures,result,verification)


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Summarize the frozen numerical grid, strict certified cases and model coverage."""
from collections import Counter
import csv

from scan import HERE, load_public, read, write
from verify import ARMS

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

INTERVAL = ('joint_interval',)*4 + ('single_original_interval','single_coherent_interval')
STYLE = (
    ('pure_joint', 'Real-sector joint', '#175983', 'o'),
    ('optical_joint', 'Optical joint', '#aa522e', 'x'),
    ('definite_joint', 'Positive definite', '#818b3c', 's'),
    ('indefinite_joint', 'Indefinite', '#533977', 'D'),
    ('single_original', 'One original beam', '#197767', 'v'),
    ('single_coherent', 'One coherent beam', '#6d7277', '^'),
)


def matches(arm, profile):
    nu = profile['proof_nu']
    return {
        'pure_joint': profile['lambda_over_mu'] == 0.,
        'optical_joint': profile['lambda_over_mu'] != 0.,
        'definite_joint': nu == 0,
        'indefinite_joint': 0 < nu < 12,
        'single_original': nu == 0,
        'single_coherent': nu == 0,
    }[arm]


def plot(sizes, fixed, caps, cases):
    fig, (left,right) = plt.subplots(1,2,figsize=(9,4.65),layout='constrained',
                                    gridspec_kw={'width_ratios':[1.28,1.]})
    plt.rcParams.update({'pdf.fonttype':42})
    for j, (arm, label, color, marker) in enumerate(STYLE):
        offsets = np.linspace(-.07,.07,len(STYLE))
        values = [sizes[str(cap)][arm]['covered_n'] for cap in caps]
        left.plot(np.arange(len(caps))+offsets[j], values, color=color, alpha=.82,
                  marker=marker, markersize=4, lw=1., label=label)
    left.set(xticks=np.arange(len(caps)),xticklabels=[f'{c:g}' for c in caps],
             yticks=range(8),xlabel=r'Tested peak cap $u$',ylabel='Certified and covered cases / 7',
             title='(a) Best selected profile in each arm',ylim=(-.25,7.3))
    left.grid(axis='y',alpha=.24)
    left.legend(frameon=False,fontsize=7.6,loc='lower left',ncol=2)
    values = [fixed[c]['joint_interval'][0]/fixed[c]['radius'] for c in cases]
    y = np.arange(len(cases))
    colors = ['#225e86' if v > 1 else '#808b95' for v in values]
    right.barh(y,values,color=colors,height=.56)
    right.axvline(1.,ls='--',lw=.9,color='#353535')
    for yi, v in enumerate(values):
        label = f'{v:.2f}' if v else '0'
        right.text(v+.21,yi,label,va='center',ha='left',fontsize=7.8)
    right.set(yticks=y,yticklabels=[c.replace('r14_new_', 'R14-') for c in cases],
              xlabel='Strict joint lower bound / CN99 radius',
              title='(b) Frozen R10 weights at cap 5',xlim=(0,1.16*max(values)))
    right.invert_yaxis()
    right.grid(axis='x',alpha=.22)
    for ax in (left,right):
        ax.set_axisbelow(True)
        ax.spines[['top','right']].set_visible(False)
    for ext in ('png','pdf'):
        fig.savefig(HERE/f'max_struct_transfer.{ext}',dpi=250)
    plt.close(fig)


def main():
    freeze = read(HERE/'frozen.json')
    grid = read(HERE/'numerical_grid.json')
    ver = read(HERE/'verified.json')
    offline = read(HERE/'offline_scalar_scores.json')
    if len(ver['receipts']) != ver['jobs_planned']:
        raise ValueError('Verification is incomplete.')
    _,_,_,_,deltas = load_public(HERE/'public_inputs.npz',freeze)
    num = {(r['case'],r['cap']):r['ranks'] for r in ver['numerical_proposals']}
    rows = []
    counts = {}
    values = ver['receipts'].values()
    fixed_geo = next(p for p in values if p['historical_r10_profile'])
    fixed = {}
    for cap in freeze['caps']:
        profiles = [r for r in values if r['cap']==cap and r['status']=='geometry_verified']
        records = [r for r in grid['values'] if r['cap']==cap]
        counts[str(cap)] = {}
        for arm, attr in zip(ARMS,INTERVAL,strict=True):
            group = []
            for case, d in zip(freeze['cases'],deltas,strict=True):
                truth = offline['rows'][case]
                candidates = []
                for p in profiles:
                    case_row = p['per_case'].get(case,{})
                    if matches(arm,p) and attr in case_row:
                        candidates.append((case_row[attr][0],p['profile'],case_row[attr],p))
                winner = max(candidates,key=lambda x:(x[0],-x[1]),default=None)
                numerical = num[case,cap][arm]
                row = dict(case=case,cap=cap,arm=arm,radius=float(d),truth_status=truth['status'],
                           truth_peak=truth['truth_peak'],
                           numeric_best_profile=numerical[0]['profile'] if numerical else None,
                           numeric_lower_best=numerical[0]['numeric'] if numerical else 0.,
                           verified_profile=None, certified_lower=None, certified_formula_upper=None,
                           verified_mu=None, verified_t=None, verified_nu=None, verified_pass=False,
                           passes_with_confirmed_coverage=False)
                if winner:
                    lower,pid,interval,record = winner
                    passes = lower>d
                    row.update(verified_profile=pid,certified_lower=lower,certified_formula_upper=interval[1],
                               verified_mu=record['mu'],verified_t=record['lambda_over_mu'],
                               verified_nu=record['proof_nu'],verified_pass=bool(passes),
                               passes_with_confirmed_coverage=bool(passes and truth['status']=='inside'))
                group.append(row)
            counts[str(cap)][arm] = dict(numeric_n=sum(r['numeric_lower_best']>r['radius'] for r in group),
                                         verified_n=sum(r['verified_pass'] for r in group),
                                         covered_n=sum(r['passes_with_confirmed_coverage'] for r in group))
            rows.extend(group)
        resolved = [r['nu'] for r in records if r['nu'] is not None]
        counts[str(cap)]['numeric_minimum_nu'] = min(resolved)
        counts[str(cap)]['numeric_inertia_distribution'] = dict(sorted(Counter(resolved).items()))
        counts[str(cap)]['numeric_unresolved_geometry_rows'] = len(records)-len(resolved)
    for case,d in zip(freeze['cases'],deltas,strict=True):
        row = fixed_geo['per_case'][case]
        covered = offline['rows'][case]['status']=='inside'
        fixed[case] = dict(radius=float(d),offline_truth_covered=covered,
                           joint_interval=row['joint_interval'],
                           one_beam_interval=row['single_original_interval'],
                           coherent_beam_interval=row['single_coherent_interval'],joint_pass=row['joint_pass'],
                           passes_with_confirmed_coverage=bool(covered and row['joint_pass']))
    with (HERE/'paired_results.csv').open('w',newline='',encoding='utf-8') as f:
        writable = csv.DictWriter(f,fieldnames=rows[0].keys())
        writable.writeheader()
        writable.writerows(rows)
    summary = dict(counts=counts, cases=freeze['cases'], independent_cases=7,
                   offline_inside_cases=[c for c in freeze['cases'] if offline['rows'][c]['status']=='inside'],
                   frozen_r10_profile=dict(profile=fixed_geo['profile'],u=fixed_geo['cap'],mu=fixed_geo['mu'],
                                           t=fixed_geo['lambda_over_mu'],proof_nu=fixed_geo['proof_nu'],
                                           proof_gap_lower=fixed_geo['proof_gap_lower'],
                                           cases=fixed,cost={k:fixed_geo[k] for k in
                                                             ('geometry_double','geometry_Arb','wall_seconds','cpu_seconds')}),
                   operations=dict(scan=grid['operations'],interval=ver['operations']),
                   cost=dict(scan_wall=grid['wall_seconds'],scan_CPU=grid['cpu_seconds'],
                             interval_wall=ver['wall_seconds'],interval_CPU=ver['cpu_seconds'],
                             offline_wall=offline['wall_seconds'],offline_CPU=offline['cpu_seconds']))
    write(HERE/'summary.json',summary)
    plot(counts,fixed,freeze['caps'],freeze['cases'])


if __name__=='__main__':
    main()

#!/usr/bin/env python3
"""Rebuild the object census and benchmark figure from sealed study receipts."""
from collections import Counter, defaultdict
from csv import DictReader, DictWriter
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
OPTICAL = STUDY.parent
RESEARCH = OPTICAL.parent
ROOT = RESEARCH.parent
ROUTING = STUDY / 'routing'
SECANTS = RESEARCH / 'z_alg_secant_tasks'
CENSUS = []
RING = OPTICAL / 'validation/max_struct'


def read(path):
    return json.loads(path.read_text())


def read_csv(path):
    with path.open(newline='') as stream:
        return list(DictReader(stream))


def identified(receipts, cohort):
    grouped = defaultdict(list)
    for row in receipts:
        if row['cohort'] == cohort:
            grouped[row['object_id']].append(row)
    return grouped


def add(object_id,family,role,kind,observations,source,grid,fit_condition,exposure,heldout,*,
        primary=False,confirmed=False,independently_counted=True):
    observations = tuple(observations)
    if not observations or len(observations)!=len(set(observations)):
        raise ValueError(f'Missing or duplicate observation IDs within {object_id}')
    CENSUS.append(dict(object_id=object_id,family=family,observed_role=role,counting_unit=kind,
        independently_counted=independently_counted,primary_native=primary,
        prospective_confirmation=confirmed,observation_views=len(observations),
        observation_ids=' | '.join(observations),physical_grid=grid,training_condition=fit_condition,
        exposure=exposure,heldout_condition=heldout,
        provenance=' | '.join(str(Path(p).relative_to(ROOT)) for p in source)))


def routing_census(primary, all_cases):
    extended = {
        r['object_id']:r['id'] for r in read(OPTICAL/'physical_revision/array/observations/manifest.json')['observations']
        if (r['sources'],r['receivers'],r['noise_pct'])==(32,32,2)
    }
    extended_manifest = OPTICAL/'physical_revision/array/observations/manifest.json'
    for r in primary:
        r_id = r['object_id']
        acq = [r['id']]
        paths = [ROOT/r['source']]
        physics = '400 MHz, 2 m square; 10 plane waves, 20 full-ring receivers; 10% noise'
        if r['cohort']=='optfresh':
            acq.append(extended[r_id])
            paths.append(extended_manifest)
            physics += '; paired 32 plane waves, 32 full-ring receivers, 2% for detailed controls'
        add(r_id,'Primary full-domain originals' if r['cohort']=='optfresh' else 'Primary full-domain transfer',
            'Born class and full-wave cap; matched solver starts; local native analysis on the archived original conditions' if r['cohort']=='optfresh' else
            'Born class and full-wave cap; matched solver starts', 'synthetic material',acq,paths,
            '96 × 96',physics,'previously exposed; route sees frozen ports only',
            'offline-only offset simulated illuminations and receivers',primary=True)
    for r_id, rs in identified(all_cases,'r14_ring7').items():
        add(r_id,'Compact annular transfer','Common cap-five transfer and same-operator routing',
            'synthetic material',[r['id'] for r in rs],[RING/'public_provenance.json',RING/'public_inputs.npz'],
            '24 × 24','12 cylindrical waves, 48 ring receivers; window 0.96 wavelengths; 1% nominal noise, CN99 radius',
            'exposed after original frozen scores','no independently acquired withheld ports')
    for r_id, rs in identified(all_cases,'loss4').items():
        width = float(rs[0]['width_m'])
        add(r_id,'Fixed-loss-angle materials','Same-phase Born and full-wave dissipative cap comparisons',
            'synthetic real-part geometry',[r['id'] for r in rs],[ROOT/r['source'] for r in rs[:1]],
            '12 × 12',f'400 MHz, {width:g} m square, 10 plane waves, 20 ring receivers; 2% noise and paired loss tangents .01, .05, .20',
            'exposed as fixed independent geometries','offline-only offset simulated illuminations and receivers')
    for r_id, rs in identified(all_cases,'v3_weak_controls').items():
        add(r_id,'Archived weak/Born control','Positive-fit or fixed-weight sensitivity control',
            'deterministic synthetic material',[r['id'] for r in rs],[ROOT/rs[0]['source']],
            '96 × 96','400 MHz, 2 m square, 10 plane waves and 20 full-ring receivers; paired 1% and 2% noise; double-disk generated with Born and full wave',
            'archived controls; repeated forward model uses the same material',
            'offline synthetic offset ports in the full-wave generating conditions')
    for r_id, rs in identified(all_cases,'new_mild_mid').items():
        add(r_id,'Weak-to-moderate routing controls','Born versus full-wave initialization; two cases also tested with CIE neighborhood rules',
            'deterministic synthetic material',[r['id'] for r in rs],[ROUTING/'mild_mid_public.npz',ROUTING/'FREEZE_ZH.md'],
            '12 × 12','400 MHz, 0.6 m square, 10 plane waves and 20 ring receivers; 2% noise',
            'predefined control shapes in this study; exploratory follow-up reused them',
            'offline-only offset simulated illuminations and receivers')
    lab = read_csv(ROUTING/'measured_port_cases.csv')
    grouped = defaultdict(list)
    for r in lab:
        grouped[(r['id'].split('/')[0],r['object_id'])].append(r)
    for (source,r_id), rs in grouped.items():
        fresnel = source=='fresnel_2005'
        family = 'Fresnel chamber' if fresnel else 'Manitoba chamber'
        acq = ('4, 6, 8 GHz; original source-specific masks with 241 real receivers per source' if fresnel else
               '3–6 GHz at 0.5 GHz spacing; original calibrated odd/even port groups and a fixed training/holdout split')
        add(f'{source}/{r_id}',family,'Calibrated Born and conditional full-wave cap comparisons',
            'measured physical object',[r['id'] for r in sorted(rs,key=lambda x:float(x['frequency_GHz']))],
            [ROUTING/'measured_port_cases.csv',OPTICAL/'physical_revision/born/ports/results.json'],
            '60 × 60' if fresnel else '40 × 40',acq,
            'historically exposed; loaded uncertainty is conditional',
            'actual acquired source and receiver holdouts used only after fitting')


def candidate_census():
    pool = STUDY/'pool'
    f = pool/'flow_results/online.json'
    for ident, receipt in read(f)['cases'].items():
        add(ident.removesuffix('__rx16'),'Archived learned-prior tests',
            'Native Flow Matching proposals and bounded exact-operator refinement','synthetic material',[ident],[f],
            '32 × 32','normalized wavelength 1; two opposing plane waves, native 16 of 32 receiver channels; clean measurement, 1% acceptance ball',
            'archived test objects of a trained model; exposed during this study',
            'other 16 saved receivers read after sealing')
    f = pool/'structured_results/online.json'
    development = STUDY/'threshold/development.json'
    for ident in read(f)['cases']:
        add(ident.split('__')[0],'Structured branch development',
            'Shifted branch gates, support ablation and selection of the fast screening rule','synthetic material',
            [ident],[f,development],'24 × 24','Window two wavelengths wide; 4 plane waves, 16 receivers; 1.25% noise, a declared 49-cell support',
            'development data for the prospective one-threshold comparison',
            'additional eight source angles and 16 receivers read offline')
    frozen = STUDY/'threshold/confirmation/public/cases.json'
    score = STUDY/'threshold/confirmation/offline.json'
    records = read(frozen)
    assert set(read(score)) == {r['id'] for r in records}
    for row in records:
        ident = row['id']
        add(ident.split('__')[0],'Prospective branch confirmation',
            'Pre-frozen one-threshold rule versus full and point-ranked refinement','synthetic material',[ident],
            [frozen,score,STUDY/'threshold/HOLDOUT_RULES.md'],'24 × 24',
            'Same 49-cell support, 4 plane waves, 16 receivers and 1.25% noise as development; independent material seed',
            'new materials generated after decision rule was fixed',
            'additional eight source angles and 16 receivers; reference opened after prediction seal',confirmed=True)


def physical_controls():
    source = STUDY/'families/results.json'
    for row in read(source)['controls']:
        ident = row['id']
        add(ident.split('__')[0],'Spatial branch controls',
            'Full native-domain range, exact physical site and background rejection tests',
            'deterministic synthetic material',[ident],[source,STUDY/'families'/f'{ident}.npz'],
            '96 × 96','400 MHz, 2 m square, native 32 plane waves and 32 full-ring receivers, 2% noise',
            'deterministic site control with independent material pattern and noise', 'none')
    box = STUDY/'boxes/results'
    a = read(box/'verified.json')['cases']
    small = next(key for key in a if key.startswith('shiftbox32'))
    disk = next(key for key in a if key.startswith('control_disk075'))
    add(small.split('__')[0],'Translated-box construction','Local exclusion on a new same-grid nonuniform material',
        'deterministic synthetic material',[small],[box/f'{small}.json',box/'verified.json'],'32 × 32',
        '400 MHz, 2 m square, 10 plane waves, 20 full-ring receivers, 2% noise',
        'new deterministic material; known generation map scored offline','none')
    original = OPTICAL/'ip_unified/numerics/feasibility/strict_endpoints'
    add(disk.split('__')[0],'Native low-peak disk control','High-peak shifted box exclusion and feasible lower-peak control',
        'deterministic synthetic material',[disk],[original/f'{disk}.json',box/'verified.json'],
        '96 × 96','400 MHz, 2 m square, 32 plane waves, 32 full-ring receivers, 2% noise',
        'archived physical disk; exposed control','none')
    r10 = 'r10_twod_0001'
    archive = SECANTS/'highdim/results'/f'{r10}.npz'
    with np.load(archive,allow_pickle=False) as saved, np.load(RING/'public_inputs.npz',allow_pickle=False) as rings:
        if any(np.array_equal(y,saved['Y_archive']) for y in rings['Y']):
            raise ValueError('The archived ring duplicates a compact-transfer observation.')
    add(r10,'Archived annular local control','Shifted cap test and older paired candidate endpoints on the same observation',
        'synthetic material',[r10],[archive,box/'verified.json'],'24 × 24',
        '12 cylindrical inputs, 48 ring receivers at k=2π; archived CN99 ball',
        'archived and previously exposed; separate from the seven compact transfer objects', 'none')
    archival = box/'ambiguity.json'
    pair = read(archival)
    ident = pair['observation_id']
    add(ident,'Historical ambiguity pair',
        'Two distinct full-wave material roots verified on one archived observation',
        'one historical observation pair',[ident],[archival],
        'Six-cell model, six independent material parameters',
        'Three illuminations and one complex receiver; shared archived data with the stated joint radius',
        'historical constructed control; both roots rechecked in this study','none',independently_counted=False)


def figure(primary, history):
    blue,green,orange = '#24607a','#1b7855','#a34f36'
    ink,muted,grid = '#202a31','#52616c','#dfe6e9'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.0,'axes.titlesize':9.5,
                         'axes.labelsize':8.0,'axes.spines.top':False,'axes.spines.right':False,
                         'text.color':ink,'axes.labelcolor':ink,'xtick.color':muted,'ytick.color':muted,
                         'pdf.fonttype':42})
    fig=plt.figure(figsize=(7.3,5.55),layout='constrained')
    layout=fig.add_gridspec(2,2,height_ratios=[1.14,1.])
    ax=fig.add_subplot(layout[0,:]); ring=fig.add_subplot(layout[1,0]); loss=fig.add_subplot(layout[1,1])
    for panel in (ax,ring,loss):
        panel.grid(axis='y',lw=.55,color=grid)
        panel.axhline(1.,color=muted,ls='--',lw=.85,zorder=2)
        panel.set_axisbelow(True)
        panel.set_ylabel('Distance lower bound / noise radius')
    x=np.arange(1,len(primary)+1)
    ax.axvspan(.55,6.5,color='#edf4f7',lw=0,zorder=0)
    for col,label,marker,color in [('l1_lower_over_radius','First Born, nonnegative','o',blue),
                               ('l2_cap05_distance_over_radius',r'Full wave, peak $\leq 0.50$','^',green),
                               ('l2_cap075_distance_over_radius',r'Full wave, peak $\leq 0.75$','v',orange)]:
        ax.scatter(x,[float(r[col]) for r in primary],marker=marker,s=18,c=color,label=label,zorder=3)
    ax.legend(frameon=False,loc='center right',fontsize=7.5)
    ax.set(title='A   Full-domain scattering, 28 independent materials',xlabel='Material (original six shaded; transfer cohort follows)',
           xlim=(.55,28.45),ylim=(-.19,6.83),xticks=[1,6,11,16,21,28])
    rr=read(RING/'summary.json')['frozen_r10_profile']['cases']
    labels=list(rr)
    ratios=[rr[k]['joint_interval'][0]/rr[k]['radius'] for k in labels]
    ring.bar(np.arange(1,8),ratios,color=[green if v>1 else '#b2bec5' for v in ratios],width=.58,zorder=3)
    ring.scatter([i+1 for i,v in enumerate(ratios) if v==0],[0 for v in ratios if v==0],
                 facecolor='white',edgecolor=muted,s=15,zorder=4)
    control=labels.index('r14_new_0018')+1
    ring.set(title=r'B   Compact annular media, peak $\leq 5$',xlabel='Material',
             xticks=range(1,8),ylim=(-.36,12.1))
    labels_at_bottom=[str(i) if i!=control else f'{i}\nFeasible peak 4.28' for i in range(1,8)]
    ring.set_xticklabels(labels_at_bottom)
    index={r['id']:float(r['distance_over_radius']) for r in history if float(r['cap_upper'])==1.}
    for n in range(1,5):
        marker='o' if n%2 else 's'
        color=blue if n<=2 else orange
        v=[index[f'loss_tm_{n:02d}__fixed_t{t}__pair2pct'] for t in ('0.01','0.05','0.2')]
        loss.plot([0,1,2],v,color=color,marker=marker,markersize=4.2,lw=1.1,zorder=3)
    loss.set_xticks([0,1,2],['0.01','0.05','0.20'])
    loss.set(title=r'C   Fixed loss angle, real peak $\leq 1$',xlabel=r'Loss tangent $t$',ylim=(-.05,13.5))
    loss.legend(handles=[Line2D([],[],color=blue,label='0.60 m domain'),Line2D([],[],color=orange,label='1.25 m domain'),
                         Line2D([],[],color=muted,marker='o',linestyle='',label='Single ellipse',markersize=4),
                         Line2D([],[],color=muted,marker='s',linestyle='',label='Ellipse and disc',markersize=4)],
                frameon=False,fontsize=6.4,ncol=2,loc='upper center',columnspacing=.9,handlelength=1.4)
    out=STUDY/'figures'
    out.mkdir(exist_ok=True)
    fig.savefig(out/'benchmark.pdf')
    fig.savefig(out/'benchmark.png',dpi=270)


def main():
    primary = read_csv(ROUTING/'principal_28.csv')
    history = read_csv(ROUTING/'historical_l2.csv')
    assert len(primary)==28 and len({r['object_id'] for r in primary})==28
    routing_census(primary,read_csv(ROUTING/'all_routing_cases.csv'))
    candidate_census()
    physical_controls()
    assert len(CENSUS)==len({r['object_id'] for r in CENSUS})
    assert sum(r['primary_native'] for r in CENSUS)==28
    assert sum(r['prospective_confirmation'] for r in CENSUS)==8
    all_ids = [x for row in CENSUS for x in row['observation_ids'].split(' | ')]
    assert len(all_ids)==len(set(all_ids)), 'A saved observation is counted twice'
    with (HERE/'objects.csv').open('w',newline='') as stream:
        writer=DictWriter(stream,fieldnames=list(CENSUS[0]))
        writer.writeheader()
        writer.writerows(CENSUS)
    figure(primary,history)
    print('Distinct materials and measured physical objects',sum(r['independently_counted'] for r in CENSUS))
    print('Historical observation controls',sum(not r['independently_counted'] for r in CENSUS))
    print('Observation views',len(all_ids))
    print(json.dumps(dict(Counter(row['family'] for row in CENSUS)),indent=2))


if __name__=='__main__': main()

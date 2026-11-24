"""Render the journal figures from saved material maps and endpoint receipts."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, StrMethodFormatter
import numpy as np

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]
OUT = HERE.parent/'figures'
INK = '#182832'
TEXT = '#293d46'
GRID = '#dce4e8'
OLD = '#526674'
NEW = '#087b72'
UPPER = '#a35536'
SECOND = '#395b9a'
LIGHT = '#788d9b'
SIZES = (10, 16, 24, 32)
LABELS = ('Slotted block', 'Four lobes', 'Double chevron', 'Crescent and plate',
          'Staggered grooves', 'Graded pinwheel')
SOURCE_ENUM = tuple(f'optfresh_{i:02d}' for i in range(1, 7))
PAPER_WIDTH = 5.88


mpl.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8, 'axes.titlesize':8.5,
    'axes.titleweight':'bold', 'axes.edgecolor':LIGHT, 'axes.labelcolor':TEXT,
    'axes.linewidth':.65, 'text.color':INK, 'xtick.color':TEXT, 'ytick.color':TEXT,
    'xtick.labelsize':7, 'ytick.labelsize':7, 'xtick.major.size':2.5, 'ytick.major.size':2.5,
    'xtick.major.pad':2.5, 'ytick.major.pad':2.5,
    'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42, 'ps.fonttype':42,
    'savefig.facecolor':'white', 'legend.frameon':False, 'legend.fontsize':7})


def load(path):
    return json.loads(path.read_text())


def study_load(relative):
    return load(STUDY/relative)


def load_numeric(relative):
    return load(HERE/relative)


def letter(ax, label, suffix):
    ax.set_title(f'({label})  {suffix}', loc='left', pad=6)


def x_grid(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='x', color=GRID, lw=.65)


def y_grid(ax):
    ax.set_axisbelow(True)
    ax.grid(axis='y', color=GRID, lw=.65)


def threshold(ax, vertical=True):
    method = ax.axvline if vertical else ax.axhline
    method(1, color=INK, lw=.9, ls=(0,(3,3)), zorder=1)


def object_name(observation):
    return observation.split('__',1)[0]


def marker(label, symbol, color, filled=True):
    return Line2D([],[], ls='', marker=symbol, color=color, mfc=color if filled else 'white',
                  label=label, markersize=5.5, markeredgewidth=1)


def save(fig, name):
    fig.savefig(OUT/f'{name}.pdf', metadata={'Creator':'plot_paper.py', 'CreationDate':None})
    fig.savefig(OUT/f'{name}.png', dpi=260)
    plt.close(fig)


def geometry():
    provenance = 'physical_revision/array/summary.json'
    summary = study_load(provenance)
    fig = plt.figure(figsize=(PAPER_WIDTH, 2.92))
    region = fig.add_gridspec(2,4,width_ratios=[1,1,1,.085], left=.093,right=.548,
                              bottom=.158,top=.824,wspace=.34,hspace=.38)
    material_meta = {}
    extent = [-1.,1.,-1.,1.]
    wrapped_labels = ('Slotted block', 'Four lobes', 'Double\nchevron', 'Crescent\nand plate',
                      'Staggered\ngrooves', 'Graded\npinwheel')
    for i, (object_id, label, displayed_label) in enumerate(zip(SOURCE_ENUM,LABELS,wrapped_labels)):
        row,col = divmod(i,3)
        ax = fig.add_subplot(region[row,col])
        path = f'synthetic/evaluation/{object_id}.npz'
        with np.load(STUDY/path,allow_pickle=False) as data:
            m = data['m_fine'].reshape(96,96)
            material_meta[label] = {'source':path, 'peak':float(m.max()), 'minimum':float(m.min())}
        image = ax.imshow(m, extent=extent, origin='lower', vmax=4, vmin=0, cmap='viridis', interpolation='nearest')
        ax.set_title(displayed_label, fontsize=7, fontweight='normal', pad=3, linespacing=1.05)
        ax.set_xticks((-1.,0.,1.) if row==1 else ())
        ax.set_yticks((-1.,0.,1.) if col==0 else ())
        ax.tick_params(labelsize=6.7, length=2, pad=2)
        ax.spines[['top','right']].set_visible(True)
        if row==1 and col==1:
            ax.set_xlabel('Horizontal position (m)', labelpad=3)
    fig.text(.028,.49,'Vertical position (m)', rotation=90, va='center', ha='center')
    color_ax = fig.add_subplot(region[:,3])
    bar = fig.colorbar(image,cax=color_ax,ticks=(0,1,2,3,4))
    bar.ax.set_title(r'$\chi$', fontsize=8.5, pad=5)
    bar.outline.set_linewidth(.45)
    counts_ax = fig.add_axes((.695,.185,.28,.645))
    data = {}
    for receivers, color, symbol, linestyle in [(20,SECOND,'o','-'),(32,NEW,'s','--')]:
        data[str(receivers)] = [summary['counts'][str(receivers)]['2'][str(tx)]['excluded_at_075'] for tx in SIZES]
        counts_ax.plot(SIZES,data[str(receivers)], marker=symbol, color=color, ls=linestyle,
                       fillstyle='none' if receivers==32 else 'full', lw=1.25, ms=4.5,
                       label=f'{receivers} receivers')
    counts_ax.set(xlim=(8.4,33.5),ylim=(0,6.3),xticks=SIZES,yticks=range(7),
                  xlabel='Number of illuminations',ylabel=r'Materials excluding $u=0.75$')
    counts_ax.legend(loc='lower right')
    y_grid(counts_ax)
    fig.text(.093,.944,'(a)  Six generating materials',fontsize=8.5,fontweight='bold')
    fig.text(.695,.944,'(b)  Illumination sweep',fontsize=8.5,fontweight='bold')
    save(fig,'geometry')
    return dict(acquisition_count_source=provenance, illuminations=list(SIZES), count_by_receivers=data,
                equal_total_fluence=True, tested_cap=.75, noise_percent=2, grid=[96,96], side_m=2,
                material_common_color_scale=[0,4], material_plots=material_meta)


def historical_cap_max():
    with (STUDY/'physical_revision/array/array_scan.csv').open(newline='') as f:
        parsed = list(csv.DictReader(f))
    data = {}
    for row in parsed:
        if (row['sources'],row['receivers'],row['noise_pct']) != ('32','32','2'):
            continue
        cap = float(row['tested_cap'])
        if row['strict_profile_distance_lower']:
            data.setdefault(cap,[]).append(float(row['strict_profile_distance_lower'])/float(row['radius']))
    return {u:max(values) for u,values in data.items()}


def full_domain():
    original = load_numeric('native/u100_best_strict.json')
    coupled = load_numeric('source/u100_best_strict_source.json')
    scalar = historical_cap_max()
    fine_grid_sources = sorted((HERE/'source').glob('cap_resolution*.json'))
    refined = {}
    for source in fine_grid_sources:
        for r in load(source)['rows']:
            if r['profile']['name'] != 'stable_scalar_optical':
                continue
            ratios = [o['distance_lower']/o['radius'] for o in r['observations']
                      if o.get('distance_lower') is not None]
            val = dict(max_lower_over_radius=max(ratios,default=None),
                       exclusions=sum(x>1 for x in ratios), available_lower_bounds=len(ratios),
                       source=str(source.relative_to(HERE)))
            previous = refined.get(r['cap'])
            if previous is None or val['available_lower_bounds'] >= previous['available_lower_bounds']:
                refined[r['cap']] = val
    joints = {object_name(r['id']):r for r in coupled['rows']}
    tight = load_numeric('native/deflated_sparse_tightening.json')
    limit = {object_name(r['id']):r for r in tight['source_upper_bounds'] if r['row_budget']==52}
    old = {object_name(r['id']):r for r in original['outward_centers']}
    fig = plt.figure(figsize=(PAPER_WIDTH, 4.62))
    a = fig.add_axes((.128, .651, .846, .288))
    b = fig.add_axes((.282, .190, .359, .276))
    c = fig.add_axes((.708, .190, .266, .276))
    cap_positions = sorted(set(scalar)|set(refined))
    positions = {value:i for i,value in enumerate(cap_positions)}
    a.scatter([positions[x] for x in scalar],list(scalar.values()),color=OLD,marker='o',s=28,
              label='Scalar-profile envelope',zorder=3)
    refined_points = [(u,val['max_lower_over_radius']) for u,val in refined.items()
                      if val['max_lower_over_radius'] is not None]
    unresolved_caps = {u for u,val in refined.items() if val['max_lower_over_radius'] is None}
    if refined_points:
        a.scatter([positions[u] for u,v in refined_points],[v for u,v in refined_points],color=NEW,marker='D',s=28,
                  label='Single fixed scalar profile',zorder=3)
    a.set(xticks=range(len(cap_positions)),
          xticklabels=[('\n' if i%2 else '')+f'{u:g}'+('†' if u in unresolved_caps else '')
                       for i,u in enumerate(cap_positions)],
          ylim=(-.075,2.37),xlabel=r'Tested contrast cap $u$',ylabel=r'Largest lower endpoint / $r$')
    a.tick_params(axis='x',labelsize=7)
    letter(a,'a','Lower bounds across caps')
    threshold(a,False)
    y_grid(a)
    a.legend(loc='upper right',borderaxespad=.35,handletextpad=.45)
    if unresolved_caps:
        a.text(.0,-.41,'† Fourier inertia gap unresolved', transform=a.transAxes,fontsize=7,color=OLD)
    y = np.arange(len(LABELS))
    per_case = []
    for j,key in enumerate(SOURCE_ENUM):
        d0 = old[key]['distance_lower']/old[key]['radius']
        d1 = joints[key]['ratio_lower']
        upper = limit[key]['ratio_upper']
        b.plot([d0,d1],[j-.13,j+.13],color=LIGHT,lw=.8,zorder=2)
        b.scatter(d0,j-.13,marker='o',s=29,facecolor='white',edgecolor=OLD,zorder=3)
        b.scatter(d1,j+.13,marker='D',s=24,facecolor=NEW,edgecolor=NEW,zorder=3)
        b.scatter(upper,j,marker='v',s=46,facecolor='white',edgecolor=UPPER,zorder=3)
        per_case.append(dict(object_id=key,name=LABELS[j],inertia_lower_over_radius=d0,
                             joint_constitutive_lower_over_radius=d1, selected_52_row_method_upper_over_radius=upper))
    b.set(yticks=y,yticklabels=LABELS,xlim=(0,1.06),ylim=(len(y)-.45,-.58),
          xlabel=r'Endpoint / $r$ at $u=1$')
    b.tick_params(axis='y',pad=4,length=0)
    b.set_xticks((0,.25,.5,.75,1.))
    threshold(b)
    x_grid(b)
    letter(b,'b','Fixed-weight bounds')
    b.legend(handles=[marker('Inertia lower','o',OLD,False),marker('Coupled lower','D',NEW),
                      marker('≤52-pixel ceiling','v',UPPER,False)],
             loc='upper left',ncol=1,bbox_to_anchor=(-.025,-.295),
             handletextpad=.45,borderaxespad=0)
    physical = []
    for j,key in enumerate(SOURCE_ENUM):
        for cap,color,symbol,dither in [(1.,SECOND,'s',-.12),(3.,UPPER,'^',.12)]:
            paths = list((HERE/'feasibility/strict_endpoints').glob(f'{key}*__u{cap:g}.json'))
            if paths:
                result = load(paths[0])
                ratio = result.get('distance_upper_over_radius')
                if ratio is not None:
                    c.scatter(ratio,j+dither,marker=symbol,s=32,facecolor='white',edgecolor=color,zorder=3)
                    physical.append(dict(object_id=key,name=LABELS[j],cap=cap,ratio=ratio,
                                         strict_candidate_source=str(paths[0].relative_to(HERE))))
    c.set(xscale='log',xlim=(.92,60),ylim=(len(y)-.45,-.58),yticks=y,yticklabels=[],
          xlabel=r'Forward residual upper / $r$')
    c.xaxis.set_major_locator(FixedLocator((1,3,10,30,60)))
    c.xaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    threshold(c)
    x_grid(c)
    letter(c,'c','Tested materials')
    c.legend(handles=[marker(r'$u=1$ physical upper','s',SECOND,False),
                      marker(r'$u=3$ physical upper','^',UPPER,False)],
             loc='upper left',bbox_to_anchor=(-.025,-.295),borderaxespad=0,
             handletextpad=.45)
    save(fig,'full-domain')
    return dict(old_lower_sources=['physical_revision/array/array_scan.csv','native/u100_best_strict.json'],
                new_lower_source='source/u100_best_strict_source.json', 
                refined_cap_sources=[str(f.relative_to(HERE)) for f in fine_grid_sources],
                historical_lower_over_radius_by_cap=scalar, fixed_scalar_profile_by_cap=refined, endpoint_u1_by_material=per_case,
                physical_trial_upper_over_radius=physical,
                limit_scope='Any set of up to 52 individually weighted constitutive rows, fixed state weights; signed source center.',
                upper_evidence='native/deflated_sparse_tightening.json')


def panel_loss_sample(document, tangent, cap):
    return next(block['rows'] for block in document['projected'] if block['tangent']==tangent and block['cap']==cap)


def state_ratio(row, name):
    return row['state'][name]['distance_lower']/row['rho']


def dissipation():
    finite = [load_numeric(f'loss/finite_{key}_verified.json') for key in (0,1)]
    scalar = load_numeric('loss/finite_0_pure_shift_verified.json')
    finite_state = load_numeric('loss/finite_0_state_verified.json')
    native = load_numeric('loss/native_verified.json')
    native_state = load_numeric('loss/native_state_verified.json')
    directional = load_numeric('loss/finite_0_directional.json')
    svd_rows = next(r['rows'] for r in directional['profiles'] if (r['tangent'],r['cap'])==(.2,1.))
    svd = {r['id']:r['methods']['right_singular']['verified']['distance_lower']/r['rho'] for r in svd_rows}
    fig = plt.figure(figsize=(PAPER_WIDTH, 5.92))
    a = fig.add_axes((.232, .741, .727, .215))
    b = fig.add_axes((.232, .405, .727, .151))
    c = fig.add_axes((.232, .115, .727, .171))
    finite_points = []
    for index,(doc,side) in enumerate(zip(finite,(.6,1.25))):
        rows = panel_loss_sample(doc,.2,1.)
        for j,row in enumerate(x for x in rows if not object_name(x['id']).endswith('_low')):
            k = 2*index+j
            r = row['rho']
            ellipsoid_lower = row['dual_lower']/r
            upper = row['primal']['primal_upper']/r
            a.plot([ellipsoid_lower,upper],[k,k],color=NEW,lw=1.4)
            a.scatter(ellipsoid_lower,k,marker='D',color=NEW,s=28,zorder=3)
            a.scatter(upper,k,marker='|',color=NEW,s=70,zorder=3)
            point = dict(electrical_domain_side_m=side,object_id=object_name(row['id']),
                         ellipse_lower=ellipsoid_lower,ellipse_upper=upper)
            if row.get('scalar'):
                scalar_upper=row['scalar']['scalar_primal_upper']/r
                a.scatter(scalar_upper,k+.12,marker='v',s=31,facecolor='white',edgecolor=SECOND,zorder=3)
                point['scalar_current_upper']=scalar_upper
            if row['id'] in svd:
                source_lower = svd[row['id']]
                a.scatter(source_lower,k+.18,marker='^',s=29,facecolor=SECOND,edgecolor=SECOND,zorder=3)
                point['right_singular_source_prices_lower'] = source_lower
            if row['id']==scalar['id']:
                shift_upper=scalar['all_tau_distance_upper']/r
                a.scatter(shift_upper,k-.14,marker='v',s=36,facecolor='white',edgecolor=UPPER,zorder=4)
                point['all_scalar_shifts_upper']=shift_upper
            finite_points.append(point)
        control = next(x for x in rows if object_name(x['id']).endswith('_low'))
        control_r = control['rho']
        low_peak = next(x for x in doc['inputs'] if x['id']==control['id'])
        physical = low_peak['forward_residual_upper']/control_r
        a.scatter(physical,4+index,marker='x',s=32,color=OLD,zorder=3)
        finite_points.append(dict(electrical_domain_side_m=side,object_id=object_name(control['id']),
                                  physical_control_upper=physical))
    a.set(yticks=range(6),yticklabels=['0.60 m · A','0.60 m · B','1.25 m · A','1.25 m · B',
           '0.60 m · low peak','1.25 m · low peak'], ylim=(5.44,-.53),xlim=(0,1.30),
          xticks=(0,.25,.5,.75,1,1.25),xlabel=r'Endpoint / $r$ at $t=0.20,\ u=1$')
    a.tick_params(axis='y',length=0,pad=4)
    threshold(a)
    x_grid(a)
    letter(a,'a','Finite TM grids')
    a.legend(handles=[marker('Receiver lower','D',NEW),marker('Scalar shift upper','v',UPPER,False),
                      marker('SVD source-price lower','^',SECOND),marker('Current cap upper','v',SECOND,False),
                      marker('Physical control','x',OLD)],
             loc='upper left',ncol=3,bbox_to_anchor=(-.025,-.33),columnspacing=.85,
             handletextpad=.35,borderaxespad=0)
    shift_comparison = next(row for row in finite_state['rows'] if row['id']==scalar['id']
                            and row['tangent']==.2 and row['cap']==1.25)
    comparisons = [('old','Passive sector',OLD),('fixed_phase','Fixed phase',SECOND),
                   ('shift','Sector + shift',NEW)]
    comp_values = {key:state_ratio(shift_comparison,key) for key,_,_ in comparisons}
    partition_source = 'loss/source_partition.json'
    if (HERE/partition_source).is_file():
        partition_run = next(row for row in load_numeric(partition_source)['runs']
                             if row['id']==shift_comparison['id'] and row['cap']==shift_comparison['cap']
                             and row['tangent']==shift_comparison['tangent'])
        partition = partition_run['frozen_selected']
        split_lower = partition['outward_selection']['distance_interval'][0]
        whole_lower, whole_upper = partition['all_sector']['total_sector_distance_interval']
        if (partition_run['radius'] != shift_comparison['rho']
                or whole_lower != shift_comparison['state']['shift']['distance_lower']
                or split_lower <= whole_upper):
            raise ValueError('Split illumination comparison requires a strict gain at the same observation and weights.')
        comparisons.append(('split','Split illuminations',NEW))
        comp_values['split'] = split_lower / partition_run['radius']
    for i,(key,label,color) in enumerate(comparisons):
        val = comp_values[key]
        b.barh(i,val,height=.64,color=color)
        b.text(val-.055,i,f'{val:.3f}',va='center',ha='right',fontsize=8,color='white')
    b.set(yticks=range(len(comparisons)),yticklabels=[name for _,name,_ in comparisons],xlim=(0,1.53),
          ylim=(len(comparisons)-.49,-.51),xlabel=r'Lower endpoint / $r$ at $t=0.20,\ u=1.25$')
    b.tick_params(axis='y',length=0,pad=4)
    threshold(b)
    x_grid(b)
    letter(b,'b','0.60 m body A')
    native_points = []
    jitter = np.linspace(-.14,.14,6)
    for j,tangent in enumerate((.01,.05,.2)):
        rows = panel_loss_sample(native,tangent,.75)
        principal = [r for r in rows if not object_name(r['id']).endswith('_low')]
        for row,offset in zip(principal,jitter):
            hi = row['primal']['primal_upper']/row['rho']
            c.scatter(hi,j+offset,marker='v',s=19,facecolor='white',edgecolor=UPPER,zorder=4,linewidth=.8)
            native_points.append(dict(object_id=object_name(row['id']), tangent=tangent,cap=.75,
                 receiver_lower=row['dual_lower']/row['rho'],receiver_upper=hi))
        control = next(r for r in rows if object_name(r['id']).endswith('_low'))
        control_ratio = control['primal']['primal_upper']/control['rho']
        c.scatter(control_ratio,j,marker='x',s=25,color=OLD,zorder=4)
        native_points.append(dict(object_id=object_name(control['id']), tangent=tangent,cap=.75,
                                  receiver_upper=control_ratio, low_peak_control=True))
    replay_points = []
    for j,case in enumerate(native_state['cases'],3):
        for row in case['rows']:
            result = row['result']
            forward = next(obs for obs in native['inputs'] if obs['id']==row['id'])
            ratio = result['distance_lower']/forward['rho']
            c.scatter(ratio,j,marker='D',color=NEW,s=28,zorder=3)
            name = LABELS[SOURCE_ENUM.index(object_name(row['id']))]
            c.text(ratio+.035,j,name,va='center',fontsize=7,color=TEXT)
            replay_points.append(dict(object_id=object_name(row['id']),arm=case['arm'],tangent=case['tangent'],
                                      cap=case['cap'],ratio_lower=ratio))
    c.set(yticks=range(5), yticklabels=['$t=0.01$','$t=0.05$','$t=0.20$','Sector + shift','Passive sector'],
          xlim=(0,1.115),ylim=(4.48,-.50),xticks=(0,.25,.50,.75,1),xlabel=r'Endpoint / $r$')
    c.tick_params(axis='y',length=0,pad=4)
    c.axhline(2.5,color=LIGHT,lw=.65)
    c.text(.72,.88,r'$u=0.75$',fontsize=7.5,ha='center',color=TEXT)
    c.text(.72,3.5,r'$u=1$, $t=0.20$',fontsize=7.5,ha='center',color=TEXT)
    threshold(c)
    x_grid(c)
    letter(c,'c','Native 96 × 96 grid')
    c.legend(handles=[marker('Receiver upper','v',UPPER,False),marker('Low-peak control','x',OLD),
                      marker('State lower','D',NEW)],loc='upper left',ncol=3,
             bbox_to_anchor=(-.025,-.43),columnspacing=.95,handletextpad=.35,borderaxespad=0)
    save(fig,'dissipation')
    return dict(source_files=['loss/finite_0_verified.json','loss/finite_1_verified.json','loss/native_verified.json',
              'loss/finite_0_state_verified.json','loss/finite_0_pure_shift_verified.json','loss/native_state_verified.json',
              'loss/finite_0_directional.json'] + ([partition_source] if 'split' in comp_values else []),
              finite_u1_t02=finite_points,finite_060m_A_u125_t02_state_lowers=comp_values,
              native_receiver_points=native_points,native_state_replays=replay_points)


def current_modes():
    needed = load_numeric('native/required_currents.json')
    physical = load_numeric('feasibility/hidden_current.json')
    if physical['projector_profile'] != 'scalar_grid_601' or needed['frozen_observation_count'] != 6:
        raise ValueError('Current comparison requires the frozen native profiles and matching geometry.')
    minima = {object_name(r['id']):r['best']['negative_current_norm_lower'] for r in needed['rows']}
    attained = {object_name(r['observation']):r for r in physical['rows']}
    fig, ax = plt.subplots(figsize=(PAPER_WIDTH, 2.70))
    fig.subplots_adjust(left=.274,right=.962,bottom=.265,top=.886)
    rows = []
    for j, key in enumerate(SOURCE_ENUM):
        lower = minima[key]
        ax.barh(j,lower,height=.57,color=NEW,zorder=2)
        ax.text(lower-.55,j,f'{np.floor(100*lower)/100:.2f}',va='center',ha='right',fontsize=7,color='white')
        row = dict(name=LABELS[j], required_to_fit_original_radius_lower=lower)
        if key in attained:
            proof = attained[key]
            pair = proof['negative_current_frobenius_interval']
            ax.plot(pair,[j,j],color=UPPER,lw=1.25,zorder=4)
            midpoint = sum(pair)/2
            ax.scatter(midpoint,j,marker='D',s=36,facecolor='white',edgecolor=UPPER,zorder=5)
            ax.text(midpoint+.72,j,f'{midpoint:.2f}',va='center',ha='left',fontsize=7,color=UPPER)
            row.update(physical_norm_interval=pair, material_peak=proof['material_peak'],
                       physical_radius_status=proof['noise_ball_status'])
        rows.append(row)
    disk = attained['control_disk075']
    disk_interval = disk['negative_current_frobenius_interval']
    ax.plot(disk_interval,[6.35,6.35],color=SECOND,lw=1.25,zorder=4)
    ax.scatter(sum(disk_interval)/2,6.35,marker='s',s=33,facecolor=SECOND,edgecolor=SECOND,zorder=5)
    ax.text(disk_interval[1]+.72,6.35,f'{disk_interval[1]:.2f}',va='center',fontsize=7,color=SECOND)
    ax.axhline(5.65,color=LIGHT,ls=(0,(3,3)),lw=.68,zorder=1)
    ax.set(yticks=[*range(6),6.35],yticklabels=[*LABELS,'Independent low-peak disk'],
           ylim=(6.88,-.53),xlim=(0,47),xticks=range(0,46,5),
           xlabel=r'Unweighted negative-space current norm $\|\Pi_-J\|_F$')
    ax.tick_params(axis='y',length=0,pad=5)
    x_grid(ax)
    ax.set_title(r'Full-domain negative-space current at $u=1$',loc='left',pad=10)
    fig.legend(handles=[marker('Lower bound for any fit','s',NEW),
                        marker('Box trial outside its ball','D',UPPER,False),
                        marker('Control in its own ball','s',SECOND)],
               loc='lower center',bbox_to_anchor=(.55,.035),ncol=3,columnspacing=1.1,
               handletextpad=.40,borderaxespad=0)
    save(fig,'current-modes')
    return dict(lower_bound_source='native/required_currents.json',physical_source='feasibility/hidden_current.json',
                contract=physical['contract'], scope=needed['bound_scope'],units=physical['current_units'],
                rows=rows, independent_low_peak_disk=dict(material_peak=disk['material_peak'],
                  physical_norm_interval=disk_interval, radius_status=disk['noise_ball_status']))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--only',nargs='+',choices=['geometry','full-domain','dissipation','current-modes'])
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    recipes = {'geometry':geometry,'full-domain':full_domain,'dissipation':dissipation,'current-modes':current_modes}
    sidecar = OUT/'figure-data.json'
    values = load(sidecar) if args.only and sidecar.exists() else {}
    for name in args.only or recipes:
        values[name] = recipes[name]()
        print(name,flush=True)
    sidecar.write_text(json.dumps(values,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()

"""Rebuild the paper figure and machine-readable complete paired table."""

import csv
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from curve_inputs import CAPS, HERE
from search import save_json

NOISE=(0.5,1,2,5,10)
LABELS=('Slotted block','Overlapping lobes','Double-folded strip','Crescent and plate',
        'Interleaved channels','Graded pinwheel')
METHODS=('global_control','global_075','local_075','fixed_075')
COLORS=dict(global_control='#64778A',global_075='#106AA2',local_075='#C35122',fixed_075='#4B815D')
DISPLAY=dict(global_control='Global real, cap 0.50',global_075='Global, cap 0.75',
             local_075='2×2 tuned, cap 0.75',fixed_075='4×4 frozen, cap 0.75')
MARKERS=dict(global_control='s',global_075='o',local_075='D',fixed_075='^')


def read(name):
    return json.loads((HERE/name).read_text())


def receipt_for_rows(verified):
    return {row['id']:row for row in verified['rows']}


def state(d):
    if d is None:
        return '', '', 'no surviving proposal'
    numeric = d.get('distance_lower','') if d['status']=='native_distance_verified' else ''
    return numeric, d.get('exact_nu',''), d['status']


def draw(rows):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.unicode_minus':False,
                         'pdf.fonttype':42,'ps.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(10.5,6.5),sharex=True,sharey=True)
    for i,ax in enumerate(axes.ravel()):
        subset=[r for r in rows if r['cap']==.75 and r['object_id']==f'optfresh_{i+1:02d}']
        levels=[r['noise_pct'] for r in subset]
        for method in METHODS:
            d=np.array([r[method]/r['radius'] if r.get(method,'')!='' else np.nan for r in subset])
            # The zero floor is meaningful for the tested profile with no
            # positive surviving source direction.
            ax.plot(levels,d,color=COLORS[method],lw=1.4,alpha=.95,marker=MARKERS[method],
                    markersize=5,markeredgewidth=1.1,markerfacecolor='white')
            occupied=(d>1.) & np.isfinite(d)
            ax.scatter(np.array(levels)[occupied],d[occupied],s=24,c=COLORS[method],
                       marker=MARKERS[method],edgecolors='white',linewidths=.35,zorder=5)
        ax.axhline(1.,c='#232830',lw=.85,ls='--',zorder=0)
        ax.set_title(f'{i+1:02d}  {LABELS[i]}',loc='left',fontweight='bold',fontsize=10.8)
        ax.set_xscale('log')
        ax.set_xlim(.46,11.5)
        ax.set_yscale('symlog',linthresh=.07,linscale=.8)
        ax.set_ylim(-.016,78)
        ax.set_xticks(NOISE,labels=['0.5','1','2','5','10'])
        ax.set_yticks([0,.1,1,10,50],labels=['0','0.1','1','10','50'])
        ax.grid(axis='y',color='#D9E1E7',lw=.6)
        ax.tick_params(which='minor',bottom=False,left=False)
        ax.spines[['left','bottom']].set_color('#A2AEB8')
    fig.suptitle('Material-peak certificates under paired full-ring noise',
                 y=.975,x=.10,ha='left',fontweight='bold',fontsize=16,color='#162535')
    fig.supxlabel('Nominal noise (% of archived numerical clean field)',y=.088,fontsize=10.5)
    fig.supylabel('Certified residual lower bound / padded noise radius',x=.045,fontsize=10.5)
    legend=[Line2D([0],[0],marker=MARKERS[k],color=COLORS[k],markerfacecolor='white',
                  markersize=5,label=DISPLAY[k]) for k in METHODS]
    legend.append(Line2D([0],[0],color='#232830',lw=.85,ls='--',label='Exclusion threshold'))
    fig.legend(handles=legend,loc='lower center',bbox_to_anchor=(.53,.011),ncol=5,
               frameon=False,fontsize=8.8)
    fig.subplots_adjust(left=.105,right=.98,top=.88,bottom=.16,hspace=.38,wspace=.18)
    for suffix in ('png','pdf'):
        fig.savefig(HERE/f'full_ring_noise.{suffix}',dpi=260,facecolor='white')
    plt.close(fig)


def main():
    raw=read('search.json')
    a,b,c=(read(x) for x in ('global.5.json','global.75.json','local.75.json'))
    control,scalar,spatial=map(receipt_for_rows,(a,b,c))
    recovery_global=read('recover_global.75.json')
    recovery_local=read('recover_local.75.json')
    for method,recovery in ((scalar,recovery_global),(spatial,recovery_local)):
        for row in recovery['rows']:
            key=row['id']
            if (method.get(key,{}).get('status')!='native_distance_verified'
                    and row['status']=='native_distance_verified'):
                method[key]=row
    if any(len(x)!=30 for x in (control,scalar,spatial)):
        raise ValueError('All 30 rows are required for each replay arm.')
    fixed=read('fixed_blocks.75.json')
    fixed_high=read('fixed_blocks_higher_caps.json')
    if len(fixed['profiles'])!=2 or len(fixed_high['cases'])!=6:
        raise ValueError('Both fixed 4×4 profiles and six high-cap probes are required.')
    by_profile={p['profile_index']:receipt_for_rows(p) for p in fixed['profiles']}
    high_4x4 = {}
    for case in fixed_high['cases']:
        high_4x4.setdefault(case['cap'],{})[case['profile_index']]=case
    rows=[]
    for r in raw['rows']:
        cid,u=r['id'],r['cap']
        y=dict(id=cid,object_id=r['object_id'],noise_pct=r['noise_pct'],
               cap=u,radius=r['radius'],sampled_global_profiles=r['n_profiles'],
               min_screened_global_nu=r['minimum_inertia_profile']['nu'])
        ds=[arm['d'] for arm in r['arms'].values() if arm is not None]
        y['best_screened_global']=max(ds,default=0.)
        for name in ('real_blocks','optical_blocks','joint_blocks'):
            proposals=[z for l in r['local'].values() for z in l['trials']
                       if z['arm']==name and z['same_inertia']]
            y['best_screened_'+name]=max((p['d'] for p in proposals),default=0.)
            y['tested_valid_'+name]=len(proposals)
        if u==.5:
            choice=control[cid]
        elif u==.75:
            choice=scalar[cid]
        else:
            choice=None
        y['strict_global_distance_lower'],y['strict_global_nu'],y['strict_global_status']=state(choice)
        y['strict_local_distance_lower'],y['strict_local_nu'],y['strict_local_status']=state(
            spatial[cid] if u==.75 else None)
        y['strict_global_exclusion']=bool(choice and choice.get('strictly_excludes'))
        y['strict_local_exclusion']=bool(u==.75 and spatial[cid].get('strictly_excludes'))
        y['global_control']=state(control[cid])[0] if u==.75 else ''
        y['global_075']=state(scalar[cid])[0] if u==.75 else ''
        y['local_075']=state(spatial[cid])[0] if u==.75 else ''
        y['fixed_075']=''
        y['fixed_4x4_exclusion']=False
        y['fixed_4x4_coherent_exclusion']=False
        for profile in sorted(by_profile):
            current=by_profile[profile][cid] if u==.75 else None
            high_case=high_4x4.get(u,{}).get(profile,{})
            y[f'fixed_{profile}_distance_lower']=current['distance_lower'] if current else ''
            y[f'fixed_{profile}_coherent_distance_upper']=current['coherent_distance_upper'] if current else ''
            y[f'fixed_{profile}_negative_count_lower']=current['nu'] if current else high_case.get(
                'native_negative_lower_count','')
            y[f'fixed_{profile}_numeric_negative_count']=current['nu'] if current else high_case.get(
                'numeric_negative_count','')
            if current:
                y['fixed_075']=max(current['distance_lower'],y['fixed_075'] or 0.)
                y['fixed_4x4_exclusion'] |= current['strictly_excludes']
                y['fixed_4x4_coherent_exclusion'] |= current['coherent_excludes']
        y['global_replay_method']=choice.get('method','') if choice else ''
        y['local_replay_method']=spatial[cid].get('method','') if u==.75 else ''
        rows.append(y)
    with (HERE/'noise_curve.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    counts=[]
    for u in CAPS:
        for pct in NOISE:
            subset=[r for r in rows if r['cap']==u and r['noise_pct']==pct]
            counts.append(dict(cap=u,noise_pct=pct, objects=6,
                               global_exclusions=sum(r['strict_global_exclusion'] for r in subset),
                               block_exclusions=sum(r['strict_local_exclusion'] for r in subset),
                               fixed_4x4_exclusions=sum(r['fixed_4x4_exclusion'] for r in subset),
                               fixed_4x4_coherent_exclusions=sum(r['fixed_4x4_coherent_exclusion'] for r in subset),
                               spatial_intersection=sum(r['strict_local_exclusion'] and r['fixed_4x4_exclusion']
                                                        for r in subset),
                               union_exclusions=sum(r['strict_global_exclusion'] or r['strict_local_exclusion']
                                                    or r['fixed_4x4_exclusion'] for r in subset)))
    save_json(HERE/'noise_curve_summary.json',dict(counts=counts,rows=rows,
              timings={name:{k:r[k] for k in ('wall_seconds','cpu_seconds')} for name,r in
                       [('public_data_proposal_search',raw),('strict_050_control',a),
                        ('strict_075_global',b),('strict_075_spatial',c),
                        ('strict_075_global_recovery',recovery_global),
                        ('strict_075_spatial_recovery',recovery_local),('fixed_4x4_cross_replay',fixed)]},
              fixed_4x4_higher_cap_evidence=fixed_high['cases']))
    draw(rows)


if __name__=='__main__':
    main()

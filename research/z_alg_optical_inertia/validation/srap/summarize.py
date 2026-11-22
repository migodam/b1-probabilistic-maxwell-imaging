"""Build a joined, portable SRAP result table and two-panel native study figure."""

import csv
from time import perf_counter, process_time

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from online import HERE, read, save_json, times


def chart(offline, base, high):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(9.6,8.5),sharey=True,gridspec_kw={'width_ratios':[1.4,1.]})
    values=list(offline.values())
    pos=np.arange(len(values))
    y=[case['object_id'].replace('s2d_','') for case in values]
    a,b=axes
    def ratio(src, field):
        return [src[row['id']][field]/row['padded_radius'] for row in values]
    a.barh(pos-.14,ratio(base,'joint_distance_lower'),.26,color='#247095',label='Joint')
    a.barh(pos+.14,ratio(base,'coherent_distance_lower'),.26,color='#9AABB8',label='Best coherent at same weights')
    b.barh(pos,ratio(high,'joint_distance_upper'),.46,color='#BC6235',
           label='Upper bound over five frozen weight choices')
    for ax in axes:
        ax.axvline(1,color='#222D3B',ls='--',lw=1.05)
        ax.grid(axis='x',color='#E0E5E9',lw=.55)
        ax.set_axisbelow(True)
        ax.spines['left'].set_color('#AAB4BD')
    a.set_yticks(pos,labels=y)
    a.invert_yaxis()
    a.set_xlim(0,3.13)
    b.set_xlim(0,1.08)
    a.set_ylabel('SRAP material ID')
    a.set_title('Cap 0.50 · Certified lower bounds',loc='left',fontweight='bold')
    b.set_title('Cap 0.75 · Certified upper bounds',loc='left',fontweight='bold')
    ha,la=a.get_legend_handles_labels()
    hb,lb=b.get_legend_handles_labels()
    fig.legend(ha+hb,la+lb,loc='lower center',bbox_to_anchor=(.53,.015),
               ncol=3,frameon=False,fontsize=8)
    fig.suptitle('Fixed multipliers on 22 additional SRAP materials',
                 x=.10,y=.973,ha='left',fontsize=15,fontweight='bold',color='#172637')
    fig.supxlabel('Certificate distance / padded noise radius (original 10% observations)',
                  y=.068,fontsize=10)
    fig.subplots_adjust(left=.10,right=.97,bottom=.125,top=.90,wspace=.23)
    for ext in ('png','pdf'):
        fig.savefig(HERE/f'srap22_transfer.{ext}',dpi=210,facecolor='white')
    plt.close(fig)


def main():
    start,cpu=perf_counter(),process_time()
    online=read(HERE/'online.json')
    offline=read(HERE/'offline.json')
    if not online['sealed'] or len(offline['rows'])!=22:
        raise ValueError('Complete original on/offline SRAP evidence before joining.')
    truth={r['id']:r for r in offline['rows']}
    cases=[]
    summary=[]
    keys=('id','object_id','stored_delta','padded_radius','exact_nu','native_retained_modes_lower',
          'joint_distance_lower','joint_distance_upper','joint_margin_lower','joint_margin_upper',
          'coherent_distance_lower','coherent_distance_upper','joint_status','coherent_status')
    base={}
    upper075={}
    for profile in online['profiles']:
        label=profile['label']
        rows=profile['rows']
        if len(rows)!=22:
            raise ValueError(f'Incomplete profile {label}.')
        counters=dict(label=label,cap=profile['cap'],positive_gap_lower=profile['positive_gap_lower'],
                      inherited_nu=profile['nu'],rows=len(rows),joint_certificates=0,joint_verified_below=0,
                      coherent_same_profile_certificates=0,coherent_same_profile_verified_below=0,
                      coherent_inertia_exhausted=0)
        for r in rows:
            t=truth[r['id']]
            cases.append(dict(profile=label,cap=profile['cap'],**{k:r[k] for k in keys},
                              true_fine_peak=t['original_fine_peak'],family=t['family'],
                              true_residual_lower=t['true_residual_lower'],
                              true_residual_upper=t['true_residual_upper'],
                              true_within_padded_radius=t['inside_padded_ball'],
                              nonvacuous_joint_exclusion=r['joint_status']=='excluded' and t['inside_padded_ball']))
            counters['joint_certificates']+=r['joint_status']=='excluded' and t['inside_padded_ball']
            counters['joint_verified_below']+=r['joint_status']=='verified_below_radius'
            counters['coherent_same_profile_certificates']+=r['coherent_status']=='excluded' and t['inside_padded_ball']
            counters['coherent_same_profile_verified_below']+=r['coherent_status']=='verified_below_radius'
            counters['coherent_inertia_exhausted']+=r['coherent_status']=='inertia_exhausts_single_source'
            if profile['cap']==.50:
                base[r['id']]=r
            elif r['joint_distance_upper']>upper075.get(r['id'],{'joint_distance_upper':-1})['joint_distance_upper']:
                upper075[r['id']]=r
        counters['best_normalized_joint_upper']=max(row['joint_distance_upper']/row['padded_radius'] for row in rows)
        counters['worst_normalized_joint_lower']=min(row['joint_distance_lower']/row['padded_radius'] for row in rows)
        summary.append(counters)
    with (HERE/'cases.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=cases[0].keys())
        writer.writeheader()
        writer.writerows(cases)
    highrow=max(upper075.values(), key=lambda x:x['joint_distance_upper']/x['padded_radius'])
    chart(truth,base,upper075)
    save_json(HERE/'summary.json',dict(profiles=summary,
              cohort_rows=len(truth),true_admissible=sum(r['inside_padded_ball'] for r in truth.values()),
              minimum_truth_margin=min(r['admissibility_margin_lower'] for r in truth.values()),
              raw_ball_truth_inside=sum(r['true_residual_upper']<=r['stored_delta'] for r in truth.values()),
              raw_ball_truth_outside=sum(r['true_residual_lower']>r['stored_delta'] for r in truth.values()),
              max_effective_relative_guard=max((r['true_residual_upper']-r['stored_delta'])
                                               /r['observation_norm_upper'] for r in truth.values()),
              true_peaks=[min(r['original_fine_peak'] for r in truth.values()),
                          max(r['original_fine_peak'] for r in truth.values())],
              lower_cap_joint_only_same_weights=[
                   dict(id=r['id'],joint_ratio_lower=r['joint_distance_lower']/r['padded_radius'],
                        coherent_ratio_upper=r['coherent_distance_upper']/r['padded_radius'])
                   for r in base.values() if r['joint_status']=='excluded' and
                   r['coherent_status']=='verified_below_radius'],
              upper_cap_most_favorable_observation=highrow['id'],
              upper_cap_most_favorable_normalized_upper=highrow['joint_distance_upper']/highrow['padded_radius'],
              online_cost={x:online[x] for x in ('wall_seconds','cpu_seconds')},
              offline_cost={x:offline[x] for x in ('wall_seconds','cpu_seconds')},
              tables_and_figure_cost=times(start,cpu)))


if __name__=='__main__':
    main()

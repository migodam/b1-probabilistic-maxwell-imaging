"""Build compact receipts and the measured-time figure from the two completed runs."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent


def load(side):
    return json.loads((HERE/f'result_{side}.json').read_text())


def main():
    rows = [load(32),load(96)]
    summary = []
    fig,axs = plt.subplots(1,2,figsize=(7.5,3.5),sharex=True,sharey=True,layout='constrained')
    points = (100,1000)
    styles = [('numerical_original','Original, floating point','#93908e','o'),
              ('numerical_moving','Moving traces, floating point','#2555a4','s'),
              ('verified_moving','Moving traces, interval proof','#a3292d','^')]
    for r,ax in zip(rows,axs,strict=True):
        side = r['side']
        prefixes = r['prefixes']
        for key,name,color,m in styles:
            vals = [prefixes[str(i)][key]['wall_seconds'] for i in points]
            ax.plot(points,vals,'-',marker=m,label=name,c=color,linewidth=1.8)
        negatives = [v['verified_state']['negative_at_least'] for v in r['boxes']]
        state_undecided = sum(v['verified_state']['unresolved'] != 0 for v in r['boxes'])
        assert all(v['update_dimension']==6 and v['update_rank']==4 for v in r['boxes'])
        assert not r['base_interval_border']['strict_exclusion_at_this_threshold']
        assert 0 < prefixes['1000']['strict_exclusions'] < 1000
        assert prefixes['1000']['numeric_full_agreement'] == 1000
        assert all(x['independent_full_counts_agree'] for x in r['interval_original_select_nodes'])
        ax.set(xscale='log',yscale='log',xticks=points,xticklabels=['100','1,000'],
               title=f'{side} × {side} grid · {r["current_dimension"]:,} current variables',
               xlabel='Distinct moving-support boxes')
        ax.tick_params(axis='both',labelsize=9)
        ax.grid(axis='y',alpha=.3,which='major')
        ax.minorticks_off()
        summary.append({k:r[k] for k in ('side','observation','current_dimension','source_channels','known_exterior_vacuum',
                                        'base_interval_pivot','base_interval_border','tail_j','threshold','prefixes',
                                        'interval_original_select_nodes_time','neumann_probe')})
        summary[-1].update(replayed_nodes=len(r['interval_original_select_nodes']),
                           verified_original_and_moving_matches=sum(x['independent_full_counts_agree'] for x in r['interval_original_select_nodes']),
                           strict_state_negative_range=[min(negatives),max(negatives)],
                           state_sign_count_undecided_boxes=state_undecided,
                           bordered_dimension_original=r['current_dimension']+r['source_channels'],
                           bordered_dimension_reduced=4+r['source_channels'],
                           timing_stages={k:r[k] for k in ('shared_operator_and_observation_load','design_batch_once',
                                                          'shared_operator_algebra_setup','moving_numeric_setup',
                                                          'moving_interval_input_and_green','moving_interval_setup',
                                                          'base_verification_time')})
    axs[0].set_ylabel('Online wall time (s)')
    handles,labels = axs[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside upper center',ncols=3,frameon=False,fontsize=8)
    fig.savefig(HERE/'moving_trace_timings.pdf')
    fig.savefig(HERE/'moving_trace_timings.png',dpi=165)
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()

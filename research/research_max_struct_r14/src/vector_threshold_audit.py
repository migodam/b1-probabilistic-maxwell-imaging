"""Existing sealed vector predictions: pair distances, chi=1 flips, tick repair."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys
import time


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();out=root/'audit/vector4_v1/threshold_precision_v1';out.mkdir(parents=True,exist_ok=False)
    sys.path[:0]=[str(root/'src'),str(root.parents[1])]
    cpu,wall=time.process_time(),time.perf_counter()
    import numpy as np
    from run_vector_factorial import dump
    from research.research_max_struct_r13.src.score_results import table
    analysis=root/'analysis/vector4_v1';rows=json.loads((analysis/'SCORES.json').read_text())
    manifest=json.loads((root/'configs/vector4_v1.json').read_text());cases=manifest['cases'];methods=manifest['methods']
    with np.load(analysis/'SCORED_MATERIAL_ROWS_AND_MASKS.npz',allow_pickle=False) as z:
        materials=z['m'].copy();truth=z['truth'].copy();support=z['truth_support'].copy();labels=z['row_labels'].copy()
    label_indices={str(label):i for i,label in enumerate(labels)}
    lookup={(r['case_id'],r['method'],r['view']):r for r in rows}
    pairs=[];flips=[];flip_masks=[];all_material=[]
    for case in cases:
        values={method:materials[label_indices[case+'__'+method+'__primary']] for method in methods}
        all_material.append(np.asarray([values[method] for method in methods]))
        index=label_indices[case+'__'+methods[0]+'__primary'];gt=truth[index];gt_support=support[index]
        amp=float(lookup[(case,methods[0],'primary')]['amplitude'])
        for first,second in itertools.combinations(methods,2):
            x,y=values[first],values[second];difference=x-y;mask=(x>1.)!=(y>1.);flip_masks.append(mask)
            ra,rb=lookup[(case,first,'primary')],lookup[(case,second,'primary')]
            n=int(mask.sum());margin=np.maximum(abs(x[mask]-1.),abs(y[mask]-1.))
            row=dict(case_id=case,method_a=first,method_b=second,common_physical_objective='original per-source-normalized vectorPhi+same prior; saved paid states, no new solve',
                L2_distance=float(np.linalg.norm(difference)),RMSE_distance=float(np.linalg.norm(difference)/np.sqrt(x.size)),
                amplitude_normalized_RMSE_distance=float(np.linalg.norm(difference)/np.sqrt(x.size)/amp),
                boxrange_normalized_RMSE_distance=float(np.linalg.norm(difference)/np.sqrt(x.size)/8),
                symmetric_relative_L2_distance=float(np.linalg.norm(difference)/max(np.linalg.norm(x),np.linalg.norm(y),1e-300)),
                Linf_distance=float(np.max(abs(difference))),Phi_a=ra['Phi'],Phi_b=rb['Phi'],
                absolute_Phi_difference=abs(ra['Phi']-rb['Phi']),
                relative_Phi_difference=abs(ra['Phi']-rb['Phi'])/max(abs(ra['Phi']),abs(rb['Phi']),1e-300),
                NRMSE_a=ra['NRMSE'],NRMSE_b=rb['NRMSE'],IoU_a=ra['IoU'],IoU_b=rb['IoU'],
                support_flips=n,support_threshold=1.,threshold_is_unchanged=True,
                largest_flipping_cell_pair_difference=float(np.max(abs(difference[mask]))) if n else None,
                max_flipping_cell_distance_to_threshold_both_sides=float(margin.max()) if n else None,
                both_values_within_fixed_diagnostic_bands={str(e):int(np.count_nonzero(margin<=e)) for e in (1e-6,1e-4,1e-3,1e-2,.1)},
                truth_support_cells=int(gt_support.sum()),truth_cells_exactly_one=int(np.count_nonzero(gt==1.)))
            pairs.append(row)
            for cell in np.flatnonzero(mask):
                flips.append(dict(case_id=case,method_a=first,method_b=second,cell=int(cell),ijk=np.unravel_index(cell,(6,6,6)),
                    chi_a=float(x[cell]),chi_b=float(y[cell]),chi_a_minus1=float(x[cell]-1),chi_b_minus1=float(y[cell]-1),
                    absolute_difference=float(abs(x[cell]-y[cell])),truth_chi=float(gt[cell]),truth_support=bool(gt_support[cell]),
                    support_a=bool(x[cell]>1),support_b=bool(y[cell]>1)))
    table(out,'PAIR_DISTANCES_AND_COMMON_PHI',pairs);table(out,'EXACT_SUPPORT_FLIPS',flips)
    np.savez_compressed(out/'PAIR_SUPPORT_ARRAYS.npz',materials=np.asarray(all_material),case_ids=np.asarray(cases),methods=np.asarray(methods),
        flip_masks=np.asarray(flip_masks),pair_labels=np.asarray([r['case_id']+'__'+r['method_a']+'__'+r['method_b'] for r in pairs]),support_threshold=np.asarray(1.))
    target=next(r for r in pairs if r['case_id']=='mixed_crossed_bracket' and r['method_a']=='pair00' and r['method_b']=='pair11')
    note=['# 3D fixed-threshold support audit','',
        'Existing sealed/scored primary predictions only. No optimizer, forward, new scoring threshold, or recovery label was introduced. Pair distances use the original inherited amplitude and fixed chi>1 masks. Phi values are the already-paid common original vector objective, not new evaluations.','',
        f"Bracket 00 versus11: amplitude-normalized pair RMSE={target['amplitude_normalized_RMSE_distance']:.12g}, symmetric relativeL2={target['symmetric_relative_L2_distance']:.12g}, Linf={target['Linf_distance']:.12g}; Phi difference={target['absolute_Phi_difference']:.12g}.",
        f"Its fixed-threshold IoU changes {target['IoU_a']:.12g}→{target['IoU_b']:.12g} through {target['support_flips']} cell flips. Maximum distance of either flipping value from chi=1 is {target['max_flipping_cell_distance_to_threshold_both_sides']}; exact cell values and truth labels are in EXACT_SUPPORT_FLIPS.",
        '', f"For this bracket comparison, {target['both_values_within_fixed_diagnostic_bands']['0.01']} flipping cells have both values within.01 ofchi1, and {target['both_values_within_fixed_diagnostic_bands']['0.1']} have both within.1. The exact margins exceed tiny floating tolerances; a pure precision explanation is unsupported here. Similar scalar NRMSE values do not imply small pairwise material distance.",
        '', 'The fixed support threshold discretizes the differing fields and may amplify their apparent support-score difference. These tables do not establish a different recovery class. No threshold/optimizer was retuned. All60 method pairs across the four conditions and every flipping voxel remain visible.']
    (out/'NOTES.md').write_text('\n'.join(note)+'\n')
    before=root/'audit/vector4_v1/figure_tick_fix_before';before.mkdir(exist_ok=False)
    name='CORRECTOR_PHI_AND_MATERIAL_SCATTER';figure_dir=analysis/'FIGURES'
    for extension in ('png','pdf'):shutil.copyfile(figure_dir/(name+'.'+extension),before/(name+'.'+extension))
    fix_cpu,fix_wall=time.process_time(),time.perf_counter()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator,FuncFormatter,NullLocator
    names={'pair00':'00 P0/P0 global','pair01':'01 P0/P1 global','pair10':'10 P1/P0 global','pair11':'11 P1/P1 global',
           'CIE_sourcewise':'CIE sourcewise','Dual_sourcewise':'Dual sourcewise'}
    colors={method:plt.get_cmap('tab10')(i) for i,method in enumerate(methods)}
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),constrained_layout=True)
    for method in methods:
        for case in cases:
            x=lookup[(case,method,'raw')];y=lookup[(case,method,'corrected')];marker='^' if 'crossgrid' in case else 'o'
            axes[0].scatter(x['Phi'],y['Phi'],color=colors[method],marker=marker,s=60,label=names[method] if case==cases[0] else None)
            axes[1].scatter(x['NRMSE'],y['NRMSE'],color=colors[method],marker=marker,s=60)
    for i,key in enumerate(('Phi','NRMSE')):
        values=[r[key] for r in rows if r['view'] in ('raw','corrected') and r.get(key) is not None]
        lower,upper=min(values),max(values);axes[i].plot([lower,upper],[lower,upper],'k--',lw=1)
        axes[i].set_xlabel('raw '+key);axes[i].set_ylabel('corrected '+key);axes[i].grid(alpha=.2)
    axes[0].set_xscale('log');axes[0].set_yscale('log');axes[0].legend(fontsize=8)
    axes[0].xaxis.set_major_locator(FixedLocator([2e-4,1e-3,5e-3]));axes[0].xaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:.0e}'))
    axes[0].xaxis.set_minor_locator(NullLocator())
    fig.suptitle('Same original physical corrector: Phi and material changes are separate; triangles repeat the crossgrid frame')
    for extension in ('png','pdf'):fig.savefig(figure_dir/(name+'.'+extension),dpi=170,bbox_inches='tight')
    plt.close(fig)
    dump(out/'FIGURE_TICK_FIX_RECEIPT.json',dict(status='COMPLETE',change='leftxaxis tick labels only: fixed2e-4/1e-3/5e-3, minorlabelsremoved; same points/scales/limits',
        old_files=str(before),new_files=[str(figure_dir/(name+'.'+extension)) for extension in ('png','pdf')],
        CPU_seconds=time.process_time()-fix_cpu,wall_seconds=time.perf_counter()-fix_wall,forward_LU=0,solver_RHS=0,optimizer_calls=0))
    dump(out/'RECEIPT.json',dict(status='COMPLETE',pairs=len(pairs),flipped_cell_pair_rows=len(flips),
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,forward_LU=0,solver_RHS=0,
        predictions_modified=False,threshold_modified=False,new_GT_files_opened=False,GT_access='previously authorized scored arrays/masks only'))
    print(json.dumps(target))


if __name__=='__main__':main()

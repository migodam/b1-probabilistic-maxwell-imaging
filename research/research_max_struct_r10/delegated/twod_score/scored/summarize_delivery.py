"""Summarize frozen scores and saved states only; no forward or GT selection."""
from pathlib import Path
import json,csv,sys,time
from collections import defaultdict
import numpy as np
P=Path(__file__).resolve().parent;ROOT=P.parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
from research.research_max_struct_r10.experiments.score_twod_comparison import LABELS,PRIMARY,csvwrite,dump
def main():
    start=time.perf_counter();cpu=time.process_time();rows=json.loads((P/'SCORES_FINAL.json').read_text());ix={(r['case_id'],r['method']):r for r in rows}
    cases=sorted({r['case_id'] for r in rows});ratios=[]
    for m in LABELS:
        z=[r for r in rows if r['method']==m];base=[ix[r['case_id'],'cie_single'] for r in z]
        ratios.append(dict(method=m,cases=len(z),all_child_CPU_seconds=sum(r['external_full_child_CPU_seconds'] for r in z),
            child_CPU_ratio_of_totals_to_CIE=sum(r['external_full_child_CPU_seconds'] for r in z)/sum(r['external_full_child_CPU_seconds'] for r in base),
            median_object_paired_CPU_ratio_to_CIE=float(np.median([a['external_full_child_CPU_seconds']/b['external_full_child_CPU_seconds'] for a,b in zip(z,base)]))))
    csvwrite(P/'FULL_CHILD_CPU_RATIOS.csv',ratios)
    failed=[];discord=[];gatechanges=[]
    for r in rows:
        if not r['material_success']:
            failed.append(dict(case_id=r['case_id'],regime=r['regime'],method=r['method'],primary_comparison=r['primary_comparison'],
                amplitude_NRMSE=r['amplitude_NRMSE'],fixed_support_IoU=r['fixed_support_IoU'],fails_NRMSE=r['amplitude_NRMSE']>.15,
                fails_fixed_IoU=r['fixed_support_IoU']<.7,fails_legality=not r['material_legal'],original_Phi=r['original_Phi'],
                corrector_stop=r['corrector_stop'],corrector_iterations=r['corrector_iterations'],external_full_child_CPU_seconds=r['external_full_child_CPU_seconds']))
        if r['material_success']!=r['legacy_success']:discord.append(r)
        if r['method'] in ('matrix_twofold','matrix_green19'):
            b=ix[r['case_id'],'cie_single']
            if r['material_success']!=b['material_success']:
                gatechanges.append(dict(case_id=r['case_id'],regime=r['regime'],method=r['method'],CIE_NRMSE=b['amplitude_NRMSE'],new_NRMSE=r['amplitude_NRMSE'],
                    NRMSE_change=r['amplitude_NRMSE']-b['amplitude_NRMSE'],CIE_fixed_IoU=b['fixed_support_IoU'],new_fixed_IoU=r['fixed_support_IoU'],
                    CIE_success=b['material_success'],new_success=r['material_success']))
    csvwrite(P/'ALL_MATERIAL_FAILURES.csv',failed);csvwrite(P/'LEGACY_PRIMARY_GATE_DISAGREEMENT.csv',discord);csvwrite(P/'NEW_VS_CIE_GATE_CHANGES.csv',gatechanges)
    with (P/'BUDGET_SCORES.csv').open() as f:b=list(csv.DictReader(f))
    group=defaultdict(list)
    for r in b:group[r['budget_axis'],r['budget'],r['method']].append(r)
    budgets=[]
    for (axis,target,m),rr in group.items():
        v=[r for r in rr if r['status']=='AVAILABLE'];budgets.append(dict(axis=axis,budget=float(target),method=m,expected_objects=len(rr),available=len(v),
            mean_amplitude_NRMSE=float(np.mean([float(r['amplitude_NRMSE']) for r in v])) if v else None,
            median_amplitude_NRMSE=float(np.median([float(r['amplitude_NRMSE']) for r in v])) if v else None,
            material_successes=sum(r['material_success']=='True' for r in v),completed_endpoint_count=sum(r['selection_kind']=='COMPLETED_PUBLIC_PHI_ENDPOINT' for r in v),
            conservative_external_snapshot_count=sum(r.get('selected_cost_semantics')=='CONSERVATIVE_EXTERNAL_CPU_UPPER_BOUND_NOT_OBSERVED_PREFIX' for r in v),
            scope='Latest saved paid Phi-incumbent; internal clock or conservative external envelope as axis states; not interpolated solver state'))
    csvwrite(P/'BUDGET_SUMMARY.csv',budgets)
    # Three cases: last native stage, original-Phi chosen corrector start, final.
    with (P/'SCORES_STAGES.csv').open() as f:stages=list(csv.DictReader(f))
    stageix={(r['case_id'],r['method'],r['state']):r for r in stages};details=[]
    for case in ('r10_twod_0001','r10_twod_0004','r10_twod_0022'):
        with np.load(ROOT/'delegated/twod_energy/gt_current_cache'/(case+'.npz'),allow_pickle=False) as z:jt=z['j_gt'].copy()
        for m in LABELS:
            r=ix[case,m];d=ROOT/'RESULTS/twod'/r['task_id'];rec=json.loads((d/'RECEIPT.json').read_text());selected=rec.get('auxiliary_selected_stage')
            for state in ('stage_03_physical','selected_before_corrector','prediction'):
                with np.load(d/(state+'.npz'),allow_pickle=False) as z:arr={k:z[k].copy() for k in z.files}
                sm=json.loads((d/(state+'.json')).read_text());met=r if state=='prediction' else stageix[case,m,state]
                dm=arr['j']-jt;pg=float(np.max(abs(arr['m']-np.clip(arr['m']-arr['gradient'],.02,20)))) if 'gradient' in arr else None
                item=dict(case_id=case,method=m,state=state,auxiliary_selected_stage=selected,amplitude_NRMSE=float(met['amplitude_NRMSE']),
                    fixed_support_IoU=float(met['fixed_support_IoU']),original_Phi=sm.get('objective'),physical_current_error_squared_energy=float(np.vdot(dm,dm).real),
                    saved_projected_gradient_inf=pg,corrector_stop=r['corrector_stop'],corrector_iterations=r['corrector_iterations'],external_full_child_CPU_seconds=r['external_full_child_CPU_seconds'])
                auxstage=3 if state=='stage_03_physical' else selected if state=='selected_before_corrector' else None
                if auxstage is not None and auxstage>=0:
                    ap=d/f'stage_{auxstage:02d}_joint.npz'
                    with np.load(ap,allow_pickle=False) as z:aux={k:z[k].copy() for k in z.files}
                    if np.array_equal(aux['m'],arr['m']):
                        ae=aux['J']-jt;gap=aux['J']-arr['j'];item.update(auxiliary_current_error_squared_energy=float(np.vdot(ae,ae).real),
                            auxiliary_minus_physical_current_squared_energy=float(np.vdot(gap,gap).real),cross_energy=2*float(np.vdot(dm,gap).real),same_m_pair=True)
                    else:item.update(same_m_pair=False)
                details.append(item)
    csvwrite(P/'THREE_CASE_NATIVE_CORRECTOR_DETAILS.csv',details)
    # Fixed budget points only; lines are visual guides, never interpolated states.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,3,figsize=(16,4.8),constrained_layout=True);colors=plt.get_cmap('tab10')
    axes=('internal_paid_CPU','external_full_child_CPU_conservative','physical_LU')
    for ax,axis in zip(axs,axes):
        for i,m in enumerate(PRIMARY):
            z=sorted([r for r in budgets if r['axis']==axis and r['method']==m],key=lambda r:r['budget'])
            ax.plot([r['budget'] for r in z],[r['mean_amplitude_NRMSE'] for r in z],':o',color=colors(i),label=LABELS[m],markersize=4)
        ax.grid(alpha=.2);ax.set_xscale('log');ax.set_ylabel('mean amplitude NRMSE over all24')
    axs[0].set_xlabel('Recorded internal paid CPU budget (nested clock)');axs[1].set_xlabel('Conservative full-child CPU budget bound');axs[2].set_xlabel('Original576-LU budget (small systems excluded)')
    axs[0].set_title('Not full child CPU');axs[1].set_title('Bounded snapshots + exact completed endpoints');axs[2].set_title('One physical counter, not equal total compute')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle('Fixed nearest-saved Phi-incumbent budget slices | lines are guides; no solver interpolation',fontsize=12)
    fp=ROOT/'FIGURES/twod_comparison/BUDGET_CURVES_SCOPE_EXPLICIT.png';fig.savefig(fp,dpi=170);fig.savefig(fp.with_suffix('.pdf'));plt.close(fig)
    dump(P/'DELIVERY_SUMMARY_RECEIPT.json',dict(status='COMPLETE_SAVED_RESULTS_ONLY',new_forward_calls=0,new_inverse_runs=0,new_sealed_truth_reads=0,
        inputs='Frozen score JSON/CSV, saved public predictions/stages, previously computed GT-current cache',
        failure_rows=len(failed),gate_disagreement_rows=len(discord),budget_cells=len(b),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start))
if __name__=='__main__':main()

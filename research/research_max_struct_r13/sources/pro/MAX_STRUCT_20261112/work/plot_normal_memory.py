"""Independent offline scoring/figures for the completed normal-memory study.

Reads completed JSON and NPZ only. Does not import or run research algorithms.
Keeps original offline_scores.json unchanged and checks agreement before adding
the integer-event branch. Truth is only used in this post-hoc scoring program.
"""
from pathlib import Path
import os
import json, math, hashlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(os.environ.get('MAX_STRUCT_RUN_ROOT', str(Path(__file__).resolve().parent.parent)))
WORK=ROOT/'work/normal_memory'
BASE=Path(os.environ.get('MAX_STRUCT_REPO', str(ROOT/'research_repo')))/'research/research_max_struct_r12/supplement_paired'
FIG=WORK/'figures';FIG.mkdir(exist_ok=True)
ORDER=['A__r10_twod_0004','B__r12_confirm24_0010','C__r8_test_0001']
SHORT={'A__r10_twod_0004':'A: slotted shell','B__r12_confirm24_0010':'B: persistent failure','C__r8_test_0001':'C: exposed NN case'}
COLORS={'ordinary':'#267c9c','event':'#e2983d','integer_event':'#754a96'}
LABELS={'ordinary':'Ordinary DR','event':'Fractional event','integer_event':'Integer event'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'figure.dpi':160,'savefig.dpi':200,'font.family':'DejaVu Sans'})

def read(p):return json.loads(p.read_text())
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(fig,name):
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight')
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

old=read(WORK/'offline_scores.json')
old_hash=digest(WORK/'offline_scores.json')
amp_map={r['case']:r['amplitude']for r in old}
old_map={(r['case'],r['arm']):r for r in old}
datasets={};scores=[];events=[];costs=[];trace_rows=[]
for case in ORDER:
    cohort,case_id=case.split('__');folder=WORK/case
    p=read(folder/'result.json');f=read(folder/'followup_result.json');g=read(folder/'integer_result.json')
    a=dict(np.load(folder/'followup_materials.npz'));b=dict(np.load(folder/'integer_materials.npz'))
    with np.load(BASE/'offline_truth'/cohort/(case_id+'.npz')) as z:truth=z['m'].copy();truth_amp=float(z['amplitude'])
    amp=amp_map[case]
    assert amp==truth_amp,(case,amp,truth_amp)
    np.testing.assert_array_equal(a['initial_material'],b['initial_material'])
    np.testing.assert_array_equal(a['AP_material'],b['AP_material'])
    arrays={'initial':a['initial_material'],'AP':a['AP_material'],'ordinary_best':a['ordinary_best_m'],
            'event_best':a['event_best_m'],'integer_event_best':b['integer_event_best_m']}
    armmeta={**f['arms'],**g['arms']}
    threshold=float(p['first_tcrit'])
    for name,m in arrays.items():
        e=m-truth
        row=dict(case=case,arm=name,amplitude_NRMSE=float(np.linalg.norm(e)/(np.sqrt(len(e))*amp)),
                 relative_L2=float(np.linalg.norm(e)/np.linalg.norm(truth)),amplitude=amp,
                 selection='Original physical Phi; GT used only for offline scoring')
        if name=='initial':
            row.update(f['initial'])
            row['selection']='Frozen original input; GT used only for offline scoring'
        elif name=='AP':
            row['residual_over_delta']=p['initial_AP_shadow_residual_over_delta']
            row['selection']='Material shadow of the shared AP pair after 3000 steps; GT used only for offline scoring'
        else:
            branch=name.removesuffix('_best')
            row.update(armmeta[branch]['best'])
            row.update(branch_CPU_seconds=armmeta[branch]['cpu_seconds'],branch_wall_seconds=armmeta[branch]['wall_seconds'],
                       first_large_box_switch_index=armmeta[branch]['first_box_jump'])
        if (case,name)in old_map:
            previous=old_map[(case,name)]
            for key in ['amplitude_NRMSE','relative_L2']:
                assert abs(row[key]-previous[key])<1e-12,(case,name,key,row[key],previous[key])
        scores.append(row)
    events.append(dict(case=case,threshold_t=threshold,ordinary_initial_t=1.,predicted_first_switch_index=int(math.floor(threshold)),
      observed_first_box_jump=f['arms']['ordinary']['first_box_jump'],
      observed_predicted_cell_jump=f['arms']['ordinary']['predicted_cell_box_jump'],
      fractional_start_t=max(1.,1.001*threshold),integer_start_t=max(1.,float(math.floor(threshold)+1)),
      AP_pair_relative_defect=p['relative_pair_error'],AP_pair_absolute_defect=p['AP_pair_error']))
    branch_cpu=sum(v['cpu_seconds']for v in armmeta.values())
    branch_wall=sum(v['wall_seconds']for v in armmeta.values())
    ns=a['initial_material'].size # overwritten below from saved state dimensions
    with np.load(folder/'states.npz') as z:ns=z['cX'].shape[1]
    # Each of the two followup processes scores the original material once before
    # starting timed branches. Those six calls are explicit in source, absent from
    # the per-branch counters, and must be added to actual physical work.
    physical_LU=int(p['physical_cost']['LU']+sum(v['physical_LU']for v in armmeta.values())+2)
    physical_RHS=int(p['physical_cost']['RHS']+sum(v['physical_RHS']for v in armmeta.values())+2*ns)
    projection=dict(p['projection_cost'])
    for v in armmeta.values():
        for k,val in v['additional_projection_cost'].items():projection[k]=projection.get(k,0)+val
    # Two followup projectors are constructed outside branch counters.
    projection['H_cholesky']+=2
    projection['port_eigh']+=2
    # Common S has one grouped receiver port per projector.
    with np.load(BASE/'inputs'/(case+'.npz')) as z:nr=z['S'].shape[0]
    projection['H_setup_rhs']+=2*nr
    costs.append(dict(case=case,probe_CPU_seconds=p['cpu_seconds'],probe_wall_seconds=p['wall_seconds'],
      probe_scope='AP 3000 + both 500-step probe continuations + probe setup/audits; not AP-only',
      continuation_branches_CPU_seconds=branch_cpu,continuation_branches_wall_sum=branch_wall,
      measured_CPU_subtotal=p['cpu_seconds']+branch_cpu,
      missing_from_CPU_subtotal='Six followup process imports, projector setup, initial physical audit, output and exit overhead; whole-process total unavailable',
      actual_physical_LU_from_receipts_and_source=physical_LU,actual_physical_RHS_from_receipts_and_source=physical_RHS,
      projection_work_from_receipts_and_source=projection))
    for x in f['trajectory']+g['trajectory']:trace_rows.append(dict(case=case,**x))
    datasets[case]=dict(probe=p,follow=f,integer=g,arrays=arrays,truth=truth,armmeta=armmeta)

assert digest(WORK/'offline_scores.json')==old_hash
summary=dict(scope='Three exposed mechanism-development cases; no corrector; independent offline scoring after frozen runs',
  existing_offline_scores_sha256=old_hash,existing_scores_unchanged=True,
  scores=scores,events=events,costs=costs,
  totals=dict(measured_CPU_subtotal_seconds=sum(c['measured_CPU_subtotal']for c in costs),
              actual_physical_LU_from_receipts_and_source=sum(c['actual_physical_LU_from_receipts_and_source']for c in costs),
              actual_physical_RHS_from_receipts_and_source=sum(c['actual_physical_RHS_from_receipts_and_source']for c in costs),
              H_cholesky=sum(c['projection_work_from_receipts_and_source']['H_cholesky']for c in costs)),
  cost_limit='Measured CPU subtotal is not whole-process or end-to-end CPU. Original probe cost includes two exploratory 500-step continuations.',
  comparison_limit='All branches restart from the same saved AP pair with different normal memory. Original published DR starts earlier from original material; its cost/quality is not matched here.')
(WORK/'combined_offline_scores.json').write_text(json.dumps(summary,indent=2)+'\n')
pd.DataFrame(scores).to_csv(WORK/'comparison.csv',index=False)
pd.DataFrame(events).to_csv(WORK/'finite_thresholds.csv',index=False)
pd.DataFrame(trace_rows).to_csv(WORK/'physical_traces.csv',index=False)

# 1. Finite ray thresholds versus actual large box switches in ordinary DR.
fig,axes=plt.subplots(2,3,figsize=(13.5,7.1),sharey=True)
for j,case in enumerate(ORDER):
    ds=datasets[case];p=ds['probe'];e=events[j]
    ray=p['ray'];t=np.array([r['t']for r in ray]);mm=np.array([r['first_cell_material']for r in ray])
    ax=axes[0,j];ax.plot(t,mm,'o-',color='#355f7b',lw=1.4,ms=4)
    ax.axvline(e['threshold_t'],color='#a14e34',ls='--',lw=1.3)
    ax.set_title(SHORT[case]+f"\nPredicted t* = {e['threshold_t']:.4f}",fontsize=11)
    ax.set_xlabel('Normal-ray distance t')
    ax.set_ylim(-.9,21.5);ax.set_yticks([.02,10,20],['0.02','10','20'])
    ax.grid(axis='y',alpha=.15)
    first=ds['follow']['arms']['ordinary']['first50']
    ax=axes[1,j];xx=np.array([r['iteration']for r in first]);yy=np.array([r['predicted_cell_m']for r in first])
    cap=max(20,e['predicted_first_switch_index']+5);mask=xx<=cap
    ax.step(xx[mask],yy[mask],where='post',color=COLORS['ordinary'],lw=1.5)
    ax.scatter(xx[mask],yy[mask],color=COLORS['ordinary'],s=14,zorder=3)
    ax.axvline(e['predicted_first_switch_index'],color='#a14e34',ls='--',lw=1.3)
    ax.set_xlabel('Ordinary DR iteration index k')
    ax.set_title(f"Predicted k = {e['predicted_first_switch_index']}; observed k = {e['observed_first_box_jump']}",fontsize=11)
    ax.set_xlim(-.5,cap+.5);ax.grid(axis='y',alpha=.15)
axes[0,0].set_ylabel('Material at predicted cell')
axes[1,0].set_ylabel('Material at the same cell')
fig.suptitle('Finite normal-memory threshold and the first large material switch',fontsize=15,y=1.01)
fig.text(.5,-.01,'Top: fixed AP normal ray. Bottom: actual ordinary DR dynamics. The AP pair is approximate (relative defect 0.09%–0.60%).\n'
                  'Index k = 0 starts at t = 1. Predicting a box switch does not predict successful reconstruction.',ha='center',fontsize=10)
fig.tight_layout(rect=[0,.05,1,.97]);save(fig,'finite_threshold_mechanism')

# 2. Physical data residual of audited candidates, not shadow gaps or incumbent minima.
fig,axes=plt.subplots(1,3,figsize=(14,4.7),sharex=True)
for j,case in enumerate(ORDER):
    ds=datasets[case];ax=axes[j]
    alltr=ds['follow']['trajectory']+ds['integer']['trajectory']
    for arm in ['ordinary','event','integer_event']:
        tr=[r for r in alltr if r['arm']==arm]
        ax.plot([r['iteration']for r in tr],[r['residual_over_delta']for r in tr],color=COLORS[arm],lw=1.55,
                label=LABELS[arm],alpha=.94)
    ax.axhline(1.,ls='--',color='#3e6e45',lw=1.25)
    ax.axhspan(.65,1,color='#c7ddcd',alpha=.25)
    ax.set_yscale('log');ax.set_ylim(.65,max(r['residual_over_delta']for r in alltr)*1.35)
    ax.set_xlim(0,3000);ax.set_xticks([0,1000,2000,3000])
    ax.set_title(SHORT[case],fontsize=12)
    ax.set_xlabel('Continuation iterations after the shared AP pair')
    ax.grid(axis='y',which='major',alpha=.16)
axes[0].set_ylabel(r'Full Maxwell residual / $\delta$')
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.035),frameon=False)
fig.text(.5,-.04,'Every curve uses full physical audits every 50 steps; the dashed line is the noise radius. All branches have 3000 steps, but actual CPU differs.\n'
                  'The shared 3000-step AP preparation and exploratory probes are additional work. No physical corrector is applied.',ha='center',fontsize=10)
fig.tight_layout(rect=[0,.04,1,.94]);save(fig,'physical_residual_trajectories')

# 3. Truth and physically selected material; identical absolute color scale.
columns=[('truth','Truth'),('initial','Original material'),('AP','AP shadow'),('ordinary_best','Ordinary DR'),
         ('event_best','Fractional event'),('integer_event_best','Integer event')]
fig,axes=plt.subplots(3,6,figsize=(14.8,8.3))
lookup={(r['case'],r['arm']):r for r in scores}
for i,case in enumerate(ORDER):
    ds=datasets[case]
    for j,(arm,title) in enumerate(columns):
        ax=axes[i,j];m=ds['truth']if arm=='truth'else ds['arrays'][arm]
        im=ax.imshow(m.reshape(24,24),origin='lower',vmin=0,vmax=20,cmap='viridis',interpolation='nearest')
        ax.set_xticks([]);ax.set_yticks([])
        if i==0:ax.set_title(title,fontsize=11,pad=10)
        if arm=='truth':sub='Offline reference'
        else:
            r=lookup[(case,arm)];sub=f"NRMSE {r['amplitude_NRMSE']:.4f}\nresidual / δ = {r['residual_over_delta']:.3f}"
            if r.get('iteration')==-1:sub+='\ninitial incumbent retained'
        ax.set_xlabel(sub,fontsize=9,labelpad=4)
        if j==0:ax.set_ylabel(SHORT[case],fontsize=11,labelpad=12)
fig.subplots_adjust(left=.065,right=.925,top=.88,bottom=.12,hspace=.50,wspace=.12)
cax=fig.add_axes([.94,.20,.014,.61]);cb=fig.colorbar(im,cax=cax);cb.set_label('Real material contrast',fontsize=11)
fig.suptitle('Exposed Maxwell mechanism experiments: physical materials selected by the original objective',fontsize=14,y=.98)
fig.text(.5,.925,'Same saved AP pair for all three continuation branches; only the initial normal memory differs. No corrector and no GT selection.',ha='center',fontsize=11)
fig.text(.5,.025,'Every continuation starts after 3000 paid AP steps. Both C event results are weaker than published DR from the original NN material (NRMSE 0.07799);\n'
                   'this figure supports a mechanism comparison, not superiority to the original DR pipeline.',ha='center',fontsize=10)
save(fig,'physical_material_comparison')

# Human-readable experiment receipt including the newly scored integer branch.
df=pd.DataFrame(scores);wide=df.pivot(index='case',columns='arm',values='amplitude_NRMSE').loc[ORDER]
lines=['# Normal-memory Maxwell mechanism experiment: offline evaluation','','## Scope','',
       'Three already exposed physical instances; exact published field/current projectors; shared states after 3000 AP steps; ordinary, fractional-event and integer-event DR continuations of 3000 steps. No corrector. No truth is used in initialization, threshold calculation, trajectories or minimum-Phi selection. This file independently scores the frozen material arrays after the runs.','',
       'The original `offline_scores.json` is preserved. Its 12 scores agree with recomputation to 1e-12. Integer-event scores are newly appended to the combined artifact.','',
       '## Material results','','|Case|Original material|AP shadow|Ordinary DR|Fractional event|Integer event|','|---|---:|---:|---:|---:|---:|']
for case in ORDER:
    r=wide.loc[case];lines.append('|'+SHORT[case]+'|'+'|'.join(f'{r[k]:.6f}'for k in ['initial','AP','ordinary_best','event_best','integer_event_best'])+'|')
lines+=['','The B continuations do not improve the original incumbent. Its best iteration is −1; therefore the selected B material is the original input, not the AP endpoint or a successful reflected reconstruction.','',
        '## Finite-event prediction','','|Case|Predicted normal t*|Predicted ordinary index|Observed large box switch|Relative AP-pair defect|','|---|---:|---:|---:|---:|']
for e in events:lines.append(f"|{SHORT[e['case']]}|{e['threshold_t']:.6f}|{e['predicted_first_switch_index']}|{e['observed_first_box_jump']}|{e['AP_pair_relative_defect']:.6g}|")
lines+=['','The ordinary branch starts at normal coordinate t=1, so the first strict switch is predicted at index floor(t*). Small interior motion before that event is not a full box switch. The three AP pairs are approximate: the prediction is checked against actual DR, not asserted as an exact trajectory theorem for these finite-precision states. The B negative example shows that a correctly predicted switch need not lead to recovery.','',
        '## Cost and limits','',f"Measured probe CPU plus all nine timed continuation branches: **{summary['totals']['measured_CPU_subtotal_seconds']:.6f} s**. This is a measured subtotal, **not** an end-to-end CPU total. The three probe times already include 3000 AP steps, both 500-step exploratory branches, setup and their physical audits. Six later process imports, projector setup, initial physical scores, serialization and process exit are outside continuation timers and cannot be reconstructed as seconds from the JSON. Do not label the saved field `shared_AP_cpu_seconds...` as AP-only cost.",'',
        f"Exact physical work derived from receipts and source: **{summary['totals']['actual_physical_LU_from_receipts_and_source']} LU / {summary['totals']['actual_physical_RHS_from_receipts_and_source']} complex RHS**, including six initial forward solves omitted from branch counters. The projector was constructed nine times, hence nine H Cholesky factorizations. Full projection-action counts are preserved per case in `combined_offline_scores.json`.",'',
        'This study does not establish an end-to-end speedup, improved generalization, or superiority over published G19/DR. It starts after a paid 3000-step AP preparation and uses only three exposed objects. In C, published DR from the original NN material reaches NRMSE 0.077988; the fractional and integer event continuations here reach 0.133181 and 0.128761. In A, ordinary and event continuations all recover part of the shell, without the published physical corrector. B remains a negative result.','',
        '## Outputs','',
        '`combined_offline_scores.json`, `comparison.csv`, `finite_thresholds.csv`, `physical_traces.csv`, and the three PNG/PDF figure pairs under `figures/`. The complete plotting/scoring program is `work/plot_normal_memory.py`.']
(WORK/'MECHANISM_EXPERIMENT_SUMMARY.md').write_text('\n'.join(lines)+'\n')
print(wide.to_string())
print(json.dumps(summary['totals'],indent=2))

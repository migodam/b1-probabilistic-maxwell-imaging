"""Post-scoring initializer survival and NN/current/readout interactions.

No private truth reads, metric recomputation, solver calls, or GT-based selection.
Requires both public seals (750+126), original scoring, and explicit root readiness.
Optional anchor scores and public projection covariates are joined, never tuned.
"""
from pathlib import Path
from collections import defaultdict
import argparse,csv,json,math,time
BEGIN_CPU=time.process_time();BEGIN_WALL=time.perf_counter()
import numpy as np
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
OUT=ROOT/'delegated/warm_mechanism_analysis'
CELLS={'cie_dual':('sourcewise','ordinary'),'global_dual':('global','ordinary'),
       'a_plus_sourcewise':('sourcewise','A+'),'a_plus':('global','A+')}
METRICS=('material_nrmse','otsu_iou','original_Phi')
COSTS=('solver_cpu_seconds','warm_deployment_cpu_seconds','actual_workflow_cpu_seconds')


def load(p):return json.loads(Path(p).read_text())
def resolve(p):
    p=Path(p);return p if p.is_absolute() else REPO/p

def finite(v):return v is not None and math.isfinite(v)
def diff(a,b):return a-b if finite(a) and finite(b) else None
def less(a,b):return a<b if finite(a) and finite(b) else None
def greater(a,b):return a>b if finite(a) and finite(b) else None
def all_known_and(*v):return None if any(x is None for x in v) else all(v)
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csvwrite(name,rows):
    with (OUT/name).open('w',newline='') as f:
        fields=list(dict.fromkeys(k for r in rows for k in r))
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)


def stamps_unchanged(seal):
    for x in seal.get('evidence',[])+[x for e in seal['entries'] for x in e['files']]:
        p=resolve(x['path']);st=p.stat()
        if st.st_size!=x['bytes'] or st.st_mtime_ns!=x['mtime_ns']:raise RuntimeError('Changed public seal file '+str(p))


def gated_inputs():
    primary=load(ROOT/'audit/WARM_ALL_PUBLIC_SEAL.json');anchor=load(ROOT/'audit/WARM_ANCHOR_PUBLIC_SEAL.json')
    for seal,status,n in [(primary,'ALL_FROZEN_WARM_TASKS_SEALED',750),(anchor,'ALL_ANCHOR_TASKS_SEALED',126)]:
        if seal.get('status')!=status or seal.get('tasks')!=n or len(seal['entries'])!=n:raise RuntimeError('Both complete 876-task public seals required')
        stamps_unchanged(seal)
    receipt=load(ROOT/'RESULTS/warm_scoring/RECEIPT.json')
    if receipt['status']!='SCORED_ALL_FROZEN_PREDICTIONS':raise RuntimeError('Original750 independent scoring must already exist')
    original=load(ROOT/'RESULTS/warm_scoring/SCORES.json')['rows']
    tt=load(ROOT/'configs/warm_all_manifest.json')['tasks']
    at=load(ROOT/'configs/warm_anchor_confirmation_manifest.json')['tasks']
    for seal,tasks in ((primary,tt),(anchor,at)):
        identities=[e['task']['task_id'] for e in seal['entries']]
        if len(identities)!=len(set(identities)) or set(identities)!={t['task_id'] for t in tasks}:
            raise ValueError('Exact frozen public seal identities required')
    if len(original)!=750 or {r['task_id'] for r in original}!={t['task_id'] for t in tt}:raise ValueError('Exact original750 scored identities required')
    cases=load(ROOT/'data_public/warm_test/test_manifest.json')['records']
    counts={c:sum(r['cohort']==c for r in cases) for c in {r['cohort'] for r in cases}}
    if len(cases)!=42 or sorted(counts.values())!=[6,12,24]:raise ValueError('Expected separate 24/12/6 cohorts')
    return original,tt,anchor,cases


def public_state(row,task,initial_fallback=None):
    d=resolve(task['output_dir']);rec=load(d/'RECEIPT.json')
    if rec['task_id']!=row['task_id'] or rec['status']!=row['status']:raise ValueError('Scored/public receipt identity or status mismatch')
    candidates=rec.get('auxiliary_candidate_scores') or []
    stage=rec.get('auxiliary_selected_stage');initial=next((c.get('objective') for c in candidates if c.get('stage')==-1),None)
    chosen=next((c for c in candidates if c.get('stage')==stage),{})
    if (d/'initial.json').exists():initial=load(d/'initial.json').get('objective',initial)
    pre=load(d/'selected_before_corrector.json').get('objective') if (d/'selected_before_corrector.json').exists() else chosen.get('objective')
    r=dict(row,initial_Phi=initial,selected_pre_corrector_Phi=pre,chosen_stage=stage,
           chosen_stage_reason=chosen.get('reason'),stage_minus1_selected=(stage==-1) if stage is not None else None,
           corrector_stop=(rec.get('corrector') or {}).get('stop'),projected_gradient_inf=row.get('projected_gradient_inf'),
           final_metric_source='Existing independently scored SCORES rows; inherited definitions, no recomputation')
    for key in ('initial_material_nrmse','initial_otsu_iou'):
        if r.get(key) is None and initial_fallback is not None:
            if task.get('initialization_file')!=initial_fallback[1].get('initialization_file'):
                raise ValueError('Cannot reuse initial material metric from a different initializer file')
            r[key]=initial_fallback[0].get(key)
            r['initial_metric_source']='Existing original same-material NN score; identical frozen initializer path'
    # Public saved state algebra only: no forward/adjoint, no test truth.
    if r['projected_gradient_inf'] is None and rec['status']=='SEALED_PUBLIC_PREDICTION':
        with np.load(d/'prediction.npz',allow_pickle=False) as z:
            if 'gradient' in z.files:
                m=z['m'].ravel();g=z['gradient'].ravel()
                r['projected_gradient_inf']=float(np.linalg.norm(m-np.clip(m-g,.02,20),np.inf))
                r['projected_gradient_source']='Saved public m/gradient, inherited [.02,20] projected-gradient algebra; 0 LU'
    r['final_minus_initial_material_nrmse']=diff(r.get('material_nrmse'),r.get('initial_material_nrmse'))
    r['pre_minus_initial_Phi']=diff(pre,initial)
    return r


def direct_comparison(row,direct):
    r=dict(row)
    if direct is None:
        r['direct_comparison_status']='MISSING_DIRECT_SCORE';return r
    r['direct_task_id']=direct['task_id']
    for key in METRICS:r['minus_NNdirect_'+key]=diff(r.get(key),direct.get(key))
    same=r['hardware']==direct['hardware'];r['direct_comparison_hardware']='SAME' if same else 'DIFFERENT_NO_ABSOLUTE_CPU_DIFFERENCE'
    for key in COSTS:r['extra_vs_NNdirect_'+key]=diff(r.get(key),direct.get(key)) if same else None
    r['selected_aux_stage_lower_pre_Phi_than_initial']=(all_known_and(r.get('chosen_stage')!=-1,less(r.get('selected_pre_corrector_Phi'),r.get('initial_Phi')))
                                                       if r.get('chosen_stage') is not None else None)
    r['selected_pre_Phi_lower_than_NNdirect_pre_Phi']=less(r.get('selected_pre_corrector_Phi'),direct.get('selected_pre_corrector_Phi'))
    r['final_Phi_worse_than_NNdirect']=greater(r.get('original_Phi'),direct.get('original_Phi'))
    r['final_material_worse_than_NNdirect']=greater(r.get('material_nrmse'),direct.get('material_nrmse'))
    r['lower_pre_Phi_but_both_final_Phi_and_material_worse_than_NNdirect']=all_known_and(
        (r.get('chosen_stage')!=-1) if r.get('chosen_stage') is not None else None,
        r['selected_pre_Phi_lower_than_NNdirect_pre_Phi'],r['final_Phi_worse_than_NNdirect'],r['final_material_worse_than_NNdirect'])
    return r


def interactions(rows,cases):
    by={(r['case_id'],r['arm'],r['initialization']):r for r in rows};out=[]
    def value(case,arm,init,metric):return by.get((case,arm,init),{}).get(metric)
    for c in cases:
        row=dict(case_id=c['case_id'],cohort=c['cohort']);case=c['case_id']
        for metric in METRICS:
            for arm in CELLS:
                row[metric+'__NN_minus_uniform__'+arm]=diff(value(case,arm,'nn',metric),value(case,arm,'uniform',metric))
            factorial={}
            for init in ('nn','uniform'):
                # (global A+ - global ordinary) - (sourcewise A+ - sourcewise ordinary).
                factorial[init]=diff(diff(value(case,'a_plus',init,metric),value(case,'global_dual',init,metric)),
                                    diff(value(case,'a_plus_sourcewise',init,metric),value(case,'cie_dual',init,metric)))
                row[metric+'__global_minus_sourcewise_Aplus_effect__'+init]=factorial[init]
            row[metric+'__NN_by_global_by_Aplus_interaction']=diff(factorial['nn'],factorial['uniform'])
        out.append(row)
    return out


def projection_join(path,cases):
    if path is None:return {},dict(status='NOT_PROVIDED; no diagnostic recomputed')
    path=resolve(path);receipt=load(path.parent/'RECEIPT.json')
    if receipt.get('truth_access')!='NONE':raise ValueError('Only existing public no-truth projection covariates may be joined')
    allowed={c['case_id'] for c in cases};data=load(path);rr=data.get('rows',data.get('records',[]));by={}
    keys=('DPCS_rank','sources','physical_current_norm','DPCS_displacement','Fourier_displacement',
          'transported_DPCS_defect','transported_Fourier_defect','signed_constitutive_cross_term',
          'full_APCS_constitutive_defect','truncated_constitutive_defect','removed_Fourier_constitutive_defect_is_helpful',
          'fixed_norm','qualification')
    for r in rr:
        if r['case_id'] in allowed:
            if r['case_id'] in by:raise ValueError('Projection covariates must have one shared record per object')
            by[r['case_id']]={'projection_'+k:r[k] for k in keys if k in r}
    return by,dict(status='JOINED_SHARED_PUBLIC_COVARIATES',path=str(path),objects=len(by),expected_objects=42,
                   interpretation='State-fixed covariates only; no online prediction, threshold search or causal attribution from associations')


def summary(rows):
    groups=defaultdict(list)
    for r in rows:groups[(r['cohort'],r['arm'],r['initialization'])].append(r)
    out=[]
    for (cohort,arm,init),rr in sorted(groups.items()):
        r=dict(cohort=cohort,arm=arm,initialization=init,all_objects=len(rr),failed=sum(x['status']=='FAILED' for x in rr))
        for key in ('stage_minus1_selected','selected_aux_stage_lower_pre_Phi_than_initial',
                    'selected_pre_Phi_lower_than_NNdirect_pre_Phi','final_Phi_worse_than_NNdirect','final_material_worse_than_NNdirect',
                    'lower_pre_Phi_but_both_final_Phi_and_material_worse_than_NNdirect'):
            vals=[x.get(key) for x in rr if x.get(key) is not None]
            r[key+'_observed']=len(vals);r[key+'_count']=sum(bool(v) for v in vals)
        # Deliberately no absolute CPU aggregation across Mac/XINAN.
        for key in METRICS:
            vals=[x[key] for x in rr if finite(x.get(key))]
            r[key+'_observed']=len(vals);r[key+'_mean']=sum(vals)/len(vals) if vals else None
        out.append(r)
    return out


def run(projection):
    original,tt,anchor_seal,cases=gated_inputs();taskmap={t['task_id']:t for t in tt}
    origby={(r['case_id'],r['arm'],r['initialization']):r for r in original}
    selected=[r for r in original if (r['arm'] in CELLS and r['initialization'] in ('nn','uniform'))
              or (r['arm']=='physical_lbfgs' and r['initialization']=='nn')
              or (r['arm']=='a_plus' and r['initialization']=='nn_joint')]
    normalized=[public_state(r,taskmap[r['task_id']]) for r in selected]
    directs={r['case_id']:r for r in normalized if r['arm']=='physical_lbfgs'}
    rows=[direct_comparison(r,directs.get(r['case_id'])) for r in normalized]
    anchor_path=ROOT/'delegated/warm_anchor_analysis/CONFIRMATION_SCORES.json'
    anchor_status='NOT_SCORED; all42 anchor slots retained as unavailable';anchor_comparisons=[]
    if anchor_path.exists():
        receipt=load(anchor_path.parent/'CONFIRMATION_RECEIPT.json')
        if receipt.get('status')!='SEPARATE126_ANCHOR_OFFLINE_SCORED':raise ValueError('Anchor score receipt incomplete')
        scores=load(anchor_path)['rows'];atasks={e['task']['task_id']:e['task'] for e in anchor_seal['entries']}
        if len(scores)!=126 or {r['task_id'] for r in scores}!=set(atasks):raise ValueError('Exact126 anchor scored identities required')
        for r in scores:
            base=origby[(r['case_id'],r['arm'],'nn')]
            ar=public_state(r,atasks[r['task_id']],(base,taskmap[base['task_id']]))
            ar=direct_comparison(ar,directs.get(r['case_id']));rows.append(ar)
            refs=[base]+([origby[(r['case_id'],'a_plus','nn_joint')]] if r['arm']=='a_plus' else [])
            for ref in refs:
                v=dict(case_id=r['case_id'],cohort=r['cohort'],arm=r['arm'],anchor_task_id=r['task_id'],reference_task_id=ref['task_id'],
                       reference_initialization=ref['initialization'],hardware=r['hardware'],reference_hardware=ref['hardware'],
                       scope='Global nn_joint controls truncated physical-current projection; other same-arm NN pairs combine initialization/affine-family changes')
                for key in METRICS:v['anchor_minus_reference_'+key]=diff(r.get(key),ref.get(key))
                for key in COSTS:v['anchor_extra_'+key]=diff(r.get(key),ref.get(key)) if r['hardware']==ref['hardware'] else None
                anchor_comparisons.append(v)
        anchor_status='126_SCORED_ANCHORS_JOINED'
    covariates,projection_status=projection_join(projection,cases)
    wide=[]
    for c in cases:
        case=c['case_id'];w=dict(case_id=case,cohort=c['cohort'],**covariates.get(case,{}))
        rr=[r for r in rows if r['case_id']==case]
        for r in rr:
            prefix=r['arm']+'__'+r['initialization']+'__'
            for key in ('status','hardware','initial_Phi','initial_material_nrmse','initial_otsu_iou','chosen_stage','stage_minus1_selected',
                        'selected_pre_corrector_Phi','original_Phi','material_nrmse','otsu_iou','corrector_stop','projected_gradient_inf',
                        'minus_NNdirect_original_Phi','minus_NNdirect_material_nrmse','extra_vs_NNdirect_solver_cpu_seconds',
                        'extra_vs_NNdirect_warm_deployment_cpu_seconds','selected_pre_Phi_lower_than_NNdirect_pre_Phi','lower_pre_Phi_but_both_final_Phi_and_material_worse_than_NNdirect'):
                w[prefix+key]=r.get(key)
        for arm in ('cie_dual','a_plus','a_plus_sourcewise'):
            if not any(r['arm']==arm and r['initialization']=='nn_affine_anchor' for r in rr):
                w[arm+'__nn_affine_anchor__status']='NOT_SCORED'
        wide.append(w)
    csvwrite('ALL42_CASE_TABLE.csv',wide);csvwrite('INITIALIZER_SURVIVAL_LONG.csv',rows)
    csvwrite('COHORT_SURVIVAL_COUNTS.csv',summary(rows));csvwrite('ALL42_NN_SOURCE_GLOBAL_INTERACTIONS.csv',interactions(rows,cases))
    csvwrite('ANCHOR_SAME_ARM_AND_GLOBAL_NN_JOINT.csv',anchor_comparisons)
    # Absolute extra CPU summaries only within each hardware and cohort.
    groups=defaultdict(list)
    for r in rows:groups[(r['cohort'],r['hardware'],r['arm'],r['initialization'])].append(r)
    costs=[]
    for (cohort,hardware,arm,init),rr in sorted(groups.items()):
        v=dict(cohort=cohort,hardware=hardware,arm=arm,initialization=init,objects=len(rr))
        for key in COSTS:
            values=[r.get('extra_vs_NNdirect_'+key) for r in rr if finite(r.get('extra_vs_NNdirect_'+key))]
            v[key+'_paired_direct_extra_observed']=len(values);v[key+'_mean_paired_direct_extra']=sum(values)/len(values) if values else None
        costs.append(v)
    csvwrite('DIRECT_EXTRA_CPU_BY_HARDWARE.csv',costs)
    dump(OUT/'MECHANISM_METADATA.json',dict(objects=42,cohorts='24 ID /12 heldout-family /6 cross-grid separately',anchor=anchor_status,
        projection=projection_status,LU=0,complex_RHS=0,truth_access='NONE: existing independently scored metrics plus public receipts/saved gradients',
        new_inverse_solves=0,threshold_search=False,metric_recomputation=False))
    (OUT/'REPORT.md').write_text('''# Post-scoring initializer survival and factorial interactions

All42 objects are retained. Cohorts24/12/6 remain separate. Material metrics and initial material metrics come exclusively from existing independent SCORES rows; no truth is opened and inherited definitions are unchanged. Failure and unavailable observations remain explicit.

Stage−1 selection means that the initial material candidate entered the corrector; it does not mean the final material stayed unchanged. InitialPhi and the selected pre-correctorPhi are public paid values. The lower-prePhi flag compares to the NNdirect branch's saved pre-corrector initialPhi and requires a selected auxiliary stage. The separate own-initialPhi flag is also retained. A lower pre-correctorPhi but worse finalPhi/material than NNdirect describes the noncommuting order of selection/correction on these paid trajectories. It is not a basin, uniqueness, barrier or causal certificate. Zero is only the sign of a saved difference; no threshold or object was tuned or selected by material quality.

The42-row wide table includes all four original sourcewise/global × ordinary/A+ arms for NN and uniform, initial/final metrics, stage, corrector stop, saved projected-gradient diagnostic and NNdirect endpoint/extra-CPU comparisons. The long table retains every relevant method record. CPU differences require identical hardware; absolute CPU summaries are stratified by cohort/hardware. Original feature/inference fields are Mac measurements; deployment prefixes are independent same-solver-hardware measurements, so those cost notions cannot be interchanged.

The per-object interaction is (global A+ − global ordinary) − (sourcewise A+ − sourcewise ordinary), computed separately for NN/uniform; their difference is the NN×global×A+ interaction. It is descriptive for this fixed finite experiment, not a newly fitted method or significance claim. Failed metrics give missing interactions, never zero substitutions.

If independently scored anchor126 outputs exist, comparisons include same-arm originalNN and the critical global A+ nn_joint control. nn_joint projects the physical initial current into the old truncated affine family; anchor keeps the full physical APCS offset. Dual/sourcewise have no corresponding truncated-physical-projection control and combine initialization/affine-space changes. Shared public projection covariates, when supplied, are state-fixed diagnostics and cannot establish causation from correlation or serve as an online success selector.

No new solver, forward, adjoint or LU call is made. Public saved-gradient algebra uses the inherited [.02,20] box only. Computational analysis cost is measured in RECEIPT.json; numerical/physical solver costs are reused, not charged again as new physical work.
''')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--authorized-post-scoring',action='store_true',help='Supply only after root explicitly announces both seals and scoring ready')
    parser.add_argument('--projection-diagnostics',type=Path,help='Existing public DIAGNOSTICS.json with sibling receipt truth_access NONE; no recomputation')
    args=parser.parse_args()
    if not args.authorized_post_scoring:parser.error('Root readiness required; analysis was not run')
    OUT.mkdir(parents=True,exist_ok=True);error=None
    try:run(args.projection_diagnostics)
    except Exception as exc:error=type(exc).__name__+': '+str(exc);raise
    finally:
        dump(OUT/'RECEIPT.json',dict(status='FAILED' if error else 'POST_SCORING_MECHANISM_ANALYSED',error=error,
             cpu_seconds=time.process_time()-BEGIN_CPU,wall_seconds=time.perf_counter()-BEGIN_WALL,
             measurement='Actual outer process including numpy import after standard-library startup',
             LU=0,complex_RHS=0,new_inverse_solves=0,truth_access='NONE: no private truth path in this analyzer'))

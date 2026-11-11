"""Freeze all independent methods/cost/selection before test dispatch or scoring."""
from pathlib import Path
from datetime import datetime,timezone
import json,subprocess


ROOT=Path(__file__).resolve().parents[1]


def main():
    stamp=datetime.now(timezone.utc).isoformat()
    source=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    protocol=dict(CPU_cap=120,proposal_fraction=.6,steps=3000,audit_every=50,
        physical_audit_CPU_interval=2.,AP_preparation=128,scale=6.,omega=1.,forecast_trace=True,
        material_box=[.02,20.],prior=dict(smooth=1e-4,tv=1e-4,epsilon=.01),
        source_commit=source,source_scope='runner plus exact deployed bundle; source Git provenance, no hash validation',
        selection='minimum original paid Phi; all native extra physical audits retained and charged',
        secondary='min original R among legal original-CN99 paid candidates, fallback primary',
        prefixes=[15,30,60,90,120],BLAS_threads=1,priority='NORMAL',P_cores=[0,2,4,6,8,10,12,14],
        corrector='unchanged physical L-BFGS-B max3000/12000; same prior, complete paid state retention',
        cadence='first complete material update after each absolute two-process-CPU-second bucket; preserve native initial/stage/final audits. Corrector/direct every paid objective evaluation.',
        truth_policy='NONE online; seal BOTH independent campaigns before any fresh material truth score',
        scope='48 independent material objects, eight hidden families, six acquisition/noise conditions; crossgrid is model-mismatch stress, not noise-qualified truth accuracy',
        parameter_policy='No fresh-result-dependent method changes, case exclusion, thresholds or budget extensions')
    methods=['cie','dual_cie','g19','g19_full_apcs','g19_full_current','dr','event_integer','random_integer','direct']
    inputs=json.loads((ROOT/'data_public/fresh48/MANIFEST.json').read_text())['case_files']
    core=[]
    for item in inputs:
        for m in methods:
            core.append(dict(task_id=item['case_id']+'__'+m,case_id=item['case_id'],method=m,algorithm=m,
                campaign='confirmation_core48',cohort='fresh48_independent',physics_kind='scalar_tm_2d',
                input_file=item['input_file'],seed=13013,observation_contract=item))
    chosen=[f'r13_fresh_{i:04}' for i in range(0,48,3)]
    extra=[]
    ordinary=['ap','anderson','a_plus','global_dual','ap_dr','online_event','online_random','nfft_equations']
    for cid in chosen:
        for m in ordinary:
            extra.append(dict(task_id=cid+'__'+m,case_id=cid,method=m,algorithm=m,campaign='confirmation_mechanism16',
                cohort='predeclared_every_third_ID16',physics_kind='scalar_tm_2d',
                input_file=f'data_public/fresh48/{cid}.npz',seed=13013))
        for m in ['dr','g19','direct']:
            extra.append(dict(task_id=cid+'__nn_'+m,case_id=cid,method='nn_'+m,algorithm=m,
                campaign='confirmation_mechanism16',cohort='predeclared_every_third_ID16',physics_kind='scalar_tm_2d',
                input_file=f'data_public/fresh16_nn/{cid}.npz',seed=13013,
                initialization='frozen R8 supervised physics NN; original observation unchanged, no new training'))
    for name,tasks in [('confirmation_core48',core),('confirmation_mechanism16',extra)]:
        data=dict(campaign=name,status='FROZEN_BEFORE_ANY_CONFIRMATION_RUN',created_at_utc=stamp,
            protocol=protocol,tasks=tasks,chosen_ID_policy='every third opaque ID, no access to material family or test truth',
            original_NFFT_qualification='nfft_equations is Eq3-6 real/prior/box/global-LBFGS adaptation, not mature author-code reproduction')
        path=ROOT/'configs'/(name+'.json')
        with path.open('x') as f:json.dump(data,f,indent=2);f.write('\n')
    record=dict(status='INDEPENDENT_ALGORITHMS_AND_ESTIMATORS_FROZEN',created_at_utc=stamp,
        source_commit=source,core_tasks=len(core),mechanism_tasks=len(extra),material_objects=48,
        selected_subset=chosen,core_methods=methods,extra_methods=ordinary+['nn_dr','nn_g19','nn_direct'],
        no_new_GT_access=True,no_confirmation_results_exist=not (ROOT/'RESULTS/confirmation_core48').exists(),
        outstanding='actual queues, acceptance, joint seals, scoring, scientific review and publication')
    (ROOT/'audit/CONFIRMATION_FREEZE.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))


if __name__=='__main__':main()

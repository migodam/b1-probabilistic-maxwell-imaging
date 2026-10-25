"""Pure saved-artifact extraction for the optional ordinary A0 polish.

No solver, Book, experiment, plotting, truth generator or top-level I/O.
"""
from pathlib import Path
import json

PARENTS=(0,3,6,9)


def read(path):return json.loads(Path(path).read_text()) if Path(path).exists() else None


def load_polish(root):
    root=Path(root);directory=root/'results/a0/adjoint_polish'
    report=read(directory/'RESULTS.json')
    output=dict(status='NOT_RUN',reason='Optional diagnostic not saved',paired_rows=[],trials=[],parent_ids=list(PARENTS),
        expected_cases=8,independent_parents=4,bootstrap='NOT_RUN_FOUR_PARENT_EXPLORATORY_DIAGNOSTIC_NO_FORMAL_CI',
        holdout_scope='EXPOSED_DEVELOPMENT_REEVALUATION_NOT_NEW_BLIND_TEST',matched_speed_claim=False,
        setup_scope='Initialization depends on all three primary arms and shared proposal bank; all preceding costs retained',
        scientific_scope='ordinary known-intrinsic optimization-error diagnostic; all original primary arms remain frozen')
    if not report:return output
    output.update(status=report.get('status','NOT_RUN'),reason=None,starting_headroom=report.get('starting_headroom'),
        budget_start_snapshot=report.get('budget_start'),budget_end_snapshot=report.get('budget_end'),
        budget_snapshot_CPU_warning='In-child snapshot may omit current wrapper CPU; final campaign ledger/job receipt remains cost authority')
    online=report.get('online',[]);by={(r['parent_id'],r['noise_repeat']):r for r in online}
    # Preserve every actual attempted trial, including failed/cached/unpaid calls.
    for p in PARENTS:
        for repeat in (0,1):
            case=directory/f'parent_{p:02d}_repeat_{repeat}';row=by.get((p,repeat),{})
            attempts=read(case/'ATTEMPTS.json') or [];start=row.get('cost_before',{}).get('rhs')
            cache={}
            for trial in attempts:
                if trial.get('objective') is not None:cache[tuple(trial['theta'])]=trial['objective']
                objective=trial.get('objective');source='SAVED_ACTUAL_COMBINED_TRIAL' if objective is not None else 'NOT_SAVED'
                if objective is None and trial.get('status')=='CACHED_COMBINED_VALUE_GRADIENT' and tuple(trial['theta']) in cache:
                    objective=cache[tuple(trial['theta'])];source='SAVED_IDENTICAL_THETA_CACHE_VALUE'
                after=trial.get('cost_after',{}).get('rhs')
                output['trials'].append(dict(parent=p,repeat=repeat,attempt=trial['attempt'],status=trial.get('status'),
                    cached=trial.get('cached'),theta=trial.get('theta'),objective=objective,objective_source=source,
                    rhs_since_case_start=after-start if after is not None and start is not None else None,
                    rhs_scope='actual reserved RHS including first native-J audit; cached/unpaid missing cost not invented'))
    terminal=read(directory/'ONLINE_COMPLETE.json')
    complete=(terminal is not None and terminal.get('status')=='ALL_ONLINE_TERMINAL_DURABLE' and
              set(by)=={(p,r) for p in PARENTS for r in (0,1)})
    output['all_online_outputs_durable']=complete
    if not complete:
        output['offline_status']='NOT_RUN_ALL_EIGHT_ONLINE_OUTPUTS_NOT_DURABLE'
        return output
    if report.get('status')!='OPTIMIZER_DIAGNOSTIC_COMPLETE_PARENT_REVIEW_REQUIRED':
        output['offline_status']='NOT_RUN_DIAGNOSTIC_SCORING_NOT_TERMINAL'
        return output
    # Only now open already-saved offline score summaries; no GT files themselves
    # are needed, and no offline information enters a choice of initialization.
    scores=read(directory/'OFFLINE_SCORES.json')
    main=read(root/'results/a0/scoring/RESULTS.json')
    if scores is None or main is None:
        output['offline_status']='NOT_RUN_SAVED_OFFLINE_SCORES_ABSENT';return output
    after={(r['parent_id'],r['noise_repeat']):r for r in scores}
    before={(r['observation_id'],r['method']):r for r in main['rows']}
    metrics=['training_noise_energy','holdout_noise_energy','pose_error_lambda0',
        'pose_squared_error_lambda0_squared','material_eps_contrast_L2','material_ell_L2',
        'training_noise_envelope_compatible','holdout_noise_envelope_compatible','both_noise_envelopes_compatible']
    for p in PARENTS:
        for repeat in (0,1):
            row=by[(p,repeat)];initial=row['saved_main_noop'];selected=row['selection'];a=after.get((p,repeat),{})
            b=before.get((row['observation_id'],initial['main_method']),{})
            paired=dict(parent=p,repeat=repeat,initial_method=initial['main_method'],status=a.get('status','NOT_RUN_MISSING_OFFLINE_SCORE'),
                fit_status=row['status'],selected_origin=selected.get('origin'),initial_main_fit_success=initial.get('main_fit_success'),
                objective_before=initial['objective'],objective_after=selected['objective'],
                objective_change=selected['objective']-initial['objective'],
                unique_joint_evaluations=row.get('unique_joint_evaluations'),attempts=row.get('attempts'),
                optimizer_nfev=(row.get('solver') or {}).get('nfev'),optimizer_njev=(row.get('solver') or {}).get('njev'),
                optimizer_success=(row.get('solver') or {}).get('success'),
                case_RHS=row['cost_after']['rhs']-row['cost_before']['rhs'],
                case_factorizations=row['cost_after']['factorizations']-row['cost_before']['factorizations'],
                nested_polish_process_cpu_seconds=row.get('process_cpu_seconds'),
                CPU_scope='measured nested polish process interval; not complete standalone campaign/setup CPU',
                native_J_audit=(row.get('gradient_audit') or {}).get('passed'),
                native_J_audit_status='SAVED' if row.get('gradient_audit') else 'NOT_RUN_BY_FIRST_CASE_ONLY_DESIGN',
                holdout_scope=output['holdout_scope'],init_depends_on_all_three_main_arms=True,matched_speed_claim=False)
            for metric in metrics:
                paired[metric+'_before']=b.get(metric);paired[metric+'_after']=a.get(metric)
                bv,av=b.get(metric),a.get(metric)
                paired[metric+'_change']=float(av)-float(bv) if av is not None and bv is not None else None
            output['paired_rows'].append(paired)
    output['offline_status']='SAVED_TERMINAL_OFFLINE_PAIR_SUMMARY'
    output['all_status_denominators']={status:sum(r['status']==status for r in output['paired_rows']) for status in sorted(set(r['status'] for r in output['paired_rows']))}
    return output

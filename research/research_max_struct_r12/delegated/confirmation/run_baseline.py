"""R12 output-owned copy of R11 one-task wrapper; scientific runner remains inherited.
Changes: output/legacy paths and observational144/96 CPU phase caps; original
iteration allowances and accepted numerical trajectories are unchanged.
"""
import os
for _name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name] = '1'
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace
ROOT = Path(__file__).resolve().parent
REPO = Path(os.environ['R12_LEGACY_REPO'] if os.environ.get('R12_LEGACY_REPO') else ROOT.parents[3]).resolve()
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO/'research/research_max_struct_r11/src'))

class ProposalCPUReached(Exception):
    """Stop auxiliary work at a native numerical prefix, then hand off."""


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest', type=Path, required=True)
    ap.add_argument('--config', type=Path, required=True)
    ap.add_argument('--task-id', required=True)
    ap.add_argument('--cpu', type=int)
    ap.add_argument('--vendor-path')
    args = ap.parse_args()
    if args.cpu is not None:
        import psutil
        psutil.Process().cpu_affinity([args.cpu])
        if os.name == 'nt': psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    if args.vendor_path: sys.path.insert(0, args.vendor_path)
    from twod_budget import WholeChildBudget, CPUBudgetReached
    from research.research_max_struct_r8.experiments import run_case as runner
    from research.research_max_struct_r8.src import benchmark_solver as solver
    from research.research_max_struct_r10.src.twod_adapter import noise_ball_rank, require_iid_raw_noise
    from threadpoolctl import threadpool_limits
    import numpy as np
    cfg = json.loads(args.config.read_text())
    manifest = json.loads(args.manifest.read_text())
    task = dict(runner.choose_task(manifest, args.task_id))
    if cfg.get('status') not in ('FROZEN_FOR_INVERSION','FROZEN_BEFORE_GENERATION_AND_CONFIRMATION_GT_ACCESS'): raise ValueError('Root has not frozen contract')
    budget = WholeChildBudget(cfg['whole_child_CPU_budget_seconds'])
    out = runner.resolve(task['output_dir']).resolve()
    if not out.is_relative_to(ROOT/'RESULTS'): raise ValueError('R12 confirmation-owned output only')
    if out.exists(): raise FileExistsError(out)
    prepared_key=task.get('engineering_attempt_id',task['task_id'])
    if not prepared_key.replace('_','').isalnum(): raise ValueError('Invalid engineering preparation key')
    prepared = ROOT/'delegated/twod/adapter_inputs'/prepared_key
    prepared.mkdir(parents=True, exist_ok=False)
    old_observer, old_dispatch, old_refine = solver.PaidStateObserver, solver.run_auxiliary, solver.refine_observed
    NativeBudgetObserver=budget.observer_class(old_observer)
    secondary_enabled=cfg.get('secondary_selection')=='CN99_MIN_SAME_PRIOR_AUDIT_ONLY_V1'
    secondary_delta=None
    if secondary_enabled and not (ROOT/'SECONDARY_SELECTION_FREEZE.json').is_file():
        raise ValueError('Secondary estimator must be registered before confirmation dispatch')
    class PhaseObserver(NativeBudgetObserver):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs)
            self.secondary_best=None;self.secondary_incumbents=[];self.secondary_audits=[]
        def record(self,*args,**kwargs):
            row=super().record(*args,**kwargs)
            row['paid_whole_child_CPU_seconds']=self.paid_states[-1]['CPU_seconds']
            if secondary_enabled:
                # Audit only an already computed physical state. No optimizer,
                # current/material iterate, primary Phi incumbent or FW changes.
                residual=float(np.linalg.norm(row['Y_prediction']-self.Y))
                prior=float(self.model.reg.evaluate(row['m'])[0])
                legal=bool(np.isfinite(row['m']).all() and np.min(row['m'])>=.02 and np.max(row['m'])<=20.)
                score=(prior,float(row['objective']))
                proposed=legal and residual<=secondary_delta and (self.secondary_best is None or score<self.secondary_best['score'])
                snapshot={k:row[k].copy() for k in ('m','j','Q','Y_prediction')} if proposed else None
                completed=time.process_time()
                deadline=144. if self.phase=='auxiliary' else budget.seconds
                feasible=residual<=secondary_delta;eligible=completed<=deadline
                rec=dict(CPU_seconds=completed,phase=self.phase,R=prior,Phi=row['objective'],
                    raw_residual_norm=residual,delta_CN99=secondary_delta,noise_feasible=feasible,legal_material=legal,
                    phase_cap_eligible=eligible,phase_deadline_CPU_seconds=deadline)
                self.secondary_audits.append(rec)
                if proposed and eligible:
                    self.secondary_best=dict(score=score,row=snapshot,metadata=rec.copy())
                    self.secondary_incumbents.append(dict(rec,m=snapshot['m']))
            return row
    solver.PaidStateObserver = PhaseObserver
    phase=dict(proposal_max_CPU_seconds=float(cfg.get('proposal_max_CPU_seconds',144.)),
               corrector_max_CPU_seconds=float(cfg.get('corrector_max_CPU_seconds',96.)),
               whole_child_max_CPU_seconds=240.,proposal_cut=False,excluded_atomic_proposal_states=0)
    if phase['proposal_max_CPU_seconds']!=144. or phase['corrector_max_CPU_seconds']!=96. or budget.seconds!=240.:
        raise ValueError('Root-frozen phase caps are144/96 with global240')
    def proposal_check():
        budget.check()
        if time.process_time()>=phase['proposal_max_CPU_seconds']:
            raise ProposalCPUReached()
    def refine_prefix(model,Y,m,*args,**kwargs):
        start=time.process_time()
        phase.update(corrector_start_whole_child_CPU_seconds=start,
            effective_corrector_deadline_CPU_seconds=min(240.,start+96.))
        # The inherited observer now applies the earlier corrector deadline to
        # completed paid-state eligibility, including atomic overrun exclusion.
        budget.seconds=phase['effective_corrector_deadline_CPU_seconds']
        return old_refine(model,Y,m,*args,**kwargs)
    solver.refine_observed=refine_prefix
    receipt = None
    try:
        with threadpool_limits(limits=1):
            budget.check()
            # The inherited public loader rejects private paths and truth-like NPZ keys.
            public = runner.load_public(task['public_files'])
            sigma = require_iid_raw_noise(public, task, cfg)
            if secondary_enabled:
                from scipy.stats import chi2
                secondary_delta=float(sigma*np.sqrt(chi2.ppf(.99,2*public['Y'].size)/2))
            decision = noise_ball_rank(public['S'], public['Y'], sigma, prepared/'RECEIVER_DPCS.npz')
            runner.dump(prepared/'DPCS_RULE.json', decision)
            if decision['rank'] is None: raise ValueError(decision['status'])
            budget.check()
            cfg['common'].pop('dpcs_rtol', None)
            cfg['common']['dpcs_rank'] = decision['rank']
            def dispatch(model, Y, m0, arm, numerical):
                numerical = dict(numerical)
                numerical['budget_check'] = proposal_check
                original_callback = numerical.get('iteration_callback')
                prefix_history=[]
                def iteration(row):
                    if original_callback: original_callback(row)
                    prefix_history.append(dict(stage=row['stage'],iteration=row['iteration'],
                        history_scope='Native accepted iteration prefix; full loss history unavailable after phase interruption'))
                    proposal_check()
                numerical['iteration_callback'] = iteration
                custom = task.get('material_readout')
                result=None
                try:
                    proposal_check()
                    if custom:
                        path = runner.resolve(custom['module_file']).resolve()
                        if path != (REPO/'research/research_max_struct_r11/src/material_readouts.py').resolve(): raise ValueError('Only inherited fixed R11 readout allowed')
                        spec = importlib.util.spec_from_file_location('r11_frozen_candidate', path)
                        mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
                        if custom.get('prepare'): numerical = getattr(mod, custom['prepare'])(model, numerical)
                        numerical['material_readout'] = getattr(mod, custom['callable'])
                        proposal_check()
                        result=solver.run_cie_fft(model, Y, m0, numerical)
                    else:result=old_dispatch(model, Y, m0, arm, numerical)
                except ProposalCPUReached as exc:
                    phase.update(proposal_cut=True,proposal_stop_whole_child_CPU_seconds=time.process_time())
                    execution=getattr(exc,'cie_execution_receipt',{})
                    execution.update(status='PROPOSAL_CPU_PREFIX_HANDOFF',
                        stopping='Native iteration boundary at root144CPU proposal allowance; no shortened mature iteration contract',
                        core_interruption_status='Historical core catches phase exception as FAILED; wrapper reclassifies intentional cutoff')
                    runner.dump(out/'AUXILIARY_EXECUTION.json',execution)
                    result=dict(status='PROPOSAL_CPU_PREFIX',history=prefix_history,
                        cost=dict(execution.get('operation_counts',{})),execution_receipt=execution)
                observer=budget.observer
                all_rows=list(observer.validation_records)
                eligible=[row for row in all_rows if row['paid_whole_child_CPU_seconds']<=144.]
                phase['excluded_atomic_proposal_states']=len(all_rows)-len(eligible)
                if not eligible:raise RuntimeError('No completed physical proposal endpoint within144CPU')
                if phase['proposal_cut'] or len(eligible)!=len(all_rows):
                    result['candidates']=[dict(m=row['m'].copy(),objective=row['objective'],
                        physical_relative_residual=row['physical_relative_residual'],stage=i-1,
                        reason='Paid true-physics endpoint before proposal144CPU') for i,row in enumerate(eligible)]
                    result['m']=min(eligible,key=lambda row:row['objective'])['m'].copy()
                observer.validation_records=eligible
                observer.best=observer.eligible_best=min(eligible,key=lambda row:row['objective'])
                observer.incumbents=[row for row in observer.incumbents if row[0]<=144.]
                for row in observer.paid_states:
                    row['eligible_within_proposal_phase']=row['CPU_seconds']<=144.
                phase['proposal_end_whole_child_CPU_seconds']=time.process_time()
                return result
            solver.run_auxiliary = dispatch
            runner.dump(prepared/'EFFECTIVE_CONFIG.json', cfg)
            runner.dump(prepared/'MANIFEST.json', dict(tasks=[task]))
            receipt, _ = runner.execute(SimpleNamespace(manifest=str(prepared/'MANIFEST.json'),
                public_list=None, task_id=args.task_id, config=str(prepared/'EFFECTIVE_CONFIG.json'), output_dir=str(out)))
    except CPUBudgetReached:
        if (out/'RECEIPT.json').exists(): receipt = json.loads((out/'RECEIPT.json').read_text())
        else: receipt = dict(task_id=args.task_id, truth_access='NONE', scoring='NOT_RUN')
        receipt['status'] = 'CPU_BUDGET_REACHED'
        auxiliary_receipt = out/'AUXILIARY_EXECUTION.json'
        if auxiliary_receipt.exists():
            partial = json.loads(auxiliary_receipt.read_text())
            if partial.get('status') == 'RUNNING':
                partial['status'] = 'INTERRUPTED_CPU_BUDGET'
                runner.dump(auxiliary_receipt, partial)
    except Exception as exc:
        if (out/'RECEIPT.json').exists(): receipt = json.loads((out/'RECEIPT.json').read_text())
        else: receipt = dict(task_id=args.task_id, truth_access='NONE', scoring='NOT_RUN')
        receipt.update(status='FAILED_R11_SETUP', error=repr(exc))
    finally:
        solver.PaidStateObserver, solver.run_auxiliary, solver.refine_observed = old_observer, old_dispatch, old_refine
        if receipt is None: receipt = dict(status='FAILED_NO_RECEIPT')
        observer = getattr(budget, 'observer', None)
        best = observer.eligible_best if observer is not None else None
        out.mkdir(parents=True, exist_ok=True)
        if secondary_enabled:
            secondary=getattr(observer,'secondary_best',None) if observer is not None else None
            info=dict(status='SELECTED_NOISE_FEASIBLE_MIN_SAME_PRIOR' if secondary is not None else 'NO_CAP_ELIGIBLE_NOISE_FEASIBLE_PAID_STATE',
                rule='raw residual<=globalCN99 radius; minimize original R(m), exact tie by original Phi',
                primary_selection='UNCHANGED_MIN_ORIGINAL_PHI',feedback_into_iterations='NONE',
                new_forward=0,new_GT=0,prior_evaluations=len(getattr(observer,'secondary_audits',[])),
                estimator_scope='Pre-registered secondary candidate selector, not a pure DR mechanism benefit',
                freeze_receipt=str(ROOT/'SECONDARY_SELECTION_FREEZE.json'))
            if secondary is not None:
                np.savez_compressed(out/'NOISE_FEASIBLE_MIN_PRIOR.npz',**secondary['row'])
                info['selection']=secondary['metadata']
            runner.dump(out/'NOISE_FEASIBLE_MIN_PRIOR.json',info)
            runner.dump(out/'NOISE_FEASIBLE_MIN_PRIOR_HISTORY.json',getattr(observer,'secondary_audits',[]))
            prefix={};prefixrows=[]
            for cap in (30,60,120,240):
                choices=[r for r in getattr(observer,'secondary_incumbents',[]) if r['CPU_seconds']<=cap]
                if choices:
                    selected=min(choices,key=lambda r:(r['R'],r['Phi']))
                    prefix['m_'+str(cap)]=selected['m']
                    prefixrows.append(dict(cap=cap,**{k:v for k,v in selected.items() if k!='m'}))
                else:prefixrows.append(dict(cap=cap,status='NO_CAP_ELIGIBLE_NOISE_FEASIBLE_PAID_STATE'))
            np.savez_compressed(out/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.npz',**prefix)
            runner.dump(out/'NOISE_FEASIBLE_MIN_PRIOR_TIMEPREFIX.json',prefixrows)
        if observer is not None:
            for label, state in [('ACTUAL_PAID_INCUMBENT',observer.best),('LAST_PAID_STATE',observer.last_paid_state)]:
                if state is not None:
                    np.savez_compressed(out/('R11_'+label+'.npz'), **{k:state[k] for k in ('m','j','Q','Y_prediction','gradient') if k in state})
                    runner.dump(out/('R11_'+label+'.json'), {k:v for k,v in state.items() if k not in ('m','j','Q','Y_prediction','gradient')})
        if best is not None:
            np.savez_compressed(out/'R11_CAP_PREDICTION.npz', **{k:best[k] for k in ('m','j','Q','Y_prediction')})
            runner.dump(out/'R11_CAP_PREDICTION.json', {k:v for k,v in best.items() if k not in ('m','j','Q','Y_prediction','gradient')})
        runner.dump(out/'R11_PAID_PHYSICAL_STATES.json', observer.paid_states if observer else [])
        points = cfg.get('CPU_budget_checkpoints_seconds', [60,120,240,480])
        checkpoints = []
        materials = []
        if observer is not None:
            for point in points:
                eligible = [r for r in observer.incumbents if r[0] <= point]
                if eligible:
                    selected = min(eligible, key=lambda r:r[2])
                    checkpoints.append(dict(CPU_cap_seconds=point, state_paid_CPU_seconds=selected[0],
                        original_Phi=selected[2], phase=selected[3], physical_model_cost=selected[4], material_index=len(materials)))
                    materials.append(selected[1])
                else: checkpoints.append(dict(CPU_cap_seconds=point, status='NO_PAID_STATE_WITHIN_CAP'))
            np.savez_compressed(out/'R11_INCUMBENT_TRAJECTORY.npz',
                CPU_seconds=np.array([r[0] for r in observer.incumbents]),
                m=np.array([r[1] for r in observer.incumbents]),
                original_Phi=np.array([r[2] for r in observer.incumbents]),
                LU=np.array([r[4].get('LU',0) for r in observer.incumbents]),
                complex_RHS=np.array([r[4].get('actual_complex_RHS',0) for r in observer.incumbents]))
        np.savez_compressed(out/'R11_CHECKPOINT_MATERIALS.npz',m=np.array(materials))
        runner.dump(out/'R11_CPU_CHECKPOINTS.json', checkpoints)
        receipt.update(study='R12', inherited_wrapper='R11 twod_run; R11 artifact prefix retained', whole_child_CPU_budget_seconds=240.,
            phase_CPU_contract=phase,effective_paid_state_CPU_deadline=budget.seconds,
            whole_child_CPU_seconds=time.process_time(), cap_stop_CPU_seconds=budget.stop_cpu,
            CPU_overrun_seconds=max(0.,time.process_time()-budget.seconds),
            cap_prediction_status='SEALED_PAID_INCUMBENT' if best is not None else 'NO_PAID_STATE_WITHIN_CAP',
            cap_physical_objective=best['objective'] if best is not None else None,
            cost_scope='Absolute child process_time includes imports, DPCS, inherited setup/aux/corrector, failures and serialization; final external GetProcessTimes authoritative. Nested timers never additive.',
            cap_policy='Only completed paid physical states at or below CPU cap; atomic overrun charged and excluded from cap prediction. No new final forward.',
            truth_access='NONE', scoring='NOT_RUN')
        runner.dump(out/'RECEIPT.json', receipt)
    print(json.dumps(dict(task_id=args.task_id,status=receipt['status'],cap_prediction_status=receipt['cap_prediction_status'])), flush=True)
    return 1 if receipt['status'].startswith('FAILED') else 0

if __name__ == '__main__': raise SystemExit(main())

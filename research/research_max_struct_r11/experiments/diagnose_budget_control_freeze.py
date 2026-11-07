"""Freeze all32 public CPU/stop-triggered controls before any new solve or GT score."""
from pathlib import Path
from datetime import datetime,timezone
import json
ROOT=Path(__file__).absolute().parents[1];REPO=ROOT.parents[1];BASE=ROOT/'delegated/diagnostics/bothblocks_confirmation';OUT=BASE/'posthoc_budget_control'
def main():
 OUT.mkdir(parents=True,exist_ok=True);path=OUT/'FROZEN_MANIFEST.json'
 if path.exists():raise FileExistsError('Rule already frozen')
 seal=json.loads((BASE/'ALL80_TERMINAL_PREDICTION_SEAL.json').read_text());ix={(e['task']['case_id'],e['task']['method']):e for e in seal['entries']};rows=[]
 for cid in sorted({e['task']['case_id'] for e in seal['entries']}):
  target=ix[cid,'fullG19']['terminal_job']['child_process_CPU_seconds']
  for method in ('cie_single','cie_dual'):
   e=ix[cid,method];task=e['task'];d=REPO/task['output_dir'];r=json.loads((d/'RECEIPT.json').read_text());old=e['terminal_job']['child_process_CPU_seconds'];stop=(r.get('corrector') or {}).get('stop');pgtol=bool((r.get('corrector') or {}).get('success') and stop and 'NORM OF PROJECTED GRADIENT' in stop.upper());trigger=bool(old<target and not pgtol)
   reason='TRIGGER_RESET_PHYSICAL_LBFGS_WITH_REMAINING_CPU' if trigger else 'ORIGINAL_CPU_AT_OR_ABOVE_G19_TARGET' if old>=target else 'ORIGINAL_PGTOL_CONVERGENCE_NO_IDLE_FILL'
   rows.append(dict(pair_id=cid+'__'+method+'__budget_to_G19',case_id=cid,method=method,classification='POSTHOC_BUDGET_CONTROL',target_G19_actual_CPU_seconds=target,original_actual_child_CPU_seconds=old,remaining_new_child_CPU_allowance_seconds=max(0.,target-old),original_stop=stop,original_stop_success=(r.get('corrector') or {}).get('success'),original_status=r['status'],original_is_PGTOL_normal_stop=pgtol,trigger=trigger,reason=reason,
    original_prediction=e['prediction_path'],original_prediction_meta=e['prediction_meta'],original_Phi=e['cap_original_Phi'],public_files=task['public_files'],output_dir=str(OUT/'runs'/(cid+'__'+method))))
 cfg=dict(status='FROZEN_PUBLIC_RULE_BEFORE_ANY_CONTINUATION',frozen_utc=datetime.now(timezone.utc).isoformat(),classification='POSTHOC_BUDGET_CONTROL_NOT_NEW_CONFIRMATION_ALGORITHM',pairs=rows,
   rule='For every16x(CIE,Dual), T=that object fullG19 exact actual child CPU. Trigger iff original exact CPU<T and not successful PGTOL stop. No GT/candidate performance read for triggers.',
   numerical='r6.refine_observed unchanged R/realbox/scale1e6/ftol1e-12/gtol.01/maxls30; maxiter1000000 only as nonbinding guard; absolute new-child CPU cap=remaining includes imports/setup/output',
   history='L-BFGS history reset from original sealed minPhi m, not seamless continuation',selection='Keep original candidate; accept lower public originalPhi only from newly paid state completed within remaining CPU. Atomic/serialization/exit overrun charged but never used to backfill eligibility.',
   cost='Actual new child wait4 CPU plus original actual child CPU; report unused and overrun. G19 target is allowed matching reference, not equal natural actual spend.',
   no_new_scoring_GT_until='All triggered continuation predictions sealed; GT already exposed but not used for rule/selection',new_threshold_scan=False,original80_and48_unchanged=True)
 path.write_text(json.dumps(cfg,indent=2)+'\n');print(json.dumps(dict(pairs=len(rows),triggered=sum(x['trigger'] for x in rows),remaining_CPU_allowance_sum=sum(x['remaining_new_child_CPU_allowance_seconds'] for x in rows if x['trigger']),trigger_ids=[x['pair_id'] for x in rows if x['trigger']]),indent=2))
if __name__=='__main__':main()

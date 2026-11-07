"""Compile completed R11 exposed DEV actual CPU receipts; no arrays or GT."""
import json,csv,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
D=ROOT/'delegated/twod'
def main():
 started=time.process_time();rows=[];controllers=[];missing=[]
 for label in ('replay','lagged'):
  p=D/('dev_'+label+'_queue/RECEIPT.json')
  if not p.exists():missing.append(str(p.relative_to(ROOT)));continue
  x=json.loads(p.read_text());controllers.append(dict(queue=label,CPU_seconds=x['controller_CPU_seconds'],status=x['status']))
  for j in x['jobs']:
   rpath=ROOT.parents[1]/j['output_dir']/'RECEIPT.json'
   r=json.loads(rpath.read_text()) if rpath.exists() else {}
   rows.append(dict(task_id=j['task_id'],role=label,status=j['status'],child_status=j.get('child_status'),external_child_CPU_seconds=j.get('child_process_CPU_seconds'),wall_seconds=j['wall_seconds'],LU=r.get('physical_model_cost',{}).get('LU'),complex_RHS=r.get('physical_model_cost',{}).get('actual_complex_RHS'),original_Phi=r.get('cap_physical_objective'),auxiliary_iterations=r.get('auxiliary_iterations'),corrector_iterations=r.get('corrector',{}).get('iterations'),corrector_stop=r.get('corrector',{}).get('stop')))
 transport=[]
 p=D/'dev_transport/TRANSPORT_COSTS.jsonl'
 if p.exists():
  all_transport=[json.loads(line) for line in p.read_text().splitlines() if line.strip()]
  transport=[]
  for row in all_transport:
   transport.append(row)
   if row['label']=='fetch_outputs':break  # Freeze DEV boundary before independent48 work.
 packaging=[]
 for p in (D/'dev_transport').glob('package_completed_*.json'):
  x=json.loads(p.read_text());packaging.append(dict(source=p.name,CPU_seconds=x['packaging_CPU_seconds']))
 extraction=[]
 for p in (D/'dev_transport').glob('extract_completed_*.json'):
  x=json.loads(p.read_text());extraction.append(dict(source=p.name,CPU_seconds=x['CPU_seconds']))
 verification=[]
 if (D/'REPLAY_COMPARISON.json').exists():
  x=json.loads((D/'REPLAY_COMPARISON.json').read_text());verification.append(dict(source='REPLAY_COMPARISON',CPU_seconds=x['CPU_seconds']))
 hook= json.loads((D/'HOOK_CHECK_RECEIPT.json').read_text()) if (D/'HOOK_CHECK_RECEIPT.json').exists() else {}
 totals=dict(hook_contract_check_internal_CPU_lower_bound_seconds=hook.get('CPU_seconds',0),inversion_child_CPU_seconds=sum(r['external_child_CPU_seconds'] or 0 for r in rows),queue_controller_CPU_seconds=sum(r['CPU_seconds'] for r in controllers),local_transport_controller_CPU_seconds=sum(r['controller_CPU_seconds'] for r in transport),local_transport_child_CPU_seconds=sum(r['child_CPU_seconds'] for r in transport),remote_packaging_CPU_seconds=sum(r['CPU_seconds'] for r in packaging),local_recorded_extraction_CPU_seconds=sum(r['CPU_seconds'] for r in extraction),read_only_replay_verification_CPU_seconds=sum(r['CPU_seconds'] for r in verification))
 result=dict(status='COMPLETED_KNOWN_CPU_LOWER_BOUND' if not missing and all(q['status']=='COMPLETED' for q in controllers) else 'PARTIAL',totals=totals,known_CPU_lower_bound_seconds=sum(totals.values()),inversions=rows,controllers=controllers,transport=transport,remote_packaging=packaging,local_extraction=extraction,verification=verification,unknown_costs=['Initial source packaging/config preparation local CPU not recorded','Initial remote extraction/vendor-copy CPU not recorded','Remote read-only status/probe process CPU not recorded','First completed_2 archive local extraction CPU not recorded','Initial pre-accounted inventory probe CPU not recorded','Final DEV remote output packaging and local full extraction CPU not recorded','Hook self-check final tail and replay verifier imports/final tail not recorded'],missing=missing,cost_rule='Disjoint actual child GetProcessTimes + controllers + separately recorded transfer/packaging/extraction/verification; internal R8/R11/aux/corrector timers overlap and are not added. Same-process latest controller CPU is one lifetime total per queue.',truth_access='NONE',GT_scoring='NOT_RUN')
 (D/'TWOD_DEV_COSTS.json').write_text(json.dumps(result,indent=2)+'\n')
 if rows:
  with (D/'TWOD_DEV_COSTS.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 lines=['# Exposed DEV actual CPU costs','','Known CPU lower bound: '+str(result['known_CPU_lower_bound_seconds'])+' seconds. Status: '+result['status']+'.','','|Task|Actual child CPU seconds|LU|Complex RHS|Stop|','|---|---:|---:|---:|---|']
 for r in rows:lines.append('|'+r['task_id']+'|'+str(r['external_child_CPU_seconds'])+'|'+str(r['LU'])+'|'+str(r['complex_RHS'])+'|'+str(r['corrector_stop'])+'|')
 lines+=['','No GT scoring. No budget padding. Child actual CPU includes setup/failure/finalization. Unknown control/setup overheads remain explicit in JSON. Nested internal timers are not charged twice.']
 (D/'TWOD_DEV_COSTS.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps(dict(status=result['status'],tasks=len(rows),totals=totals,compiler_CPU_seconds=time.process_time()-started)))
if __name__=='__main__':main()

"""Reconcile already-terminal supplementary64 public artifacts; never rerun science/GT."""
import json,tarfile,time,subprocess,resource
from pathlib import Path
import twod_dev_transport as t
from twod_postprocess_cost import attach_external_cpu
D=t.ROOT/'delegated/twod/crossgrid8_paired';P=t.ROOT/'delegated/diagnostics/crossgrid8_protocol'
def main():
 t.OUT=D;archive=D/'PUBLIC_OUTPUTS.tgz'
 if not archive.exists():raise RuntimeError('Existing completed public download required; no repack or solver dispatch')
 c=time.process_time()
 with tarfile.open(archive) as f:
  if any('sealed_truth' in n or 'data_private' in n for n in f.getnames()):raise RuntimeError('Private transport forbidden')
  f.extractall(t.REPO,filter='data')
 extraction=time.process_time()-c;q=json.loads((D/'queue/RECEIPT.json').read_text())
 if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS') or q['running'] or q['pending'] or len(q['jobs'])!=64:raise RuntimeError('All64 terminal required')
 out=D/'public_seal';cost=D/'PUBLIC_SEAL_PROCESS_COST.json';runtime=json.loads((t.ROOT/'configs/twod_local_runtime.json').read_text())['python']
 if not out.exists():
  before=resource.getrusage(resource.RUSAGE_CHILDREN);c=time.process_time();wall=time.perf_counter();p=subprocess.run([runtime,str(P/'seal_paired_grid8.py'),'--queue',str(D/'queue/RECEIPT.json'),'--output',str(out)],capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
  cost.write_text(json.dumps(dict(returncode=p.returncode,whole_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-wall,stdout=p.stdout,stderr=p.stderr,truth_access='NONE'),indent=2)+'\n')
  if p.returncode:raise RuntimeError(p.stderr)
  attach_external_cpu(out,cost)
 elif json.loads((out/'RECEIPT.json').read_text())['status']!='ALL64_PUBLIC_PAIRED_GRID_CHOICES_SEALED_BEFORE_GT':raise RuntimeError('Preserve invalid seal attempt before retry')
 report=dict(status='ALL64_TERMINAL_PUBLIC_SEALED_GT_HELD',scientific_terminal=64,scientific_failed=sum(j['status']!='FINISHED' for j in q['jobs']),independent_materials=8,child_GetProcessTimes_CPU_seconds=sum(j['child_process_CPU_seconds'] for j in q['jobs']),queue_controller_CPU_seconds=q['controller_CPU_seconds'],local_extraction_CPU_seconds=extraction,remote_packaging_CPU_seconds=None,unknown_costs=['Initial controller remote packaging CPU did not reach final receipt if Python3.9 extraction failed; transport ledger remains billed'],source_transport_ledger=str(D/'TRANSPORT_COSTS.jsonl'),public_seal_external_cost=str(cost),GT_scoring='NOT_RUN',truth_access='NONE',solver_reruns=0)
 transport=[json.loads(line) for line in (D/'TRANSPORT_COSTS.jsonl').read_text().splitlines()];prepared=json.loads((D/'PREPARED.json').read_text());seal=json.loads(cost.read_text())
 totals=dict(inverse_whole_child_GetProcessTimes_CPU_seconds=report['child_GetProcessTimes_CPU_seconds'],queue_controller_CPU_lower_bound_seconds=report['queue_controller_CPU_seconds'],local_SSH_scp_and_controller_CPU_seconds=sum(r['child_CPU_seconds']+r['controller_CPU_seconds'] for r in transport),initial_public_extraction_remote_CPU_seconds=prepared['remote_extraction_CPU_seconds'],local_terminal_extraction_CPU_seconds=extraction,public_seal_external_child_CPU_seconds=seal['whole_child_CPU_seconds'],public_seal_controller_CPU_seconds=seal['controller_CPU_seconds'])
 report.update(totals=totals,known_campaign_CPU_lower_bound_seconds=sum(totals.values()),exact_campaign_CPU_seconds=None,excluded_already_root_accounted_producers=[str(P/'GENERATION_V2_EXTERNAL_COST.json'),str(P/'PUBLIC_PACKAGE_RECEIPT.json')],cost_deduplication='This ledger includes only this inverse queue, transport and public seal/postprocessing. Generation/public-package are excluded because root ledger charges those producers. Shared Windows post-run environment metadata charged once in primary48 ledger, not here.',unknown_costs=report['unknown_costs']+['Local controller import/final tails and remote metadata probe CPU unmeasured; not substituted0'])
 (D/'COSTS.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()

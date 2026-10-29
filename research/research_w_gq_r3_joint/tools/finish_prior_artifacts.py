"""Pure I/O completion from cached prior endpoint; never calls a physical solver."""
from pathlib import Path
import json,time,csv
ROOT=Path(__file__).resolve().parents[1]
def main():
    records=[json.loads(l) for l in (ROOT/'PHYSICS_LEDGER.jsonl').read_text().splitlines() if l]
    for case in ['blind_normal','blind_cavity','blind_strong']:
        out=ROOT/'results/online'/f'{case}_r0'/'PRIOR_ONLY'
        if (out/'ONLINE_SEALED.json').exists():raise RuntimeError('Prior output already sealed')
        if not (out/'endpoint.npz').exists():raise FileNotFoundError('Cached forward required')
        construction=json.loads((out/'reconstruction.json').read_text());summary=construction['summary'];source=ROOT/'results/jobs'/f'control_{case}_r0_prior'/'RECEIPT.json';receipt=json.loads(source.read_text());pid=receipt['pid'];rr=[r for r in records if r.get('pid')==pid and r.get('kind')=='reservation'];cost=dict(full_fe=sum(r.get('factorizations',0) for r in rr),rhs=sum(r.get('rhs',0) for r in rr))
        trace=dict(time_unix=(out/'endpoint.npz').stat().st_mtime,candidate=1,file='endpoint.npz',objective=summary['objective'],stage='prior_parameter_mean',kind='prior_cached_forward')
        (out/'TRACE.jsonl').write_text(json.dumps(trace)+'\n')
        with (out/'optimization_trace.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(trace));w.writeheader();w.writerow(trace)
        (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
        (out/'ONLINE_SEALED.json').write_text(json.dumps(dict(time_unix=time.time(),selected_by='independent parameter prior means; no measurement optimization',truth_access=False,parameter_names=construction['parameter_names'],summary=summary,artifact_recovery='cached forward; no PDE rerun'),indent=2)+'\n')
        runtime=dict(pid=pid,method='PRIOR_ONLY',process_cpu_seconds=receipt['process_cpu_seconds'],wall_seconds=receipt['wall_seconds'],peak_rss_bytes=None,polled_peak_RSS_bytes=receipt['polled_peak_RSS_bytes'],job_physical_cost=cost,missing_metrics_status='CPU_RECORDED_RSS_POLLED_ONLY',original_exit_code=1,artifact_recovery=True,CPU_scope='original forward and failed output job; shared artifact-only recovery recorded separately as finish_prior_artifacts_once in JOB_LEDGER, not allocated among three parents',source_receipt=str(source.relative_to(ROOT)))
        (out/'runtime_cost.json').write_text(json.dumps(runtime,indent=2)+'\n');(out/'ARTIFACT_IO_REPAIR.json').write_text(json.dumps(dict(reason='frozen save_endpoint expected TRACE for optimized arms; prior representative has no optimization',forward_reused=True,new_physics=False,original_failure_log=f'results/jobs/control_{case}_r0_prior/OUTPUT.log',shared_repair_job='finish_prior_artifacts_once',original_cost=cost),indent=2)+'\n')
    print('Three cached prior endpoints sealed without new PDE; original failures/cost retained')
if __name__=='__main__':main()

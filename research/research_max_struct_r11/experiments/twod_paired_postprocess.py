"""Explicit paired8-only authorized scoring/plots with disjoint external CPU ledger."""
import json,subprocess,resource,time,argparse
from pathlib import Path
from twod_postprocess_cost import attach_external_cpu
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'delegated/twod/crossgrid8_paired';P=ROOT/'delegated/diagnostics/crossgrid8_protocol'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['score','plots','common-b'],required=True);a=ap.parse_args();auth=ROOT/'audit/PAIRED_GRID8_ROOT_POSTSEAL_AUTHORIZATION.json'
 if json.loads(auth.read_text())['status']!='ROOT_EXPLICIT_POST64SEAL_SCORE':raise PermissionError('Paired-specific root authorization required')
 runtime=json.loads((ROOT/'configs/twod_local_runtime.json').read_text())['python'];out=D/('scored' if a.stage=='score' else 'common_b_reindex' if a.stage=='common-b' else 'plots')
 if a.stage=='score':cmd=[runtime,str(P/'score_paired_grid8.py'),'--public-seal',str(D/'public_seal'),'--output',str(out),'--authorize','ROOT_EXPLICIT_POST64SEAL_SCORE']
 elif a.stage=='common-b':cmd=[runtime,str(ROOT/'experiments/twod_paired_common_b.py'),'--output',str(out)]
 else:cmd=[runtime,str(ROOT/'experiments/twod_paired_plots.py'),'--scores',str(D/'scored'),'--output',str(out),'--authorization',str(auth)]
 b=resource.getrusage(resource.RUSAGE_CHILDREN);c=time.process_time();w=time.perf_counter();process=subprocess.run(cmd,capture_output=True,text=True);e=resource.getrusage(resource.RUSAGE_CHILDREN);ledger=D/'POSTPROCESS_EXECUTIONS';ledger.mkdir(exist_ok=True);cost=ledger/(a.stage+'_'+str(time.time_ns())+'.json');x=dict(stage=a.stage,command=cmd,returncode=process.returncode,whole_child_CPU_seconds=e.ru_utime+e.ru_stime-b.ru_utime-b.ru_stime,controller_CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,stdout=process.stdout,stderr=process.stderr,authorization=str(auth),scope='Paired8 only; external whole child CPU replaces nested self timestamps',new_FW=0);cost.write_text(json.dumps(x,indent=2)+'\n')
 if process.returncode:print(process.stderr);raise SystemExit(process.returncode)
 attach_external_cpu(out,cost);campaign=json.loads((D/'COSTS.json').read_text());rows=[json.loads(p.read_text()) for p in sorted(ledger.glob('*.json'))];campaign['totals']['authorized_paired_postprocess_external_child_CPU_seconds']=sum(r['whole_child_CPU_seconds'] for r in rows);campaign['totals']['authorized_paired_postprocess_controller_CPU_seconds']=sum(r['controller_CPU_seconds'] for r in rows);campaign['known_campaign_CPU_lower_bound_seconds']=sum(campaign['totals'].values());campaign['GT_scoring']=json.loads((D/'scored/RECEIPT.json').read_text())['status'];campaign['status']='ALL64_SCORED_POST_PUBLIC_SEAL' if a.stage=='score' else 'ALL64_SCORED_AND_PLOTTED_POST_PUBLIC_SEAL';campaign['authorization']=str(auth);(D/'COSTS.json').write_text(json.dumps(campaign,indent=2)+'\n');print(process.stdout);print(json.dumps({k:x[k] for k in ['stage','returncode','whole_child_CPU_seconds','controller_CPU_seconds']}))
if __name__=='__main__':main()

"""Root-invoked external CPU wrapper; forwards authorization unchanged, never auto-scores."""
import argparse,json,subprocess,resource,time
from pathlib import Path
from twod_postprocess_cost import attach_external_cpu
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'delegated/twod/independent48'
SCRIPTS={'primary':'twod_independent_score.py','secondary':'twod_secondary_score.py','plots':'twod_independent_plots.py','fairness-seal':'twod_actual_cpu_fairness_seal.py','endpoint-bookkeeping':'../delegated/diagnostics/endpoint_bookkeeping/endpoint_bookkeeping.py','plot-margin-audit':'twod_plot_margin_audit.py'}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--kind',choices=SCRIPTS,required=True);ap.add_argument('child_args',nargs=argparse.REMAINDER);a=ap.parse_args();args=a.child_args[1:] if a.child_args[:1]==['--'] else a.child_args;runtime=json.loads((ROOT/'configs/twod_local_runtime.json').read_text())['python'];out=Path(args[args.index('--output')+1]) if '--output' in args else (D/'scored' if a.kind=='primary' else D/'actual_cpu_fairness/public_seal' if a.kind=='fairness-seal' else D/'secondary/scored')
 c=time.process_time();wall=time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN);p=subprocess.run([runtime,str((ROOT/'experiments'/SCRIPTS[a.kind]).resolve()),*args],capture_output=True,text=True,errors='replace');after=resource.getrusage(resource.RUSAGE_CHILDREN);duration=time.perf_counter()-wall
 ledger=D/'POSTPROCESS_EXECUTIONS';ledger.mkdir(exist_ok=True);cost=ledger/(a.kind+'_'+str(time.time_ns())+'.json');child_receipt=json.loads((out/'RECEIPT.json').read_text()) if p.returncode==0 and (out/'RECEIPT.json').exists() else {}
 r=dict(kind=a.kind,returncode=p.returncode,child_status=child_receipt.get('status','PROCESS_FAILED_SEE_TRANSCRIPT'),output=str(out),whole_child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-c,wall_seconds=duration,scope='Disjoint external whole child; replaces rather than adds child selfCPU',stdout=p.stdout,stderr=p.stderr,GT_scoring_status=child_receipt.get('status','FAILED_PROCESS_GT_ACCESS_NOT_ASSUMED') if a.kind=='primary' else 'NO_ADDITIONAL_GT_SCORING_IN_WRAPPER')
 cost.write_text(json.dumps(r,indent=2)+'\n')
 if p.returncode==0:attach_external_cpu(out,cost)
 print(p.stdout,end='');print(p.stderr,end='');print(json.dumps({k:r[k] for k in ('kind','returncode','whole_child_CPU_seconds','controller_CPU_seconds','child_status')}))
 return p.returncode
if __name__=='__main__':raise SystemExit(main())

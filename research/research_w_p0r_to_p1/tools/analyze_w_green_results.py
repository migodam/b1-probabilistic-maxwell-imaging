"""Read completed engineering JSON only; never physics, GT, Book or arrays."""
import csv,json,time,statistics as st
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
 cpu,wall=time.process_time(),time.perf_counter()
 base=ROOT/'results/P1B/green_reuse_control'
 s=json.loads((base/'SUMMARY.json').read_text());r=json.loads((base/'RESULTS.json').read_text())
 if s['status']!='ENGINEERING_CONTROL_COMPLETE' or len(r)!=112:raise RuntimeError('All registered16x7 results required')
 rows=[]
 for m in sorted({x['method'] for x in r}):
  a=[x for x in r if x['method']==m]
  rows.append(dict(method=m,n=len(a),response_error_white_max=max(x['response_error_white'] for x in a),
   mean_warm_CPU=st.mean(x['warm_deployment_cost']['cpu_seconds'] for x in a),
   mean_warm_wall=st.mean(x['warm_deployment_cost']['wall_seconds'] for x in a),
   paid_warm_RHS=sum(x['warm_deployment_cost']['new_rhs'] for x in a),
   successful_responses=sum(x['response_error_white']<=.05 and x['status']!='GMRES_NOT_CONVERGED' for x in a)))
 a17=[x for x in r if x['method'].startswith('A17')]
 ev=[json.loads(l) for l in (base/'EVENTS.jsonl').read_text().splitlines() if l.strip()]
 kernel=[x for x in ev if x.get('event')=='engineering_operation' and x.get('label')=='KRYLOV_base_LU_apply_1RHS']
 out=ROOT/'tables/green_control';out.mkdir(exist_ok=False)
 with (out/'METHODS.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 result=dict(status='DESCRIPTIVE_COMPLETED_ENGINEERING_READBACK',methods=rows,
   same_current_space_relative_difference_max=max(x['same_current_space_relative_difference'] for x in a17),
   same_current_space_response_difference_white_max=max(x['same_current_space_response_difference_white'] for x in a17),
   exact_selector_active_rows=s['setup']['exact_selector_support']['active_row_count'],
   paid_exact_setup_RHS=s['setup']['exact_R0U']['cost']['new_rhs'],
   pooled_stage_scope_cost={k:v for k,v in s['cost'].items() if k!='reserve_ids'},
   root_process_receipt='results/jobs/p1b-full-Woodbury-ROM-GMRES-control/RECEIPT.json',
   kernels_only=dict(Krylov_base_actions=len(kernel),CPU=sum(x['operation_cpu_seconds'] for x in kernel),wall=sum(x['operation_wall_seconds'] for x in kernel)),
   setup_components_overlap=True,standalone_operator_preparation_paid_in_Full_pool=True,
   Full_fallback_policy_evaluated=False,posterior_deltaPhi_or_certification_claim=False,
   source_files=['results/P1B/green_reuse_control/'+x for x in ['RESULTS.json','SUMMARY.json','SETUP.json','EVENTS.jsonl']],
   new_RHS=0,new_LU=0,CPU=time.process_time()-cpu,wall=time.perf_counter()-wall)
 (out/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k in ['status','new_RHS','new_LU','CPU','wall']}))
if __name__=='__main__':main()

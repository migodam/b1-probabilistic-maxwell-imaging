"""Postseal oracle boundary over already scored paid Phi-incumbents. Not an algorithm."""
from pathlib import Path
import json,csv,time
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1]
def main():
 wall,cpu=time.perf_counter(),time.process_time()
 d=ROOT/'delegated/diagnostics/bothblocks_confirmation'
 final={(x['case_id'],x['method']):x for x in json.loads((d/'SINGLE_ARM_SCORES.json').read_text())}
 scores={}
 with (d/'PAID_INCUMBENT_FRONTIER.csv').open() as f:
  for x in csv.DictReader(f):
   k=(x['case_id'],x['method']);v=float(x['amplitude_NRMSE'])
   if k not in scores:scores[k]=dict(points=0,best_NRMSE=v,best_index=int(x['source_index']),best_CPU=float(x['CPU_seconds']),any_material_success=False,first_success_index=None)
   s=scores[k];s['points']+=1
   if v<s['best_NRMSE']:s.update(best_NRMSE=v,best_index=int(x['source_index']),best_CPU=float(x['CPU_seconds']))
   if x['material_success']=='True':
    s['any_material_success']=True
    if s['first_success_index'] is None:s['first_success_index']=int(x['source_index'])
 rows=[];groups=defaultdict(list)
 for (cid,method),s in sorted(scores.items()):
  f=final[cid,method]
  row=dict(case_id=cid,method=method,scope='POSTHOC_ORACLE_BOUNDARY_ONLY; no online GT selection',**s,
    final_NRMSE=f['amplitude_NRMSE'],final_material_success=f['material_success'],
    final_minus_oracle_NRMSE=f['amplitude_NRMSE']-s['best_NRMSE'],
    recovery_visited_then_lost=bool(s['any_material_success'] and not f['material_success']))
  rows.append(row);groups[method].append(row)
 summary=[dict(method=m,objects=len(a),final_successes=sum(x['final_material_success'] for x in a),
    any_paid_incumbent_successes=sum(x['any_material_success'] for x in a),
    visited_success_then_lost=sum(x['recovery_visited_then_lost'] for x in a),
    mean_final_minus_oracle_NRMSE=sum(x['final_minus_oracle_NRMSE'] for x in a)/len(a)) for m,a in groups.items()]
 out=d/'oracle_trace_boundary';out.mkdir(exist_ok=False)
 (out/'ROWS.json').write_text(json.dumps(rows,indent=2)+'\n')
 (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
 (out/'RECEIPT.json').write_text(json.dumps(dict(status='COMPLETE',scope='Postseal scored CSV only; no new solver, no new GT array access, no online selection',physical_forward=0,physical_adjoint=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,warning='Best NRMSE and existence of a passing state are separate GT oracles, not a deployable common selector; only saved physical-Phi incumbents considered'),indent=2)+'\n')
 print(json.dumps(summary,indent=2))
if __name__=='__main__':main()

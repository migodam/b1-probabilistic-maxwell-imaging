"""Run frozen direct inverse using public geometry and observations only."""
import sys,json,time
from pathlib import Path
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model
from algebraic_solver import reconstruct_identifiable
from validate import CountingBook,load_case,validate_endpoint,write

def main():
 manifest=json.loads((P/'experiments/public/manifest.json').read_text())
 for row in manifest['cases']:
  case=row['case'];out=P/'results'/f'blind_{case}_receiver_svd.json'
  if out.exists():raise RuntimeError('Frozen result already exists; preserve original')
  _,g,y=load_case(case);book=CountingBook();model=Model(book,**g);t=time.perf_counter();cpu=time.process_time()
  answer=reconstruct_identifiable(model,y,rank_rtol=1e-12,node_threshold=1e-20,constitutive_rtol=1e-8,box=(.5,4.))
  r={'case':case,'method':'receiver_svd','status':answer['status'],'algebraic_wall_seconds':time.perf_counter()-t,'algebraic_CPU_seconds':time.process_time()-cpu,'receiver_rank':answer['receiver_rank'],'algebraic_diagnostics':{k:v for k,v in answer.items() if k not in ('candidate','reconstructed_current','node_field_power')}}
  if answer['candidate'] is not None:
   r.update(m=answer['candidate'],validation=validate_endpoint(model,answer['candidate'],y))
  r.update(accounting=book.summary(model.ns,model.nr),model_counts=dict(model.counts),model_timings=dict(model.timings));write(out,r)
  print(case,r['status'],r.get('validation',{}).get('response_relative'))
if __name__=='__main__':main()

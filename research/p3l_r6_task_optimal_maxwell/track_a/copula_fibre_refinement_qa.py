"""Independent GL64/128/256 prior-fibre QA; AST only, no physics imports."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import ast,json,csv,time
from pathlib import Path
import numpy as np
from scipy.special import ndtr,ndtri,betaincinv
from numpy.polynomial.legendre import leggauss
TRACK=Path(__file__).resolve().parent;SRC=TRACK/'copula_counterexample.py'
allowed={'LO','HI','Z0','ZEDGE','PH','ULO','UHI','RLO','RHI','ILO','IHI','AREA'};body=[]
for node in ast.parse(SRC.read_text()).body:
 if isinstance(node,ast.Assign) and {x.id for target in node.targets for x in ast.walk(target) if isinstance(x,ast.Name)}<=allowed:body.append(node)
 elif isinstance(node,ast.FunctionDef) and node.name in ['H','conditional_fibre']:body.append(node)
ns=dict(np=np,ndtr=ndtr,ndtri=ndtri,betaincinv=betaincinv,leggauss=leggauss)
exec(compile(ast.Module(body=body,type_ignores=[]),str(SRC),'exec'),ns)
def main():
 tick=time.process_time();records=[];comparison=[]
 for case in ['REGULAR','PURE_OVERLAP']:
  with np.load(TRACK/f'raw/copula_counterexample/{case}.npz') as z:v=complex(*z['saved_material_mode'])
  rows=[]
  for n in [64,128,256]:
   f=ns['conditional_fibre'](v,n);row=dict(case=case,GL=n,V_real=v.real,V_imag=v.imag,density=f['density'],pure_probability=f['pure_probability'],E_Q=f['moments']['1'],E_Q2=f['moments']['2'],E_invQ=f['moments']['-1'],E_invQ2=f['moments']['-2'],CS_excess=f['CS_excess'],Cov_Q_deltaR=f['conditional_Cov_Q_deltaR'],product_mean_bias_real=f['exact_product_mean_bias'][0],product_mean_bias_imag=f['exact_product_mean_bias'][1],product_material_complex_variance=f['exact_product_complex_material_variance']);records.append(row);rows.append(row)
  ref=rows[-1]
  for row in rows[:-1]:comparison.append(dict(case=case,GL=row['GL'],reference_GL=256,max_abs_moment_difference=max(abs(row[k]-ref[k]) for k in ['E_Q','E_Q2','E_invQ','E_invQ2']),density_relative_difference=abs(row['density']/ref['density']-1),product_bias_real_abs_difference=abs(row['product_mean_bias_real']-ref['product_mean_bias_real']),variance_abs_difference=abs(row['product_material_complex_variance']-ref['product_material_complex_variance']),qualification='observed float64 quadrature refinement only; no interval or globalposterior certificate'))
 result=dict(status='INDEPENDENT_ZERO_RHS_FIBRE_REFINEMENT_QA',records=records,comparisons=comparison,scope='fixedmaterialmode conditional prior fibre only; ASTextraction doesnotimport sharedphysics/Book/Target',physics_RHS=0,training_jobs=0,CPU_threads=1,CPU_seconds=time.process_time()-tick,old_raw_or_protocol_modified=False,A_total_actual_RHS=158144)
 (TRACK/'COPULA_FIBRE_REFINEMENT_QA.json').write_text(json.dumps(result,indent=2)+'\n')
 with (TRACK/'COPULA_FIBRE_REFINEMENT_QA.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(records[0]));w.writeheader();w.writerows(records)
 print(json.dumps(dict(CPU_seconds=result['CPU_seconds'],comparisons=comparison),indent=2))
if __name__=='__main__':main()

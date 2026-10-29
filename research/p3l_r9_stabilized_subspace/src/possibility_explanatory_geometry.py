"""Post-seal, post-scoring explanation only; no forward, fitting or online changes."""
from pathlib import Path
import os
for k in ['OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS']:os.environ[k]='1'
import numpy as np,json,time,csv
P=Path(__file__).resolve().parents[1];V=P/'results/possibility_volume';O=V/'posthoc_geometry'
assert (V/'offline_statistics/SCORING_COMPLETE.json').exists()
O.mkdir(exist_ok=False)
(O/'DESIGN.json').write_text(json.dumps({'utc':time.time(),'role':'POSTHOC_OFFLINE_EXPLANATION_AFTER_FROZEN_PRIMARY_SCORE','not_primary':True,'no_New_Maxwell_RHS':0,'no_model_or_candidate_changes':True,'questions':['GT material error projection onto initial weak8 vs full material tangent','GT error projection onto span of actual qualified candidate increments','GT noise eligibility and supplied latent distance to anchor (not a gauge-invariant feasibility certificate)'],'projection_warning':'Affine-span projection may be negative weighted, decoder illegal, data inconsistent and GT-selected; it is only a representation diagnostic, never an imaging algorithm.'},indent=2))
def pack(z):return np.r_[z.real,z.imag]/24

def proj_error(e,A):
 if not A.size:return float(e@e),0
 u,s,_=np.linalg.svd(A,full_matrices=False);r=int(np.sum(s>max(float(s[0])*1e-8,1e-12)))
 ee=e-u[:,:r]@(u[:,:r].T@e)
 return float(ee@ee),r
rows=[]
for pd in sorted((V/'parents').glob('R9V_FRESH_*')):
 with np.load(P/'private/possibility_volume/truth'/f'{pd.name}.npz') as t, np.load(pd/'anchor/ANCHOR.npz') as a,np.load(pd/'ANCHOR_FRAME.npz') as f:
  e=pack(t['material']-a['material']);base=float(e@e);wf,wr=proj_error(e,f['M']@f['Qw']);tf,tr=proj_error(e,f['M'])
  common={'parent':pd.name,'anchor_MSE':base,'truth_Phi':float(np.sum(np.abs(t['noise'])**2)/float(t['sigma'])**2),'truth_latent_anchor_L2':float(np.linalg.norm(t['x']-a['x'])),'weak8_linear_residual_MSE':wf,'weak8_explained_fraction':1-wf/base,'material_tangent_linear_residual_MSE':tf,'tangent_explained_fraction':1-tf/base,'rank_weak':wr,'rank_tangent':tr}
  for arm in ['FIXED_Q8_RAY','DYNAMIC_Q8_RAY','RANDOM66_FULL']:
   materials=[]
   for rp in sorted((pd/arm).glob('*/RESULT.json')):
    r=json.loads(rp.read_text())
    if r.get('primary_qualified') and r.get('selected_path'):
     with np.load(r['selected_path']) as z:materials.append(pack(z['material']-a['material']))
   A=np.stack(materials,axis=1) if materials else np.empty((1152,0));ce,cr=proj_error(e,A)
   rows.append(dict(common,method=arm,candidate_increment_span_rank=cr,candidate_affine_oracle_residual_MSE=ce,candidate_affine_explained_fraction=1-ce/base))
with (O/'PARENT_DIAGNOSTICS.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={'role':'POSTHOC_EXPLANATORY_ONLY','parents':32,'new_Maxwell_RHS':0,'mean':{}}
for arm in sorted(set(r['method'] for r in rows)):
 rs=[r for r in rows if r['method']==arm]
 summary['mean'][arm]={k:float(np.mean([r[k] for r in rs])) for k in rs[0] if k not in ['parent','method']}
summary['truth_noise_qualified_parents']=sum(r['truth_Phi']<=282.8737921622988 for r in rows if r['method']=='FIXED_Q8_RAY')
summary['latent_distance_over4_parents']=sum(r['truth_latent_anchor_L2']>4 for r in rows if r['method']=='FIXED_Q8_RAY')
summary['warning']='latent preimage distance is not same-material quotient distance; affine GT projection is not legal or deployable.'
(O/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))

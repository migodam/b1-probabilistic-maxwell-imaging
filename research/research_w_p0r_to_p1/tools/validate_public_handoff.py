"""Read the exported native hypothesis package only, never GT/physics/Book."""
from pathlib import Path
import json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def main():
 cpu,wall=time.process_time(),time.perf_counter();p=ROOT/'MaterialHypothesisPackage'
 d=json.loads((p/'EXAMPLE.json').read_text());manifest=json.loads((p/'MANIFEST.json').read_text())
 if d['schema']!='W_PROBABILITY_HANDOFF_V1' or d['observation']!=0 or d['truth_access'] is not False:raise ValueError('Native example identity changed')
 for r in manifest['files']:
  f=p/r['public_path']
  if not f.resolve().is_relative_to(p.resolve()) or not f.is_file() or f.stat().st_size!=r['bytes']:raise ValueError('Missing/incomplete exported evidence')
 e=d['event_partition'];events=e['events']
 if len(events)!=2 or e['topK_renormalization'] is not False or e['represented_mass']!=1 or not e['exhaustive_within_controlled_two_wall_family']:raise ValueError('Event partition semantics')
 if abs(sum(x['posterior_mass'] for x in events)-1)>1e-12 or sum(x['prior_mass'] for x in events)!=1:raise ValueError('Mass conservation')
 if d['outside_family_mass'] is not None or any(x['point_probability'] is not None for x in events):raise ValueError('Point mass/outside-family semantics')
 downstream=d['downstream']
 with np.load(p/downstream['finite_joint_weights']['path'],allow_pickle=False) as a:
  w=a['posterior_weights'].reshape(650,-1);branch=a['material_branches'];aa=a['amplitude_axis'];pp=a['phase_axis']
 if not np.isfinite(w).all() or (w<0).any() or abs(w.sum()-1)>1e-12:raise ValueError('Joint weights')
 finite_mass=[float(w[branch==b].sum()) for b in [0,1]]
 discrepancy=max(abs(finite_mass[b]-events[b]['posterior_mass']) for b in [0,1])
 if discrepancy>.02:raise ValueError('Finite/continuous mass discrepancy above declared limit')
 with np.load(p/downstream['public_material_response']['path'],allow_pickle=False) as a:F=a['F_full'].copy()
 a,ph=np.meshgrid(aa,pp,indexing='ij');gamma=(a*np.exp(1j*ph)).ravel()
 expected=np.tensordot(w@gamma,F,axes=(0,0))
 with np.load(p/downstream['correlation_factors']['path'],allow_pickle=False) as a:
  saved=a['posterior_predictive_approx'];rep=a['Fraw'].copy()
 error=float(np.max(abs(expected-saved)))
 if error>1e-12:raise ValueError('Joint predictive correlation lost')
 traced=[]
 for row in downstream['weighted_samples']:
  with np.load(p/row['path'],allow_pickle=False) as a:
   idx=a[row['final_particle_indexes_key']];weights=a[row['final_weights_key']]
   history=a[row['all_repeated_states_key']]
   if idx.shape!=weights.shape or (idx<0).any() or (idx>=w.size).any() or not np.isfinite(weights).all() or (weights<0).any() or abs(weights.sum()-1)>1e-12:raise ValueError('Invalid final particles')
   if (history<0).any() or (history>=w.size).any():raise ValueError('Invalid repeated trajectory state')
   traced.append(dict(schedule=row['schedule'],seed=row['seed'],final_particles=len(idx),history_shape=list(history.shape)))
 result=dict(status='PUBLIC_PACKAGE_ARRAY_AND_SEMANTIC_READBACK_PASS',native_schema_preserved=True,files=len(manifest['files']),
  posterior_masses_continuous=[x['posterior_mass'] for x in events],finite_masses=finite_mass,
  max_finite_vs_continuous_event_error=discrepancy,posterior_predictive_max_raw_readback_difference=error,
  final_trajectories=traced,represented_mass=1,outside_family_mass=None,
  not_a_general_unknown_wall_posterior=True,no_new_probability_qualification_inferred=True,
  CPU=time.process_time()-cpu,wall=time.perf_counter()-wall,new_RHS=0,new_LU=0,GT_read=False)
 out=ROOT/'AUDIT/public_handoff';out.mkdir(exist_ok=False)
 (out/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k in ['status','files','max_finite_vs_continuous_event_error','posterior_predictive_max_raw_readback_difference','CPU','wall','new_RHS','GT_read']}))
if __name__=='__main__':main()

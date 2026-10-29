"""Post-freeze truth scorer and independent total-field endpoint validator."""
import argparse,csv,json,time,shutil,os
from pathlib import Path
import numpy as np
from scipy.linalg import lu_factor,lu_solve
from scipy.linalg.lapack import get_lapack_funcs
from scipy.ndimage import uniform_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
METHODS=['fullj_native','fullj_reciprocal','globaldual']
LABELS={'fullj_native':'Full-J native','fullj_reciprocal':'Full-J reciprocal','globaldual':'Global dual + Full-J'}
def dump(p,x):
 p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def metrics(m,t,mask):
 d=m-t;foreground=m>=5.;union=np.logical_or(foreground,mask).sum()
 # Standard local SSIM: 7x7 uniform window, unbiased sample covariance,
 # range 19.5, K1=.01/K2=.03; crop 3-pixel borders before averaging.
 x=m.reshape(24,24);y=t.reshape(24,24);ux=uniform_filter(x,7);uy=uniform_filter(y,7)
 vx=49/48*(uniform_filter(x*x,7)-ux*ux);vy=49/48*(uniform_filter(y*y,7)-uy*uy);vxy=49/48*(uniform_filter(x*y,7)-ux*uy)
 c1=(.01*19.5)**2;c2=(.03*19.5)**2;s=((2*ux*uy+c1)*(2*vxy+c2))/((ux*ux+uy*uy+c1)*(vx+vy+c2))
 return dict(RMSE=float(np.sqrt(np.mean(d*d))),MAE=float(np.mean(abs(d))),foreground_RMS=float(np.sqrt(np.mean(d[mask]**2))),background_RMS=float(np.sqrt(np.mean(d[~mask]**2))),IoU=float(np.logical_and(foreground,mask).sum()/max(union,1)),SSIM=float(s[3:-3,3:-3].mean()))
def independent(m,ops,Y):
 start=time.perf_counter();cpu=time.process_time();G,H,U0=ops['G'],ops['H'],ops['U0'];n=len(m)
 # Distinct total-field system, versus solver's current-system LU.
 A=np.eye(n,dtype=complex)-G*m[None,:];factor,pivots=lu_factor(A,check_finite=False);E=lu_solve((factor,pivots),U0,check_finite=False);j=m[:,None]*E;y=H@j
 normA1=float(np.linalg.norm(A,1));condition_start=time.perf_counter();condition_cpu=time.process_time()
 rcond,condition_info=get_lapack_funcs('gecon',(factor,))(factor,normA1,norm='1')
 condition_seconds=time.perf_counter()-condition_start;condition_cpu_seconds=time.process_time()-condition_cpu
 # Current equation is evaluated, never solved a second time.
 residual=j-m[:,None]*(U0+G@j)
 return dict(independent_validator_version='total-field-LU-zgecon-v2',independent_state_relative_residual=float(np.linalg.norm(A@E-U0)/np.linalg.norm(U0)),independent_normwise_backward_error_1=float(np.linalg.norm(A@E-U0,1)/(normA1*np.linalg.norm(E,1)+np.linalg.norm(U0,1))),independent_current_relative_residual=float(np.linalg.norm(residual)/max(np.linalg.norm(j),1e-300)),independent_data_relative_residual=float(np.linalg.norm(y-Y)/np.linalg.norm(Y)),independent_field_Frobenius_norm=float(np.linalg.norm(E)),independent_field_max_abs=float(abs(E).max()),independent_field_incident_norm_ratio=float(np.linalg.norm(E)/np.linalg.norm(U0)),independent_matrix_norm_1=normA1,independent_rcond_1_estimate=float(rcond),independent_condition_1_estimate=float(1/rcond) if rcond>0 else None,independent_zgecon_info=int(condition_info),independent_zgecon_calls=1,independent_zgecon_seconds=condition_seconds,independent_zgecon_cpu_seconds=condition_cpu_seconds,independent_zgecon_internal_RHS='OPAQUE_LAPACK_ESTIMATOR_WORK_NOT_EXPOSED_NOT_ZERO',independent_seconds=time.perf_counter()-start,independent_cpu_seconds=time.process_time()-cpu,independent_complex_RHS=U0.shape[1],independent_LU=1),y
def find_output(split,method,ident):
 candidates=[ROOT/'results'/split/method/f'{ident}.json',ROOT/'results'/split/f'{ident}_{method}.json',ROOT/'results'/split/f'{method}_{ident}.json']
 for p in candidates:
  if p.exists():return p
 raise FileNotFoundError(f'No output {split}/{method}/{ident}')
def freeze(split):
 input_split='dev' if split.startswith('dev') else split
 ids=sorted(p.stem for p in (ROOT/'public_data'/input_split).glob('*.npz'));expected=2 if split.startswith('dev') else 9
 if len(ids)!=expected:raise ValueError(f'Expected {expected} cases, got {len(ids)}')
 dest=ROOT/'experiments/frozen_predictions'/split
 if dest.exists():raise ValueError('Refusing overwrite prediction freeze')
 files=[]
 for method in METHODS:
  for ident in ids:
   p=find_output(split,method,ident);r=json.loads(p.read_text())
   if 'm' not in r and 'candidates' not in r:raise ValueError(f'Missing prediction {p}')
   q=dest/method/f'{ident}.json';q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
   files.append(dict(method=method,observation_id=ident,source=str(p.relative_to(ROOT)),snapshot=str(q.relative_to(ROOT)),bytes=q.stat().st_size))
 receipt=dict(split=split,cases=ids,methods=METHODS,files=files,predictions=len(files),frozen_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),truth_read_before_freeze=split=='dev_v2',scope='Known DEV tuning reuse after v1 truth release' if split=='dev_v2' else 'First prediction set; scorer reads only snapshots')
 dump(ROOT/f'audit/EXPERIMENT_PREDICTIONS_FREEZE_{split.upper()}.json',receipt);print(json.dumps(receipt))
def candidates(r):
 cs=r.get('candidates')
 if cs is None:return [dict(m=r['m'],objective=r.get('objective',r.get('score')),candidate_index=0)]
 return cs
def score(split):
 gate=ROOT/f'audit/EXPERIMENT_PREDICTIONS_FREEZE_{split.upper()}.json'
 if not gate.exists():raise ValueError('Predictions must freeze before truth scoring')
 input_split='dev' if split.startswith('dev') else split
 freeze_record=json.loads(gate.read_text());ops=dict(np.load(ROOT/'public_data/operators.npz'));mapping=json.loads((ROOT/'private_truth'/input_split/'mapping_PRIVATE.json').read_text());lookup={r['observation_id']:r for r in mapping['records']}
 rows=[];allrows=[];curves=[];runs={};images={};truth_diagnostics=[];oracle_validations=[];validation_cost=0;score_start=time.perf_counter();truthout=ROOT/'truth_after_freeze'/split
 truthout.mkdir(parents=True,exist_ok=True)
 for ident in freeze_record['cases']:
  truth=dict(np.load(ROOT/'private_truth'/input_split/f'{ident}.npz'));t=truth['m'];mask=truth['mask'].astype(bool);digit=lookup[ident]['digit'];Y=np.load(ROOT/'public_data'/input_split/f'{ident}.npz')['Y']
  np.savez_compressed(truthout/f'{ident}.npz',m=t,mask=mask)
  images[ident]={'GT':t,'digit':digit}
  gt_validation,gt_y=independent(t,ops,Y);validation_cost+=gt_validation['independent_seconds']
  born=ops['H']@(t[:,None]*ops['U0'])
  truth_diagnostics.append(dict(observation_id=ident,digit=digit,**gt_validation,generation_independent_forward_relative_disagreement=float(np.linalg.norm(gt_y-truth['Y_clean'])/np.linalg.norm(truth['Y_clean'])),first_Born_relative_disagreement=float(np.linalg.norm(born-truth['Y_clean'])/np.linalg.norm(truth['Y_clean'])),actual_noise_relative_norm=float(np.linalg.norm(Y-truth['Y_clean'])/np.linalg.norm(truth['Y_clean']))))
  for method in METHODS:
   r=json.loads((ROOT/'experiments/frozen_predictions'/split/method/f'{ident}.json').read_text());cs=candidates(r)
   runs[(ident,method)]=r
   selected=int(r.get('selected_index',r.get('top1_index',0)))
   # Prefer explicit solver-selected top1 material. Candidate scores are never GT selected.
   material=r.get('m') if r.get('m') is not None else (cs[selected]['m'] if cs else None)
   top=np.asarray(material,float).reshape(-1) if material is not None else np.array([])
   if top.size!=576 or not np.isfinite(top).all():
    images[ident][method]=np.full(576,np.nan)
    rows.append(dict(observation_id=ident,digit=digit,method=method,status=r.get('status','FAILED_NO_FINITE_ENDPOINT'),selection='physical_objective_top1',K=len(cs),RMSE=None,MAE=None,foreground_RMS=None,background_RMS=None,IoU=None,SSIM=None,independent_data_relative_residual=None,oracle_best_at_K_RMSE=None,within_bounds=False,finite_endpoint=False));continue
   met=metrics(top,t,mask);val,y=independent(top,ops,Y);validation_cost+=val['independent_seconds'];images[ident][method]=top
   row=dict(observation_id=ident,digit=digit,method=method,status=r.get('status','UNKNOWN'),selection='physical_objective_top1',K=len(cs),**met,**val)
   row['material_min']=float(top.min());row['material_max']=float(top.max());row['within_bounds']=bool(top.min()>=.5-1e-9 and top.max()<=20+1e-9)
   row['finite_endpoint']=True
   finalstages=[s for s in r.get('stages',[]) if s.get('parameter_dimension')==576];finalstage=finalstages[-1] if finalstages else {};fc=finalstage.get('counts',{})
   row.update(final576_evaluated=bool(finalstage.get('has_evaluated_candidate')),final576_fullj_builds=fc.get('FullJ',0),final576_smallgram_gradients=fc.get('smallgram_gradients',0),final576_Gram_factorizations=fc.get('Gram_factorizations',0),final576_smallmatrix_Gram_RHS=fc.get('Gram_RHS',0),final576_dual_accepted_iterations=fc.get('dual_optimizer_accepted_iterations',None),final576_gn_iterations=sum(e.get('kind')=='GN_ITERATION' and e.get('stage')==24 for e in r.get('events',[])),final576_accepted_gn_steps=sum(e.get('kind')=='GN_ITERATION' and e.get('stage')==24 and e.get('accepted',False) for e in r.get('events',[])))
   row['solver_reported_data_relative_residual']=r.get('data_relative_residual',r.get('relative_residual',(r.get('top1') or {}).get('relative_residual')))
   if row['solver_reported_data_relative_residual'] is not None:row['independent_residual_disagreement']=abs(float(row['solver_reported_data_relative_residual'])-val['independent_data_relative_residual'])
   candrows=[]
   for i,c in enumerate(cs):
    cm=np.asarray(c['m'],float).reshape(-1);cr=dict(observation_id=ident,digit=digit,method=method,candidate_index=i,reported_objective=c.get('objective',c.get('score')),**metrics(cm,t,mask));candrows.append(cr);allrows.append(cr)
   best=min(candrows,key=lambda z:z['RMSE']);row['oracle_best_at_K_RMSE']=best['RMSE'];row['oracle_best_at_K_index']=best['candidate_index']
   archive=np.array([c['m'] for c in cs],float);distances=np.sqrt(np.mean((archive[:,None,:]-archive[None,:,:])**2,axis=2));pairs=distances[np.triu_indices(len(cs),1)]
   representatives=[]
   for i in range(len(cs)):
    if not representatives or min(distances[i,j] for j in representatives)>=1.:representatives.append(i)
   row.update(archive_pairwise_material_RMS_min=float(pairs.min()) if pairs.size else 0.,archive_pairwise_material_RMS_median=float(np.median(pairs)) if pairs.size else 0.,archive_pairwise_material_RMS_max=float(pairs.max()) if pairs.size else 0.,archive_distinct_greedy_RMS1=len(representatives),archive_interpretation='optimizer iterates; not certified branches or independent samples')
   objective_values=[c.get('objective',c.get('score')) for c in cs];row['archive_sorted_by_physical_objective']=all(objective_values[i]<=objective_values[i+1]+1e-15 for i in range(len(cs)-1))
   for k in range(1,len(candrows)+1):
    kb=min(candrows[:k],key=lambda z:z['RMSE']);curves.append(dict(observation_id=ident,digit=digit,method=method,K=k,oracle_best_at_K_RMSE=kb['RMSE'],oracle_best_at_K_index=kb['candidate_index'],oracle_best_at_K_IoU=kb['IoU'],candidate_order='ascending frozen physical objective'))
   if not np.array_equal(archive[best['candidate_index']],top):
    ov,oy=independent(archive[best['candidate_index']],ops,Y);validation_cost+=ov['independent_seconds'];oracle_validations.append(dict(observation_id=ident,digit=digit,method=method,oracle_best_at_K_index=best['candidate_index'],**ov));row['oracle_endpoint_independently_validated']=True
   else:row['oracle_endpoint_independently_validated']='same as top1'
   rows.append(row)
 publicmap=[dict(observation_id=r['observation_id'],digit=r['digit'],parameters=r['parameters']) for r in mapping['records']]
 dump(truthout/'mapping_seedless.json',dict(records=publicmap,private_entropy_omitted=True,predictions_freeze=str(gate.relative_to(ROOT))))
 summary={}
 for method in METHODS:
  mr=[r for r in rows if r['method']==method];valid=[r for r in mr if r.get('finite_endpoint')]
  summary[method]={key:float(np.mean([r[key] for r in valid])) if valid else None for key in ['RMSE','MAE','foreground_RMS','background_RMS','IoU','SSIM','independent_data_relative_residual','oracle_best_at_K_RMSE']}
  summary[method].update(finite_endpoint_count=len(valid),total_count=len(mr),aggregate_scope='finite endpoints only; failures retained separately')
 coarse_audit=[]
 for ident in freeze_record['cases']:
  for stage in [8,12,24]:
   inits={method:next((s.get('initial_material') for s in runs[(ident,method)].get('stages',[]) if s.get('stage')==stage),None) for method in METHODS}
   available={method:np.array(v,float) for method,v in inits.items() if v is not None};pair_rms=[float(np.sqrt(np.mean((a-b)**2))) for i,a in enumerate(available.values()) for b in list(available.values())[i+1:]]
   coarse_audit.append(dict(observation_id=ident,stage=stage,available_methods=list(available),all3_initial_material_exact_match=len(available)==3 and all(np.array_equal(list(available.values())[0],v) for v in available.values()),maximum_pairwise_initial_material_RMS=max(pair_rms) if pair_rms else None))
 dump(ROOT/f'results/{split}_coverage.json',dict(split=split,coarse_initialization_comparison=coarse_audit,rows=[{key:r.get(key) for key in ['observation_id','method','status','finite_endpoint','final576_evaluated','final576_fullj_builds','final576_smallgram_gradients','final576_Gram_factorizations','final576_smallmatrix_Gram_RHS','final576_dual_accepted_iterations','final576_gn_iterations','final576_accepted_gn_steps','archive_distinct_greedy_RMS1','archive_pairwise_material_RMS_min','archive_pairwise_material_RMS_median','archive_pairwise_material_RMS_max']} for r in rows]))
 dump(ROOT/f'results/{split}_scores.json',dict(split=split,rows=rows,aggregate_means=summary,truth_physics_diagnostics=truth_diagnostics,metric_definition=dict(material_range=[.5,20],foreground_threshold=5,SSIM='standard local 7x7 uniform window, sample covariance, data_range19.5, K1=.01 K2=.03, cropped valid interior',archive='retained optimizer iterates sorted by physical objective; not certified branches or independent samples',distinct='greedy ascending-physical-rank representatives at material RMS>=1; descriptive, no acceptance gate',oracle='GT-selected best RMSE among first K retained physical-ranked optimizer iterates; zero-based candidate index'),all_candidate_metrics=allrows,best_at_K_curves=curves))
 with (ROOT/f'results/{split}_scores.csv').open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in rows))));writer.writeheader();writer.writerows(rows)
 validrows=[r for r in rows if r.get('finite_endpoint')]
 validation_records=validrows+truth_diagnostics+oracle_validations
 dump(ROOT/f'audit/EXPERIMENT_VALIDATION_{split.upper()}.json',dict(split=split,validator_version='total-field-LU-zgecon-v2',endpoint_count=len(validrows),failed_endpoint_count=len(rows)-len(validrows),truth_count=len(truth_diagnostics),additional_oracle_endpoint_count=len(oracle_validations),independent_total_field_LU=len(validation_records),independent_explicit_source_RHS=len(validation_records)*12,independent_complex_RHS=len(validation_records)*12,condition_estimator_calls=len(validation_records),condition_estimator_seconds=sum(v['independent_zgecon_seconds'] for v in validation_records),condition_estimator_internal_RHS='OPAQUE_LAPACK_ZGECON_WORK_NOT_EXPOSED; additional work beyond explicit source RHS, not zero',minimum_estimated_rcond_1=min(v['independent_rcond_1_estimate'] for v in validation_records),maximum_estimated_condition_1=max(v['independent_condition_1_estimate'] for v in validation_records if v['independent_condition_1_estimate'] is not None),validation_seconds=validation_cost,total_scoring_seconds=time.perf_counter()-score_start,all_bounds_pass=all(r['within_bounds'] for r in rows),maximum_state_relative_residual=max(v['independent_state_relative_residual'] for v in validation_records),truth_diagnostics=truth_diagnostics,oracle_endpoint_validations=oracle_validations,scope='same discrete scalar TM VIE; independent total-field system; no cross-model certification'))
 order=sorted(images,key=lambda ident:images[ident]['digit']);fig,axes=plt.subplots(len(order),4,figsize=(9,2.05*len(order)),squeeze=False,layout='constrained')
 im=None
 for i,ident in enumerate(order):
  for col,key in enumerate(['GT']+METHODS):
   ax=axes[i,col];im=ax.imshow(images[ident][key].reshape(24,24),vmin=.5,vmax=20,cmap='viridis',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
   if i==0:ax.set_title('Ground truth' if key=='GT' else LABELS[key],fontsize=11)
   if col==0:ax.set_ylabel(f'Digit {images[ident]["digit"]}',fontsize=10)
   else:
    row=next(r for r in rows if r['observation_id']==ident and r['method']==key);ax.set_xlabel(f'RMSE {row["RMSE"]:.2f} | IoU {row["IoU"]:.2f}' if row.get('finite_endpoint') else 'FAILED: no finite endpoint',fontsize=8)
 fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.5,label='Real TM contrast m; shared [0.5,20]')
 fig.suptitle(f'R3-DIGITS {split.upper()}: frozen physical-objective top1',fontsize=13)
 (ROOT/'figures').mkdir(exist_ok=True);fig.savefig(ROOT/f'figures/{split}_all_digits_main.png',dpi=170);fig.savefig(ROOT/f'figures/{split}_all_digits_main.pdf');plt.close(fig)
 fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
 for method in METHODS:
  ks=sorted(set(c['K'] for c in curves if c['method']==method));means=[];low=[];high=[]
  for k in ks:
   vals=[c['oracle_best_at_K_RMSE'] for c in curves if c['method']==method and c['K']==k];means.append(np.mean(vals));low.append(min(vals));high.append(max(vals))
  line=ax.plot(ks,means,label=LABELS[method],linewidth=2)[0];ax.fill_between(ks,low,high,color=line.get_color(),alpha=.12)
 ax.set_xlabel('K retained optimizer iterates in ascending physical-objective order');ax.set_ylabel('Oracle best@K material RMSE');ax.set_title(f'{split.upper()}: oracle archive curve; mean with case min–max');ax.legend();ax.grid(alpha=.2)
 fig.savefig(ROOT/f'figures/{split}_oracle_best_at_K.png',dpi=170);plt.close(fig)
 print(json.dumps(dict(split=split,aggregates=summary,validation_seconds=validation_cost)))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','score']);ap.add_argument('split',choices=['dev','dev_v2','holdout']);a=ap.parse_args();freeze(a.split) if a.action=='freeze' else score(a.split)
if __name__=='__main__':main()

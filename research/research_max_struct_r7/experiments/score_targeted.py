"""Independent targeted scorer; truth only after explicit TARGETED_SEALED authorization."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,json,time
import numpy as np
from score_primary import load,write,csvwrite,resolve,stamp,unchanged,material_metrics,noise_metrics,fixed_prior,public_observation
from analyse_primary import aggregate,bootstrap_pair,pairwise_frontiers
ROOT=Path(__file__).resolve().parents[1]
ARMS=('a_plus','cie_dual','global_dual','cie_single')
COMPARISONS=[('a_plus','cie_dual','TARGETED_SECONDARY_APLUS_VS_CIE_DUAL'),('a_plus','global_dual','TARGETED_SECONDARY_APLUS_VS_GLOBAL_DUAL'),('a_plus','cie_single','TARGETED_AMENDED_APLUS_VS_CIE_SINGLE')]

def plot(rows,truths,pictures,folder,frontiers):
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 folder.mkdir(parents=True,exist_ok=True)
 for group in range(4):
  ids=sorted(truths)[group*5:(group+1)*5];fig,axes=plt.subplots(5,5,figsize=(15,14),layout='constrained')
  for i,ident in enumerate(ids):
   entries=[r for r in rows if r['case_id']==ident];vmax=max([float(np.max(truths[ident]))]+[float(np.max(pictures[(ident,a)])) for a in ARMS if (ident,a) in pictures])
   axes[i,0].imshow(truths[ident].reshape(24,24),vmin=.02,vmax=vmax);axes[i,0].set_title(ident+' GT\n'+entries[0]['family'],fontsize=9)
   for j,arm in enumerate(ARMS,1):
    r=next(r for r in entries if r['arm']==arm);m=pictures.get((ident,arm))
    if m is not None:axes[i,j].imshow(m.reshape(24,24),vmin=.02,vmax=vmax);title=f"{arm}\nN={r['material_nrmse']:.3f} IoU={r['otsu_iou']:.2f}\nr={r['relative_residual']:.3g}"
    else:title=arm+' '+r['status']
    axes[i,j].set_title(title,fontsize=9)
   for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
  fig.suptitle('Secondary selected-family confirmation: all20 x4, group%d/4; shared row color range [.02,max(GT,all predictions)]'%(group+1));fig.savefig(folder/f'all20_group_{group+1:02d}.png',dpi=150);plt.close(fig)
 fig,axes=plt.subplots(3,3,figsize=(14,12),layout='constrained')
 for i,comp in enumerate(COMPARISONS):
  for j,axis in enumerate(('CPU','LU','complex_RHS')):
   ax=axes[i,j];g=[r for r in frontiers if r['comparison']==comp[2] and r['cost_axis']==axis and r.get('paired_objects',0)]
   x=[r['window_fraction'] for r in g];y=[r['mean_NRMSE_difference'] for r in g];lo=[r['mean_NRMSE_difference_CI95'][0] for r in g];hi=[r['mean_NRMSE_difference_CI95'][1] for r in g]
   ax.plot(x,y,'.-');ax.fill_between(x,lo,hi,alpha=.2);ax.axhline(0,color='k',lw=.7);ax.set_title(comp[0]+' minus '+comp[1]+' | '+axis);ax.set_xlabel('Fraction of observed paired cost window');ax.set_ylabel('Mean NRMSE difference');ax.grid(alpha=.2)
 fig.suptitle('Targeted20 paired actual-cost frontiers; 95% paired-object bootstrap CI; no extrapolation');fig.savefig(folder/'paired_cost_frontiers.png',dpi=160);plt.close(fig)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--authorization',required=True);a=ap.parse_args()
 if a.authorization!='TARGETED_SEALED':raise RuntimeError('Explicit parent TARGETED_SEALED required before any truth access')
 t,c=time.perf_counter(),time.process_time();sealpath=ROOT/'audit/targeted/PUBLIC_SEAL.json';sealed=load(sealpath)
 assert sealed['status']=='ALL_TARGETED80_PUBLIC_OUTPUTS_SEALED' and sealed['tasks_count']==80 and sealed['objects']==20
 for s in sealed['evidence']+[s for e in sealed['tasks'] for s in e['files']+e['public_files']]:
  if not unchanged(s):raise RuntimeError('Sealed public artifact changed '+s['path'])
 out=ROOT/'RESULTS/targeted_analysis'
 if out.exists():raise FileExistsError(out)
 cases={r['case_id']:r for r in load(ROOT/'data_public/targeted/cases.json')['records']}
 # First private opens are below this authorization and complete public seal verification.
 mapping={r['case_id']:r for r in load(ROOT/'sealed_truth/targeted/mapping_PRIVATE.json')['records']}
 rows=[];checkpoints=[];truths={};supports={};pictures={};witness=[];budget_audit=[]
 for entry in sealed['tasks']:
  task=entry['task'];ident=task['case_id'];case=cases[ident];directory=resolve(task['output_dir']);rec=load(directory/'RECEIPT.json');cost=rec.get('physical_model_cost') or {};corr=rec.get('corrector') or {}
  if ident not in truths:
   with np.load(resolve(mapping[ident]['truth_path']),allow_pickle=False) as z:
    truths[ident]=z['m'].ravel().copy();supports[ident]=z['mask'].ravel().astype(bool);clean=z['Y_clean'].copy()
   public=public_observation(task);rn=float(np.linalg.norm(clean-public['Y']));rr=rn/float(np.linalg.norm(public['Y']));prior=fixed_prior(truths[ident]);witness.append(dict(case_id=ident,GT_Phi=.5*rr**2+prior,GT_prior=prior,GT_relative_residual=rr,GT_legal=bool(truths[ident].min()>=.02 and truths[ident].max()<=20),new_physical_calls=0))
  row=dict(task_id=task['task_id'],case_id=ident,object_id=ident,arm=task['arm'],status=entry['status'],condition=0 if case['family']=='curved_bar' else 1,family=case['family'],inverse_crime=True,denominator_included=True,material_image_success=False,same_grid_material_noise_success=False,material_nrmse=None,otsu_iou=None,relative_residual=None,original_Phi=rec.get('physical_objective'),cpu_seconds=rec.get('cpu_seconds'),wall_seconds=rec.get('wall_seconds'),auxiliary_cpu_seconds=rec.get('auxiliary_cpu_seconds'),LU=cost.get('LU'),actual_complex_RHS=cost.get('actual_complex_RHS'),corrector_stop=corr.get('stop'),corrector_success=corr.get('success'))
  if entry['status']=='SEALED_PUBLIC_PREDICTION':
   with np.load(directory/'prediction.npz',allow_pickle=False) as z:pred={k:z[k].copy() for k in z.files}
   m=pred['m'].ravel();public=public_observation(task);row.update(material_metrics(m,truths[ident],supports[ident],case['amplitude']));row.update(noise_metrics(pred['Y_prediction'],public,True));row['same_grid_material_noise_success']=bool(row['material_image_success'] and row['noise_CN99_pass']);row['independent_Phi']=.5*row['relative_residual']**2+fixed_prior(m);row['Phi_receipt_difference']=row['independent_Phi']-rec['physical_objective'];pictures[(ident,task['arm'])]=m
   with np.load(directory/'BUDGET_STATES.npz',allow_pickle=False) as z:columns=[str(v) for v in z['columns']];costs=z['costs'].copy();materials=z['m'].copy()
   assert columns==['CPU_seconds_including_process_prefix','LU','complex_RHS','original_Phi'] and materials.shape==(len(costs),576)
   assert not np.any(np.diff(costs[:,:3],axis=0)<-1e-9)
   for idx,(cc,mat) in enumerate(zip(costs,materials)):
    checkpoints.append(dict(case_id=ident,arm=task['arm'],index=idx,CPU=float(cc[0]),LU=float(cc[1]),complex_RHS=float(cc[2]),original_Phi=float(cc[3]),**material_metrics(mat,truths[ident],supports[ident],case['amplitude'])))
   budget_audit.append(dict(case_id=ident,arm=task['arm'],checkpoints=len(costs),final_CPU_unobserved_tail=float(rec['cpu_seconds']-costs[-1,0]),final_material_float32_max_difference=float(np.max(abs(materials[-1].astype(float)-m)))))
  rows.append(row)
 rng=np.random.default_rng(2026110429);summary=aggregate(rows,['arm']);families=aggregate(rows,['family','arm']);pairs=[bootstrap_pair(rows,*comp,rng) for comp in COMPARISONS]
 for p in pairs:p['multiple_comparison_adjustment']='NONE; secondary selected-family confirmation; original two contrasts declared before targeted draws; cie_single added based only on disclosed primary30 before any targeted truth scoring'
 family_pairs=[dict(family=f,**bootstrap_pair([r for r in rows if r['family']==f],*comp,rng)) for f in ('curved_bar','slit_shell') for comp in COMPARISONS]
 frontiers,frontsummary=pairwise_frontiers(rows,checkpoints,rng,comparisons=COMPARISONS)
 out.mkdir();write(out/'SCORES.json',dict(rows=rows,objects=20,arms=ARMS,new_physical_calls=0));csvwrite(out/'SCORES.csv',rows);write(out/'CHECKPOINT_SCORES.json',dict(rows=checkpoints,new_physical_calls=0));write(out/'ARM_SUMMARY.json',summary);write(out/'FAMILY_SUMMARY.json',families);write(out/'PAIRED_BOOTSTRAP.json',pairs);write(out/'FAMILY_PAIRED_BOOTSTRAP.json',family_pairs);write(out/'PAIRWISE_COST_FRONTIERS.json',frontiers);write(out/'PAIRWISE_COST_FRONTIER_SUMMARY.json',frontsummary);write(out/'BUDGET_AUDIT.json',budget_audit);write(out/'GT_FEASIBLE_WITNESS.json',dict(records=witness,new_physical_calls=0,scope='Offline diagnostic only, no uniqueness or global-optimum certificate'))
 plot(rows,truths,pictures,ROOT/'FIGURES/targeted',frontsummary)
 lines=['# R7 targeted selected-family secondary confirmation','','Twenty fresh independent materials, all retained. Ten curved_bar at contrast12/six full sources; ten slit_shell at contrast18/12 half sources; relative complex noise .001. This conditional study follows selection of two successful family/acquisition strata in primary30; results are never pooled and cannot establish universal recovery. Primary material success requires amplitude-normalized NRMSE<=.15, Otsu IoU>=.7, and legal material box. CN99 is separate.','','|Arm|Material successes /20|Material+CN99 /20|Mean NRMSE|Median NRMSE|CPU total s|LU total|Complex RHS total|','|---|---:|---:|---:|---:|---:|---:|---:|']
 for r in summary:lines.append('|%s|%d|%d|%.5g|%.5g|%.2f|%.0f|%.0f|'%(r['arm'],r['material_successes'],r['same_grid_joint_successes'],r['material_nrmse']['mean'],r['material_nrmse']['median'],r['cpu_seconds']['total'],r['LU']['total'],r['actual_complex_RHS']['total']))
 lines+=['','|Family|Arm|Material successes /10|Mean NRMSE|Median NRMSE|','|---|---|---:|---:|---:|']
 for r in families:lines.append('|%s|%s|%d|%.5g|%.5g|'%(r['family'],r['arm'],r['material_successes'],r['material_nrmse']['mean'],r['material_nrmse']['median']))
 lines+=['','|Comparison|Metric|Mean left-minus-right|95% paired-object bootstrap CI|','|---|---|---:|---|']
 for p in pairs:
  for field in ('material_nrmse','material_image_success','cpu_seconds','LU','actual_complex_RHS'):
   q=p['metrics'].get(field)
   if q:lines.append('|%s vs %s|%s|%.6g|[%.6g, %.6g]|'%(p['left'],p['right'],field,q['mean_difference_left_minus_right'],*q['mean_difference_CI95']))
 lines+=['','Ten-thousand paired20-object bootstrap resamples, seed2026110429, no multiplicity adjustment. Family10 tables are conditional descriptive evidence. All20 material realizations are independent; arm/acquisition rows are not new objects. Failures remain in success denominators; continuous metrics disclose finite pairs.','','Actual CPU/LU/RHS comparisons use each object/pair observed common window [max(first saved cost),min(last saved cost)] at11 fixed fractions, taking the last paid original-Phi incumbent at or below budget. No interpolation, extrapolation or GT iterate choice. Sparse auxiliary and every-ten-observation checkpoint resolution and float32 storage can leave unused cost; exact budgets and tails are recorded. CPU includes setup and auxiliary work; equal LU/RHS is not equal CPU. CN99 reflects same-grid circular complex Gaussian noise. No new forward, LU or RHS calls occur in scoring.']
 (out/'ANALYSIS.md').write_text('\n'.join(lines)+'\n')
 write(out/'SCORE_RECEIPT.json',dict(status='COMPLETE',authorization='TARGETED_SEALED explicitly received',created_utc=datetime.now(timezone.utc).isoformat(),objects=20,tasks=80,truth_access='POST_COMPLETE80_TASK_PUBLIC_SEAL_ONLY',new_physical_calls=0,new_physical_LU=0,new_physical_RHS=0,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,bootstrap_replicates=10000,bootstrap_seed=2026110429,scope='SECONDARY_SELECTED_FAMILIES_NOT_POOLED'))
 print(json.dumps(dict(status='TARGETED_SCORED',rows=len(rows),objects=20,out=str(out))))
if __name__=='__main__':main()

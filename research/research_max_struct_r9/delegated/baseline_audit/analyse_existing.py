"""Read-only R4--R8 evidence audit; no solver imports, forward solves or GT reads."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import csv,json,time,collections
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
OUT=Path(__file__).resolve().parent
REPO=OUT.parents[3]
OLD=REPO.parent/'r8_worktree'
R8=REPO/'research/research_max_struct_r8'
OLD8=OLD/'research/research_max_struct_r8'

def read(p):return json.loads(p.read_text())
def dump(name,v):(OUT/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def loadcsv(p):return list(csv.DictReader(p.open()))
def table(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def f(r,k):return float(r[k])
def yes(r,k):return str(r.get(k,'')).lower()=='true'
def fields(r,keys):return {k:r.get(k) for k in keys}

def main():
 cpu=time.process_time();wall=time.perf_counter()
 rows=loadcsv(R8/'delegated/warm_analysis/PRIMARY_630_object_rows.csv');assert len(rows)==630
 manifest=read(R8/'configs/warm_primary_manifest.json');tasks={t['task_id']:t for t in manifest['tasks']}
 uni=[r for r in rows if r['initialization']=='uniform'];assert len(uni)==294
 keys=['case_id','cohort','arm','initialization','hardware','amplitude','sources','aperture_fraction','relative_complex_RMS_noise','generation_grid','inverse_grid','material_nrmse','otsu_iou','material_image_success','noise_CN99_pass','relative_residual','original_Phi','initial_material_nrmse','initial_original_Phi','selected_auxiliary_stage','auxiliary_iterations','corrector_iterations','corrector_evaluations','corrector_success','corrector_stop','projected_gradient_inf','solver_cpu_seconds','auxiliary_cpu_seconds','LU','complex_RHS']
 table('R8_ENDPOINTS_630.csv',[fields(r,keys) for r in rows])
 strata=[]
 for kind in ['ALL','AMPLITUDE','NOISE','HIGH_VS_LOW','COHORT','HARDWARE']:
  def group(r):
   return {'ALL':'all42','AMPLITUDE':r['amplitude'],'NOISE':r['relative_complex_RMS_noise'],'HIGH_VS_LOW':'high>=12' if f(r,'amplitude')>=12 else 'low<12','COHORT':r['cohort'],'HARDWARE':r['hardware']}[kind]
  for g in sorted(set(map(group,uni))):
   for arm in sorted(set(r['arm'] for r in uni)):
    rr=[r for r in uni if group(r)==g and r['arm']==arm]
    strata.append(dict(stratum=kind,value=g,arm=arm,n=len(rr),successes=sum(yes(r,'material_image_success') for r in rr),mean_nrmse=float(np.mean([f(r,'material_nrmse') for r in rr])),median_nrmse=float(np.median([f(r,'material_nrmse') for r in rr])),corrector_CAP_or_failure=sum(not yes(r,'corrector_success') for r in rr),selected_stage_minus1=sum(int(r['selected_auxiliary_stage'])==-1 for r in rr),mean_CPU=float(np.mean([f(r,'solver_cpu_seconds') for r in rr])),hardware_scope=','.join(sorted(set(r['hardware'] for r in rr))),scope='CPU mixedhardware descriptiveonly unless HARDWARE stratum; all DEV now'))
 table('R8_UNIFORM_STRATIFIED.csv',strata)
 by={(r['case_id'],r['arm']):r for r in uni};pairs=[]
 for a,b,label in [('a_plus','a_material','Aplus_minus_A_same_global'),('a_plus','global_dual','Aplus_minus_Global_same_global'),('a_plus_sourcewise','cie_dual','AplusSource_minus_Dual_same_sourcewise'),('a_plus','cie_dual','AplusGlobal_minus_Dual_CONFOUNDED_CURRENT')]:
  for case in sorted(set(r['case_id'] for r in uni)):
   x,y=by[case,a],by[case,b];same=x['hardware']==y['hardware'];pairs.append(dict(case_id=case,comparison=label,amplitude=f(x,'amplitude'),noise=f(x,'relative_complex_RMS_noise'),cohort=x['cohort'],arm_a=a,arm_b=b,nrmse_a=f(x,'material_nrmse'),nrmse_b=f(y,'material_nrmse'),nrmse_difference=f(x,'material_nrmse')-f(y,'material_nrmse'),success_a=yes(x,'material_image_success'),success_b=yes(y,'material_image_success'),a_only_success=yes(x,'material_image_success') and not yes(y,'material_image_success'),b_only_success=yes(y,'material_image_success') and not yes(x,'material_image_success'),same_hardware=same,CPU_difference=f(x,'solver_cpu_seconds')-f(y,'solver_cpu_seconds') if same else None,aux_CPU_difference=f(x,'auxiliary_cpu_seconds')-f(y,'auxiliary_cpu_seconds') if same else None,final_Phi_difference=f(x,'original_Phi')-f(y,'original_Phi')))
 table('R8_MATCHED_PAIRS.csv',pairs)
 ps=[]
 for label in sorted(set(r['comparison'] for r in pairs)):
  for g in ['all42','high>=12','low<12']:
   rr=[r for r in pairs if r['comparison']==label and (g=='all42' or (r['amplitude']>=12)==(g=='high>=12'))];ps.append(dict(comparison=label,stratum=g,n=len(rr),mean_nrmse_difference=float(np.mean([r['nrmse_difference'] for r in rr])),median_nrmse_difference=float(np.median([r['nrmse_difference'] for r in rr])),a_successes=sum(r['success_a'] for r in rr),b_successes=sum(r['success_b'] for r in rr),a_only=sum(r['a_only_success'] for r in rr),b_only=sum(r['b_only_success'] for r in rr)))
 table('R8_PAIR_SUMMARY.csv',ps)
 selections=[];candidates=[];missing=[]
 for row in rows:
  t=tasks[row['task_id']];p=OLD/t['output_dir']/'RECEIPT.json'
  if not p.exists():missing.append(str(p));continue
  r=read(p);cc=r.get('auxiliary_candidate_scores',[]);selected=r.get('auxiliary_selected_stage');chosen=next((x for x in cc if x['stage']==selected),None);last=cc[-1] if cc else None
  for c in cc:candidates.append(dict(task_id=row['task_id'],case_id=row['case_id'],arm=row['arm'],initialization=row['initialization'],amplitude=f(row,'amplitude'),**c,selected=c['stage']==selected))
  selections.append(dict(**fields(row,keys),selected_pre_Phi=chosen['objective'] if chosen else None,last_auxiliary_stage=last['stage'] if last else None,last_auxiliary_Phi=last['objective'] if last else None,initial_candidate_Phi=cc[0]['objective'] if cc else None,corrector_Phi_reduction=(chosen['objective']-r['physical_objective']) if chosen else None,selected_is_last_auxiliary_stage=bool(chosen and last and chosen['stage']==last['stage']),receipt_path=str(p),sourcewise_cg=r.get('auxiliary_execution_receipt',{}).get('sourcewise_cg')))
 table('R8_SELECTION_AND_CORRECTOR_630.csv',selections);table('R8_AUXILIARY_CANDIDATE_SCORES.csv',candidates)
 # Spectral diagnostic only: SVD of public receiver matrix; no Maxwell solve,
 # clean observation or material truth loaded. CN noise expectation is conditional
 # on fixed S and the declared independent circular complex sigma.
 spectra={};noise=[]
 case_rows={r['case_id']:r for r in uni}
 with threadpool_limits(limits=1):
  for case,row in sorted(case_rows.items()):
   task=tasks[row['task_id']];opref=task['public_files'][0]
   if opref not in spectra:
    with np.load(OLD/opref) as z:S=z['S'].copy()
    U,s,Vh=svd(S,full_matrices=False,check_finite=False);rank=int(np.sum(s>.001*s[0]));spectra[opref]=(U,s,Vh,rank)
   U,s,Vh,rank=spectra[opref]
   obsref=task['public_files'][1]
   with np.load(OLD/obsref) as z:Y=z['Y'].copy();sigma=float(z['sigma'])
   coef=(U[:,:rank].conj().T@Y)/s[:rank,None];JDnorm=float(np.linalg.norm(coef));expected=Y.shape[1]*sigma*sigma*float(np.sum(s[:rank]**-2));noise.append(dict(case_id=case,operator=opref,amplitude=f(row,'amplitude'),relative_complex_RMS_noise=f(row,'relative_complex_RMS_noise'),sigma=sigma,sources=Y.shape[1],receiver_rank=rank,retained_condition_number=float(s[0]/s[rank-1]),s_min_retained_relative=float(s[rank-1]/s[0]),expected_JD_noise_Frobenius_norm=float(np.sqrt(expected)),observed_JD_Frobenius_norm=JDnorm,expected_noise_norm_over_observed_JD=float(np.sqrt(expected)/JDnorm),qualification='Expected RMS norm not realized error; conditional on declared iidCN noise; JD denominator itself noisy; no GT, forward or inversion calls'))
 spectrum_rows=[dict(operator=k,rank=v[3],s_max=float(v[1][0]),s_min_retained=float(v[1][v[3]-1]),s_min_retained_relative=float(v[1][v[3]-1]/v[1][0]),first_excluded_relative=float(v[1][v[3]]/v[1][0]) if v[3]<len(v[1]) else None,retained_condition_number=float(v[1][0]/v[1][v[3]-1]),singular_values=v[1].tolist()) for k,v in spectra.items()]
 dump('DPCS_PUBLIC_OPERATOR_SPECTRA.json',spectrum_rows);table('DPCS_CASE_NOISE_SEVERITY.csv',noise)
 r7=loadcsv(REPO/'research/research_max_struct_r7/RESULTS/scoring/primary/SCORES.csv')
 table('R7_ENDPOINTS_180.csv',[fields(r,['case_id','arm','family','amplitude','sources','aperture','noise_relative_complex_RMS','material_nrmse','otsu_iou','material_image_success','relative_residual','original_Phi','cpu_seconds','auxiliary_cpu_seconds','corrector_iterations','corrector_success','corrector_stop']) for r in r7])
 hist=read(R8/'RESULTS/history_warm/SCORES.json')['rows'];table('HISTORICAL_CONTROL_SCORES.csv',[fields(r,['case_id','method','material_nrmse','otsu_iou','material_image_success','relative_residual','original_Phi','cpu_seconds','corrector_iterations','corrector_stop','prediction_source','receipt_source']) for r in hist])
 r3=read(REPO/'research/research_max_struct_r3/results/diagnostics/NOISE_FLOOR_DISTRIBUTION.json');r3['qualification_added']='Fixed published GT anchor only; not final material error lowerbound; clean/release controls still fail.';dump('R3_HARDLOCK_EVIDENCE.json',r3)
 table('R3_DEVELOPMENT_CONTROL_SCORES.csv',loadcsv(REPO/'research/research_max_struct_r3/results/DEVELOPMENT_SCORES.csv'))
 wide=loadcsv(R8/'delegated/warm_mechanism_analysis/ALL42_CASE_TABLE.csv')
 survival=[]
 for row in wide:
  base=case_rows[row['case_id']]
  for arm in ['cie_dual','global_dual','a_plus','a_plus_sourcewise']:
   for init in ['uniform','nn']:
    p=arm+'__'+init+'__';survival.append(dict(case_id=row['case_id'],amplitude=f(base,'amplitude'),noise=f(base,'relative_complex_RMS_noise'),arm=arm,initialization=init,chosen_stage=row.get(p+'chosen_stage'),initial_Phi=row.get(p+'initial_Phi'),selected_pre_Phi=row.get(p+'selected_pre_corrector_Phi'),final_Phi=row.get(p+'original_Phi'),final_nrmse=row.get(p+'material_nrmse'),initial_nrmse=row.get(p+'initial_material_nrmse'),lower_pre_Phi_but_final_Phi_and_material_worse_NNdirect=row.get(p+'lower_pre_Phi_but_both_final_Phi_and_material_worse_than_NNdirect'),DPCS_rank=row.get('projection_DPCS_rank'),current_DPCS_displacement=row.get('projection_DPCS_displacement'),current_Fourier_displacement=row.get('projection_Fourier_displacement')))
 table('R8_INITIALIZER_SELECTION_SURVIVAL.csv',survival)
 dump('ANALYSIS_RECEIPT.json',dict(status='COMPLETED_READ_ONLY_EVIDENCE_ANALYSIS',R8_test_status='PUBLIC_DEVELOPMENT',input_R8_score_rows=len(rows),receipt_rows_read=len(selections),missing_receipts=missing,operator_SVDs=len(spectra),physical_forward_calls=0,physical_adjoint_calls=0,inverse_runs=0,training_runs=0,new_GT_files_opened=0,online_choices=0,SHA_checks=0,scope='Reads existing public S/Y/sigma and sealed paid scores/receipts; writes only baseline_audit; no solver import; no GT-based new online choices',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
 print(json.dumps(dict(rows=len(rows),missing=len(missing),operators=len(spectra),noise_ratios=[min(r['expected_noise_norm_over_observed_JD'] for r in noise),max(r['expected_noise_norm_over_observed_JD'] for r in noise)],pair_summary=ps)))

if __name__=='__main__':main()

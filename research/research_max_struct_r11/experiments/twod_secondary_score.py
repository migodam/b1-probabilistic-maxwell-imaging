"""Score pre-GT sealed conditional views by reindexing authorized primary scores only."""
import argparse,csv,json,time
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def csvrows(p):
 with Path(p).open(newline='') as f:return list(csv.DictReader(f))
def finite(v):
 try:return float(v) if v not in ('',None) else None
 except (ValueError,TypeError):return None
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--public-seal',type=Path,required=True);ap.add_argument('--primary-scores',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--authorize-postscore',required=True);ap.add_argument('--global-appendix-seal',type=Path,default=ROOT/'delegated/twod/independent48/secondary/global_appendix/public_seal');a=ap.parse_args();wall=time.perf_counter()
 if a.authorize_postscore!='ROOT_APPROVED_SECONDARY_POSTSCORE':raise PermissionError('Root authorization required; no GT is read by this script')
 seal=read(a.public_seal/'RECEIPT.json');primary=read(a.primary_scores/'RECEIPT.json');gate=read(a.primary_scores/'PUBLIC_GATE_BEFORE_TRUTH.json')
 if seal['status']!='ALL96_PUBLIC_DEPLOYMENTS_SEALED_BEFORE48_GT' or primary['status']!='SCORED_POST_ALL288_TERMINAL' or gate['status']!='PASS_ALL288_TERMINAL':raise ValueError('Pre-GT96 seal and completed authorized primary scores required')
 if 'created_utc' not in gate or datetime.fromisoformat(seal['seal_finished_utc'])>datetime.fromisoformat(gate['created_utc']):raise ValueError('Public secondary seal must precede primaryGT gate timestamp')
 appendix=read(a.global_appendix_seal/'RECEIPT.json')
 if appendix['status']!='ADDITIONAL48_GLOBAL_PUBLIC_CHOICES_SEALED_BEFORE48_GT' or datetime.fromisoformat(appendix['seal_finished_utc'])>datetime.fromisoformat(gate['created_utc']):raise ValueError('Append-only48 Global seal must also precede primaryGT gate')
 deployments=read(a.public_seal/'DEPLOYMENT_SELECTIONS.json')['rows']+read(a.global_appendix_seal/'DEPLOYMENT_SELECTIONS.json')['rows'];scored=[r for r in csvrows(a.primary_scores/'SCORES_AND_FIXED_BUDGETS.csv') if r['view']=='FINAL_CAP_ELIGIBLE'];lookup={(r['case_id'],r['method']):r for r in scored}
 if len(deployments)!=144 or len(lookup)!=288 or len({r['case_id'] for r in deployments})!=48:raise ValueError('48objects,144deployments and288primary source rows required')
 a.output.mkdir(parents=True,exist_ok=False);rows=[]
 metric_keys=['NRMSE','RMSE','fixed_support_IoU','material_legal','material_success','material_min','material_max','true_low_cells','RMSE_true_low_region','true_high_cells','RMSE_true_high_region','true_background_cells','RMSE_true_background_region','data_relative_residual','data_residual_Frobenius_norm','data_residual_over_sigma','CN99_radius_over_sigma','CN99_pass','complex_observations','public_sigma']
 for deployment in deployments:
  case=deployment['case_id'];baseline=lookup[(case,'cie_single')];chosen=deployment['selected'];row=dict(case_id=case,method=deployment['method'],family=baseline['family'],regime=baseline['regime'],high_chi=baseline['high_chi'],trigger=deployment['trigger'],selected_task_id=chosen['task_id'] if chosen else None,selected_source_method=chosen['method'] if chosen else None,selected_original_Phi=chosen['Phi'] if chosen else None,selection='Frozen public originalPhi, exact ties CIE; no GT',deployment_scientific_prefix_child_CPU_seconds=deployment['deployment_full_child_CPU_seconds'],deployment_end_to_end_CPU_exact=deployment.get('whole_deployment_end_to_end_CPU_seconds'),deployment_end_to_end_CPU_lower_bound=deployment.get('whole_deployment_CPU_lower_bound'),deployment_end_to_end_CPU_upper_observed_bound=deployment.get('whole_deployment_CPU_upper_observed_bound'),exact_end_to_end_cost_qualified=deployment.get('exact_end_to_end_cost_qualified',True),engineering_affected_object=deployment.get('engineering_affected_object',False),deployment_physical_LU=deployment['deployment_physical_LU'],deployment_physical_complex_RHS=deployment['deployment_physical_complex_RHS'],deployment_cost_additive=False,available=chosen is not None)
  if chosen:
   source=lookup[(case,chosen['method'])]
   if finite(source.get('original_Phi'))!=float(chosen['Phi']):raise ValueError('Selected publicPhi differs from immutable primary score source')
   for key in metric_keys:row[key]=source.get(key)
   row['NRMSE_difference_from_CIE']=finite(source.get('NRMSE'))-finite(baseline.get('NRMSE')) if finite(source.get('NRMSE')) is not None and finite(baseline.get('NRMSE')) is not None else None
  else:
   for key in metric_keys:row[key]=None
   row['material_success']=False;row['NRMSE_difference_from_CIE']=None
  rows.append(row)
 keys=list(rows[0])
 with (a.output/'SECONDARY_DEPLOYMENT_SCORES.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 summary=[]
 for name in ('CIE_then_G19_PGP','CIE_then_DualCIE','CIE_then_GlobalDual'):
  rr=[r for r in rows if r['method']==name];nr=sorted(finite(r['NRMSE']) for r in rr if finite(r['NRMSE']) is not None);cpu=[r['deployment_scientific_prefix_child_CPU_seconds'] for r in rr];middle=len(nr)//2;median=(nr[middle] if len(nr)%2 else (nr[middle-1]+nr[middle])/2) if nr else None
  summary.append(dict(method=name,objects=48,finite_material_scores=len(nr),missing_material_scores=48-len(nr),material_success=sum(r['material_success'] is True or r['material_success']=='True' for r in rr),triggered=sum(r['trigger'] for r in rr),mean_NRMSE=sum(nr)/len(nr) if nr else None,median_NRMSE=median,deployment_view_scientific_prefix_child_CPU_seconds=sum(cpu) if all(c is not None for c in cpu) else None,source_choices={m:sum(r['selected_source_method']==m for r in rr) for m in ('cie_single','fullG19','cie_dual','global_dual')},cost_view_not_additive=True))
 receipt=dict(status='SECONDARY_SCORED_BY_FROZEN_SELECTION_REINDEX',summary=summary,GT_access='NO_ADDITIONAL_GT_READ; copies metrics from already-authorized primary scored CSV',new_physical_calls=0,new_inverse_runs=0,new_candidate_choice=0,policy_registered_utc=seal['registered_utc'],public_seal_finished_utc=seal['seal_finished_utc'],primary_GT_gate_created_utc=gate['created_utc'],registration_context=seal['registration_context'],global_appendix_registered_utc=appendix['registered_utc'],global_appendix_seal_finished_utc=appendix['seal_finished_utc'],original96_unchanged=True,total_deployment_choices=144,budget_scope='Each conditional deployment can pay two480CPU attempt allowances. This is not the same actual CPU or budget as one CIE trial; no same-cost superiority inferred.',CPU_seconds=time.process_time(),CPU_scope='Absolute whole postscore reindex process CPU including standard imports',wall_seconds=time.perf_counter()-wall)
 dump(a.output/'RECEIPT.json',receipt);dump(a.output/'SUMMARY.json',summary);print(json.dumps(dict(status=receipt['status'],deployments=len(rows))))
if __name__=='__main__':main()

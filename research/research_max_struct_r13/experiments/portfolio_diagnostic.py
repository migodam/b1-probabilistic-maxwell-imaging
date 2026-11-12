"""Three fixed POSTHOC_COUNTERFACTUAL_FROM_PAID_RUNS; no new inverse runs.

Phase1 reads only public primary metadata and whole-job exit costs. It writes
selection references BEFORE importing already-scored truth metrics offline.
"""
import argparse,csv,json,math,platform,statistics,sys,time
from pathlib import Path
PARENTS=('cie','dual_cie','g19')
STATUS='POSTHOC_COUNTERFACTUAL_FROM_PAID_RUNS'
def read(path):return json.loads(path.read_text())
def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def csvfile(path,rows):
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with path.open('w',newline='') as stream:w=csv.DictWriter(stream,fieldnames=keys);w.writeheader();w.writerows(rows)
def addcost(a,b):return a+b if a is not None and b is not None else None
def main():
 cpu,wall=time.process_time(),time.perf_counter();p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();root=a.root.resolve();auth=read(a.authorization)
 if not auth.get('authorized_joint_offline_scoring') or not auth.get('authorized_aggregation'):raise PermissionError('Parent joint authorization required')
 sys.path.insert(0,str(root/'src'));from seal_results import validate_seal
 for campaign in auth['required_campaigns']:
  seal=read(root/auth['seals'][campaign]);validate_seal(seal,require_authorization=True)
  if Path(seal['root']).resolve()!=root:raise ValueError('Different sealed root')
 out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);sources=[str(a.authorization.resolve())];campaign='confirmation_core48';queue=root/'queues'/campaign/'RECEIPT.json';sources.append(str(queue));q=read(queue);jobs={(j['case_id'],j['method']):j for j in q['jobs']};cases=[f'r13_fresh_{i:04d}' for i in range(48)];public={}
 # No score/GT tables are read in this phase. No NPZ is ever opened.
 for cid in cases:
  for method in PARENTS+('dr',):
   job=jobs.get((cid,method));folder=root/'RESULTS'/campaign/(job['task_id'] if job else cid+'__'+method);path=folder/'primary.json';sources.append(str(path));value=read(path) if path.is_file() else{};public[(cid,method)]={k:value.get(k) for k in ['candidate','Phi','legal','eligible','residual_over_delta']};public[(cid,method)].update(task_id=job['task_id'] if job else None,external_CPU_seconds=(job.get('external') or{}).get('child_process_CPU_seconds') if job else None)
 selections=[]
 for parent in PARENTS:
  for cid in cases:
   pp=public[(cid,parent)];dr=public[(cid,'dr')];ratio=pp.get('residual_over_delta');trigger=ratio>1 if isinstance(ratio,(int,float)) and math.isfinite(ratio) else None;options=[]
   for method,row in [(parent,pp)]+([('dr',dr)] if trigger else[]):
    if row.get('legal') is True and row.get('eligible') is True and isinstance(row.get('Phi'),(int,float)) and math.isfinite(row['Phi']) and row.get('candidate') is not None:options.append((row['Phi'],int(row['candidate']),0 if method==parent else 1,method,row))
   selected=min(options) if options and trigger is not None else None
   selections.append(dict(status=STATUS,case_id=cid,strategy=parent+'_to_uniform_DR',parent_method=parent,triggered=trigger,bypass=trigger is False,parent_candidate=pp.get('candidate'),parent_Phi=pp.get('Phi'),parent_residual_over_delta=ratio,uniform_DR_candidate=dr.get('candidate') if trigger else None,uniform_DR_Phi=dr.get('Phi') if trigger else None,selected_method=selected[3] if selected else None,selected_task_id=selected[4]['task_id'] if selected else None,selected_candidate=selected[4]['candidate'] if selected else None,selected_Phi=selected[4]['Phi'] if selected else None,parent_whole_job_CPU_seconds=pp['external_CPU_seconds'],trigger_uniform_DR_whole_job_CPU_seconds=dr['external_CPU_seconds'] if trigger else 0 if trigger is False else None,counterfactual_CPU_seconds=addcost(pp['external_CPU_seconds'],dr['external_CPU_seconds'] if trigger else 0 if trigger is False else None),host='i9_remote',selection_rule='trigger parent residual/delta>1; minimum original Phi among legal eligible parent and triggered uniform-start DR primary; ties candidate index then parent; no GT',available=selected is not None))
 write(out/'PUBLIC_SELECTIONS_BEFORE_SCORE_ATTACHMENT.json',dict(status=STATUS,truth_fields_read=False,NPZ_read=False,selections=selections));csvfile(out/'PUBLIC_SELECTIONS.csv',selections)
 # Phase2 attaches existing offline metrics without changing any public choice.
 scorepath=root/'analysis/combined/ALL_ESTIMATOR_ROWS.json';deploypath=root/'analysis/combined/CONDITIONAL_DEPLOYMENT_ALL48.json';costpath=root/'analysis/combined/CONDITIONAL_DEPLOYMENT_COSTS.json';sources +=[str(scorepath),str(deploypath),str(costpath)];scored={(r['case_id'],r['method']):r for r in read(scorepath) if r['campaign']==campaign and r['view']=='primary'};registered=[r for r in read(deploypath) if r['method'] in ['dr_rescue','event_rescue'] and r['view']=='primary'];costs={(r['case_id'],r['method']):r for r in read(costpath)};metric_keys=['NRMSE','RMSE','IoU','boxrange_NRMSE','relativeL2','D','R','Phi','residual_over_delta','noise_feasible'];rows=[]
 for selection in selections:
  r=dict(selection);s=scored.get((r['case_id'],r['selected_method']),{})
  if s and s.get('selection_candidate')!=r['selected_candidate']:raise ValueError('Scored primary differs from frozen public selection')
  r.update({k:s.get(k) for k in metric_keys});r['offline_metric_available']=bool(s.get('available'));rows.append(r)
 for s in registered:
  c=costs.get((s['case_id'],s['method']),{});r=dict(status='REGISTERED_ALL48_BAD_ENDPOINT_RESCUE_COMPARISON',strategy='registered_g19_to_'+('DR' if s['method']=='dr_rescue' else 'event'),case_id=s['case_id'],host=s['host'],triggered=s.get('triggered'),bypass=s.get('bypass'),available=s.get('available'),offline_metric_available=s.get('available'),counterfactual_CPU_seconds=c.get('deployment_external_CPU_seconds'),parent_whole_job_CPU_seconds=c.get('parent_external_CPU_seconds'),trigger_uniform_DR_whole_job_CPU_seconds=None,registered_rescue_CPU_seconds=c.get('rescue_external_CPU_seconds'));r.update({k:s.get(k) for k in metric_keys});rows.append(r)
 csvfile(out/'ALL48_PORTFOLIO_AND_REGISTERED_ROWS.csv',rows);write(out/'ALL48_PORTFOLIO_AND_REGISTERED_ROWS.json',rows);summary=[]
 for strategy in sorted({r['strategy'] for r in rows}):
  rr=[r for r in rows if r['strategy']==strategy];item=dict(strategy=strategy,status=rr[0]['status'],expected_objects=48,rows=len(rr),metric_available_objects=sum(r['offline_metric_available'] for r in rr),triggered=sum(r.get('triggered') is True for r in rr),bypass=sum(r.get('bypass') is True for r in rr),CPU_known_objects=sum(r.get('counterfactual_CPU_seconds') is not None for r in rr),total_counterfactual_CPU_seconds=sum(r['counterfactual_CPU_seconds'] for r in rr) if all(r.get('counterfactual_CPU_seconds') is not None for r in rr) else None,host='i9_remote',measured_new_conditional_wall_seconds=None)
  for metric in ['NRMSE','RMSE','IoU']:
   vals=[r[metric] for r in rr if isinstance(r.get(metric),(int,float)) and math.isfinite(r[metric])];item[metric+'_count']=len(vals);item[metric+'_mean']=statistics.mean(vals) if vals else None;item[metric+'_median']=statistics.median(vals) if vals else None
  summary.append(item)
 csvfile(out/'SUMMARY.csv',summary);write(out/'SUMMARY.json',summary)
 comparisons=[];look={(r['case_id'],r['strategy']):r for r in rows}
 for parent in PARENTS:
  for reference in ['registered_g19_to_DR','registered_g19_to_event']:
   for cid in cases:
    left=look[(cid,parent+'_to_uniform_DR')];right=look[(cid,reference)];comparisons.append(dict(case_id=cid,strategy=parent+'_to_uniform_DR',reference=reference,portfolio_NRMSE=left.get('NRMSE'),registered_NRMSE=right.get('NRMSE'),delta_NRMSE=left['NRMSE']-right['NRMSE'] if left.get('NRMSE') is not None and right.get('NRMSE') is not None else None,portfolio_CPU_seconds=left.get('counterfactual_CPU_seconds'),registered_CPU_seconds=right.get('counterfactual_CPU_seconds')))
 csvfile(out/'PAIRED_TO_REGISTERED_ALL48.csv',comparisons)
 (out/'REPORT.md').write_text('# POSTHOC_COUNTERFACTUAL_FROM_PAID_RUNS\n\nExactly three fixed parents CIE/Dual-CIE/G19 on all48 original uniform-start cases. Parent public residual/delta>1 triggers the already-frozen core DR run, also from uniform0.1; choose lower originalPhi among legal eligible primaries, ties candidate index then parent. No GT fields entered selection; public references were serialized before offline score attachment. No primary indices, policy files or seals changed.\n\nThis is a posthoc counterfactual portfolio assembled from paid independent runs, not a new blind confirmation or an actually executed conditional queue. It has no measured new end-to-end wall. Every48 parent job is fully charged; triggered DR gets its full measured whole-job CPU including setup, without proration. Strategies describe alternative deployments and are not additive experimental consumption.\n\nRegistered G19→DR/event comparisons restart from a bad endpoint. Core DR instead starts at the original uniform material: the starts and dispatch experiments differ. No NN prior is used. Missing rows remain in48 denominators with unavailable metrics. No new optimization, candidate, forward, parameter search, success threshold, GT-selected winner or scientific conclusion is introduced. Parent owns interpretation.\n')
 write(out/'RECEIPT.json',dict(status=STATUS,host='Apple_local',platform=platform.platform(),CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,public_selection_rows=len(selections),offline_rows=len(rows),paired_rows=len(comparisons),input_source_files=sorted(set(sources)),new_Maxwell_LU=0,new_RHS=0,new_optimization=0,NPZ_or_truth_pixels_read=False,selection_uses_GT=False,measured_conditional_end_to_end_wall=None,timer_scope='Single diagnostic public selection + offline attachment + serialization process; source solver costs remain references, not new work'))
 print(json.dumps(dict(status=STATUS,public_rows=len(selections),all_strategy_rows=len(rows),paired_rows=len(comparisons),new_Maxwell_LU=0)))
if __name__=='__main__':main()

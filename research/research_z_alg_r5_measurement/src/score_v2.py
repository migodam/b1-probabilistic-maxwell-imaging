"""First permitted material GT reader after all inference/scoring inputs freeze."""
from common import *
import csv,datetime

def decode(x):
 if isinstance(x,dict) and set(x)=={'real','imag'}:return complex(x['real'],x['imag'])
 if isinstance(x,dict):return {k:decode(v) for k,v in x.items()}
 if isinstance(x,list):return [decode(v) for v in x]
 return x
def read(f):return decode(json.loads(Path(f).read_text()))
def rms(m,t):return float(np.linalg.norm(np.asarray(m)-t)/np.sqrt(len(t)))
assert (ROOT/'inputs/prediction_freeze.json').exists()
freeze=read(ROOT/'inputs/prediction_freeze.json')
for path,sha in freeze['files'].items():assert digest(ROOT/path)==sha,path
book=Ledger('GT_scoring_only');truth=np.load(ROOT/'private/truth.npz');bank=np.load(ROOT/'private/response_bank.npz');output={};table=[];root_table=[];source_table=[];costs=[]
qual=json.loads((ROOT/'results/qualification.json').read_text())
for case in ['T0','T1','T2','T3','T4','T5']:
 t=truth[case];output[case]=dict(truth=t,arms={});s=len(t)//2
 for arm in ('clean','noisy'):
  path=ROOT/'results'/case/arm;ar=dict(methods={},strategies={});design=read(path/'design_freeze.json')
  for method in ('Q','TRF64','reciprocal64'):
   d=read(path/(method+'.json'));roots=d['accepted_roots'];errs=[rms(m,t) for m in roots];certs=d['certificates']
   ar['methods'][method]=dict(root_count=len(roots),root_errors=errs,GT_covered_within_1e_5=any(v<1e-5 for v in errs),GT_covered_within_0_1=any(v<=.1 for v in errs),discovery_cpu=d['discovery_CPU_total'],starts_completed=d.get('starts_completed'),note='Set-coverage diagnostic only, not oracle reconstruction')
   for i,m in enumerate(roots):
    cert=next((r.get('certificate',{}) for r in certs if np.linalg.norm(np.array(r['m'])-m)<1e-8),{})
    root_table.append(dict(case=case,arm=arm,method=method,root=i,material=json.dumps(m),GT_RMS=errs[i],certificate=cert.get('status','UNKNOWN')))
   costs.append(dict(case=case,arm=arm,method=method,CPU_seconds=d['discovery_CPU_total'],wall_seconds=d['wall_seconds'],**d['operation_counts']))
  q=read(path/'Q.json');roots=q['accepted_roots']
  for strategy in ('fixed','random','maxmin','fisher','TRF64_active','reciprocal64_active'):
   tr=read(path/(strategy+'_trace.json'));m=tr['corrected']['best'];a=np.array(tr['source_a']);gap=abs((bank[case+'_shifted']-bank[case+'_clean'])@a)
   rr=dict(corrected_RMS=rms(m,t),corrected_success=rms(m,t)<=.1,new_measurements=1,source_label=tr['source_label'],model_error_actual=float(gap),declared_model_allowance_ok=bool(gap<=MODEL_ALLOWANCE),selected_material=m)
   if 'hard_retained' in tr:
    hard=tr['hard_retained'];rr.update(hard_retained=hard,hard_selected_RMS=rms(roots[hard[0]],t) if len(hard)==1 else None,hard_any_close=any(rms(roots[i],t)<=.1 for i in hard))
    tubes=tr['tubes'];kept=[i for i,x in enumerate(tubes) if x['compatible_found']];rr.update(tube_retained=kept,tube_fit_errors=[rms(x['fit'],t) for x in tubes],tube_any_close=any(rms(tubes[i]['fit'],t)<=.1 for i in kept))
   ar['strategies'][strategy]=rr
   source_table.append(dict(case=case,arm=arm,strategy=strategy,source=rr['source_label'],RMS=rr['corrected_RMS'],model_error=gap,model_bound_ok=rr['declared_model_allowance_ok'],new_measurements=1,hard_retained=str(rr.get('hard_retained',[])),tube_retained=str(rr.get('tube_retained',[]))))
  direct=read(path/'direct.json');joint=read(path/'same_active_joint_TRF64.json');ar['direct']=dict(raw_complex_RMS=rms(direct['raw_material'],t),projected_RMS=rms(direct['projected_material'],t),joint_RMS=rms(direct['corrected']['best'],t),new_measurements=direct['new_measurements'],condition=direct['condition'],rank=direct['rank'])
  ar['same_active_joint_TRF64']=dict(RMS=rms(joint['best'],t),success=rms(joint['best'],t)<=.1)
  delta=abs(bank[case+'_shifted'][:s]-bank[case+'_clean'][:s]);ar['old_model_error']=delta;ar['old_model_allowance_ok']=bool(np.all(delta<=MODEL_ALLOWANCE));ar['qualification']=qual[case][arm]['status'];output[case]['arms'][arm]=ar
  table.append(dict(case=case,arm=arm,Q_roots=ar['methods']['Q']['root_count'],TRF_roots=ar['methods']['TRF64']['root_count'],reciprocal_roots=ar['methods']['reciprocal64']['root_count'],qualification=ar['qualification'],maxmin_RMS=ar['strategies']['maxmin']['corrected_RMS'],fixed_RMS=ar['strategies']['fixed']['corrected_RMS'],fisher_RMS=ar['strategies']['fisher']['corrected_RMS'],TRF_active_RMS=ar['strategies']['TRF64_active']['corrected_RMS'],reciprocal_active_RMS=ar['strategies']['reciprocal64_active']['corrected_RMS'],same_active_TRF64_RMS=ar['same_active_joint_TRF64']['RMS'],direct_raw_RMS=ar['direct']['raw_complex_RMS'],direct_joint_RMS=ar['direct']['joint_RMS'],direct_new_measurements=ar['direct']['new_measurements'],old_model_bound_ok=ar['old_model_allowance_ok']))
for name,rows in [('score_table',table),('candidate_registry',root_table),('source_design_ledger',source_table),('discovery_costs',costs)]:
 with (ROOT/'results'/(name+'.csv')).open('x') as f:
  keys=list(dict.fromkeys(k for row in rows for k in row));writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
save(ROOT/'results/scores.json',output);save(ROOT/'results/scoring_receipt.json',book.receipt())
for r in table:print(r)

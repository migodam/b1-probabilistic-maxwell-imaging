"""Join already-authorized main48 CSV to exact-match endpoint JSON; no GT/FW."""
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import csv,json,math,time,itertools
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
BASE=ROOT/'delegated/twod/independent48';METHODS=('cie_single','cie_dual','global_dual','fullG19','freezeA','contractedT');AFFECTED='r11_independent_0022';TOL=.01/1e6
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def boolean(x):
 if x in (True,'True'):return True
 if x in (False,'False'):return False
 if x in ('',None):return None
 raise ValueError('Invalid boolean '+repr(x))
def number(x):
 if x in ('',None):return None
 y=float(x)
 if not math.isfinite(y):raise ValueError('Nonfinite field')
 return y
def quantile(values,q):
 values=sorted(values)
 if not values:return None
 k=(len(values)-1)*q;i=int(k);j=min(i+1,len(values)-1)
 return values[i]+(values[j]-values[i])*(k-i)
def stats(rows,key):
 vals=[number(r.get(key)) for r in rows];v=[x for x in vals if x is not None]
 return dict(known=len(v),sum=sum(v) if len(v)==len(rows) else None,known_subtotal=sum(v),mean=sum(v)/len(v) if v else None,median=quantile(v,.5),min=quantile(v,0),q10=quantile(v,.1),q90=quantile(v,.9),max=quantile(v,1))
def write(p,rows):
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
def main():
 wall=time.perf_counter();src=BASE/'scored/SCORES_AND_FIXED_BUDGETS.csv';ep=BASE/'endpoint_bookkeeping_corrected_path/ENDPOINT_BOOKKEEPING.json'
 with src.open(newline='') as f:primary=[r for r in csv.DictReader(f) if r['view']=='FINAL_CAP_ELIGIBLE']
 end=read(ep);left={(r['case_id'],r['method']):r for r in primary};right={(r['case_id'],r['method']):r for r in end};cases=sorted({k[0] for k in left})
 if len(primary)!=288 or len(left)!=288 or len(end)!=288 or len(right)!=288 or len(cases)!=48 or set(left)!=set(right):raise ValueError('Strict288-to288 unique case/method join required')
 if set(left)!={(c,m) for c in cases for m in METHODS}:raise ValueError('Balanced48x6 original identities required')
 joined=[]
 for key in sorted(left):
  r=dict(left[key]);e=right[key]
  if e['case_id']!=r['case_id'] or e['method']!=r['method']:raise ValueError('Identity mismatch')
  r.update({'endpoint_'+k:v for k,v in e.items() if k not in ('case_id','method')});pg=number(e.get('projected_gradient_inf')) if e['available'] else None
  r.update(PG_stationarity_threshold=TOL,PG_satisfies_1e_minus8=pg<=TOL if pg is not None else None,stationarity_scope='Descriptive original solvergtol/scale only; not new material success gate or local-minimum/information certificate',exclude_from_47_cost_sensitivity=r['case_id']==AFFECTED)
  joined.append(r)
 write(HERE/'JOINED_FINAL288.csv',joined);summaries=[];cross=[];stops=[];grouped={}
 for cohort in ('ALL48','SENSITIVITY47_EXCLUDE_CASE0022'):
  subset=[r for r in joined if cohort=='ALL48' or r['case_id']!=AFFECTED]
  for m in METHODS:
   rr=[r for r in subset if r['method']==m];n=len(rr);pgs=stats(rr,'endpoint_projected_gradient_inf');c=stats(rr,'actual_child_CPU_seconds');lu=stats(rr,'final_physical_LU');rhs=stats(rr,'final_complex_RHS');exact=stats(rr,'exact_end_to_end_CPU');lb=stats(rr,'end_to_end_CPU_lower_bound');ub=stats(rr,'end_to_end_CPU_upper_bound')
   row=dict(cohort=cohort,method=m,materials=n,material_success=sum(boolean(r['material_success']) is True for r in rr),CN99_pass=sum(boolean(r['CN99_pass']) is True for r in rr),PG_available=pgs['known'],PG_satisfies_1e_minus8=sum(r['PG_satisfies_1e_minus8'] is True for r in rr),PG_exceeds_1e_minus8=sum(r['PG_satisfies_1e_minus8'] is False for r in rr),PG_missing=sum(r['PG_satisfies_1e_minus8'] is None for r in rr),actual_valid_child_CPU_total=c['sum'],actual_valid_child_CPU_mean=c['mean'],actual_valid_child_CPU_median=c['median'],exact_end_to_end_CPU_total=exact['sum'],exact_end_to_end_known=exact['known'],end_to_end_CPU_lower_bound_total=lb['sum'],end_to_end_CPU_upper_observed_bound_total=ub['sum'],uncertain_end_to_end_tasks=sum(boolean(r.get('exact_end_to_end_cost_qualified')) is not True for r in rr),physical_LU_total=lu['sum'],physical_LU_mean=lu['mean'],physical_LU_median=lu['median'],complex_RHS_total=rhs['sum'],complex_RHS_mean=rhs['mean'],complex_RHS_median=rhs['median'],mean_NRMSE=stats(rr,'NRMSE')['mean'],median_NRMSE=stats(rr,'NRMSE')['median'],PG_min=pgs['min'],PG_q10=pgs['q10'],PG_median=pgs['median'],PG_q90=pgs['q90'],PG_max=pgs['max'],corrector_stop_counts=dict(Counter(r['corrector_stop'] or 'UNRECORDED' for r in rr)),CPU_scope='Valid final scientific childOSCPU only; original interrupted/setup failedattempt CPU not included here, end-to-end remainsnull when unknown; pair views nonadditive')
   summaries.append(row)
   for reason,count in row['corrector_stop_counts'].items():stops.append(dict(cohort=cohort,method=m,corrector_stop_reason=reason,count=count))
   for cn,pg in itertools.product((True,False,None),(True,False,None)):
    v=[r for r in rr if boolean(r['material_success']) is False and boolean(r['CN99_pass']) is cn and r['PG_satisfies_1e_minus8'] is pg]
    cross.append(dict(cohort=cohort,method=m,material_status='FAIL',CN99='IN' if cn is True else 'OUT' if cn is False else 'MISSING',PG='SATISFIES' if pg is True else 'EXCEEDS' if pg is False else 'MISSING',count=len(v),case_ids=[r['case_id'] for r in v],interpretation='Descriptive cross-count only, not information nonidentifiability or local-minimum classification'))
 write(HERE/'METHOD_SUMMARY.csv',summaries);dump(HERE/'METHOD_SUMMARY.json',summaries);write(HERE/'FAIL_CN99_PG_CROSSCOUNTS.csv',cross);write(HERE/'CORRECTOR_STOP_COUNTS.csv',stops)
 differences=[]
 for candidate,baseline in itertools.product(('fullG19','contractedT'),('global_dual','cie_single','cie_dual')):
  for cid in cases:
   a,b=next(r for r in joined if r['case_id']==cid and r['method']==candidate),next(r for r in joined if r['case_id']==cid and r['method']==baseline)
   rr=dict(case_id=cid,family=a['family'],regime=a['regime'],high_chi=a['high_chi'],candidate=candidate,baseline=baseline,candidate_success=boolean(a['material_success']),baseline_success=boolean(b['material_success']),candidate_only_success=boolean(a['material_success']) and not boolean(b['material_success']),baseline_only_success=boolean(b['material_success']) and not boolean(a['material_success']),exclude_from_47_cost_sensitivity=cid==AFFECTED,precise_end_to_end_CPU_pair=boolean(a.get('exact_end_to_end_cost_qualified')) is True and boolean(b.get('exact_end_to_end_cost_qualified')) is True)
   for field in ('NRMSE','fixed_support_IoU','original_Phi','data_residual_over_sigma','actual_child_CPU_seconds','final_physical_LU','final_complex_RHS','endpoint_projected_gradient_inf'):
    av,bv=number(a.get(field)),number(b.get(field));rr.update({candidate+'_'+field:av,baseline+'_'+field:bv,'delta_'+field:av-bv if av is not None and bv is not None else None})
   for side,r in [('candidate',a),('baseline',b)]:rr[side+'_CN99_pass']=boolean(r['CN99_pass']);rr[side+'_PG_satisfies']=r['PG_satisfies_1e_minus8'];rr[side+'_corrector_stop']=r['corrector_stop'];rr[side+'_exact_end_to_end_CPU']=number(r.get('exact_end_to_end_CPU'))
   differences.append(rr)
 write(HERE/'G19_CONTRACTED_VS_BASELINES_BY_OBJECT.csv',differences);write(HERE/'G19_CONTRACTED_VS_BASELINES_47_SENSITIVITY.csv',[r for r in differences if not r['exclude_from_47_cost_sensitivity']])
 sets=[];equal=[]
 for cohort in ('ALL48','SENSITIVITY47_EXCLUDE_CASE0022'):
  subset=[r for r in joined if cohort=='ALL48' or r['case_id']!=AFFECTED]
  for family in ['ALL']+sorted({r['family'] for r in subset}):
   group=[r for r in subset if family=='ALL' or r['family']==family];ms={m:sorted(r['case_id'] for r in group if r['method']==m and boolean(r['material_success']) is True) for m in METHODS}
   for m,ss in ms.items():sets.append(dict(cohort=cohort,family=family,method=m,success_count=len(ss),success_case_ids=ss))
   for l,r in itertools.combinations(METHODS,2):equal.append(dict(cohort=cohort,family=family,left=l,right=r,success_sets_identical=ms[l]==ms[r],left_only=sorted(set(ms[l])-set(ms[r])),right_only=sorted(set(ms[r])-set(ms[l]))))
 write(HERE/'SUCCESS_SETS_BY_FAMILY.csv',sets);dump(HERE/'SUCCESS_SETS_BY_FAMILY.json',sets);write(HERE/'SUCCESS_SET_EQUALITY.csv',equal)
 mainrows=[r for r in summaries if r['cohort']=='ALL48'];text=['# Main48 endpoint table integration','', 'All288 authorized final cap rows strictly joined1:1 by(case_id,method) with the corrected exact-material-match endpoint table. No scoring/bootstrap/solver rerun; no GT arrays read.','', '|Method|MaterialPASS/48|CN99/48|PG<=1e-8/48|PGmedian|Scientific childCPU total/mean/median s|LU/RHS|','|---|---:|---:|---:|---:|---:|---:|']
 for r in mainrows:text.append(f"|{r['method']}|{r['material_success']}|{r['CN99_pass']}|{r['PG_satisfies_1e_minus8']}|{r['PG_median']:.3g}|{r['actual_valid_child_CPU_total']:.3f} / {r['actual_valid_child_CPU_mean']:.3f} / {r['actual_valid_child_CPU_median']:.3f}|{int(r['physical_LU_total'])} / {int(r['complex_RHS_total'])}|")
 text+=['','PGcriterion is the descriptive originalgtol=.01/scale1e6=1e-8 applied to the saved unscaled originalPhi gradient including prior. It is not a new material success gate, a second-order local-minimum certificate, or information nonidentifiability evidence. Strict bound-sign counts and near1e-8 counts remain in the joined table.','', 'case0022 remains in the main48 table. Its four interrupted/retry arms have unknown original full end-to-endCPU; valid successful-childOSCPU is known, but cannot be relabeled complete task expenditure. Exact end-to-end sums are null when affected, and source lower/observed-upper bounds remain separate. The47 sensitivity excludes the entire case0022 for all methods, preserving matched object sets. No retry is a new material.','', 'See METHOD_SUMMARY/STOP_COUNTS, FAIL_CN99_PG_CROSSCOUNTS, full288 object differences and their47-case sensitivity, SUCCESS_SETS_BY_FAMILY and SUCCESS_SET_EQUALITY. All successes, failures and missing values retained.']
 (HERE/'SUMMARY.md').write_text('\n'.join(text)+'\n')
 dump(HERE/'RECEIPT.json',dict(status='COMPLETE_MECHANICAL_POSTSCORE_INTEGRATION',created_utc=datetime.now(timezone.utc).isoformat(),source_primary=str(src),source_endpoint=str(ep),strict_join_rows=288,independent_materials=48,endpoint_matches=sum(boolean(r['endpoint_available']) is True for r in joined),stationarity_threshold=TOL,criterion_is_new_success_gate=False,new_FW=0,new_adjoint=0,new_solver=0,new_scoring=0,new_bootstrap=0,GT_arrays_read=0,fresh3D_or_paired8_GT_access='NONE',CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall))
 print(json.dumps(mainrows))
if __name__=='__main__':main()

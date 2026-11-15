"""Re-extract named saved scalar tables and source lines; never import solvers/load arrays."""
import csv,json,time
from pathlib import Path
T=time.perf_counter();C=time.process_time()
O=Path(__file__).resolve().parent
R=Path(__file__).resolve().parents[3];W=R.parent.parent
reads=[]
def table(p):
 reads.append(str(p));return [dict(x,source_table=str(p),source_csv_line=i) for i,x in enumerate(csv.DictReader(p.open()),2)]
def save(name,rows):
 fields=list(dict.fromkeys(k for x in rows for k in x))
 with (O/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 return len(rows)
def js(p):reads.append(str(p));return json.loads(p.read_text())
local=table(R/'audit/local_closed_loop_v1/TABLES.csv');s=js(R/'audit/local_closed_loop_v1/SUMMARY.json')
save('LOCAL_DIRECTIONAL.csv',[x for x in local if x['row_type']=='directional'])
save('LOCAL_ADDITIVITY.csv',[x for x in local if x['row_type']=='additivity'])
save('LOCAL_IDEAL_MAPS.csv',[dict(x,source_table=str(R/'audit/local_closed_loop_v1/SUMMARY.json'),source_json_key='ideal_maps') for x in s['ideal_maps']])
new=table(R/'analysis/independent24_v1/ALL_VIEWS.csv');cost=table(R/'audit/cost_ledger_solver_complete_v1/SOLVER_TASK_COSTS.csv');ci={x['task_id']:x for x in cost if x['campaign']=='independent24_v1'}
methods=['g19','cie','dual_cie','pair01','pair11','dr'];views=['primary','raw','corrected']
def costs(x):
 y=dict(x);c=ci.get(x['task_id'],{})
 for k in ['actual_child_OS_CPU_seconds','OS_CPU_null_reason','physical_LU','physical_complex_RHS','physical_matrix_dimension','host','cost_scope']:
  y['whole_child_'+k]=c.get(k,'')
 y['whole_child_cost_source']=c.get('source_table','');y['whole_child_cost_source_line']=c.get('source_csv_line','')
 return y
fields=['case_id','solver','pair','view','campaign','task_id','run_status','initialization','provenance','source_file','source_view','secondary_fallback','available','availability_reason','NRMSE','RMSE','relativeL2','IoU','Phi','D','R','CPU_seconds','wall_seconds','residual_over_delta','noise_feasible','noise_qualification','physical_cost_at_candidate','phase','raw_original_source','raw_original_candidate','raw_original_Phi','initial_Phi','raw_initial_fallback','raw_audit_comparable','source_table','source_csv_line']
filtered=[costs({k:x.get(k,'') for k in fields}) for x in new if x['solver'] in methods and x['view'] in views]
save('INDEPENDENT24_SIX_METHODS_THREE_VIEWS.csv',filtered)
new21=[costs({k:x.get(k,'') for k in fields}) for x in new if x['case_id']=='r14_new_0021' and x['view'] in views];save('R14_NEW_0021_ALL11.csv',new21)
old=table(W/'research/research_max_struct_r13/supplement_aplus48/analysis/stack_comparison/SCOREBOARD_COMPACT.csv')
old21=[{k:v for k,v in x.items() if k!='family'} for x in old if x['case_id']=='r13_fresh_0021' and x['initialization']=='uniform' and x['view'] in views and x['solver'] in ['g19','cie','dual_cie','dr','a_plus_mm_r13','a_plus_arc_r7']];save('R13_FRESH_0021_SIX_METHODS.csv',old21)
winners=[]
for cohort,rows in [('r14_new_0021',new21),('r13_fresh_0021',old21)]:
 for view in views:
  for metric,direction in [('NRMSE','min'),('IoU','max'),('Phi','min')]:
   a=[x for x in rows if x['view']==view and x.get(metric) not in ['',None]];v=(min if direction=='min' else max)(float(x[metric]) for x in a)
   winners.append(dict(case_id=cohort,view=view,metric=metric,direction=direction,value=v,tied_solvers=';'.join(x['solver'] for x in a if float(x[metric])==v),scope='all11' if cohort.startswith('r14') else 'saved_uniform_six_historical'))
save('CASE0021_METRIC_WINNERS.csv',winners)
spaces=[('material h/dm/ec transported/em','real','64','material_error_norm/current_error_material_norm/idealgap_norm','144:161'),('J,dJ,ec','complex','64x6','current_error_norm Frobenius; 384 complex entries','144:160'),('native full U','complex','64x64 with rank48','redundant native current columns; not independent64','79:80'),('independent U,K','complex','64x48','48 independent spatial current directions per source','79:80'),('D observation restricted','complex','16x48','16 receivers; six sources','79:80'),('X_s','complex','48x64 (six blocks)','material direction to independent current coefficients','87:96'),('B,T','real','64x64','material Gram and ideal closed-loop material map','97:112'),('main G19 material','real','576','full material variables; Green strong rank19 denominator','R11 METHODS:19')]
save('SPACE_LABELS.csv',[dict(object=a,scalar=b,shape=c,label=d,source_lines=e) for a,b,c,d,e in spaces])
sources=[('R14/src/local_closed_loop.py','local_closed_loop.py','79:114','Independent current basis; current normal and material maps; rho vs common-B00 norm'),('R14/src/local_closed_loop.py','local_closed_loop.py','142:176','Current error transport, material residual, exact residual closure and finite additivity'),('R13/src/full_current.py','full_current.py','17;137:146','Imports prepare_full and uses it for actual g19 readout'),('R11/src/material_readouts.py','material_readouts.py','97:101','prepare_full calls prepare_green_strong; alias full_material_readout'),('R10/src/lowrank_denominator.py','lowrank_denominator.py','16:34','G/S right singular vectors, threshold .1; green_strong=Vg'),('R10/src/lowrank_denominator.py','lowrank_denominator.py','37:79','Galerkin K=P_G G P_G action; full denominator material gradient'),('R10/src/lowrank_denominator.py','lowrank_denominator.py','86:128','Full-dimensional bounded LBFGS maxiter8, seed/fallback; selected green_strong'),('R10/METHODS_AND_ABLATIONS.md','METHODS_AND_ABLATIONS.md','12:15;31','Saved naming provenance: Green19 and19/30 denominator ranks, not material dimension'),('R11/METHODS_AND_ABLATIONS.md','METHODS_AND_ABLATIONS.md','3:19','Experimental name inherited R10, full derivative vs full Green, material576')]
sm=[]
for label,fn,lines,claim in sources:
 cohort=label.split('/')[0][1:];path=W/f'research/research_max_struct_r{cohort}'/('/'.join(label.split('/')[1:]));reads.append(str(path));assert path.exists();sm.append(dict(path=str(path),lines=lines,evidence=claim))
save('SOURCE_MAP.csv',sm)
assert len(filtered)==432 and len(new21)==33 and len(old21)==18
receipt=dict(status='COMPLETE',scope='saved scalar CSV/JSON and source only; no GT arrays, no scorer, no physics, no remote',rows=dict(local_directional=36,local_additivity=12,new_six_three_views=len(filtered),new0021_all11=len(new21),old0021_six=len(old21)),read_paths=sorted(set(reads)),CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,prior_interactive_inspection_CPU_seconds=None,prior_interactive_inspection_CPU_null_reason='not instrumented; not claimed zero',historical_compute='reused; candidate timestamps and wholechild costs separate, never summed')
(O/'EXTRACTION_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt['rows']));print('CPU',receipt['CPU_seconds'],'wall',receipt['wall_seconds'])

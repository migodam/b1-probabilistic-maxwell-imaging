"""Read-only finite evidence inventory. No numpy/scientific module imports."""
import ast,csv,json,struct,time,zipfile
from pathlib import Path
from collections import Counter,defaultdict
from datetime import datetime,timezone
W=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
R=W/'research/research_max_struct_r14';S=R/'supplement_observable19'
OLD=W.parent/'r13_worktree/research/research_max_struct_r13'
cpu,wall=time.process_time(),time.perf_counter()
files=[];errors=[]
def csvrows(p):
 with p.open(newline='',encoding='utf8') as f:return list(csv.DictReader(f))
def dump(name,v):(OUT/name).write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n')
def table(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (OUT/name).open('w',newline='',encoding='utf8') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
def inspect(p,role,cohort='',case='',method='',view='',headers=True):
 row=dict(path=str(p),role=role,cohort=cohort,case_id=case,method=method,view=view,exists=p.is_file(),bytes=p.stat().st_size if p.is_file() else None)
 if p.is_file() and p.suffix=='.npz' and headers:
  try:
   with zipfile.ZipFile(p) as z:
    row['keys']=[x[:-4] if x.endswith('.npy') else x for x in z.namelist()]
    if role=='public input' or view in ('primary','FINAL_AUXILIARY_STATE'):
     meta={}
     for name in z.namelist():
      if not name.endswith('.npy'):continue
      with z.open(name) as f:
       if f.read(6)!=b'\x93NUMPY':raise ValueError('Invalid NPY signature')
       major,minor=f.read(2);n=struct.unpack('<H' if major==1 else '<I',f.read(2 if major==1 else 4))[0]
       info=ast.literal_eval(f.read(n).decode('latin1' if major<3 else 'utf8'))
       meta[name[:-4]]=dict(shape=info['shape'],dtype=info['descr'])
     row['NPY_headers_only']=meta
  except Exception as exc:row['header_error']=repr(exc);errors.append(dict(path=str(p),error=repr(exc)))
 files.append(row);return row
# Supplement finite cohort, existing final scalar tables.
allrows=csvrows(S/'analysis/observable19_v1/ALL_VIEWS.csv')
summary=csvrows(S/'analysis/observable19_v1/SUMMARY.csv')
selected=[r for r in allrows if r['view'] in ('primary','raw','corrected')]
compact=[{k:r[k] for k in ['case_id','solver','view','source_file','available','NRMSE','IoU','Phi','relative_residual','noise_feasible','noise_certificate_eligible','CPU_seconds','run_status','archived_selector_consistent']} for r in selected]
table('OBSERVABLE19_PRIMARY_RAW_CORRECTED.csv',compact)
table('OBSERVABLE19_SUMMARY_VIEWS.csv',[r for r in summary if r['view'] in ('primary','raw','corrected')])
config=json.loads((S/'configs/observable19_exposed24_v1.json').read_text())
cases=sorted({r['case_id'] for r in allrows});methods=sorted({r['solver'] for r in allrows})
truth_by={r['case_id']:r for r in csvrows(S/'analysis/observable19_v1/TRUTH_SOURCES.csv')}
cohort=[]
for c in cases:
 task=next(t for t in config['tasks'] if t['case_id']==c)
 p=(S/task['input_file']).resolve();inspect(p,'public input','observable19 exposed24',c)
 truth=R/'sealed_truth/independent24'/(c+'.npz')
 cohort.append(dict(case_id=c,exposure_status=config['exposure_status'],public_input=str(p),public_input_exists=p.is_file(),
  truth_recorded_scoring_path=truth_by[c]['path'],truth_key_recorded=truth_by[c]['key'],local_original_generation_truth_path=str(truth),local_truth_exists=truth.is_file(),truth_arrays_not_read=True,
  inverse_side=task['observation_contract']['inverse_side'],forward_side=task['observation_contract']['forward_side'],
  model_mismatch_status=task['observation_contract']['model_mismatch_status'],supplement_local_RESULTS_exists=(S/'RESULTS/observable19_exposed24_v1'/task['task_id']).exists()))
table('COHORT_PUBLIC_TRUTH_SEPARATION.csv',cohort)
# Every local historical new24 primary in four relevant existing baseline families.
for c in cases:
 for method in ['dr','global_dual','g19','direct']:
  task=R/'RESULTS/independent24_v1'/f'{c}__{method}'
  for v in ['primary','raw','corrected','FINAL_AUXILIARY_STATE']:
   inspect(task/(v+'.npz'),'R14 historical endpoint','independent24 originally, now exposed',c,method,v)
# Explicit old exposed48 campaign only: do not roam other campaigns/cohorts.
oldroot=OLD/'RESULTS/confirmation_core48'
if oldroot.exists():
 for task in sorted(oldroot.iterdir()):
  if not task.is_dir():continue
  case,_,method=task.name.partition('__')
  if method not in ['dr','g19_full_current','g19_full_apcs']:continue
  for view in ['primary','raw','corrected','FINAL_AUXILIARY_STATE']:
   inspect(task/(view+'.npz'),'R13 historical endpoint','r13 exposed48 (different from r14 new24)',case,method,view)
# Old exposed48 public/truth separation comes from the actual old scoring receipt.
oldcohort=[]
for truth in csvrows(OLD/'analysis/confirmation_core48/TRUTH_SOURCES.csv'):
 c=truth['case_id']; receipt=json.loads((oldroot/f'{c}__dr/RECEIPT.json').read_text())
 t=receipt['task']; public=OLD/t['input_file'];inspect(public,'public input','r13 exposed48',c)
 contract=t['observation_contract']; oldcohort.append(dict(case_id=c,public_input=str(public),public_input_exists=public.is_file(),
   recorded_truth_path=truth['path'],recorded_truth_path_exists=Path(truth['path']).is_file(),truth_arrays_not_read=True,
   inverse_side=contract['inverse_side'],sources=contract['sources'],receivers=contract['receivers'],receiver_aperture=contract['receiver_aperture'],
   mismatch=contract['model_mismatch_status'],qualification=truth['qualification']))
table('OLD48_PUBLIC_TRUTH_SEPARATION.csv',oldcohort)
# Compact checkpoint subset already specifically curated, no raw scans.
for c in ['0012','0021']:
 for pair in ['01','11']:
  for k in [128,512,1500]:inspect(R/f'RESULTS/independent24_v1/r14_new_{c}__pair{pair}/CHECKPOINT_{k:05}.npz','curated complete checkpoint','r14 new24 exposed',f'r14_new_{c}','pair'+pair,f'call{k}')
# Already saved GN-Green operator spectra/derived errors; headers only, never array bodies.
for name in ['fixed_spectra.npz','error_arrays.npz','cumulative_mode_energies.npz']:
 inspect(W/'research/research_gn_green_error_20261106/results'/name,'saved spectrum or truth-derived diagnostic; not prospective public input')
# Remote fullstate inventory: metadata only. Do not fetch any release asset.
planpath=S/'archives/raw_shards_v1/PACKAGE_PLAN.json';plan=json.loads(planpath.read_text())
manifest=[]
for package in plan['packages']:
 for item in package.get('files',[]):
  p=item['path']
  if p.endswith('.npz') and any('/'+v+'.npz' in p for v in ['primary','raw','corrected','FINAL_AUXILIARY_STATE']):
   manifest.append(dict(path=p,bytes=item['bytes'],local_exists=(S/p).is_file(),evidence='remote archive PACKAGE_PLAN only; NPZ fields not examined'))
table('SUPPLEMENT_ARCHIVED_ENDPOINTS.csv',manifest)
# Source evidence definitions: transcription/references, not independent proof review.
evidence=[]
def ev(topic,path,line,claim,limit):evidence.append(dict(topic=topic,path=str(W/path),line=line,evidence=claim,scope_limit=limit))
ev('GN projection contract','research/research_gn_nn_preservation_r4/MATHEMATICAL_THEORY.md',15,'VG: G relative singular value >=.1, rank19. VS: all48 nonzero S directions. W=orth[VG,VS]. PW=I-VGVG*; P0=I-VSVS*; PD=I-WW*. Ranks557/528/509, bright union67.','Old fixed geometry only; these auxiliary current spaces are not material-Jacobian nullspaces. S visible48 is not S strong9.')
ev('Absolute energy convention','research/research_gn_nn_preservation_r4/MATHEMATICAL_THEORY.md',5,'Projected squared Frobenius error changes use 2Re<Pe,Pu>+||Pu||²; absolute error changes are primary.','Not electromagnetic stored energy/Poynting; updates in strong directions do not alone imply lost accuracy.')
ev('Implementation/two orders','research/research_gn_nn_preservation_r4/src/energy.py',24,'Frozen tau=.1, numerical cutoff1e-12; numerical intersection tolerance1e-9; I_to_O and O_to_I retained separately, plus intersection remainder.','Do not call sequential products orthogonal intersection projectors.')
ev('Green strong weak complement','research/research_gn_green_error_20261106/MATHEMATICAL_DEFINITIONS.md',28,'strong sigma/sigmamax>=.1; weak 1e-12<relative sigma<.1; complement remainder including structural right nullspace. Old G:19/557/0; S:9/39/528.','Ranks belong to original G576x576/S48x576 geometry; not inferred for new24.')
ev('Sequential energy','research/research_gn_green_error_20261106/MATHEMATICAL_DEFINITIONS.md',47,'Xab=Pob,b Pin,a X; reversed order separate; sum squared norms=total though components need not be mutually orthogonal.','Not strict intersections; propagated components require signed cross terms.')
ev('Numerical intersections','research/research_gn_green_error_20261106/MATHEMATICAL_DEFINITIONS.md',65,'Fixed sin-angle1e-9; weak/weak20 and weak/complement509; remaining47 must be preserved.','Numerical intersection, not symbolic certification; sensitivity1e-10/1e-8 saved.')
ev('Projected spans and actual chart','research/research_gn_green_error_20261106/MATHEMATICAL_DEFINITIONS.md',80,'Zab=Pob,b Vin,a, then orthogonalized span. Actual chart=[Vo,strong,orth((I-Po,strong)Vi,strong)].','Projected spans can overlap and energies exceed100%; projecting separately into weak and null with independent coefficients changes chart.')
ev('Literal chart implementation','research/research_gn_green_error_20261106/src/tsom_chart.py',15,'Receiver and internal strong threshold .1; residual projected basis absolute rank cutoff1e-10; orthogonal receiver9+secondfold19 chart.','Shared internal coefficients in weak+complement projection, not separately parameterized union of those two projected spans.')
ev('R10 union not intersection','research/research_max_struct_r10/src/lowrank_denominator.py',16,'Vg,Vo strong threshold.1. Residual=(I-VoVo*)Vg. extra SVD columns>1e-10. U=[Vo,extra].','This chart spans sum/union of receiver/internal strong spaces, not their intersection. Code only inspected, not executed.')
ev('R10 geometry-specific ranks','research/research_max_struct_r10/THEORY_DERIVATIONS.md',39,'Historical geometry G19/S9/twofold28; fresh24 G19/S11/twofold30; early28 fresh label was corrected without changing basis.','R10 fresh24 are r10_twod IDs, not R14 new24. Never carry ranks across cohorts without actual operator evidence.')
ev('R10 distinct diagnostics','research/research_max_struct_r10/THEORY_DERIVATIONS.md',102,'Sequential, numerical intersection, projected span and actual FFT/APCS admissible space are recorded separately; absolute energies required.','No equivalence or material-identifiability conclusion inferred.')
ev('R14 supplement fixed19','research/research_max_struct_r14/supplement_observable19/configs/observable19_exposed24_v1.json',1,str(config['basis_contract']),'S/O0 first19 is not the historical all48 visible space or .1-threshold S-strong space; O0 is current-to-port transport, not material Jacobian.')
ev('C2 physical endpoint schema','research/research_max_struct_r14/runtime/equal_budget_C2_runner.py',268,'Saves endpoint m,j,Q,Y_prediction; FINAL_AUXILIARY_STATE saved separately only if final_states exists.','Local selected-m exports may omit j/Q even though original runner wrote them; endpoint physical fields and independent DR iterator are different objects.')
ev('Full-current existing source','research/research_max_struct_r13/src/full_current.py',198,'mode full_current uses empty actual_vd (no DPCS projection); restricted/full_apcs keep original Vd.','Source-level description only; no method executed. Full-current freedom is not supplied by m-only export.')
table('SOURCE_EVIDENCE.csv',evidence);table('FILES_AND_NPZ_HEADERS.csv',files)
counts=[]
for key,group in __import__('itertools').groupby(sorted(files,key=lambda r:(r['role'],r['cohort'],r['method'],r['view'])),key=lambda r:(r['role'],r['cohort'],r['method'],r['view'])):
 group=list(group);counts.append(dict(role=key[0],cohort=key[1],method=key[2],view=key[3],expected=len(group),exists=sum(x['exists'] for x in group),key_patterns=dict(Counter(json.dumps(x.get('keys',[])) for x in group if x['exists']))))
table('AVAILABILITY_COUNTS.csv',counts)
result=dict(status='COMPLETE_MECHANICAL_INVENTORY',cases=cases,case_count=len(cases),arms=methods,selected_view_rows=len(selected),selected_views=dict(Counter(r['view'] for r in selected)),view_counts=dict(Counter(r['view'] for r in allrows)),
 local_files_count=len(files),local_existing_count=sum(r['exists'] for r in files),NPZ_header_errors=errors,
 remote_archived_endpoint_entries=len(manifest),remote_endpoint_local_present=sum(r['local_exists'] for r in manifest),basis_contract=config['basis_contract'],
 source_root=str(W),old_R13_root=str(OLD),truth_arrays_read=False,physics_runs=0,scorer_runs=0,new_SVD=0,remote_access=False,hash_operations=False,
 CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,created_utc=datetime.now(timezone.utc).isoformat())
dump('RECEIPT.json',result)
print(json.dumps({k:result[k] for k in ['case_count','selected_view_rows','local_files_count','local_existing_count','remote_archived_endpoint_entries','remote_endpoint_local_present','CPU_seconds','wall_seconds']},indent=2))

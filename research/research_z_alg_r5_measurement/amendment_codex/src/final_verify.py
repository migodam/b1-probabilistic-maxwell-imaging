from complete_quotient import *
import ast,subprocess,csv

def main():
 start=time.process_time();checks={};errors=[]
 for fn in ['solver_freeze.json','prediction_freeze.json']:
  o=json.loads((BASE/'inputs'/fn).read_text());bad=[p for p,h in o['files'].items() if not (BASE/p).exists() or digest(BASE/p)!=h];checks[fn+'_unchanged']=not bad;errors+=bad
  if 'historical_files' in o:
   bad=[p for p,h in o['historical_files'].items() if not (BASE.parent/p).exists() or digest(BASE.parent/p)!=h];checks['historical_imports_unchanged']=not bad;errors+=bad
 for path in list((BASE/'src').glob('*.py'))+list((OUT/'src').glob('*.py')):ast.parse(path.read_text(),filename=str(path))
 checks['all_research_python_syntax']=True
 source_checks=[];start_counts=[]
 for case in [f'T{i}' for i in range(6)]:
  g=geometry(case);pos,Q=dictionary(case,g)
  for arm in ['clean','noisy']:
   path=BASE/'results'/case/arm;d=json.loads((path/'design_freeze.json').read_text());source_checks.append(len(d['menu'])<=48)
   for row in d['menu']:
    a=np.array([complex(v['real'],v['imag']) for v in row['a']]);source_checks.append(bool(abs(np.vdot(a,Q@a)-1)<1e-10 and np.linalg.norm(a)<=4+1e-10))
   for alg in ['TRF64','reciprocal64']:
    o=json.loads((path/(alg+'.json')).read_text());start_counts.append(dict(case=case,arm=arm,algorithm=alg,completed=o['starts_completed'],requested=o['starts_requested']));assert o['starts_requested']==64 and o['starts_completed']==64
 checks['all_24_fixed_64_start_arms_completed']=True;checks['original_sources_meet_numeric_constraints']=all(source_checks)
 for fn in ['constrained_source_design_freeze.json']:
  d=json.loads((OUT/'results'/fn).read_text());g=geometry('T0');pos,Q=dictionary('T0',g)
  checks['amendment_sources_meet_numeric_constraints']=all(abs(v['power']-1)<1e-10 and v['feed_norm']<=4+1e-10 for v in d['designs'].values())
 for fn in ['constrained_source_predictions_frozen.json']:
  d=json.loads((OUT/'results'/fn).read_text());checks[fn+'_unchanged']=all(digest(OUT/p)==h for p,h in d['files'].items())
 checks['all_12_central_material_box_certificates']=all(r['material_box_complete'] for r in json.loads((OUT/'results/complete_batch.json').read_text())['rows'])
 checks['saved_184_complex_boxes_real_form_recertified']=json.loads((OUT/'results/saved_certificate_audit.json').read_text())['passed']
 checks['original_formula_tests']=json.loads((BASE/'results/formula_tests.json').read_text())['passed'];checks['amendment_formula_tests']=json.loads((OUT/'results/amendment_formula_tests.json').read_text())['passed']
 checks['critical_point_interval_IVT']=json.loads((OUT/'results/certified_fold.json').read_text())['certified_critical_value_inside_old_noise_ball']
 git=subprocess.run(['git','diff','--name-only'],cwd=BASE.parents[1],text=True,capture_output=True,check=True).stdout.splitlines();checks['only_existing_repository_file_modified_is_authorized_AGENTS']=git==['AGENTS.md']
 versions=dict(python=sys.version,numpy=np.__version__,scipy=__import__('scipy').__version__,mpmath=mp.__version__,matplotlib=__import__('matplotlib').__version__)
 save(BASE/'results/FINAL_VERIFICATION.json',dict(passed=all(checks.values()),checks=checks,hash_errors=errors,starts=start_counts,versions=versions,CPU_seconds=time.process_time()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,visual_review='Original five figures inspected; figure02 legend corrected. Amendments06/07 inspected;07 annotation overlap corrected and re-inspected.',scope='Internal scientific/code review, not independent external/formal verification. Reproduction wrappers syntax/path inspected, whole batch not rerun just for wrapper.'))
 print(checks);assert all(checks.values())
if __name__=='__main__':main()

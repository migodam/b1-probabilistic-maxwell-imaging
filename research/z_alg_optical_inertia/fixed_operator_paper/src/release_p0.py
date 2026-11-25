"""Aggregate P0 from independently rebuilt mechanism checks and a native 128-bit replay."""
from pathlib import Path
import json
p=Path(__file__).resolve().parents[1]/'results'
def read(path):return json.loads((p/path).read_text())
a=read('mechanisms/fixed_space_checks.json'); b=read('mechanisms/active_tm_results.json'); c=read('mechanisms/active_clean_results.json')
d=read('mechanisms/active_ambiguity_closed_loop.json');e=read('mechanisms/background_compensation_exact.json'); n=read('native_replay.json'); f=read('p0_core.json')
passed=(len(read('mechanisms/runs.json'))==6 and all(not r['returncode'] for r in read('mechanisms/runs.json')) and a['exact_inertia_cases_passed']==45 and b['future_ellipse_inertia_mismatches']==b['no_information_mismatches']==0 and c['nu']==c['kold'] and c['true_response_quad']<c['beta_clean_no_added_t'] and d['certified_noisy_separation_lower']>0 and d['exact_ambiguous_alternative_inside_second_box'] and e['proofs'][0]['polynomial_sign_change'] and f['status']=='PASS' and len(n)>=1 and all(x['witness']['fits_noise_ball'] and x['cap']['strictly_excludes'] for x in n))
if not passed:raise RuntimeError('P0 discrepancy')
doc=dict(status='PASS',passed=True,proof='Arb128 P0 full/reduced and singular, cover and retention controls plus native disk replay',positive_native_replays=[x['id'] for x in n],core='p0_core.json',independent_sources=['mechanisms/runs.json','native_replay.json'])
tmp=p/'p0.tmp';tmp.write_text(json.dumps(doc,indent=2));tmp.replace(p/'p0.json')

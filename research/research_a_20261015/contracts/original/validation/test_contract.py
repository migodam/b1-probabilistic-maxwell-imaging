"""Synthetic schema/semantic tests only: no fixture represents real data."""
from pathlib import Path
import copy,json,sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from validate_package import validate
schema=json.loads((ROOT/'MaterialHypothesisPackage.schema.json').read_text())

def sample(s):
    if '$ref' in s: return sample(schema['$defs'][s['$ref'].split('/')[-1]])
    if 'const' in s: return s['const']
    if 'enum' in s: return s['enum'][0]
    if 'anyOf' in s: return None if any(v.get('type')=='null' for v in s['anyOf']) else sample(s['anyOf'][0])
    typ=s.get('type')
    if typ=='object': return {k:sample(v) for k,v in s.get('properties',{}).items()}
    if typ=='array': return []
    if typ=='string': return '0'*64 if 'sha256' in str(s) or s.get('pattern')=='^[a-f0-9]{64}$' else 'synthetic-unit-fixture'
    if typ in ('integer','number'): return max(s.get('minimum',0),s.get('exclusiveMinimum',-1)+1)
    if typ=='boolean': return False
    if typ=='null': return None
    raise ValueError(str(s))

fixture=sample(schema)
fixture.update(status='handoff_ready',requested_k=1,valid_k=1,weight_reliability='candidate_only',weights=None)
fixture['candidates']=[sample(schema['$defs']['candidate'])]
fixture['candidates'][0]['posterior_weight']=None
fixture['next_measurement_interface']['posterior_information_gain_allowed']=False
results={'scope':'Synthetic format tests, not a validated experimental hypothesis package','positive':validate(fixture,ROOT),'negative':{}}
mutations={
 'candidate_only_fake_probability':lambda d:d['candidates'][0].update(posterior_weight=.5),
 'wrong_target_identity':lambda d:d['candidates'][0]['validation'].update(physics_target_hash='wrong'),
 'invalid_residual':lambda d:d['candidates'][0]['fit'].update(whitened_residual_l2=2.,squared_residual=4.),
 'wrong_k':lambda d:d.update(valid_k=2),
 'truth_leak':lambda d:d['candidates'][0]['compensation']['edit_parameters'].update(truth_mae=.1),
 'fake_bayesian_handoff':lambda d:d['next_measurement_interface'].update(posterior_information_gain_allowed=True)
}
for name,mutation in mutations.items():
    d=copy.deepcopy(fixture); mutation(d)
    try:
        validate(d,ROOT)
    except Exception:
        results['negative'][name]='REJECTED_AS_EXPECTED'
    else: raise AssertionError('Bad fixture accepted: '+name)
(ROOT/'validation'/'contract_test_results.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))

#!/usr/bin/env python3
"""Validate schema, cross-field science semantics, and optional asset hashes.
Requires jsonschema. Example:
python validate_package.py path/to/manifest.json --check-assets
This is a format/consistency validator, not a Maxwell solver or a certificate.
"""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path
import jsonschema


def walk(value):
    if isinstance(value,dict):
        yield value
        for sub in value.values(): yield from walk(sub)
    elif isinstance(value,list):
        for sub in value: yield from walk(sub)


def validate(data, root, check_assets=False):
    schema=json.loads((Path(__file__).parent/'MaterialHypothesisPackage.schema.json').read_text())
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.Draft202012Validator(schema).validate(data)
    def require(ok,msg):
        if not ok: raise ValueError(msg)
    candidates=data['candidates']; ids=[c['id'] for c in candidates]
    require(len(ids)==len(set(ids)), 'Candidate ids must be unique.')
    require(data['valid_k']==len(candidates), 'valid_k differs from candidate count.')
    require(data['valid_k']<=data['requested_k'], 'Export is larger than requested K.')
    for node in walk(data):
        require(not (set(node)&{'truth','ground_truth','truth_material','truth_mae','best_of_k_truth'}), 'Truth-derived fields must not enter the online handoff.')
        for value in node.values():
            if isinstance(value,float): require(math.isfinite(value), 'Non-finite floating point value.')
    for c in candidates:
        f=c['fit']
        require(f['whitened_residual_l2']<=f['threshold_l2'], 'A retained candidate fails its data-fit gate.')
        require(math.isclose(f['squared_residual'], f['whitened_residual_l2']**2,rel_tol=1e-7,abs_tol=1e-12), 'Residual norm and squared residual disagree.')
        require(c['validation']['physics_target_hash']==data['physics_target']['hash'], 'Candidate target mismatch.')
        require(c['validation']['material_law_id']==data['physics_target']['material_law_id'], 'Candidate material-law mismatch.')
        require(f['noise_model_id']==data['noise']['id'], 'Candidate noise-model mismatch.')
        require(f['nuisance_policy']!='bounded_profile' or f['likelihood_semantics']!='proper_marginal', 'Profile likelihood is not a proper marginal likelihood.')
    seen=set()
    for pair in data['pairs']:
        a,b=pair['candidate_i'],pair['candidate_j']; key=tuple(sorted((a,b)))
        require(a in ids and b in ids and a!=b, 'Pair refers to invalid candidates.')
        require(key not in seen, 'Duplicate pair.'); seen.add(key)
        require(pair['connection_status']!='verified_path' or pair['validated_path_asset'] is not None,'A verified connection needs a saved validated path.')
    weights=data['weights']
    if weights is not None:
        require(len(weights['values'])==len(ids), 'Wrong number of weights.')
        require(all(c['posterior_weight'] is not None for c in candidates),'Weighted package needs candidate weights.')
        require(all(math.isclose(c['posterior_weight'],w,rel_tol=1e-8,abs_tol=1e-12) for c,w in zip(candidates,weights['values'])),'Candidate/global weight mismatch.')
        mass=sum(weights['values'])
        if weights['normalization_domain']=='candidate_pool_only':
            require(math.isclose(mass,1.,rel_tol=1e-7,abs_tol=1e-9),'Pool-normalized weights must sum to one.')
        else:
            require(mass<=1.+1e-8,'Region/sample mass exceeds one.')
        require(data['weight_reliability'] in ('reference','corrected'), 'Posterior weights require declared reliable semantics.')
    if check_assets:
        root=root.resolve()
        for node in walk(data):
            if {'path','sha256','format','shape','dtype'}<=set(node):
                path=(root/node['path']).resolve()
                require(path.is_relative_to(root), 'Asset path escapes package directory.')
                require(path.is_file(), 'Missing file asset: '+str(path))
                h=hashlib.sha256()
                with path.open('rb') as f:
                    for block in iter(lambda:f.read(1048576),b''): h.update(block)
                require(h.hexdigest()==node['sha256'], 'Asset hash mismatch: '+str(path))
    return {'schema':'PASS','cross_field_checks':'PASS','asset_hashes':'PASS' if check_assets else 'NOT_CHECKED','science_scope':'No Maxwell equations or continuous-model error certified by this validator.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('manifest',type=Path); parser.add_argument('--check-assets',action='store_true')
    args=parser.parse_args()
    try:
        data=json.loads(args.manifest.read_text(encoding='utf-8'))
        print(json.dumps(validate(data,args.manifest.parent,args.check_assets),indent=2))
    except (OSError,ValueError,jsonschema.ValidationError) as exc:
        parser.exit(1,'VALIDATION FAILED: '+str(exc)+'\n')

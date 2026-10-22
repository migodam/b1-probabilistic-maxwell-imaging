"""Separate immutable-metadata contract for the second256-parent mechanism study.

No GT reader is defined here. File identity checks use declared size/mtime and
explicit protocol/source versions, never hashes or the original256 contract.
"""
from pathlib import Path
import json,re
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
VERSION=1
KIND='p3r2_mechanism_confirmation_v1'
FAMILIES=('full_geometry','local_anchor','capacity','permuted','sensitivity_only','backprojection_only','decoder_random_sensitivity')
SEEDS=(11,22,33)
PARENTS=[f'R2_locked_{i:05d}' for i in range(256,512)]
PRIMARY=(('full_geometry','backprojection_only'),('full_geometry','sensitivity_only'),('sensitivity_only','decoder_random_sensitivity'))
PRIMARY_CONFIDENCE=1-.05/3
MODEL_COUNT=len(FAMILIES)*len(SEEDS)
KINDS=dict(full_geometry='local_geometry_cnn',local_anchor='local_anchor_cnn',capacity='local_capacity_cnn',permuted='local_permuted_geometry_cnn',sensitivity_only='mechanism_sensitivity_only',backprojection_only='mechanism_backprojection_only',decoder_random_sensitivity='mechanism_decoder_random_sensitivity')
REQUIRED_SOURCE_FILES=('src/p3r2_mechanism/__init__.py','src/p3r2_mechanism/contract.py','src/p3r2_mechanism/masked_risk.py','src/p3r2_mechanism/random_sensitivity.py',
 'src/p3r2/residual_risk.py','src/p3r2/risk_models.py','src/p3r2/risk_training.py','src/p1_material/prior.py',
 'tools/train_p3r2_mechanism_controls.py','tools/train_p3r2_random_sensitivity.py','tools/infer_p3r2_mechanism_confirmation.py','tools/evaluate_p3r2_mechanism_confirmation.py',
 'tools/run_p3r2_anchors.py','tools/generate_p3r2_data.py')
ONLINE_FIELDS=('observation','residual','prediction','q0','chi0','J')


def resolve(value,out):
    p=Path(value)
    if not p.is_absolute(): p=(ROOT if p.parts and p.parts[0] in ('src','tools') else out)/p
    return p.resolve()


def file_record(path):
    path=Path(path);s=path.stat();return dict(path=str(path.resolve()),size_bytes=s.st_size,mtime_ns=s.st_mtime_ns)


def check_record(record,out):
    path=resolve(record['path'],out);s=path.stat()
    if s.st_size!=int(record['size_bytes']) or s.st_mtime_ns!=int(record['mtime_ns']): raise ValueError('frozen file metadata changed: '+str(path))
    return path


def validate_freeze(out,manifest_path):
    out=Path(out).resolve();manifest_path=Path(manifest_path).resolve();m=json.loads(manifest_path.read_text())
    if m.get('version')!=VERSION or m.get('protocol_kind')!=KIND or m.get('status')!='FROZEN': raise ValueError('new explicit mechanism freeze version/status required')
    if not isinstance(m.get('source_revision'),str) or not m['source_revision']: raise ValueError('source revision required')
    if m.get('estimator_action','raw_unprojected_no_blend')!='raw_unprojected_no_blend': raise ValueError('new confirmation requires raw unprojected CNN actions without blend')
    if m.get('primary_metric','complex_fullgrid_MSE')!='complex_fullgrid_MSE': raise ValueError('primary mechanism risk is complex full-grid MSE')
    if m.get('expected_count')!=256 or m.get('parents')!=PARENTS: raise ValueError('exact new cohort R2_locked_00256..00511 required')
    frozen=float(m['frozen_unix'])
    if frozen<=0 or manifest_path.stat().st_mtime<frozen: raise ValueError('invalid method freeze time')
    sources=m.get('frozen_source_files')
    if not isinstance(sources,list) or not set(REQUIRED_SOURCE_FILES).issubset(sources): raise ValueError('new protocol source whitelist incomplete')
    if any(not isinstance(s,str) or not (s.startswith('src/') or s.startswith('tools/')) or '..' in Path(s).parts for s in sources): raise ValueError('invalid source path')
    models=m.get('models',[])
    if len(models)!=MODEL_COUNT or {(x.get('family'),x.get('seed')) for x in models}!={(f,s) for f in FAMILIES for s in SEEDS}: raise ValueError('all seven families × three frozen seeds required')
    names=[x.get('name','') for x in models]
    if len(set(names))!=MODEL_COUNT or any(not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,79}',n) for n in names): raise ValueError('unique safe method names required')
    for entry in models:
        expected='random_sensitivity' if entry['family']=='decoder_random_sensitivity' else 'mechanism' if entry['family'] in ('sensitivity_only','backprojection_only') else 'original'
        if entry.get('loader')!=expected: raise ValueError('explicit family-aware model loader mismatch')
    if m.get('primary_comparisons')!=[list(p) for p in PRIMARY]: raise ValueError('three preregistered primary comparisons cannot change')
    if m.get('bootstrap_replicates',2000)!=2000 or m.get('bootstrap_seed',2026102393)!=2026102393 or m.get('primary_CI')!=PRIMARY_CONFIDENCE: raise ValueError('fixed parent bootstrap/CI contract')
    required=[resolve(s,out) for s in sources]+[resolve(x['checkpoint'],out) for x in models]
    target=resolve(m.get('target_config','configs/TARGET.json'),out);geometry=resolve(m.get('geometry_config',ROOT/'research/p1_material_target_20261018/inputs/geometry.json'),out)
    required += [target,geometry]
    records=m.get('frozen_files',[])
    if not isinstance(records,list): raise ValueError('explicit frozen metadata records required')
    registered={resolve(r['path'],out):r for r in records}
    if len(registered)!=len(records): raise ValueError('duplicate frozen file metadata')
    for path in required:
        if path not in registered: raise ValueError('missing frozen file registration: '+str(path))
    for record in records:
        path=check_record(record,out)
        if path.stat().st_mtime>frozen: raise ValueError('frozen artifact/source postdates freeze: '+str(path))
    return m,dict(target=target,geometry=geometry,models={x['name']:resolve(x['checkpoint'],out) for x in models}),file_record(manifest_path)


def validate_online_inputs(out,m):
    frozen=float(m['frozen_unix'])
    for parent in PARENTS:
        for path in (out/'inputs/online'/f'{parent}.npz',out/'L/contexts'/f'{parent}.npz',out/'L/contexts'/f'{parent}.json'):
            if not path.is_file() or path.stat().st_mtime<frozen: raise ValueError('new post-freeze context/observation required: '+str(path))
        side=json.loads((out/'L/contexts'/f'{parent}.json').read_text())
        if side.get('parent')!=parent or side.get('native_cap')!=96 or side.get('probit_cap')!=64 or side.get('GT_read') is not False: raise ValueError('truth-free96+64 context required: '+parent)


def read_online(path):
    with np.load(path,allow_pickle=False) as z: d={k:z[k].copy() for k in ONLINE_FIELDS}
    shapes=dict(observation=(32,8),residual=(32,8),prediction=(32,8),q0=(66,),chi0=(576,),J=(512,66))
    for key,shape in shapes.items():
        if d[key].shape!=shape or not np.isfinite(d[key]).all(): raise ValueError('invalid online context '+key)
    return {k:v[None] for k,v in d.items()}


def completed_jobs(out,m):
    generator=resolve(m.get('generator_closed_json','L/jobs/data_locked_256_256/CLOSED.json'),out)
    data_manifest=generator.parent/'DATA_MANIFEST.json'
    if json.loads(generator.read_text()).get('status')!='COMPLETE' or generator.stat().st_mtime<float(m['frozen_unix']) or data_manifest.stat().st_mtime<float(m['frozen_unix']): raise ValueError('post-freeze COMPLETE data subjob and data manifest required')
    data=json.loads(data_manifest.read_text())
    if data.get('group')!='locked' or data.get('parents')!=PARENTS: raise ValueError('new generated cohort manifest mismatch')
    closures=m.get('context_closed_jsons',['L/jobs/anchors_locked_256_256_96_64/CLOSED.json'])
    if not closures: raise ValueError('new anchor subjob closures required')
    for value in closures:
        path=resolve(value,out)
        if json.loads(path.read_text()).get('status')!='COMPLETE' or path.stat().st_mtime<float(m['frozen_unix']): raise ValueError('post-freeze completed anchor subjob required')
    online=resolve(m['online_inference_closure'],out)
    if json.loads(online.read_text()).get('status')!='COMPLETE' or online.stat().st_mtime<float(m['frozen_unix']):
        raise ValueError('post-freeze online inference job must close before GT')
    return [generator]+[resolve(v,out) for v in closures]+[online]


def checkpoint_metadata(path,entry):
    import torch
    from .masked_risk import FORMAT,mask_contract
    p=torch.load(path,map_location='cpu',weights_only=False)
    if p.get('kind')!=KINDS[entry['family']] or int(p.get('seed',-1))!=entry['seed']: raise ValueError('checkpoint family/seed mismatch: '+entry['name'])
    expected=5266 if entry['family']=='local_anchor' else 5842
    count=sum(v.numel() for k,v in p['state_dict'].items() if k.startswith('local.') or k.startswith('head.'))
    if count!=expected: raise ValueError('frozen CNN parameter count mismatch')
    if entry['loader']=='mechanism':
        if p.get('format')!=FORMAT or p.get('mask')!=mask_contract(entry['family']) or p.get('objective')!='SUPERVISED_RISK': raise ValueError('mask-aware checkpoint version mismatch')
    elif entry['loader']=='random_sensitivity':
        from .random_sensitivity import validate_payload
        validate_payload(p)
    elif p.get('format')!='p3r2_anchor_material_risk_v1': raise ValueError('original checkpoint format mismatch')
    return dict(name=entry['name'],family=entry['family'],seed=entry['seed'],kind=p['kind'],format=p['format'],
        parameter_count=count,best_update=int(p['best_update']),mask=p.get('mask'),feature_contract=p.get('feature_contract'),objective='SUPERVISED_RISK',
        ROI_cells=len(p['geometry']['roi']),width=p['width'],dimensions=p['dimensions'],residual_scale=p['residual_scale'])

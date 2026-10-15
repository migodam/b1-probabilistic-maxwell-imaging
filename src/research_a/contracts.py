"""Zero-physics Research A native-profile export and versioned common exchange.

Never assigns candidate probabilities, computes digests, or invokes a solver.
The original A schema requires digests; the explicitly named execution profile
permits null digests in accordance with the human's no-hash policy.
"""
from __future__ import annotations

import copy
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

CONTRACT_DIR = Path(__file__).resolve().parents[2] / 'research/research_a_20261015/contracts'
PROFILE = 'research-a-native-1.0.0-execution-no-digest-v1'
EXCHANGE_VERSION = 'research-a-common-exchange/1.0.0'
TRUTH_KEYS = {'truth', 'ground_truth', 'truth_material', 'truth_mae', 'best_of_k_truth'}


def utcnow() -> str:
    return datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z')


def _clean(x: Any) -> Any:
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [_clean(v) for v in x]
    return x


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(value), indent=2, allow_nan=False) + '\n')


def _walk(value: Any):
    if isinstance(value, dict):
        yield value
        for sub in value.values():
            yield from _walk(sub)
    elif isinstance(value, list):
        for sub in value:
            yield from _walk(sub)


def _require(ok: bool, message: str) -> None:
    if not ok:
        raise ValueError(message)


def validate_schema_subset(value: Any, schema: dict, root: dict | None = None, path: str = '$') -> None:
    """Dependency-free validator for all keywords used by our saved schemas.

    This is deliberately scoped, rather than an implementation of arbitrary JSON
    Schema. It also rejects non-finite JSON numbers, which standard JSON excludes.
    """
    root = schema if root is None else root
    if '$ref' in schema:
        target = root
        for part in schema['$ref'][2:].split('/'):
            target = target[part]
        validate_schema_subset(value, target, root, path)
    if 'anyOf' in schema:
        errors = []
        for choice in schema['anyOf']:
            try:
                validate_schema_subset(value, choice, root, path)
                break
            except ValueError as exc:
                errors.append(str(exc))
        else:
            raise ValueError(f'{path}: no anyOf alternative: {errors}')
    for clause in schema.get('allOf', []):
        validate_schema_subset(value, clause, root, path)
    if 'if' in schema:
        try:
            validate_schema_subset(value, schema['if'], root, path)
        except ValueError:
            if 'else' in schema:
                validate_schema_subset(value, schema['else'], root, path)
        else:
            if 'then' in schema:
                validate_schema_subset(value, schema['then'], root, path)
    if 'const' in schema:
        _require(value == schema['const'] and (not isinstance(schema['const'], bool) or isinstance(value, bool)), f'{path}: const mismatch')
    if 'enum' in schema:
        _require(value in schema['enum'], f'{path}: invalid enum')
    typ = schema.get('type')
    types = {
        'object': isinstance(value, dict), 'array': isinstance(value, list),
        'string': isinstance(value, str), 'number': isinstance(value, (int, float)) and not isinstance(value, bool),
        'integer': isinstance(value, int) and not isinstance(value, bool),
        'boolean': isinstance(value, bool), 'null': value is None,
    }
    if typ:
        _require(types[typ], f'{path}: expected {typ}')
    if isinstance(value, (float, int)) and not isinstance(value, bool):
        _require(math.isfinite(value), f'{path}: non-finite number')
        if 'minimum' in schema:
            _require(value >= schema['minimum'], f'{path}: below minimum')
        if 'exclusiveMinimum' in schema:
            _require(value > schema['exclusiveMinimum'], f'{path}: below exclusive minimum')
    if isinstance(value, str):
        _require(len(value) >= schema.get('minLength', 0), f'{path}: short string')
        if 'pattern' in schema:
            _require(re.search(schema['pattern'], value) is not None, f'{path}: pattern mismatch')
    if isinstance(value, list):
        _require(len(value) >= schema.get('minItems', 0), f'{path}: insufficient items')
        for i, item in enumerate(value):
            if 'items' in schema:
                validate_schema_subset(item, schema['items'], root, f'{path}[{i}]')
    if isinstance(value, dict):
        _require(set(schema.get('required', [])) <= value.keys(), f'{path}: missing keys {set(schema.get("required", []))-value.keys()}')
        props = schema.get('properties', {})
        if schema.get('additionalProperties') is False:
            _require(value.keys() <= props.keys(), f'{path}: unexpected keys {value.keys()-props.keys()}')
        for key, sub in props.items():
            if key in value:
                validate_schema_subset(value[key], sub, root, f'{path}.{key}')


def _asset(root: Path, name: str, value: Any) -> dict:
    rel = f'assets/{name}'
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(value, (dict, list)):
        rel += '.json'
        _write_json(root / rel, value)
        fmt, shape, dtype = 'json', [], 'json'
    else:
        array = np.asarray(value)
        _require(array.dtype.kind != 'O', f'{name}: object arrays are forbidden')
        _require(np.isfinite(array).all(), f'{name}: non-finite array')
        rel += '.npy'
        np.save(root / rel, array, allow_pickle=False)
        fmt, shape, dtype = 'npy', list(array.shape), str(array.dtype)
    return {'path': rel, 'sha256': None, 'format': fmt, 'key': None, 'shape': shape, 'dtype': dtype}


def zero_cost(**updates) -> dict:
    """Only use zeros for work which really did not occur; callers supply receipts."""
    out = dict(new_full_state_frequency_attempts=0, vector_rhs=0, forward_rhs=0,
               adjoint_rhs=0, tangent_rhs=0, cache_hits=0, proposal_count=0,
               rejected_count=0, failed_count=0, wall_seconds=0.,
               peak_memory_bytes=None, setup_allocation_seconds=0., historical_cost_ledger=None)
    out.update(updates)
    return out


def error_record(value=None, status='unknown', method='not quantified', norm='relative_l2') -> dict:
    return dict(value=value, norm=norm, status=status, method=method)


def _cost(value: dict) -> dict:
    schema = json.loads((CONTRACT_DIR/'MaterialHypothesisPackage.execution.schema.json').read_text())
    out = zero_cost(**value)
    validate_schema_subset(out, schema['$defs']['cost'], schema)
    return out


def _asset_validate(asset: dict, root: Path) -> np.ndarray | Any:
    path = (root/asset['path']).resolve()
    _require(path.is_relative_to(root.resolve()), 'Asset path escapes package root')
    _require(path.is_file(), f'Missing asset {asset["path"]}')
    if asset['format'] == 'json':
        return json.loads(path.read_text())
    _require(asset['format'] == 'npy', 'Only JSON and NPY are materialized by this adapter')
    array = np.load(path, allow_pickle=False)
    _require(list(array.shape) == asset['shape'], f'Array shape mismatch: {asset["path"]}')
    _require(str(array.dtype) == asset['dtype'], f'Array dtype mismatch: {asset["path"]}')
    _require(np.isfinite(array).all(), f'Non-finite asset {asset["path"]}')
    return array


def validate_native(manifest: dict, root: str | Path, *, check_assets=True) -> dict:
    root = Path(root)
    json.dumps(manifest, allow_nan=False)
    schema = json.loads((CONTRACT_DIR/'MaterialHypothesisPackage.execution.schema.json').read_text())
    validate_schema_subset(manifest, schema)
    ids = [c['id'] for c in manifest['candidates']]
    _require(len(ids) == len(set(ids)), 'Duplicate candidate ids')
    _require(manifest['valid_k'] == len(ids) <= manifest['requested_k'], 'Invalid K')
    _require(manifest['weight_reliability'] == 'candidate_only' and manifest['weights'] is None,
             'This adapter only exports candidate-only semantics')
    for node in _walk(manifest):
        _require(not TRUTH_KEYS.intersection(node), 'Truth-derived fields in online handoff')
        if check_assets and {'path', 'sha256', 'format', 'shape', 'dtype'} <= node.keys():
            _asset_validate(node, root)
    for candidate in manifest['candidates']:
        f, v = candidate['fit'], candidate['validation']
        _require(f['whitened_residual_l2'] <= f['threshold_l2'], 'Candidate fails fit acceptance')
        _require(math.isclose(f['squared_residual'], f['whitened_residual_l2']**2, rel_tol=1e-7, abs_tol=1e-12), 'Inconsistent residual square')
        _require(v['physics_target_hash'] == manifest['physics_target']['hash'], 'Candidate target mismatch')
        _require(v['material_law_id'] == manifest['physics_target']['material_law_id'], 'Material-law mismatch')
        _require(v['solver_id'] == manifest['physics_target']['solver_id'], 'Solver identity mismatch')
        _require((v['constitutive_relative_defect_status'] == 'measured') == (v['constitutive_relative_defect'] is not None), 'An unexposed internal constitutive defect cannot be assigned a fabricated zero')
        _require(f['noise_model_id'] == manifest['noise']['id'], 'Noise identity mismatch')
        _require(not(f['nuisance_policy'] == 'bounded_profile' and f['likelihood_semantics'] == 'proper_marginal'), 'Profile residual is not marginal evidence')
        _require(candidate['posterior_weight'] is None, 'Candidate-only has no posterior probability')
    seen = set()
    for pair in manifest['pairs']:
        a, b = pair['candidate_i'], pair['candidate_j']
        _require(a in ids and b in ids and a != b, 'Invalid pair references')
        key = tuple(sorted((a, b)))
        _require(key not in seen, 'Duplicate pair')
        seen.add(key)
        _require(pair['connection_status'] != 'verified_path' or pair['validated_path_asset'] is not None,
                 'Verified path must have saved full-validated path asset')
    return dict(execution_profile=PROFILE, schema='PASS', semantic_checks='PASS',
                asset_arrays='PASS' if check_assets else 'NOT_CHECKED', digest_computation='NOT_RUN_BY_USER_POLICY',
                original_strict_schema='NOT_APPLICABLE: required digest fields are explicitly null under execution overlay',
                scope='Consistency only; validator does not invoke or independently certify Maxwell physics')


def write_packages(output_dir: str | Path, context: dict, candidates: list[dict], *, base: dict | None = None,
                   pair_distances: list[dict] | None = None) -> dict:
    """Write native-profile manifest, arrays and independent common exchange.

    No solver calls occur here. Every input candidate must already be accepted by
    the parent's full-Maxwell verifier. See INPUT_GUIDE.md for exact required data.
    All stored predictions must be effective (gain-applied) old-channel predictions.
    """
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    _require(bool(candidates), 'No actual accepted candidates; no handoff-ready package may be fabricated')
    for c in candidates:
        _require(c.get('validated_full_maxwell') is True, f'{c.get("id")}: missing full-Maxwell validation flag')
        _require(c.get('physics_target_id', context['target_id']) == context['target_id'], 'Candidate target differs')
    base = candidates[0] if base is None else base
    obs = np.asarray(context['observation'], dtype=float).reshape(-1)
    white = np.asarray(context['whitener'], dtype=float)
    cov = np.asarray(context['noise_covariance'], dtype=float)
    _require(white.shape == cov.shape == (len(obs), len(obs)), 'Whitening/covariance dimensions disagree')
    _require(np.allclose(white@cov@white.T, np.eye(len(obs)), atol=1e-7, rtol=1e-7), 'Whitening is inconsistent with covariance')
    now = utcnow()
    target, noise_id = context['target_id'], context['noise_id']
    material_law = context['material_law_id']
    acq = dict(context['acquisition'])
    acq.setdefault('n_real_measurements', len(obs))
    acq['manifest'] = _asset(root, 'acquisition', context['acquisition'])
    _require(acq['n_real_measurements'] == len(obs), 'Acquisition measurement count mismatch')
    assets = {}
    native_candidates = []
    bmat = np.asarray(base['material'])
    bpred = np.asarray(base['prediction'], dtype=float).reshape(-1)
    _require(bpred.shape == obs.shape, 'Base prediction shape mismatch')
    requested = int(context.get('requested_k', len(candidates)))
    for i, c in enumerate(candidates):
        cid = str(c['id'])
        _require(re.fullmatch(r'[A-Za-z0-9_.-]+', cid) is not None, 'Unsafe candidate asset id')
        mat = np.asarray(c['material'])
        pred = np.asarray(c['prediction'], dtype=float).reshape(-1)
        _require(mat.shape == bmat.shape and pred.shape == obs.shape, 'Material/prediction dimensions disagree')
        residual = float(np.linalg.norm(white@(pred-obs)))
        _require(math.isclose(residual, float(c['whitened_residual_l2']), rel_tol=1e-6, abs_tol=1e-8), f'{cid}: prediction does not replay declared residual')
        threshold = float(c.get('fit_threshold', context['fit_threshold']))
        _require(residual <= threshold, f'{cid}: rejected fit cannot enter archive')
        val = c['validation']
        _require(val.get('solver_relative_residual') is not None, f'{cid}: native validation requires measured solver_relative_residual')
        _require('constitutive_relative_defect' in val, f'{cid}: explicitly declare the constitutive defect or null')
        defect_status = val.get('constitutive_relative_defect_status', 'measured' if val['constitutive_relative_defect'] is not None else None)
        _require(defect_status in ('measured', 'NOT_EXPOSED_BY_ANALYTIC_MIE_INTERFACE'), f'{cid}: missing constitutive defect status')
        _require((defect_status == 'measured') == (val['constitutive_relative_defect'] is not None), f'{cid}: unexposed constitutive defect must remain null, never fabricated zero')
        _require(val.get('material_constraints_passed') is True, f'{cid}: parent must explicitly confirm physical material constraints')
        _require(val.get('all_illuminations_same_material') is True, f'{cid}: parent must explicitly confirm shared material')
        material_asset = _asset(root, f'{cid}_material', mat)
        prediction_asset = _asset(root, f'{cid}_prediction', pred)
        structure = dict(c['structure'])
        descriptor = structure.pop('descriptor', {'support_volume': structure['support_volume']})
        structure['descriptor_asset'] = _asset(root, f'{cid}_structure', descriptor)
        comp = dict(c['compensation'])
        comp['active_material_coordinates'] = _asset(root, f'{cid}_coordinates', comp.pop('coordinates'))
        comp.setdefault('path_asset', None)
        comp.setdefault('parent_candidate_id', None)
        comp.setdefault('protected_structure_constraints', [])
        cp = _cost(c['cost'])
        nc = dict(id=cid, material_field=material_asset,
                  material_parameterization_id=c['material_parameterization_id'], structure=structure,
                  full_predicted_measurement=prediction_asset,
                  fit=dict(whitened_residual_l2=residual, squared_residual=residual**2, threshold_l2=threshold,
                           noise_model_id=noise_id, nuisance_policy=c.get('nuisance_policy', 'bounded_profile'),
                           nuisance_parameters=c.get('nuisance_parameters', {}), log_likelihood=None,
                           likelihood_semantics='profile_not_marginal' if c.get('nuisance_policy', 'bounded_profile') == 'bounded_profile' else 'not_assigned',
                           likelihood_numerical_error=error_record(), fit_acceptance=True),
                  validation=dict(passed=True, scope='full_discrete_maxwell', solver_id=context['solver_id'],
                                  physics_target_hash=target, all_illuminations_same_material=True,
                                  constitutive_relative_defect=float(val['constitutive_relative_defect']) if val['constitutive_relative_defect'] is not None else None,
                                  constitutive_relative_defect_status=defect_status,
                                  solver_relative_residual=float(val['solver_relative_residual']),
                                  material_constraints_passed=True, material_law_id=material_law,
                                  numerical_error=val.get('numerical_error', error_record(status='solver_tolerance_only', method='parent supplied full numerical solver')),
                                  independent_mesh_checked=bool(val.get('independent_mesh_checked', False)),
                                  independent_mesh_result=val.get('independent_mesh_result'), timestamp_utc=c.get('timestamp_utc', now)),
                  compensation=comp, material_distance_to_base=float(c['material_distance_to_base']),
                  measurement_distance_to_base=float(np.linalg.norm(white@(pred-bpred))),
                  distance_definition_id=context['distance_definition_id'], proper_log_prior=c.get('proper_log_prior'),
                  prior_score_is_normalized_density=bool(c.get('prior_score_is_normalized_density', False)),
                  proposal_log_density=None, posterior_weight=None, conditional_samples_asset=None, cost=cp,
                  rom=c.get('rom', dict(used=False, engine='full_solver', anchor_id=None, rank=None,
                                       modified_current_rank=None, refresh_count=0, validity='not_applicable',
                                       numerical_error=error_record(), uniform_certificate_claimed=False)))
        native_candidates.append(nc)
        assets[cid] = dict(material=material_asset, prediction=prediction_asset)
    pairs = []
    for p in pair_distances or []:
        q = dict(p)
        q.setdefault('nuisance_policy', 'bounded_profile')
        q.setdefault('comparison_uncertainty', error_record())
        q.setdefault('hard_to_distinguish', None)
        q.setdefault('decision_rule_id', None)
        q.setdefault('connection_status', 'not_tested')
        q.setdefault('validated_path_asset', None)
        pairs.append(q)
    basset = assets[base['id']]['material'] if base['id'] in assets else _asset(root, 'base_material', bmat)
    interface = context['new_action_interface']
    native_interface = dict(evaluator_entrypoint=interface['evaluator_entrypoint'],
                            environment_manifest=_asset(root, 'environment', interface['environment']),
                            allowed_actions=_asset(root, 'allowed_actions', interface['allowed_actions']),
                            response_cache_manifest=_asset(root, 'cache_manifest', interface.get('cache_manifest', {})),
                            new_action_calls_require_separate_budget=True, posterior_information_gain_allowed=False,
                            candidate_only_design_semantics='Candidate discrimination/minimax or geometric design only; no posterior information gain without a separately justified probability target.')
    native = dict(version='1.0.0', package_id=context['package_id'], created_utc=now,
                  status='handoff_ready', observation=_asset(root, 'observation', obs),
                  noise=dict(id=noise_id, realification_order=context.get('realification_order', 'all-real-then-all-imaginary'),
                             covariance=_asset(root, 'noise_covariance', cov), whitener=_asset(root, 'noise_whitener', white),
                             units=context.get('noise_units', 'dimensionless scattering amplitude'),
                             model_discrepancy_policy=context['model_discrepancy_policy']), acquisition=acq,
                  physics_target=dict(hash=target, material_law_id=material_law, time_convention=context['time_convention'],
                                      grid=_asset(root, 'grid', context['grid']), backend_commit=context['backend_commit'],
                                      solver_id=context['solver_id'], self_term_convention=context['self_term_convention'],
                                      polarizability_map=context['polarizability_map'], boundary_conditions=context['boundary_conditions'],
                                      material_support=_asset(root, 'material_support', context['material_support']),
                                      nuisance_support=_asset(root, 'nuisance_support', context['nuisance_support']),
                                      prior_manifest=_asset(root, 'prior_manifest', context['prior_manifest']) if context.get('prior_manifest') is not None else None,
                                      dispersion_model=context['dispersion_model']),
                  base_reconstruction=dict(material_field=basset, nuisance_parameters=base.get('nuisance_parameters', {}),
                                           method=base.get('method', 'parent baseline reconstruction'), full_validation_passed=base['validated_full_maxwell'], cost=_cost(base['cost'])),
                  requested_k=requested, valid_k=len(candidates), candidate_ordering_rule=context['candidate_ordering_rule'],
                  candidates=native_candidates, pairs=pairs, weight_reliability='candidate_only', weights=None,
                  unexplored_uncertainty=dict(within_support='partially_audited', outside_support='unquantified',
                                             unexplored_regions=context['unexplored_regions'], modeled_residual_status='Not assigned a posterior mass; finite search archive only',
                                             modeled_residual_asset=None, known_search_restrictions=context['search_restrictions'], coverage_guarantee='none'),
                  next_measurement_interface=native_interface, total_cost=_cost(context['total_cost']),
                  upstream_manifest=_asset(root, 'upstream_manifest', context['upstream_manifest']),
                  failed_attempts_ledger=_asset(root, 'failed_attempts', context['failed_attempts']))
    _write_json(root/'manifest.json', native)
    audit = validate_native(native, root)
    exchange = native_to_exchange(native, root, histories={c['id']: c.get('history', []) for c in candidates},
                                  action_contract=interface)
    _write_json(root/'common_exchange_v1.json', exchange)
    _write_json(root/'CONTRACT_AUDIT.json', audit)
    for schema_name in ['MaterialHypothesisPackage.execution.schema.json', 'CommonMaterialExchange.v1.schema.json']:
        (root/schema_name).write_text((CONTRACT_DIR/schema_name).read_text())
    (root/'Original_A_MaterialHypothesisPackage.schema.json').write_text((CONTRACT_DIR/'original/MaterialHypothesisPackage.schema.json').read_text())
    _write_json(root/'EXECUTION_PROFILE.json', dict(profile=PROFILE, native_layout_version='1.0.0',
                original_strict_compatibility=False, reasons=['Human no-digest policy: asset.sha256=null',
                'Analytic Mie exposes boundary residual, not internal constitutive residual; explicit null/status supported'],
                candidate_only=True, posterior_probabilities=None, adapter_physics_calls=0))
    return dict(native_manifest=str(root/'manifest.json'), exchange_manifest=str(root/'common_exchange_v1.json'),
                audit=audit, manifest=native, exchange=exchange)


def native_to_exchange(native: dict, root: str | Path, *, histories: dict | None = None,
                       action_contract: dict | None = None) -> dict:
    """A-profile → independent exchange; does NOT claim compatibility with B-native schema."""
    root = Path(root)
    validate_native(native, root)
    histories = histories or {}
    action_contract = action_contract or {}
    out = dict(exchange_version=EXCHANGE_VERSION, source_profile=PROFILE,
               source_native_version=native['version'], package_id=native['package_id'], created_utc=native['created_utc'],
               target_identity=native['physics_target']['hash'], acquisition_identity=native['acquisition']['id'],
               noise_identity=native['noise']['id'], physics=native['physics_target'], acquisition=native['acquisition'],
               observation=native['observation'], noise=native['noise'], material_worlds=[], pair_relations=native['pairs'],
               probability_semantics=dict(status='candidate_only', weights=None, represented_mass=None,
                                          unshown_modeled_mass=None, unexplored_mass='unquantified', outside_model_support_mass=None,
                                          likelihood_profile_is_posterior=False),
               shared_calibration=dict(rule='One declared nuisance value per material world, shared by all specified frequencies/sources; never independent per channel',
                                       support=native['physics_target']['nuisance_support']),
               unexplored_structures=native['unexplored_uncertainty'],
               measurement_evaluation=dict(entrypoint=native['next_measurement_interface']['evaluator_entrypoint'],
                                           input_fields=['material_world_id', 'action_id', 'budget_authorization'],
                                           output_fields=['target_identity', 'action_id', 'realification_order', 'prediction', 'cost_receipt', 'validation_status'],
                                           allowed_actions_asset=native['next_measurement_interface']['allowed_actions'],
                                           allowed_action_ids=[a['id'] for a in action_contract.get('allowed_actions', [])],
                                           requires_new_explicit_budget=True, auto_evaluation=False, evaluator_calls_made_by_adapter=0,
                                           probability_design='not available: candidate-only archive',
                                           callable_contract='evaluate(candidate_material_array, action_dict, budget_authorization) -> prediction + cost_receipt; caller owns authorization and physics accounting'),
               total_cost=native['total_cost'], upstream_manifest=native['upstream_manifest'], failed_attempts=native['failed_attempts_ledger'],
               provenance=dict(native_manifest='manifest.json', original_schema='Research A 1.0.0 (preserved in repository)',
                               execution_overlay=PROFILE, digest_policy='not computed or checked; null digests are explicit',
                               b_native_schema_compatibility='not asserted; Research B must explicitly import this versioned common exchange'))
    for c in native['candidates']:
        out['material_worlds'].append(dict(id=c['id'], full_material_array=c['material_field'],
                                           material_parameterization_id=c['material_parameterization_id'], structure=c['structure'],
                                           old_measurement_prediction=c['full_predicted_measurement'], fit=c['fit'],
                                           full_maxwell_validation=c['validation'], compensation=c['compensation'],
                                           material_distance_to_base=c['material_distance_to_base'], measurement_distance_to_base=c['measurement_distance_to_base'],
                                           distance_definition_id=c['distance_definition_id'], generation_history=histories.get(c['id'], []),
                                           posterior_weight=None, physical_legality=True, cost=c['cost']))
    validate_exchange(out, root)
    return out


def validate_exchange(exchange: dict, root: str | Path) -> dict:
    root = Path(root)
    json.dumps(exchange, allow_nan=False)
    schema = json.loads((CONTRACT_DIR/'CommonMaterialExchange.v1.schema.json').read_text())
    validate_schema_subset(exchange, schema)
    _require(exchange['probability_semantics']['weights'] is None, 'Candidate-only cannot have weights')
    ids = [w['id'] for w in exchange['material_worlds']]
    _require(bool(ids) and len(ids) == len(set(ids)), 'Actual unique material worlds required')
    for node in _walk(exchange):
        _require(not TRUTH_KEYS.intersection(node), 'Truth leakage in exchange')
        if {'path', 'format', 'shape', 'dtype'} <= node.keys():
            _asset_validate(node, root)
    obs = _asset_validate(exchange['observation'], root)
    white = _asset_validate(exchange['noise']['whitener'], root)
    for c in exchange['material_worlds']:
        _require(c['posterior_weight'] is None, 'Invented posterior probability')
        v = c['full_maxwell_validation']
        _require(v['passed'] and v['physics_target_hash'] == exchange['target_identity'], 'Unvalidated or wrong-target candidate')
        _require(v['all_illuminations_same_material'] and v['material_constraints_passed'], 'Shared physical material required')
        _require((v['constitutive_relative_defect_status'] == 'measured') == (v['constitutive_relative_defect'] is not None), 'Unexposed constitutive defect must remain null')
        prediction = _asset_validate(c['old_measurement_prediction'], root)
        _require(prediction.shape == obs.shape, 'Old acquisition prediction shape mismatch')
        residual = np.linalg.norm(white@(prediction-obs))
        _require(math.isclose(residual, c['fit']['whitened_residual_l2'], rel_tol=1e-6, abs_tol=1e-8), 'Exchange prediction residual fails replay')
        _require(residual <= c['fit']['threshold_l2'], 'Exchange candidate is not measurement-compatible')
    _require(exchange['measurement_evaluation']['requires_new_explicit_budget'] is True,
             'New acquisitions require explicit budget authorization')
    return {'schema': 'PASS', 'semantics': 'PASS', 'physics_evaluations': 0, 'B_native_compatibility': 'NOT_ASSERTED'}


def evaluate_new_action(exchange: dict, package_root: str | Path, material_world_id: str, action_id: str,
                        budget_authorization: dict, evaluator: Callable[[np.ndarray, dict, dict], dict]) -> dict:
    """Explicit downstream hook, never run by export/validation.

    A real caller-supplied evaluator is mandatory. This routine does not create a
    budget, infer permission from a cached action, or silently run a solver.
    """
    validate_exchange(exchange, package_root)
    _require(budget_authorization.get('authorized') is True, 'New action has no explicit budget authorization')
    _require(budget_authorization.get('remaining_attempts', 0) > 0, 'New action budget exhausted')
    _require(bool(budget_authorization.get('ledger_id')), 'New action needs named accounting ledger')
    _require(budget_authorization.get('target_identity') == exchange['target_identity'], 'Budget target mismatch')
    actions = _asset_validate(exchange['measurement_evaluation']['allowed_actions_asset'], Path(package_root))
    action = next((a for a in actions if a['id'] == action_id), None)
    _require(action is not None, 'Unknown measurement action')
    world = next((w for w in exchange['material_worlds'] if w['id'] == material_world_id), None)
    _require(world is not None, 'Unknown material-world id')
    material = _asset_validate(world['full_material_array'], Path(package_root))
    result = evaluator(material, copy.deepcopy(action), copy.deepcopy(budget_authorization))
    _require(result.get('target_identity') == exchange['target_identity'], 'Evaluator changed target')
    _require(result.get('action_id') == action_id, 'Evaluator returned wrong action')
    _require('cost_receipt' in result and 'prediction' in result, 'Evaluator must return actual cost and prediction')
    _require(bool(result.get('validation_status')), 'Evaluator must state physical/numerical validation status')
    _require(result.get('realification_order') == exchange['noise']['realification_order'], 'New action changed declared real packing')
    prediction = np.asarray(result['prediction'])
    _require(prediction.ndim == 1 and prediction.dtype.kind in 'fiu' and np.isfinite(prediction).all(), 'New action prediction must be a finite real-packed vector')
    if 'n_real_measurements' in action:
        _require(len(prediction) == action['n_real_measurements'], 'New action prediction length differs from action specification')
    receipt = result['cost_receipt']
    _require(isinstance(receipt, dict) and bool(receipt), 'New action needs an actual accounting receipt')
    if 'new_full_state_frequency_attempts' in receipt:
        _require(0 <= receipt['new_full_state_frequency_attempts'] <= budget_authorization['remaining_attempts'], 'New action receipt exceeds its authorization')
    return result


def make_mie_context(*, package_id: str, target_id: str, backend_commit: str,
                     observation: np.ndarray, sigma: float | np.ndarray,
                     fit_threshold: float, acquisition: dict, grid: dict,
                     material_support: dict, nuisance_support: dict,
                     total_cost: dict, upstream_manifest: dict, failed_attempts: list | dict,
                     new_action_interface: dict, distance_definition_id: str,
                     noise_id: str, requested_k: int, **overrides) -> dict:
    """Small helper for the already-frozen concentric B1 Mie target.

    sigma is per-real-channel SD fixed from the observation protocol, never a
    candidate-dependent RMS. No solver, cache, ledger, prior or permission is
    created. Geometry grid must contain the complete radial/layer convention.
    """
    obs = np.asarray(observation, dtype=float).reshape(-1)
    std = np.broadcast_to(np.asarray(sigma, dtype=float), obs.shape).copy()
    _require(np.all(np.isfinite(std)) and np.all(std > 0), 'Positive finite per-real-channel noise SD required')
    out = dict(package_id=package_id, target_id=target_id, backend_commit=backend_commit,
               observation=obs, noise_covariance=np.diag(std**2), whitener=np.diag(1/std),
               fit_threshold=fit_threshold, acquisition=acquisition, grid=grid,
               material_support=material_support, nuisance_support=nuisance_support,
               total_cost=total_cost, upstream_manifest=upstream_manifest, failed_attempts=failed_attempts,
               new_action_interface=new_action_interface, distance_definition_id=distance_definition_id,
               noise_id=noise_id, requested_k=requested_k,
               solver_id='B1-vector-Mie-tangential-boundary-matching',
               material_law_id='passive-isotropic-relative-permittivity-mu1',
               time_convention='exp(-i omega t)', realification_order='all-real-then-all-imaginary',
               model_discrepancy_policy='Numerical full-Mie target only; model discrepancy outside concentric finite-band support is unquantified',
               self_term_convention='not applicable: analytic spherical boundary matching',
               polarizability_map='direct relative-epsilon layers; no voxel polarizability approximation',
               boundary_conditions='regular at origin, tangential electric/magnetic continuity, outgoing radiation',
               dispersion_model='shared frequency-independent complex epsilon over the declared finite band',
               prior_manifest=None, candidate_ordering_rule='discovery order under material-space deduplication; no truth used',
               unexplored_regions=['Undiscovered radial material worlds within declared support', 'Nonconcentric and outside-support materials are unquantified'],
               search_restrictions=['Finite declared radial representation and bounded CPU search budget'],
               noise_units='dimensionless scattering amplitude')
    out.update(overrides)
    return out


def _main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description='Validate Research A package consistency; zero solver calls, no digest computation.')
    parser.add_argument('manifest', type=Path)
    args = parser.parse_args()
    data = json.loads(args.manifest.read_text())
    if 'exchange_version' in data:
        report = validate_exchange(data, args.manifest.parent)
    else:
        report = validate_native(data, args.manifest.parent)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    _main()

#!/usr/bin/env python3
"""Freeze -> observation-only inference -> separate-PID generator scoring.

freeze --methods SPEC.json --out FROZEN.json [--target TARGET.json --geometry geometry.json]
infer --freeze FROZEN.json --contexts inputs/online --out ONLINE --device cpu
score --freeze FROZEN.json --online ONLINE --truth inputs/offline_only --out SCORES

No subcommand launches training, physics, candidate filtering, or selection.
CFM numeric NLL requires a pre-frozen density_audit with explicit parents,
nfes=[16,32,64] and an agreement_tolerance. Otherwise its NLL is NOT_RUN.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
from scipy.special import ndtr

from p3r4.probability_eval import *
from p3r4.evaluation import paired_bootstrap


def load_json(path):
    return json.loads(Path(path).read_text())


def resources(wall, cpu):
    return dict(wall_seconds=time.monotonic()-wall, cpu_seconds=time.process_time()-cpu,
                peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))


def variants(frozen):
    result = []
    for method in frozen['methods']:
        for nfe in ([16, 32, 64] if method['kind'] == 'cfm' else [None]):
            result.append(dict(id=method['id']+(f'__nfe{nfe}' if nfe else ''),
                               family=method['family']+(f'__nfe{nfe}' if nfe else ''),
                               seed=method['seed'], nfe=nfe, method=method))
    return result


def load_predictor(method, device):
    if method.get('loader') == 'anchor':
        from p3r4.anchor_density import load_anchor_density
        return load_anchor_density(method['checkpoint'], device)
    from p3r4.density import load_density
    return load_density(method['checkpoint'], device)


def sample_log_prob(predictor, x, y, condition, batch_size):
    values = []
    for start in range(0, len(x), batch_size):
        count = len(x[start:start+batch_size])
        c = {key: np.repeat(value, count, axis=0) for key, value in condition.items()}
        values.append(predictor.log_prob(x[start:start+batch_size], np.repeat(y, count, axis=0), batch_size=batch_size, **c))
    return np.concatenate(values)


def freeze(args):
    import torch
    spec = load_json(args.methods)
    methods, parents = spec['methods'], spec['parents']
    if not methods or not parents or len(parents) != len(set(parents)) or not all(map(safe_id, parents)):
        raise ValueError('nonempty unique safe parents and methods required')
    ids = [m['id'] for m in methods]
    if len(ids) != len(set(ids)) or not all(map(safe_id, ids)):
        raise ValueError('unique safe method IDs required')
    target = load_json(args.target)
    decoder = decoder_from_geometry(target, load_json(args.geometry)['points'])
    for method in methods:
        if not safe_id(method['family']) or not isinstance(method['seed'], int):
            raise ValueError('each method needs safe family and integer training seed')
        identity = checkpoint_identity(method['checkpoint'])
        payload = torch.load(identity['path'], map_location='cpu', weights_only=False)
        if payload.get('format') not in ('p3r4_full66_density_v1', 'p3r4_anchor_full66_density_v1'):
            raise ValueError('full66 R4 checkpoint required')
        if method['seed'] != payload['seed']:
            raise ValueError('declared training seed does not match checkpoint')
        method.update(checkpoint=identity['path'], checkpoint_identity=identity, kind=payload['kind'],
                      loader='anchor' if payload['format'] == 'p3r4_anchor_full66_density_v1' else 'raw_y',
                      j_summary=payload.get('anchor_contract', {}).get('j_summary', 'none'),
                      coordinate_preconditioner=payload.get('anchor_contract', {}).get('coordinate_preconditioner', 'none'),
                      anchor_base=payload.get('contract', {}).get('cfm_base') == 'anchor')
        audit = method.get('density_audit')
        if audit is not None:
            if method['kind'] != 'cfm' or audit.get('nfes') != [16, 32, 64]:
                raise ValueError('CFM density audit requires explicit nfes=[16,32,64]')
            selected = audit.get('parents', [])
            if not 0 < len(selected) <= 16 or len(set(selected)) != len(selected) or not set(selected) <= set(parents):
                raise ValueError('predeclare1..16 unique density-audit parents')
            if not np.isfinite(audit['agreement_tolerance']) or audit['agreement_tolerance'] <= 0:
                raise ValueError('positive numeric NLL agreement tolerance required')
    spec.update(format='p3r4_probability_evaluation_v1', status='FROZEN_BEFORE_TRUTH',
                freeze_pid=os.getpid(), freeze_time_ns=time.time_ns(), truth_read=False,
                sample_count=128, fixed_sample_index=0,
                sample_seed=int(spec.get('sample_seed', 20261025031)),
                decoder=decoder, projections=fixed_projections(decoder),
                coverage_levels=list(COVERAGE_LEVELS), phase_thresholds=list(PHASE_THRESHOLDS),
                bootstrap=dict(draws=4000, seed=BOOTSTRAP_SEED, level=CI_LEVEL, unit='independent parent; seed losses averaged before resampling'),
                online_allowlist=['observation', 'parent identity', 'target_id', 'optional y-only base_mean/base_chol',
                                  'explicit anchor-model context_x and optional context_J/J'],
                samples_selected_by_GT=False, physical_likelihood_status='NOT_RUN_AWAITING_PARENT_PAID_PREDICTIONS',
                comparisons=spec.get('comparisons', []))
    known_families = {v['family'] for v in variants(spec)}
    for comparison in spec['comparisons']:
        if comparison['candidate_family'] not in known_families or comparison['reference_family'] not in known_families:
            raise ValueError('comparison family not declared (CFM families have __nfe16/32/64 suffix)')
        if comparison.get('metric') not in ('raw_sample_mean_full_mse', 'material_energy_score', 'full66_nll'):
            raise ValueError('predeclare comparison metric raw mean MSE/material energy/full66 NLL')
    write_new_json(args.out, spec)
    return dict(status='FROZEN', path=str(args.out), methods=len(methods), parents=len(parents))


def infer(args, frozen, wall, cpu):
    if args.truth is not None:
        raise ValueError('truth path prohibited during online inference')
    if not 1 <= args.batch_size <= 16:
        raise ValueError('online inference batch must be<=16')
    import torch
    torch.set_num_threads(1)
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError('online destination must be new/empty; no implicit replay')
    args.out.mkdir(parents=True, exist_ok=True)
    artifacts, failures, records = [], [], []
    decoder = frozen['decoder']
    for method in frozen['methods']:
        if checkpoint_identity(method['checkpoint']) != method['checkpoint_identity']:
            raise ValueError('checkpoint changed since freeze')
        try:
            predictor = load_predictor(method, args.device)
            load_error = None
        except Exception as exc:
            predictor, load_error = None, repr(exc)
        for index, parent in enumerate(frozen['parents']):
            path = args.contexts/f'{parent}.npz' if args.contexts.is_dir() else args.contexts
            for variant in [v for v in variants(frozen) if v['method']['id'] == method['id']]:
                record = dict(parent=parent, variant=variant['id'], status='FAILED')
                try:
                    if load_error:
                        raise RuntimeError(load_error)
                    y, base = read_online_observation(path, parent, need_base=method['anchor_base'] and method.get('loader') != 'anchor',
                                                       target_id=decoder['config']['target_id'])
                    if method.get('loader') == 'anchor':
                        base = read_online_anchor_context(path, parent, j_summary=method['j_summary'],
                            coordinate_preconditioner=method.get('coordinate_preconditioner', 'none'))
                    # Stable parent seed, identical across method families and
                    # CFM NFEs; sample index0 is fixed before any truth access.
                    seed = int(np.random.SeedSequence([frozen['sample_seed'], index]).generate_state(1)[0])
                    kwargs = dict(nfe=variant['nfe']) if variant['nfe'] else {}
                    if method.get('loader') == 'anchor':
                        x = predictor.sample(y, n_samples=128, seed=seed, batch_size=args.batch_size, **base, **kwargs)[0]
                    else:
                        x = predictor.sample(y, 128, seed, args.batch_size, **base, **kwargs)[0]
                    chi, fraction, boundary_count = decode_gaussianized(x, decoder)
                    if not np.isfinite(chi).all():
                        raise FloatingPointError('nonfinite decoded world')
                    mean = chi.mean(0)
                    values = dict(parent=np.array(parent), observation_features=y, x_samples=x,
                                  q_samples=np.column_stack((x[:, :64], ndtr(x[:, 64:]))),
                                  chi_samples=chi, material_mean=mean,
                                  material_variance=np.mean(abs(chi-mean)**2, axis=0),
                                  material_variance_real=chi.real.var(0), material_variance_imag=chi.imag.var(0),
                                  sample_index=np.arange(128), fixed_sample_index=np.array(0),
                                  phase_event_probability=(fraction[:, decoder['roi'], None] >= np.array(PHASE_THRESHOLDS)).mean(0),
                                  probit_floating_boundary_count=np.array(boundary_count))
                    values.update(base)
                    if method['kind'] != 'cfm':
                        values['log_q_x'] = sample_log_prob(predictor, x, y, base, args.batch_size)
                        if not np.isfinite(values['log_q_x']).all():
                            raise FloatingPointError('nonfinite explicit density for samples')
                    dest = args.out/variant['id']/f'{parent}.npz'
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    with dest.open('xb') as file:
                        np.savez_compressed(file, **values)
                    identity = checkpoint_identity(dest)
                    artifacts.append(dict(relative_path=str(dest.relative_to(args.out)), identity=identity))
                    record.update(status='COMPLETE', relative_path=str(dest.relative_to(args.out)),
                                  seed=seed, sample_index=0, nfe=variant['nfe'],
                                  density='ALGEBRAIC_FULL66' if method['kind'] != 'cfm' else 'NOT_SAVED_NUMERICAL_CFM',
                                  probit_floating_boundary_count=boundary_count)
                except Exception as exc:
                    record['error'] = repr(exc)
                    failures.append(record.copy())
                records.append(record)
        del predictor
    gate = dict(format='p3r4_probability_evaluation_v1', status='ONLINE_COMPLETE', pid=os.getpid(),
                completed_time_ns=time.time_ns(), freeze=frozen, truth_read=False, physics_calls=0,
                optimizer_updates=0, records=records, artifacts=artifacts, failures=failures,
                no_fallback_law=True, no_candidate_selection=True, fixed_sample_index=0,
                K=128, maximum_network_batch=args.batch_size, **resources(wall, cpu))
    write_new_json(args.out/'ONLINE_COMPLETE.json', gate)
    return dict(status='ONLINE_COMPLETE', parents=len(frozen['parents']), variants=len(variants(frozen)), failures=len(failures))


def score(args, frozen, wall, cpu):
    if 'offline_only' not in args.truth.resolve().parts:
        raise ValueError('GT must come from a declared offline_only directory')
    gate = load_json(args.online/'ONLINE_COMPLETE.json')
    check_online_gate(frozen, gate, args.online)  # Precedes every truth read.
    expected_records = {(v['id'], p) for v in variants(frozen) for p in frozen['parents']}
    actual_records = [(r['variant'], r['parent']) for r in gate['records']]
    if len(actual_records) != len(expected_records) or set(actual_records) != expected_records:
        raise ValueError('online record set does not cover every frozen parent/method/NFE')
    if args.out.exists() and any(args.out.iterdir()):
        raise FileExistsError('scoring destination must be new/empty')
    args.out.mkdir(parents=True, exist_ok=True)
    import torch
    torch.set_num_threads(1)
    parents, decoder = frozen['parents'], frozen['decoder']
    lookup = {(r['variant'], r['parent']): r for r in gate['records']}
    all_rows, all_summary, family, density_audit_cache = {}, {}, {}, {}
    for variant in variants(frozen):
        method = variant['method']
        predictor = load_predictor(method, args.device)
        rows = []
        numeric_audit = []
        for index, parent in enumerate(parents):
            record = lookup[(variant['id'], parent)]
            if record['status'] != 'COMPLETE':
                rows.append(None)
                continue
            # First truth access occurs only after all inference artifacts pass
            # the gate above, and this process has a different PID.
            x_truth, chi_truth = read_truth_after_gate(args.truth/f'{parent}.npz')
            with np.load(args.online/record['relative_path'], allow_pickle=False) as z:
                if str(z['parent'].item()) != parent or int(z['fixed_sample_index']) != 0 or not np.array_equal(z['sample_index'], np.arange(128)):
                    raise ValueError('sealed sample identity mismatch')
                x, chi = z['x_samples'], z['chi_samples']
                y = z['observation_features']
                if method.get('loader') == 'anchor':
                    base = {key: z[key] for key in ('context_x', 'context_J', 'base_mean', 'base_chol') if key in z}
                else:
                    base = {key: z[key] for key in ('base_mean', 'base_chol')} if method['anchor_base'] else {}
                if not np.allclose(z['material_mean'], chi.mean(0), atol=1e-12, rtol=0):
                    raise ValueError('sample mean mismatch')
            decoded, fraction, _ = decode_gaussianized(x, decoder)
            true_decoded, true_fraction, _ = decode_gaussianized(x_truth[None], decoder)
            if not np.allclose(decoded, chi, atol=1e-12, rtol=0) or not np.allclose(true_decoded[0], chi_truth, atol=1e-10, rtol=0):
                raise ValueError('original decoder/truth material mismatch')
            seed = np.random.SeedSequence([BOOTSTRAP_SEED, index])
            metrics = parent_scores(x, chi, fraction, x_truth, chi_truth, true_fraction[0],
                                    frozen['projections'], decoder['roi'], seed)
            if method['kind'] != 'cfm':
                metrics['full66_nll'] = float(-predictor.log_prob(x_truth[None], y, batch_size=1, **base)[0])
            else:
                audit = method.get('density_audit')
                if audit and parent in audit['parents']:
                    cache_key = (method['id'], parent)
                    if cache_key not in density_audit_cache:
                        density_audit_cache[cache_key] = [-float(predictor.log_prob(x_truth[None], y, batch_size=1, allow_approximate=True, nfe=nfe, **base)[0])
                                                        for nfe in audit['nfes']]
                    nll = density_audit_cache[cache_key]
                    numeric_audit.append(dict(parent=parent, nfes=audit['nfes'], numeric_continuous_nll=nll,
                                              difference_32_64=abs(nll[1]-nll[2]), agreement_tolerance=audit['agreement_tolerance'],
                                              within_tolerance=abs(nll[1]-nll[2]) <= audit['agreement_tolerance'],
                                              exact_discrete_map_density=False))
            rows.append(metrics)
        keys = set().union(*(r.keys() for r in rows if r is not None))
        arrays = {'parent': np.array(parents), 'success': np.array([r is not None for r in rows])}
        for key in sorted(keys):
            shape = np.shape(next(r[key] for r in rows if r is not None))
            arrays[key] = np.stack([np.asarray(r[key]) if r is not None else np.full(shape, np.nan) for r in rows])
        np.savez_compressed(args.out/(variant['id']+'_per_parent.npz'), **arrays)
        summaries = {}
        targets = ['material_energy_score', 'raw_sample_mean_full_mse', 'raw_sample_mean_full_mae', 'full66_nll']
        targets += [key for key in keys if 'coverage' in key or key.endswith('crps') or key.endswith('brier_per_threshold') or key.endswith('jeffreys_per_threshold')]
        for key in targets:
            if key in arrays:
                value = arrays[key].reshape(len(parents), -1).mean(1)
                summaries[key] = bootstrap_mean(value)
                all_rows.setdefault(variant['id'], {})[key] = value
        for key in ('phase_brier_per_threshold', 'phase_log_score_jeffreys_per_threshold'):
            if key in arrays:
                summaries[key+'_by_threshold'] = [bootstrap_mean(arrays[key][:, j]) for j in range(3)]
        if 'projection_crps' in arrays:
            summaries['projection_crps_by_frozen_projection'] = [bootstrap_mean(arrays['projection_crps'][:, j]) for j in range(20)]
        all_summary[variant['id']] = dict(family=variant['family'], training_seed=method['seed'], kind=method['kind'], nfe=variant['nfe'],
            planned_parents=len(parents), successful_parents=int(arrays['success'].sum()), failures=len(parents)-int(arrays['success'].sum()),
            metrics=summaries, cfm_numeric_density_audit=numeric_audit,
            full66_nll_status='EVALUATED_ON_ALL_PARENTS' if method['kind'] != 'cfm' and arrays['success'].all() else
                              ('NOT_EVALUATED_INCOMPLETE' if method['kind'] != 'cfm' else 'NOT_EXACT_CFM_NUMERIC_AUDIT_ONLY' if numeric_audit else 'NOT_RUN_CFM_AUDIT_NOT_FROZEN'))
        family.setdefault(variant['family'], []).append(variant)
        del predictor
    families = {}
    family_values = {}
    for name, group in family.items():
        seeds = [v['seed'] for v in group]
        if len(set(seeds)) != len(seeds):
            raise ValueError('duplicate training seed within family/NFE')
        keys = set.intersection(*(set(all_rows.get(v['id'], {})) for v in group))
        families[name] = dict(training_seeds=seeds, metrics={}, aggregation='average per-parent per-training-seed losses; no image/sample ensemble')
        family_values[name] = {}
        for key in keys:
            value = np.stack([all_rows[v['id']][key] for v in group]).mean(0)
            family_values[name][key] = value
            families[name]['metrics'][key] = bootstrap_mean(value)
    comparisons = []
    for comparison in frozen['comparisons']:
        a = family_values[comparison['candidate_family']].get(comparison['metric'])
        b = family_values[comparison['reference_family']].get(comparison['metric'])
        if a is None or b is None or not np.isfinite(a).all() or not np.isfinite(b).all():
            comparisons.append(dict(**comparison, status='NOT_EVALUATED_INCOMPLETE_OR_DENSITY_UNAVAILABLE'))
        else:
            comparisons.append(dict(**comparison, status='EVALUATED', **paired_bootstrap(a, b, seed=BOOTSTRAP_SEED)))
    result = dict(format='p3r4_probability_evaluation_v1', status='OFFLINE_SCORE_COMPLETE', scoring_pid=os.getpid(), online_pid=gate['pid'],
        n_parent=len(parents), sample_count=128, fixed_sample_index=0, methods=all_summary, families=families, comparisons=comparisons,
        all_online_failures=gate['failures'], no_candidate_selection=True, physics_calls=0, optimizer_updates=0,
        full_Maxwell_likelihood='NOT_RUN_AWAITING_PARENT_PAID_PREDICTIONS',
        Gamma256_residual_law='NOT_RUN_AWAITING_PARENT_PAID_PREDICTIONS',
        data_dependent_SBC='NOT_RUN_AWAITING_PARENT_PAID_PREDICTIONS',
        latent_marginal_rank_diagnostic='EVALUATED_WHERE_ONLINE_SUCCEEDED; not data-dependent SBC or posterior qualification',
        coverage_interpretation='marginal coordinate/projection coverage only; randomized inclusive ranks handle atoms separately from ordinary empirical intervals',
        CI_interpretation='pointwise intervals at the frozen protocol level; no simultaneous coverage guarantee over all diagnostic metrics',
        phase_score_interpretation='raw sample event probability Brier; Jeffreys finite logscore uses(count+.5)/(K+1); raw zero-support counts/logscores saved separately',
        material_variance_interpretation='empirical ddof0 per-pixel complex E|chi-Echi|^2 plus real/imag marginal variances; no full covariance claimed',
        truth_access='separate PID after checking complete sealed prediction set', posterior_qualification='NOT_ESTABLISHED',
        **resources(wall, cpu))
    write_new_json(args.out/'RESULTS.json', result)
    return dict(status='OFFLINE_SCORE_COMPLETE', parents=len(parents), variants=len(all_summary), comparisons=len(comparisons))


def main(argv=None):
    stage = ROOT/'research/p3l_r4_conditional_transport'
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('freeze', 'infer', 'score'))
    parser.add_argument('--methods', type=Path)
    parser.add_argument('--freeze', type=Path)
    parser.add_argument('--contexts', type=Path)
    parser.add_argument('--online', type=Path)
    parser.add_argument('--truth', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--target', type=Path, default=stage/'configs/TARGET.json')
    parser.add_argument('--geometry', type=Path, default=ROOT/'research/p1_material_target_20261018/inputs/geometry.json')
    parser.add_argument('--device', choices=('cpu', 'mps'), default='cpu')
    parser.add_argument('--batch-size', type=int, default=16)
    args = parser.parse_args(argv)
    wall, cpu = time.monotonic(), time.process_time()
    if args.command == 'freeze':
        result = freeze(args)
    else:
        frozen = load_json(args.freeze)
        if frozen.get('status') != 'FROZEN_BEFORE_TRUTH' or frozen.get('sample_count') != 128 or frozen.get('fixed_sample_index') != 0:
            raise ValueError('invalid frozen probability evaluation contract')
        result = infer(args, frozen, wall, cpu) if args.command == 'infer' else score(args, frozen, wall, cpu)
    print(json.dumps(result, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()

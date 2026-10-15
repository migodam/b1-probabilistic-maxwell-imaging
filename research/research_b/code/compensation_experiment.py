"""Bounded nonlinear compensation witnesses; no solver imports or execution.

Only ``run`` may request paid validation, through the supplied callback. All
optimization uses the explicitly approximate Response and keeps geometry fixed.
Distances are profile residuals, not evidences or posterior probabilities.
"""
from __future__ import annotations

import csv
import json
import os
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from src.b1_joint.material import BOUNDS, from_geometry, geometry
try:
    from .sensing import field_l1
except ImportError:
    from sensing import field_l1

OLD_FREQUENCY = 1.6
MAX_CALLBACK_ATTEMPTS = 30  # Includes reference, endpoint and failed callbacks.
MAX_NFEV = 80
GEOMETRIES = (('small_cavity', .08, 0.), ('large_cavity', .25, 0.),
              ('annular_cavity', 0., .08))
STARTS = ((2.5, .1, 1.), (2.25, .05, 1.15), (2.75, .15, .85))


def _plain(obj):
    if isinstance(obj, dict):
        return {str(k): _plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_plain(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return _plain(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _dump(path, obj):
    path = Path(path)
    path.write_text(json.dumps(_plain(obj), ensure_ascii=False, indent=2,
                               allow_nan=False) + '\n')


def _append(path, obj):
    with Path(path).open('a') as stream:
        stream.write(json.dumps(_plain(obj), ensure_ascii=False, allow_nan=False) + '\n')
        stream.flush()
        os.fsync(stream.fileno())


def _field(value):
    if isinstance(value, tuple):
        value = value[0]
    value = np.asarray(value)
    if np.iscomplexobj(value):
        raise ValueError('Complex callback output requires explicit Re16/Im16 packing')
    value = np.asarray(value, dtype=float)
    if value.shape == (1, 32):
        value = value[0]
    if value.shape != (32,) or not np.isfinite(value).all():
        raise ValueError('Expected one finite 32-real-channel Re16/Im16 field')
    return value


class _PaidVerifier:
    """A local resumable cap; global billing remains the callback's authority."""

    def __init__(self, callback, path):
        self.callback, self.path = callback, Path(path)
        self.attempts, self.cache, self.unavailable = 0, {}, {}
        if self.path.exists():
            for line in self.path.read_text().splitlines():
                row = json.loads(line)
                key = row['key']
                if row['event'] == 'BEGIN':
                    self.attempts += 1
                    self.unavailable[key] = 'previous attempt has no successful receipt'
                elif row.get('status') == 'PASS':
                    self.cache[key] = row
                    self.unavailable.pop(key, None)
                else:
                    self.unavailable[key] = row.get('error', 'previous callback failure')
        if self.attempts > MAX_CALLBACK_ATTEMPTS:
            raise RuntimeError('Existing compensation callback log exceeds fixed cap')

    def call(self, x, frequency, kind):
        x = np.asarray(x, dtype=float)
        key = json.dumps([x.tolist(), float(frequency)], separators=(',', ':'))
        if key in self.cache:
            return _field(self.cache[key]['field']), 'RESTORED_OR_SHARED_FULL_CACHE'
        if key in self.unavailable:
            return None, self.unavailable[key]
        if self.attempts >= MAX_CALLBACK_ATTEMPTS:
            return None, 'LOCAL_CALLBACK_ATTEMPT_CAP_EXHAUSTED'
        self.attempts += 1
        info = dict(key=key, callback_attempt=self.attempts, latent=x,
                    frequency=float(frequency), kind=kind)
        _append(self.path, dict(event='BEGIN', **info))
        start = time.perf_counter()
        try:
            field = _field(self.callback(x.copy(), float(frequency), kind=kind))
        except Exception as error:
            reason = type(error).__name__ + ': ' + str(error)
            _append(self.path, dict(event='END', status='FAIL', error=reason,
                                   wall_seconds=time.perf_counter()-start, **info))
            self.unavailable[key] = reason
            # Do not turn the orchestrator's global stop into more callback attempts.
            if type(error).__name__ in ('BudgetExceeded', 'Limit'):
                raise
            return None, reason
        row = dict(event='END', status='PASS', field=field,
                   wall_seconds=time.perf_counter()-start, **info)
        _append(self.path, row)
        self.cache[key] = _plain(row)
        self.unavailable.pop(key, None)
        return field, 'FULL_CALLBACK_SUCCESS'


def _gain(field, target, sigma):
    a, b = field/sigma, target/sigma
    return float(np.clip(np.dot(a, b)/max(float(np.dot(a, a)), 1e-300), .7, 1.3))


def _distances(fields, targets, sigmas, gains, calibration_sd=None):
    parts = [float(np.linalg.norm((g*f-t)/s))
             for f, t, s, g in zip(fields, targets, sigmas, gains)]
    cal = None if calibration_sd is None else float((gains[0]-1)/calibration_sd)
    return dict(per_frequency=parts, field_distance=float(np.linalg.norm(parts)),
                calibration_residual=cal,
                total_distance=float(np.linalg.norm(parts + ([] if cal is None else [cal]))))


def _fit(response, x, frequencies, targets, sigmas, label, trial_path,
         independent=False, calibration_sd=None):
    """Three bounded starts; all attempts retained; epsilon stays shared."""
    gain_count = len(frequencies) if independent else 1
    lower = np.r_[BOUNDS[2:, 0], np.repeat(.7, gain_count)]
    upper = np.r_[BOUNDS[2:, 1], np.repeat(1.3, gain_count)]

    def prediction(parameters):
        material = x.copy()
        material[2:] = parameters[:2]
        assert np.array_equal(material[:2], x[:2]), 'Geometry must not change'
        fields = response.predict(material, frequencies)[0].reshape(len(frequencies), 32)
        if not np.isfinite(fields).all():
            raise FloatingPointError('Nonfinite proxy response')
        return material, fields

    def residual(parameters):
        _, fields = prediction(parameters)
        gains = parameters[2:] if independent else np.repeat(parameters[2], len(frequencies))
        out = [(g*f-target)/sd for f, target, sd, g in zip(fields, targets, sigmas, gains)]
        if calibration_sd is not None:
            out.append(np.array([(parameters[2]-1)/calibration_sd]))
        return np.concatenate(out)

    candidates, failures = [], []
    for start_id, (er, ei, gain) in enumerate(STARTS):
        initial = np.r_[er, ei, np.repeat(gain, gain_count)]
        start = time.perf_counter()
        try:
            opt = least_squares(residual, initial, bounds=(lower, upper),
                                max_nfev=MAX_NFEV, xtol=1e-11, ftol=1e-11, gtol=1e-11)
            if not np.isfinite(opt.x).all() or not np.isfinite(opt.fun).all():
                raise FloatingPointError('Nonfinite optimizer endpoint')
            trial = dict(label=label, start_id=start_id, initial=initial,
                         endpoint=opt.x, success=bool(opt.success), status=int(opt.status),
                         message=str(opt.message), nfev=int(opt.nfev),
                         njev=None if opt.njev is None else int(opt.njev),
                         objective_squared=float(np.dot(opt.fun, opt.fun)),
                         optimality=float(opt.optimality), active_mask=opt.active_mask,
                         wall_seconds=time.perf_counter()-start)
            candidates.append((trial, opt.x))
        except Exception as error:
            trial = dict(label=label, start_id=start_id, initial=initial,
                         success=False, status='EXCEPTION', error=type(error).__name__+': '+str(error),
                         wall_seconds=time.perf_counter()-start)
            failures.append(trial)
        _append(trial_path, trial)
    if not candidates:
        return dict(status='ALL_PROXY_OPTIMIZER_STARTS_FAILED', failures=failures)
    trial, parameters = min(candidates, key=lambda item: item[0]['objective_squared'])
    material, fields = prediction(parameters)
    gains = parameters[2:] if independent else np.repeat(parameters[2], len(frequencies))
    return dict(status='PROXY_BEST_FOUND' if trial['success'] else 'PROXY_NONCONVERGED_FEASIBLE_ENDPOINT',
                target='proxy candidate response against noiseless full-reference fields',
                minimum_status='bounded best-found endpoint; no global-minimum guarantee',
                latent=material, epsilon=parameters[:2], gains=gains,
                geometry=geometry(material), geometry_frozen=True,
                frequencies=list(frequencies), proxy_unscaled_fields=fields,
                distances=_distances(fields, targets, sigmas, gains, calibration_sd),
                material_L1=float(field_l1(from_geometry(0., 0.), material)[0, 0]),
                active_bounds=[bool(a) for a in (np.minimum((parameters-lower)/(upper-lower),
                                                         (upper-parameters)/(upper-lower)) < 1e-6)],
                parameter_names=['epsilon_real', 'epsilon_imag'] +
                    ([f'gain_{k:g}' for k in frequencies] if independent else ['shared_gain']),
                fraction_to_lower_bound=(parameters-lower)/(upper-lower),
                fraction_to_upper_bound=(upper-parameters)/(upper-lower),
                selected_start=trial['start_id'], optimizer_success=trial['success'],
                failed_start_count=len(failures),
                nonconverged_start_count=sum(not a[0]['success'] for a in candidates),
                gain_model='ILLEGAL independent-frequency gain negative control' if independent else 'one shared gain per candidate world',
                new_measurement_blocks=int(any(k != OLD_FREQUENCY for k in frequencies)) + int(calibration_sd is not None))


def _write_outputs(directory, result, verifier, final=False):
    _dump(directory/('COMPENSATION_MECHANISMS.json' if final else 'COMPENSATION_PARTIAL.json'), result)
    fields = ['geometry', 'new_frequency', 'status', 'raw_proxy_joint_distance',
              'old_only_proxy_distance', 'new_only_proxy_distance', 'joint_shared_proxy_distance',
              'independent_gain_proxy_distance', 'calibration_joint_proxy_distance',
              'joint_shared_full_distance', 'full_gain_conflict_squared', 'material_L1']
    with (directory/'COMPENSATION_MECHANISMS.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in result['records']:
            fit_distance = lambda key: row.get(key, {}).get('distances', {}).get('total_distance')
            writer.writerow(dict(geometry=row['geometry'], new_frequency=row['new_frequency'],
                status=row['status'], raw_proxy_joint_distance=row.get('raw', {}).get('total_distance'),
                old_only_proxy_distance=fit_distance('old_only'), new_only_proxy_distance=fit_distance('new_only'),
                joint_shared_proxy_distance=fit_distance('joint_shared'),
                independent_gain_proxy_distance=fit_distance('independent_gain_negative_control'),
                calibration_joint_proxy_distance=fit_distance('calibration_plus_joint'),
                joint_shared_full_distance=row.get('full_validation', {}).get('shared', {}).get('total_distance'),
                full_gain_conflict_squared=row.get('full_validation', {}).get('gain_conflict_squared'),
                material_L1=row.get('joint_shared', {}).get('material_L1')))
    paid = list(verifier.cache.values())
    proxy = []
    for row in result['records']:
        for name, fit in row.items():
            if not isinstance(fit, dict) or 'proxy_unscaled_fields' not in fit:
                continue
            for frequency, field in zip(fit['frequencies'], fit['proxy_unscaled_fields']):
                proxy.append((fit['latent'], frequency, field, row['geometry']+'/'+name))
    np.savez_compressed(directory/'COMPENSATION_RAW_PREDICTIONS.npz',
        full_latent=np.array([r['latent'] for r in paid]).reshape(-1, 4),
        full_frequency=np.array([r['frequency'] for r in paid]),
        full_field=np.array([r['field'] for r in paid]).reshape(-1, 32),
        proxy_latent=np.array([r[0] for r in proxy]).reshape(-1, 4),
        proxy_frequency=np.array([r[1] for r in proxy]),
        proxy_field=np.array([r[2] for r in proxy]).reshape(-1, 32),
        proxy_label=np.array([r[3] for r in proxy], dtype=str))


def run(response, paid_predict_callback, output_dir):
    """Run only inside the orchestrator's charged E1 job, returning a table.

    Callback signature: ``callback(latent, frequency, kind=...) -> real32``;
    the legacy forward tuple is also accepted. Exactly 29 callbacks are planned
    (5 reference + 12 endpoints * 2 frequencies), at most 30 including failures.
    Callback cache hits may be cheaper, but are never assumed free here.
    """
    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)
    final = directory/'COMPENSATION_MECHANISMS.json'
    if final.exists():
        return json.loads(final.read_text())
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    start_proxy_evaluations = int(response.evaluations)
    if not callable(paid_predict_callback):
        raise TypeError('A globally metered paid_predict_callback is required')
    cfg = json.loads((Path(__file__).resolve().parents[1]/'configs/E1_FROZEN.json').read_text())
    new_frequencies = tuple(float(k) for k in cfg['new_frequencies'])
    if new_frequencies != (1.2, 2., 2.8, 3.2):
        raise ValueError('Compensation endpoint allocation was frozen for four new frequencies')
    frequencies = (OLD_FREQUENCY, *new_frequencies)
    if any(k not in response.coeff for k in frequencies):
        raise ValueError('All frozen frequencies must exist before compensation optimization')
    if not np.array_equal(np.asarray(response.bounds)[2:], BOUNDS[2:]):
        raise ValueError('Response material bounds differ from the frozen epsilon support')
    reference = from_geometry(0., 0., 2.5, .1)
    verifier = _PaidVerifier(paid_predict_callback, directory/'COMPENSATION_CALLBACK_ATTEMPTS.jsonl')
    reference_fields, sigma, reference_status = {}, {}, {}
    for k in frequencies:
        field, status = verifier.call(reference, k, 'RESEARCH_B_COMPENSATION_FIXED_REFERENCE')
        reference_status[k] = status
        if field is not None:
            sd = .12*np.sqrt(np.mean(field*field))
            if sd <= 0 or not np.isfinite(sd):
                raise ValueError('Canonical fixed-noise anchor must have positive RMS')
            reference_fields[k], sigma[k] = field, float(sd)
    result = dict(schema='research-b-compensation-v1',
        target='finite-frequency passive concentric vector-Mie calibration benchmark',
        hypotheses='new frequency may expose shared compensation conflict; no universal benefit assumed',
        prior_bounds=BOUNDS, gain_bounds=[.7, 1.3], reference_latent=reference, reference_gain=1.,
        reference_status=reference_status, sigma_per_real_channel=sigma,
        noise_rule='.12 canonical solid full-response real RMS per frequency',
        dispersion_assumption=cfg['dispersion'],
        gain_scope='fixed reference-world gain1 is a witness; candidate worlds may have different lawful gains',
        profile_is_not_posterior=True, global_minimum_established=False,
        selected_full_endpoint_policy='all 3 fixed geometries x 4 frozen new frequencies, joint-shared only',
        callback_attempt_cap=MAX_CALLBACK_ATTEMPTS, planned_callback_attempts=29,
        frozen_config_id=cfg['id'], response_class=type(response).__module__+'.'+type(response).__name__,
        response_frequencies=sorted(float(k) for k in response.coeff),
        nonjoint_endpoints_full_validated=False, records=[])
    trial_path = directory/'COMPENSATION_OPTIMIZER_TRIALS.jsonl'
    for name, r, h in GEOMETRIES:
        nominal = from_geometry(r, h, 2.5, .1)
        old_fit, calibration_old_fit = None, None
        for k in new_frequencies:
            row = dict(geometry=name, fixed_radius=r, fixed_annular_halfwidth=h,
                       new_frequency=k, status='NOT_AVAILABLE_REFERENCE_CALLBACK_FAILURE')
            if OLD_FREQUENCY not in reference_fields or k not in reference_fields:
                result['records'].append(row)
                continue
            ks = (OLD_FREQUENCY, k)
            targets = [reference_fields[a] for a in ks]
            sigmas = [sigma[a] for a in ks]
            nominal_fields = response.predict(nominal, ks)[0].reshape(2, 32)
            row['raw'] = dict(**_distances(nominal_fields, targets, sigmas, [1., 1.]),
                target='PROXY_NOMINAL_VS_FULL_REFERENCE; not full raw validation',
                latent=nominal, frequencies=ks, proxy_unscaled_fields=nominal_fields,
                material_L1=float(field_l1(reference, nominal)[0, 0]))
            if old_fit is None:
                old_fit = _fit(response, nominal, [OLD_FREQUENCY], targets[:1], sigmas[:1],
                               name+'/old_only', trial_path)
            row['old_only'] = old_fit
            row['new_only'] = _fit(response, nominal, [k], targets[1:], sigmas[1:], name+f'/{k:g}/new_only', trial_path)
            row['joint_shared'] = _fit(response, nominal, ks, targets, sigmas, name+f'/{k:g}/joint_shared', trial_path)
            row['independent_gain_negative_control'] = _fit(response, nominal, ks, targets, sigmas,
                name+f'/{k:g}/independent_gain', trial_path, independent=True)
            cal_sd = sigma[OLD_FREQUENCY]/np.linalg.norm(reference_fields[OLD_FREQUENCY])
            if calibration_old_fit is None:
                calibration_old_fit = _fit(response, nominal, [OLD_FREQUENCY], targets[:1], sigmas[:1],
                    name+'/calibration_old', trial_path, calibration_sd=cal_sd)
            row['calibration_plus_old'] = calibration_old_fit
            row['calibration_plus_joint'] = _fit(response, nominal, ks, targets, sigmas,
                name+f'/{k:g}/calibration_joint', trial_path, calibration_sd=cal_sd)
            row['calibration'] = dict(gain_statistic=1., gain_statistic_sd=float(cal_sd),
                direct_material_sensitivity=0., new_blocks=1,
                interpretation='independent known-reference acquisition sharing the same candidate-world calibration model')
            fit = row['joint_shared']
            if 'latent' in fit:
                full, statuses = [], []
                for frequency in ks:
                    field, status = verifier.call(fit['latent'], frequency, 'RESEARCH_B_COMPENSATION_JOINT_ENDPOINT')
                    full.append(field)
                    statuses.append(status)
                if all(f is not None for f in full):
                    gains = [_gain(f, target, sd) for f, target, sd in zip(full, targets, sigmas)]
                    shared = _distances(full, targets, sigmas, fit['gains'])
                    separate = _distances(full, targets, sigmas, gains)
                    weighted_f = np.concatenate([f/sd for f, sd in zip(full, sigmas)])
                    weighted_t = np.concatenate([t/sd for t, sd in zip(targets, sigmas)])
                    gshared = float(np.clip(np.dot(weighted_f, weighted_t)/
                                           max(float(np.dot(weighted_f, weighted_f)), 1e-300), .7, 1.3))
                    refit_shared = _distances(full, targets, sigmas, [gshared, gshared])
                    gcal = float(np.clip((np.dot(weighted_f, weighted_t)+1/cal_sd**2)/
                                        (np.dot(weighted_f, weighted_f)+1/cal_sd**2), .7, 1.3))
                    calibration_at_same_material = _distances(full, targets, sigmas, [gcal, gcal], cal_sd)
                    row['full_validation'] = dict(status='FULL_ENDPOINT_VALIDATED', callback_statuses=statuses,
                        latent=fit['latent'], shared=shared, gain=fit['gains'][0],
                        required_gain_old=gains[0], required_gain_new=gains[1],
                        required_gain_difference=gains[1]-gains[0],
                        shared_gain_refit_at_same_epsilon=refit_shared,
                        refitted_shared_gain=gshared,
                        independent_gain_at_same_epsilon=separate,
                        gain_conflict_squared=max(0., refit_shared['field_distance']**2-separate['field_distance']**2),
                        gain_conflict_scope='exact bounded gain profiling at one fixed full-validated epsilon; not joint epsilon global optimization',
                        calibration_at_same_epsilon=calibration_at_same_material,
                        calibration_gain_at_same_epsilon=gcal,
                        calibration_scope='same joint-shared epsilon endpoint; not reoptimized full epsilon minimum',
                        material_L1=fit['material_L1'],
                        minimum_status='full evaluation of proxy-guided feasible endpoint; not global minimum')
                    row['status'] = 'FULL_JOINT_WITNESS_WITH_PROXY_CONTROLS'
                else:
                    row['status'] = 'FULL_VALIDATION_FAILED_PROXY_RESULTS_RETAINED'
                    row['full_validation'] = dict(status='FAILED', callback_statuses=statuses)
            else:
                row['status'] = 'JOINT_PROXY_OPTIMIZATION_FAILED'
            if 'gains' in row['old_only'] and 'gains' in row['new_only']:
                row['separate_fit_gain_difference'] = float(row['new_only']['gains'][0]-row['old_only']['gains'][0])
                row['separate_fit_gain_scope'] = 'different fitted epsilon values; diagnostic, not isolated gain conflict'
            result['records'].append(row)
            result['callback_attempts_including_restored'] = verifier.attempts
            _write_outputs(directory, result, verifier)
    result['callback_attempts_including_restored'] = verifier.attempts
    result['failed_or_unresolved_callback_keys'] = verifier.unavailable
    result['full_validated_rows'] = sum(row['status'] == 'FULL_JOINT_WITNESS_WITH_PROXY_CONTROLS'
                                        for row in result['records'])
    failures = bool(verifier.unavailable) or result['full_validated_rows'] != len(result['records'])
    result['execution_status'] = 'COMPLETED_WITH_DECLARED_FAILURES' if failures else 'COMPLETED'
    result['body_wall_seconds_before_final_export'] = time.perf_counter()-start_wall
    result['process_cpu_seconds_before_final_export'] = time.process_time()-start_cpu
    result['proxy_material_frequency_evaluations_this_invocation'] = int(response.evaluations)-start_proxy_evaluations
    trials = [json.loads(line) for line in trial_path.read_text().splitlines()] if trial_path.exists() else []
    result['optimizer_trials_including_prior_invocations'] = len(trials)
    result['nonconverged_or_failed_optimizer_trials_including_prior_invocations'] = sum(not row['success'] for row in trials)
    result['cost_scope'] = 'body timer excludes final export; parent charged job includes it. Callback attempt count is not additional physics beyond the global callback ledger.'
    _write_outputs(directory, result, verifier, final=True)
    return _plain(result)

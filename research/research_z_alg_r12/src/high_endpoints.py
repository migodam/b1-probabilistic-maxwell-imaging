"""Physical continuation of archived scalar-TM endpoints, with scoring in a separate script."""
from argparse import ArgumentParser
from dataclasses import dataclass
from io import BytesIO
import json
import os
from pathlib import Path
import sys
import time
from zipfile import ZipFile

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ.setdefault(variable, '1')

import numpy as np
from scipy.linalg import svdvals
from scipy.optimize import minimize

R12 = Path(__file__).resolve().parents[1]
_LOCAL_HISTORY = R12.parents[3] / 'max_struct_r1' / 'r14_worktree'
HISTORY = Path(os.environ.get('ZALG_HISTORY', _LOCAL_HISTORY if _LOCAL_HISTORY.is_dir() else R12.parents[1])).resolve()
sys.path.insert(0, str(HISTORY))
from research.research_z_alg_r4_diversity.src.physics import Model, Regularizer

R14 = HISTORY / 'research/research_max_struct_r14'
_LOCAL_R13 = HISTORY.parent / 'r13_worktree'
R13ROOT = Path(os.environ.get('ZALG_R13_HISTORY', _LOCAL_R13 if _LOCAL_R13.is_dir() else HISTORY)).resolve()
R13 = R13ROOT / 'research/research_max_struct_r13'
PAIRED = HISTORY / 'research/research_max_struct_r12/supplement_paired'
R11 = HISTORY / 'research/research_max_struct_r11/RESULTS/exposed_endpoint_paths'


@dataclass(frozen=True)
class Endpoint:
    label: str
    case_id: str
    method: str
    input_path: Path
    material_path: Path | None = None
    archive_member: str | None = None

    def material(self, public):
        if self.material_path is None:
            return public['m_init'].copy()
        if self.archive_member is None:
            with np.load(self.material_path, allow_pickle=False) as data:
                return data['m'].copy().ravel()
        with ZipFile(self.material_path) as archive:
            return np.load(BytesIO(archive.read(self.archive_member)), allow_pickle=False)['m'].copy().ravel()


def paired(label, index, method, repeat=None):
    case_id = f'r10_twod_{index:04}'
    prefix = f'A__{case_id}'
    output = PAIRED / 'RESULTS' / f'{prefix}__{repeat}' / 'primary.npz' if repeat else None
    return Endpoint(label, case_id, method, PAIRED / 'inputs' / f'{prefix}.npz', output)


def independent(label, index, method):
    case_id = f'r14_new_{index:04}'
    input_path = R14 / 'data_public/independent24' / f'{case_id}.npz'
    member = f'RESULTS/independent24_v1/{case_id}__{method}/primary.npz'
    archive = R14 / 'archives/independent24_PRIMARY_MATERIALS.zip'
    if (R14 / member).is_file():
        return Endpoint(label, case_id, method, input_path, R14 / member)
    return Endpoint(label, case_id, method, input_path, archive, member)


def r13(label, index, method):
    case_id = f'r13_fresh_{index:04}'
    return Endpoint(label, case_id, method, R13 / 'data_public/fresh48' / f'{case_id}.npz',
                    R13 / 'RESULTS/confirmation_core48' / f'{case_id}__{method}' / 'primary.npz')


def historical_path(path):
    for base in (HISTORY, R13ROOT):
        if path.is_relative_to(base):
            return str(path.relative_to(base))
    return str(path)


ENDPOINTS = {
    'ring': paired('ring', 1, 'CIE-single'),
    'shell': paired('shell', 4, 'matrix-G19'),
    'bridge': paired('bridge', 22, 'matrix-G19'),
    'ring_ap': paired('ring_ap', 1, 'A+', 'ap'),
    'shell_dr': paired('shell_dr', 4, 'DR', 'dr'),
    'bridge_dr': paired('bridge_dr', 22, 'DR', 'dr'),
    'disconnected': independent('disconnected', 7, 'pair00'),
    'islands': independent('islands', 19, 'g19'),
    'porous': independent('porous', 23, 'dr'),
    'patchwork': independent('patchwork', 21, 'dr'),
    'bridge_ap': independent('bridge_ap', 20, 'pair01'),
    'continuous_crossgrid': independent('continuous_crossgrid', 4, 'pair01'),
    'continuous_mixed': independent('continuous_mixed', 14, 'cie'),
    'continuous_dr31': r13('continuous_dr31', 31, 'dr'),
    'continuous_dr32': r13('continuous_dr32', 32, 'dr'),
}


def load_endpoint(label):
    endpoint = ENDPOINTS[label]
    with np.load(endpoint.input_path, allow_pickle=False) as source:
        public = {key: source[key] for key in source.files}
    material = endpoint.material(public).astype(float)
    return endpoint, public, material


def bounds(public):
    return tuple(map(float, public.get('material_bounds', [public.get('lower', .02), public.get('upper', 20.)])))


def noise_radius(public):
    return float(public['delta'] if 'delta' in public else public['sigma'] * public['noise_radius'])


def physics(public):
    return Model(public), Regularizer(public['grid_shape'], smooth=1e-4, tv=1e-4, epsilon=.01)


def evaluate(model, reg, m, target, with_gradient=True):
    if model.state is None or not np.array_equal(model.state['m'], m):
        model.forward(m)
    state = model.state
    norm = float(np.linalg.norm(state['y'] - target))
    ynorm = float(np.linalg.norm(target))
    prior, g_prior, _ = reg.evaluate(m)
    data_loss = .5 * (norm / ynorm) ** 2
    result = dict(data_loss=data_loss, prior=prior, phi=data_loss + prior,
                  data_residual=norm, relative_data_residual=norm / ynorm,
                  state_residual=state['state_relative_residual'])
    if with_gradient:
        return result, model.data_gradient(state, target, ynorm**2), g_prior
    return result


def gradient_mapping(material, gradient, lo, hi):
    return material - np.clip(material - gradient, lo, hi)


def around_endpoint(model, reg, m, target, lo, hi, delta, with_jacobian=True):
    result, gd, gr = evaluate(model, reg, m, target)
    active_low = m <= lo + 1e-8
    active_high = m >= hi - 1e-8
    free = ~(active_low | active_high)
    projections = {name: gradient_mapping(m, gradient, lo, hi)
                   for name, gradient in [('phi', gd + gr), ('data', gd)]}
    result.update(residual_over_noise_radius=result['data_residual'] / delta,
                  lower_active=int(active_low.sum()), upper_active=int(active_high.sum()), free=int(free.sum()),
                  phi_gradient_inf=float(np.linalg.norm(gd + gr, np.inf)),
                  data_gradient_inf=float(np.linalg.norm(gd, np.inf)),
                  phi_projected_inf=float(np.linalg.norm(projections['phi'], np.inf)),
                  data_projected_inf=float(np.linalg.norm(projections['data'], np.inf)))
    if with_jacobian:
        jac = model.jacobian(model.state)
        singular = svdvals(jac, check_finite=False)
        result['jacobian'] = singular_metrics(singular, delta, jac.shape)
        result['jacobian_gradient_relative_gap'] = float(np.linalg.norm(
            jac.T @ model.pack(model.state['y'] - target) / np.linalg.norm(target)**2 - gd
        ) / max(np.linalg.norm(gd), np.finfo(float).tiny))
        if free.any():
            free_jacobian = jac[:, free]
            free_singular = svdvals(free_jacobian, check_finite=False)
            result['jacobian_free_box_face'] = singular_metrics(free_singular, delta, free_jacobian.shape)
    return result


def singular_metrics(singular, delta, shape):
    indices = [0, len(singular)//4, len(singular)//2, 3*len(singular)//4, len(singular) - 1]
    return dict(shape=shape, samples={str(i+1): float(singular[i]) for i in indices},
                numerical_rank=int(np.count_nonzero(singular > singular[0] * max(shape) * np.finfo(float).eps)),
                singulars_above_raw_noise_radius=int(np.count_nonzero(singular > delta)),
                raw_noise_comparison='One unit Euclidean material perturbation per right singular direction')


def continue_lbfgs(model, reg, m, target, lo, hi, delta, arm, max_evaluations):
    history = []
    accepted_iterations = 0
    best = (float('inf'), m.copy())
    calls0 = model.counts.copy()
    wall0, cpu0 = time.perf_counter(), time.process_time()

    def objective(x):
        nonlocal best
        if len(history) >= max_evaluations:
            raise StopIteration
        values, gd, gr = evaluate(model, reg, x, target)
        score, gradient = (values['phi'], gd + gr) if arm == 'phi' else (values['data_loss'], gd)
        current = dict(evaluation=len(history)+1, value=score, phi=values['phi'], data_loss=values['data_loss'],
                       residual_over_noise_radius=values['data_residual'] / delta)
        history.append(current)
        if score < best[0]:
            best = score, x.copy()
        return score * 1e6, gradient * 1e6

    def callback(_):
        nonlocal accepted_iterations
        accepted_iterations += 1

    try:
        opt = minimize(objective, m.copy(), method='L-BFGS-B', jac=True, bounds=[(lo, hi)] * m.size, callback=callback,
                       options=dict(maxiter=max_evaluations, maxfun=max_evaluations, maxls=30,
                                    ftol=1e-15, gtol=1e-5, maxcor=20))
        status, iterations = str(opt.message), int(opt.nit)
    except StopIteration:
        status, iterations = 'PHYSICAL_EVALUATION_LIMIT', accepted_iterations
    cost = dict(wall_seconds=time.perf_counter()-wall0, cpu_seconds=time.process_time()-cpu0,
                counts=dict(model.counts - calls0))
    output = best[1]
    calls1 = model.counts.copy()
    wall1, cpu1 = time.perf_counter(), time.process_time()
    diagnostics = around_endpoint(model, reg, output, target, lo, hi, delta, with_jacobian=False)
    diagnosis_cost = dict(wall_seconds=time.perf_counter()-wall1, cpu_seconds=time.process_time()-cpu1,
                          counts=dict(model.counts - calls1))
    return output, dict(status=status, accepted_iterations=iterations, evaluated_materials=len(history),
                        cost=cost, final_diagnostic_cost=diagnosis_cost, diagnostics=diagnostics, trajectory=history,
                        max_material_change=float(np.linalg.norm(output-m, np.inf)),
                        material_change_L2=float(np.linalg.norm(output-m)))


def run(label, output, limit):
    endpoint, public, start = load_endpoint(label)
    lo, hi = bounds(public)
    target, delta = public['Y'], noise_radius(public)
    model, reg = physics(public)
    start = np.clip(start, lo, hi)
    t0, c0 = time.perf_counter(), time.process_time()
    initial = around_endpoint(model, reg, start, target, lo, hi, delta)
    initial.update(cost=dict(wall_seconds=time.perf_counter()-t0, cpu_seconds=time.process_time()-c0,
                             counts=dict(model.counts)))
    arrays, results = {'initial': start}, {'initial': initial}
    for arm in ('phi', 'data'):
        arrays[arm], results[arm] = continue_lbfgs(model, reg, start, target, lo, hi, delta, arm, limit)
    geometry = dict(grid_shape=public['grid_shape'].tolist(), k=float(public['k']),
                    cell_radius=float(public['radius']), source_count=model.ns, receiver_count=model.nr,
                    sources=public['sources'].tolist(), receivers=public['receivers'].tolist(),
                    material_bounds=[lo, hi], inverse_side=int(public['grid_shape'][0]),
                    forward_side=int(public.get('forward_side', public['grid_shape'][0])))
    meta = dict(label=label, case_id=endpoint.case_id, method=endpoint.method,
                status='exposed historical development material', geometry=geometry,
                operator_file=historical_path(endpoint.input_path),
                archived_material=historical_path(endpoint.material_path or endpoint.input_path),
                material_member=endpoint.archive_member, sigma=float(public['sigma']), noise_radius=delta,
                model_noise_scope='matched' if geometry['forward_side']==geometry['inverse_side'] else 'crossgrid noise-only threshold',
                reg=dict(smooth=reg.smooth, tv=reg.tv, epsilon=reg.epsilon), max_evaluations=limit,
                arms=results, total_wall_seconds=time.perf_counter()-t0, total_cpu_seconds=time.process_time()-c0,
                total_counts=dict(model.counts))
    np.savez_compressed(output / f'{label}.npz', **arrays)
    (output / f'{label}.json').write_text(json.dumps(meta, indent=2, allow_nan=False) + '\n')
    print(label, f"{initial['residual_over_noise_radius']:.3f}",
          *(f"{results[arm]['diagnostics']['residual_over_noise_radius']:.3f}" for arm in ('phi','data')),
          f"{initial['phi_projected_inf']:.3g}", f"{initial['data_projected_inf']:.3g}",
          f"{meta['total_wall_seconds']:.2f}s", flush=True)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=list(ENDPOINTS), nargs='+')
    parser.add_argument('--limit', type=int, default=300)
    parser.add_argument('--output', type=Path, default=R12 / 'results/high_endpoints')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for label in args.case:
        run(label, args.output, args.limit)

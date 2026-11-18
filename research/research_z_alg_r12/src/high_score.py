"""Score already saved high-contrast runs and historical comparators."""
from argparse import ArgumentParser
from csv import DictWriter
from io import BytesIO
import json
from pathlib import Path
import time
from zipfile import ZipFile

import numpy as np

from high_endpoints import ENDPOINTS, HISTORY, PAIRED, R11, R12, R13, R14, bounds, historical_path
from high_endpoints import independent, load_endpoint, noise_radius, physics, r13


def json_file(path):
    return json.loads(path.read_text())


def zip_file(path, member):
    with ZipFile(path) as archive:
        return archive.read(member)


def truth_for(case_id):
    if case_id.startswith('r10_'):
        path = PAIRED/'offline_truth/A'/f'{case_id}.npz'
        with np.load(path, allow_pickle=False) as z:
            return z['m'].ravel(), z['support'].ravel().astype(bool), float(z['amplitude']), historical_path(path)
    if case_id.startswith('r13_'):
        path = R13/'sealed_truth/fresh48'/f'{case_id}.npz'
        with np.load(path, allow_pickle=False) as z:
            truth = z['material_coarse'].ravel()
        source = historical_path(path)
    else:
        path = R14/'archives/offline_truth_v1/R14_SYNTHETIC_OFFLINE_TRUTH_v1.zip'
        member = f'research/research_max_struct_r14/sealed_truth/independent24/{case_id}.npz'
        extracted = HISTORY/member
        if extracted.is_file():
            with np.load(extracted, allow_pickle=False) as z:
                truth = z['material_coarse'].ravel()
            source = historical_path(extracted)
        else:
            with np.load(BytesIO(zip_file(path, member)), allow_pickle=False) as z:
                truth = z['material_coarse'].ravel()
            source = f'{historical_path(path)}::{member}'
    return truth, truth > 1, float(truth.max()-.1), source


def material_scores(material, truth, support, amplitude):
    error = np.linalg.norm(material-truth)
    seg = material > 1.
    union = np.count_nonzero(seg | support)
    return dict(RMSE=float(error/np.sqrt(truth.size)),
                NRMSE=float(error/np.sqrt(truth.size)/amplitude),
                relative_L2=float(error/np.linalg.norm(truth)),
                IoU=float(np.count_nonzero(seg & support)/union) if union else 1.,
                IoU_threshold=1., amplitude=amplitude)


def historical_comparators(endpoint, public):
    case_id = endpoint.case_id
    if case_id.startswith('r10_'):
        for method in ('ap', 'dr', 'g19', 'direct'):
            folder = PAIRED/'RESULTS'/f'A__{case_id}__{method}'
            path = folder/'primary.npz'
            with np.load(path, allow_pickle=False) as src:
                m = src['m'].ravel()
            yield f'paired_archive_{method}', m, json_file(folder/'RECEIPT.json'), historical_path(path), 'historical same-bad-start continuation'
        method = 'matrix_green19' if case_id.endswith('0001') else 'cie_single'
        material_file = R11/f'{case_id}.npz'
        with np.load(material_file, allow_pickle=False) as src:
            material = src['m'][-1].ravel()
        r10 = HISTORY.parent/'r10_worktree' if (HISTORY.parent/'r10_worktree').exists() else HISTORY
        receipt = r10/f'research/research_max_struct_r10/RESULTS/twod/{case_id}__{method}__noise_ball_99/RECEIPT.json'
        cost = json_file(receipt) if receipt.is_file() else {}
        yield f'r10_alternate_{method}', material, cost, historical_path(material_file)+'::m[-1]', 'historical separately initialized comparator'
    elif case_id.startswith('r13_'):
        for method in ('cie', 'g19', 'dr', 'direct'):
            spec = r13(endpoint.label, int(case_id[-4:]), method)
            path = spec.material_path
            yield f'r13_archive_{method}', spec.material(public), json_file(path.parent/'RECEIPT.json'), historical_path(path), 'historical separately initialized comparator'
    else:
        review_path = R14/'archives/independent24_v1_REVIEW.zip'
        for method in ('pair01', 'cie', 'g19', 'dr', 'direct'):
            spec = independent(endpoint.label, int(case_id[-4:]), method)
            member = spec.archive_member or str(spec.material_path.relative_to(R14))
            receipt_member = str(Path(member).with_name('RECEIPT.json'))
            receipt = R14/receipt_member
            try:
                cost = json_file(receipt) if receipt.is_file() else json.loads(zip_file(review_path, receipt_member))
            except (KeyError, FileNotFoundError):
                cost = {}
            yield f'r14_archive_{method}', spec.material(public), cost, (f'{historical_path(spec.material_path)}::{spec.archive_member}' if spec.archive_member else historical_path(spec.material_path)), 'historical separately initialized comparator'


def cost_fields(receipt):
    counters = receipt.get('physical_cost', receipt.get('physics_cost', receipt.get('physical_model_cost', receipt.get('counts', {}))))
    return dict(wall_seconds=receipt.get('wall_seconds'),
                cpu_seconds=next((receipt[k] for k in ('cpu_seconds','CPU_seconds','whole_child_CPU_seconds','whole_process_CPU_seconds') if k in receipt), None),
                physical_forward=counters.get('FULL_FORWARD'), complex_RHS=counters.get('actual_complex_RHS'),
                auxiliary_cost=next((receipt[k] for k in ('auxiliary_physical_model_cost','auxiliary_cost','projection_cost','projector_cost') if k in receipt), None))


def run(labels, output):
    saved = R12/'results/high_endpoints'
    folded = R12/'results/high_fold'
    scores, materials, ordering = [], {}, {}
    for label in labels:
        _, public, _ = load_endpoint(label)
        spec = ENDPOINTS[label]
        with np.load(saved/f'{label}.npz', allow_pickle=False) as f:
            initial_materials = {k:f[k].ravel() for k in f.files}
        continuation = json_file(saved/f'{label}.json')
        if (folded/f'{label}.npz').is_file():
            with np.load(folded/f'{label}.npz', allow_pickle=False) as f:
                fold_materials = {k:f[k].ravel() for k in f.files}
            fold = json_file(folded/f'{label}.json')
        else:
            fold_materials, fold = {}, {'arms': {}}
        long_label = R12/'results/high_long'/label
        if long_label.with_suffix('.npz').is_file():
            with np.load(long_label.with_suffix('.npz'), allow_pickle=False) as f:
                long_materials = {k:f[k].ravel() for k in ('phi','data')}
            long_run = json_file(long_label.with_suffix('.json'))
        else:
            long_materials, long_run = {}, {'arms': {}}
        archived = list(historical_comparators(spec, public))
        truth, support, amplitude, gt_location = truth_for(spec.case_id)
        materials[label], ordering[label] = {'truth': truth}, []

        def save(method, value, data, work, source, role):
            row = dict(case=label, case_id=spec.case_id, archived_start_method=spec.method, method=method, role=role,
                       model_noise_scope=continuation['model_noise_scope'], gt_source=gt_location, material_source=source,
                       residual_over_noise=float(data['residual_over_noise_radius']), Phi=float(data['phi']),
                       **material_scores(value, truth, support, amplitude), **work)
            ordering[label].append(row)
            scores.append(row)
            materials[label][method] = value

        for name, m in initial_materials.items():
            arm = continuation['arms'][name]
            data = arm['diagnostics'] if name != 'initial' else arm
            work = cost_fields(arm['cost']) if name != 'initial' else dict(wall_seconds=None,cpu_seconds=None,physical_forward=None,complex_RHS=None,auxiliary_cost=None)
            save(name, m, data, work, f'high_endpoints/{label}.npz::{name}', 'archived starting material' if name=='initial' else 'same archived start L-BFGS-B continuation')
        for name, m in long_materials.items():
            arm = long_run['arms'][name]
            save('long_'+name, m, arm['diagnostics'], cost_fields(arm['cost']),
                 f'high_long/{label}.npz::{name}', 'same archived start extended L-BFGS-B continuation')
        for name, arm in fold['arms'].items():
            work = dict(wall_seconds=arm.get('matched_search_wall_seconds', arm['wall_seconds']),
                        cpu_seconds=arm.get('matched_search_cpu_seconds', arm['cpu_seconds']),
                        physical_forward=arm['matched_search_forwards'], complex_RHS=arm['matched_search_rhs'], auxiliary_cost=None)
            save(name, fold_materials[name], arm['result'], work, f'high_fold/{label}.npz::{name}', 'same archived start and declared full physical search ceiling')
        model, reg = physics(public)
        lo, hi = bounds(public)
        t0 = time.perf_counter()
        for name, m, receipt, source, role in archived:
            y = model.forward(m)['y']
            residual = np.linalg.norm(y-public['Y'])
            prior = reg.evaluate(m)[0]
            data = dict(residual_over_noise_radius=residual/noise_radius(public), phi=.5*(residual/np.linalg.norm(public['Y']))**2+prior)
            save(name, m, data, cost_fields(receipt), source, role)
        validation = dict(wall_seconds=time.perf_counter()-t0, counts=dict(model.counts))
        for row in ordering[label]:
            row['historical_comparator_new_validation_cost'] = validation

    output.mkdir(parents=True, exist_ok=True)
    (output/'scores.json').write_text(json.dumps(scores, indent=2, allow_nan=False)+'\n')
    with (output/'scores.csv').open('w') as f:
        table = DictWriter(f, fieldnames=list(scores[0]))
        table.writeheader()
        table.writerows(scores)
    figures = R12/'figures'
    figures.mkdir(exist_ok=True)
    core = [label for label in ('ring','shell','bridge','disconnected','islands','patchwork','continuous_dr31','continuous_dr32') if label in materials]
    if core:
        plot(figures, core, materials, ordering)
    return scores


def plot(figures, labels, materials, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    axes_names = ('Truth', 'Archived start', 'Physical continuation', 'Best historical comparator')
    ncols = len(axes_names)
    fig, axes = plt.subplots(len(labels), ncols, figsize=(8.7, 2.03*len(labels)), constrained_layout=True, squeeze=False)
    for i, label in enumerate(labels):
        row = {r['method']:r for r in rows[label]}
        last = min((r for r in rows[label] if r['role'].startswith('historical')), key=lambda r:r['Phi'])
        fresh = row.get('continue_fullj', row['phi'])
        methods = ('truth','initial',fresh['method'],last['method'])
        for j, (ax, method) in enumerate(zip(axes[i],methods)):
            im = ax.imshow(materials[label][method].reshape(24,24), cmap='magma', vmin=0, vmax=20, interpolation='nearest')
            ax.set_xticks([]); ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(label.replace('_',' '), fontweight='bold')
            if i == 0:
                ax.set_title(axes_names[j], fontsize=9)
            if method != 'truth':
                x = row[method]
                ax.set_xlabel(f"{x['residual_over_noise']:.2f}× noise, IoU {x['IoU']:.2f}", fontsize=8)
    fig.colorbar(im, ax=axes, shrink=.5, label='Real contrast χ')
    fig.savefig(figures/'high_materials.pdf')
    fig.savefig(figures/'high_materials.png', dpi=190)
    plt.close(fig)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('case', choices=list(ENDPOINTS), nargs='+')
    parser.add_argument('--output', type=Path, default=R12/'results/high_scores')
    args = parser.parse_args()
    run(args.case, args.output)

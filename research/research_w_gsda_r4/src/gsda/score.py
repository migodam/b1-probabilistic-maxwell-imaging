"""Post-seal array scoring and fixed-scale plots; never calls forward physics."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import label

SUPPORT_THRESHOLD = 1.0
CONTRAST_LIMITS = (0., 8.)
LOSS_LIMITS = (0., .25)
CROSS4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=int)
FULL8 = np.ones((3, 3), dtype=int)


def _rms(a):
    return float(np.sqrt(np.mean(np.abs(a) ** 2)))


def _ratio(a, b):
    b = float(b)
    return float(a / b) if b > 0 else None


def topology(support):
    support = np.asarray(support, bool)
    fg4, n4 = label(support, structure=CROSS4)
    fg8, n8 = label(support, structure=FULL8)
    bg4, nb = label(~support, structure=CROSS4)
    edges = set(np.r_[bg4[0], bg4[-1], bg4[:, 0], bg4[:, -1]]) - {0}
    holes = [k for k in range(1, nb + 1) if k not in edges]
    sizes = np.bincount(fg8.ravel())[1:]
    return dict(components4=int(n4), components8=int(n8), holes_bg4=len(holes),
                connected4=bool(n4 == 1), connected8=bool(n8 == 1),
                largest_component_fraction8=_ratio(sizes.max() if len(sizes) else 0, support.sum()),
                euler_fg8_bg4=int(n8 - len(holes))), (bg4, edges, holes)


def array_metrics(z, ztruth, x_edges, y_edges):
    z, ztruth = np.asarray(z, float), np.asarray(ztruth, float)
    if z.shape != (110,) or ztruth.shape != (110,):
        raise ValueError('Expected exactly 110 coordinates')
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(ztruth)):
        raise ValueError('Non-finite endpoint or truth')
    chi, loss = z[12:61].reshape(7, 7), z[61:110].reshape(7, 7)
    ct, lt = ztruth[12:61].reshape(7, 7), ztruth[61:110].reshape(7, 7)
    support, actual = chi >= SUPPORT_THRESHOLD, ct >= SUPPORT_THRESHOLD
    tp, fp = int(np.sum(support & actual)), int(np.sum(support & ~actual))
    fn = int(np.sum(~support & actual))
    ttop, (tlab, _, tholes) = topology(actual)
    ptop, (plab, pedges, _) = topology(support)
    representatives = []
    for k in tholes:
        points = np.argwhere(tlab == k)
        representatives.append(points[np.argmin(np.sum((points - points.mean(axis=0)) ** 2, axis=1))])
    mapped = [int(plab[tuple(p)]) for p in representatives]
    enclosed = [k for k in mapped if k > 0 and k not in pedges]
    xx, yy = np.meshgrid((np.asarray(x_edges)[:-1] + np.asarray(x_edges)[1:]) / 2,
                         (np.asarray(y_edges)[:-1] + np.asarray(y_edges)[1:]) / 2)
    def centroid(a):
        a = np.maximum(a, 0.)
        return np.array([(a * xx).sum(), (a * yy).sum()]) / a.sum() if a.sum() > 0 else None
    c, tc = centroid(chi), centroid(ct)
    delta = (chi - ct) + 1j * (loss - lt)
    pchi = np.maximum(chi, 0)
    return dict(support_threshold=SUPPORT_THRESHOLD,
        support_IoU=_ratio(tp, tp + fp + fn), support_true_cells=int(actual.sum()),
        support_predicted_cells=int(support.sum()), true_positive_cells=tp,
        false_positive_cells=fp, false_negative_cells=fn,
        false_positive_ROI_fraction=float(fp / 49),
        false_positive_background_fraction=_ratio(fp, np.sum(~actual)),
        false_positive_discovery_fraction=_ratio(fp, support.sum()),
        offsupport_positive_contrast_mass_fraction=_ratio(pchi[~actual].sum(), pchi.sum()),
        target_complex_contrast_relative_L2_f1=_ratio(np.linalg.norm(delta), np.linalg.norm(ct + 1j * lt)),
        target_real_contrast_RMSE=_rms(chi - ct), target_loss_RMSE=_rms(loss - lt),
        contrast_centroid_error=None if c is None or tc is None else float(np.linalg.norm(c - tc)),
        predicted_centroid=None if c is None else c.tolist(), truth_centroid=None if tc is None else tc.tolist(),
        wall_epsilon_RMSE=_rms(z[:6] - ztruth[:6]), wall_loss_RMSE=_rms(z[6:12] - ztruth[6:12]),
        true_topology=ttop, predicted_topology=ptop,
        components4_match=ttop['components4'] == ptop['components4'],
        components8_match=ttop['components8'] == ptop['components8'],
        holes_bg4_match=ttop['holes_bg4'] == ptop['holes_bg4'],
        true_hole_representative_components=mapped,
        all_true_holes_distinct_enclosed=(len(enclosed) == len(tholes) and len(set(enclosed)) == len(tholes))
            if tholes else None,
        topology_convention='Report both foreground4/8; holes use background4; no topology PASS gate')


def prediction_metrics(endpoint, truth, observation, sigma):
    def prediction(key_options):
        return next((np.asarray(endpoint[k]) for k in key_options if k in endpoint), None)
    data = prediction(('data', 'prediction', 'predicted_data', 'pred'))
    heldout = prediction(('heldout_data', 'heldout_prediction', 'heldout_pred'))
    out = dict(prediction_train_status='NOT_SAVED', prediction_heldout_status='NOT_RUN')
    if data is not None:
        y = np.asarray(observation)
        if data.shape != y.shape: raise ValueError('Saved train prediction shape mismatch')
        if not np.all(np.isfinite(data)): raise ValueError('Nonfinite saved prediction')
        out.update(prediction_train_status='SAVED_ACTUAL_ARRAY',
                   observation_white_RMS=_rms((data - y) / sigma) / np.sqrt(2),
                   observation_white_norm=float(np.linalg.norm((data - y) / sigma)))
        if 'clean' in truth:
            clean = np.asarray(truth['clean'])
            if data.shape != clean.shape: raise ValueError('Train clean shape mismatch')
            out['prediction_train_clean_white_RMS'] = _rms((data - clean) / sigma) / np.sqrt(2)
    if heldout is not None and 'heldout_clean' in truth:
        clean = np.asarray(truth['heldout_clean'])
        if heldout.shape != clean.shape: raise ValueError('Saved heldout prediction shape mismatch')
        if not np.all(np.isfinite(heldout)): raise ValueError('Nonfinite heldout prediction')
        out.update(prediction_heldout_status='SAVED_ACTUAL_ARRAY',
                   prediction_heldout_clean_white_RMS=_rms((heldout - clean) / sigma) / np.sqrt(2),
                   prediction_heldout_relative_L2=_ratio(np.linalg.norm(heldout - clean), np.linalg.norm(clean)))
    elif heldout is not None:
        out['prediction_heldout_status'] = 'PRIVATE_HELDOUT_TRUTH_UNAVAILABLE'
    return out


def require_method_freeze(root, case_id):
    if case_id in {f'c{i:03d}' for i in range(1, 10)}:
        root = Path(root)
        candidates = [root / name for name in ('METHOD_FREEZE', 'METHOD_FREEZE.json', 'METHOD_FREEZE.md')]
        if not any(p.is_file() for p in candidates):
            raise RuntimeError('METHOD_FREEZE required before digit truth/scoring access')


def mechanism_metrics(seal, endpoint_dir):
    """Count canonical actions once; seal/trace aliases do not create new actions."""
    events = seal.get('absorption_events')
    source = 'ONLINE_SEALED.absorption_events' if events is not None else 'NOT_RECORDED'
    if events is None:
        trace_path = Path(endpoint_dir) / 'TRACE.jsonl'
        trace_events = [json.loads(line) for line in trace_path.read_text().splitlines()] if trace_path.exists() else []
        if any(e.get('kind') in ('ABSORB', 'EXACT_ABSORB', 'REOPEN') for e in trace_events):
            events, source = trace_events, 'TRACE.jsonl fallback'
    risk_controlled = 'EMPIRICAL_RISK' in str(seal.get('method', '')) or str(seal.get('baseline_class', '')).startswith('B4_')
    actions = {}
    for index, event in enumerate(events or []):
        raw = event.get('kind')
        if raw not in ('ABSORB', 'EXACT_ABSORB', 'REOPEN'): continue
        canonical = 'ABSORB' if raw in ('ABSORB', 'EXACT_ABSORB') else 'REOPEN'
        start, iteration = event.get('start'), event.get('iteration')
        # Frozen CSI emits one operation of each kind at a given start/iteration.
        identity = (canonical, start, iteration) if start is not None or iteration is not None else (canonical, event.get('event_id', index))
        if identity in actions:
            actions[identity]['raw_kind_aliases'] = sorted(set(actions[identity]['raw_kind_aliases'] + [raw]))
            continue
        if canonical == 'ABSORB':
            mode = 'EMPIRICAL_SCREEN_TRIGGERED' if risk_controlled else 'SCHEDULED'
        elif event.get('trigger') or risk_controlled:
            mode = 'AUTOMATIC_CONFLICT_TRIGGERED'
        else:
            mode = 'SCHEDULED_FORCED'
        actions[identity] = dict(kind=canonical, raw_kind_aliases=[raw], start=start, iteration=iteration,
            action_mode=mode, trigger=event.get('trigger'), selection=event.get('selection'),
            selected_estimated_wall_cells=event.get('selected_estimated_wall_cells'))
    action_list = list(actions.values())
    def count(kind, mode=None):
        return None if events is None else sum(e['kind']==kind and (mode is None or e['action_mode']==mode) for e in action_list)
    return dict(completed_iterations=seal.get('completed_iterations'), action_counts_source=source,
        ABSORB_count=count('ABSORB'), REOPEN_count=count('REOPEN'),
        ABSORB_scheduled_count=count('ABSORB','SCHEDULED'),
        ABSORB_empirical_screen_count=count('ABSORB','EMPIRICAL_SCREEN_TRIGGERED'),
        REOPEN_scheduled_count=count('REOPEN','SCHEDULED_FORCED'),
        REOPEN_automatic_conflict_count=count('REOPEN','AUTOMATIC_CONFLICT_TRIGGERED'),
        action_events=action_list, action_alias_deduplication='canonical action + start + iteration; use one source only',
        target_risk_status='COMPUTED' if seal.get('paid_target_risk_computed') else ('NOT_COMPUTED' if 'paid_target_risk_computed' in seal else 'NOT_RECORDED'),
        paid_target_risk_screens=None if events is None else sum(bool(e.get('additional_full_J_computed')) for e in events))


def score_endpoint(root, endpoint_dir, case_id):
    """Truth is first opened after checking the saved endpoint and online seal."""
    root, endpoint_dir = Path(root), Path(endpoint_dir)
    require_method_freeze(root, case_id)
    seal_path, ep = endpoint_dir / 'ONLINE_SEALED.json', endpoint_dir / 'endpoint.npz'
    if not seal_path.exists() or not ep.exists(): raise RuntimeError('No post-seal endpoint available')
    seal = json.loads(seal_path.read_text())
    if ep.stat().st_mtime > seal_path.stat().st_mtime + .001:
        raise RuntimeError('Endpoint changed after its online seal')
    # No physics imports or solver calls occur below this private-data boundary.
    with np.load(ep, allow_pickle=False) as a:
        endpoint = {k: a[k] for k in a.files}
    with np.load(root / 'inputs' / 'private' / case_id / 'truth.npz', allow_pickle=False) as a:
        truth = {k: a[k] for k in a.files}
    if 'frequency_indices' in endpoint and list(np.asarray(endpoint['frequency_indices'], int)) != [0, 1, 2]:
        raise RuntimeError('Endpoint did not use all three primary frequencies')
    if seal.get('full_objective_comparable') is False:
        raise RuntimeError('Subset/interrupted objective is not comparable with full-band primary objective')
    heldout_file = endpoint_dir / 'heldout_prediction.npz'
    if heldout_file.exists():
        if not (endpoint_dir / 'HELDOUT_RECEIPT.json').exists():
            raise RuntimeError('Saved heldout prediction requires HELDOUT_RECEIPT.json')
        with np.load(heldout_file, allow_pickle=False) as h:
            endpoint['heldout_data'] = next((h[k] for k in ('heldout_data', 'data', 'y') if k in h.files), None)
        if endpoint['heldout_data'] is None: raise ValueError('No heldout prediction array in heldout file')
    identity = json.loads((root / 'inputs' / 'private' / case_id / 'identity.json').read_text())
    config = json.loads((root / 'inputs' / 'online' / case_id / 'config.json').read_text())
    with np.load(root / 'inputs' / 'online' / case_id / 'observation.npz', allow_pickle=False) as a:
        y = a['y']
    sigma = float(config.get('sigma_real_imag', config.get('sigma', .01)))
    if sigma <= 0: raise ValueError('Non-positive sigma')
    current_selection = {}
    currents_path = endpoint_dir / 'best_currents.npz'
    if currents_path.exists():
        with np.load(currents_path, allow_pickle=False) as currents:
            for key in ('selected_start', 'selected_iteration'):
                if key in currents: current_selection[key] = int(currents[key])
        current_selection['selected_initial_checkpoint'] = current_selection.get('selected_iteration') == 0
    target = endpoint['z'][12:110]
    current_selection['target_contrast_uniform'] = bool(np.ptp(target[:49]) == 0)
    current_selection['target_loss_uniform'] = bool(np.ptp(target[49:]) == 0)
    if current_selection['target_contrast_uniform']: current_selection['uniform_contrast'] = float(target[0])
    if current_selection['target_loss_uniform']: current_selection['uniform_loss'] = float(target[49])
    mechanism = mechanism_metrics(seal, endpoint_dir)
    metrics = array_metrics(endpoint['z'], truth['z'], truth['target_x_edges'], truth['target_y_edges'])
    metrics.update(prediction_metrics(endpoint, truth, y, sigma))
    guard_path = endpoint_dir / 'FILE_ACCESS_SUMMARY.json'
    if not guard_path.exists() and (endpoint_dir.parent / 'BATCH_RECEIPT.json').exists() and (endpoint_dir.parent / 'FILE_ACCESS_SUMMARY.json').exists():
        guard_path = endpoint_dir.parent / 'FILE_ACCESS_SUMMARY.json'
    guard = json.loads(guard_path.read_text()) if guard_path.exists() else None
    costs = json.loads((endpoint_dir / 'RECEIPT.json').read_text()) if (endpoint_dir / 'RECEIPT.json').exists() else {}
    end_to_end_path = endpoint_dir / 'END_TO_END_CPU.json'
    end_to_end_cost = json.loads(end_to_end_path.read_text()) if end_to_end_path.exists() else None
    additional_cost = None
    cost_basis = 'warm_case_partition' if costs.get('warm_case_only') else 'standalone_wrapper'
    if seal.get('not_cold_start') and seal.get('cohort') == 'DATA_EXPOSED_MEASUREMENT_ONLY_REFINEMENT':
        additional_cost = dict(costs)
        upstream = costs.get('b0_cost', seal.get('b0_cost'))
        if not upstream: raise RuntimeError('GN warm refinement requires recorded upstream B0 cost')
        costs = dict(costs, derived_cumulative_B0_plus_GN=True, upstream_cost=upstream,
            additional_cost=additional_cost, accounting='Display trajectory cost only; root receipts remain additional native jobs')
        for key in ('full_fe', 'rhs', 'wrapper_process_cpu_seconds', 'wall_seconds'):
            costs[key] = upstream[key] + additional_cost[key]
        cost_basis = 'cumulative_B0_plus_GN_DATA_EXPOSED'
    heldout_cost = json.loads((endpoint_dir / 'HELDOUT_RECEIPT.json').read_text()) if (endpoint_dir / 'HELDOUT_RECEIPT.json').exists() else {}
    relative = endpoint_dir.relative_to(root)
    campaign = relative.parts[1] if len(relative.parts) > 2 and relative.parts[0] == 'results' else endpoint_dir.parent.name
    return dict(case_id=case_id, **{k: identity[k] for k in ('cohort', 'shape', 'digit')},
        endpoint_dir=str(relative), method=seal.get('method', endpoint_dir.parent.name), campaign=campaign,
        status='POST_SEAL_ARRAY_SCORED', metrics=metrics, objective=seal.get('objective'),
        physical_data_objective=float(endpoint['data_objective']) if 'data_objective' in endpoint else .5 * metrics['observation_white_norm'] ** 2,
        stop=seal.get('stop'), converged=seal.get('converged'), online_status=seal.get('status'),
        starts=seal.get('starts'), current_selection=current_selection, mechanism=mechanism, error=seal.get('error'), cost=costs, additional_cost=additional_cost, end_to_end_cost=end_to_end_cost,
        cost_basis=cost_basis, scope=seal.get('cohort'), not_cold_start=seal.get('not_cold_start'), independent_blind_claim=seal.get('independent_blind_claim'), guard=guard,
        guard_source=str(guard_path.relative_to(root)) if guard_path.exists() else None,
        selected_scaled_projected_gradient_inf=seal.get('selected_scaled_projected_gradient_inf'),
        full_objective_comparable=seal.get('full_objective_comparable'),
        selection=seal.get('selection'), total_objective_calls=seal.get('total_objective_calls'), heldout_cost=heldout_cost, physics_calls_in_scoring=0)


def plot_digit_rows(root, records, method, output):
    """Nine actual truth/reconstruction/error rows. Missing endpoints stay blank."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root, output = Path(root), Path(output)
    digit_rows = {r['digit']: r for r in records if r.get('cohort') == 'digit' and r.get('method') == method}
    fig, axes = plt.subplots(9, 6, figsize=(15, 23), layout='constrained')
    extent = None
    nactual = 0
    for d in range(1, 10):
        row = digit_rows.get(d)
        if row is None:
            for ax in axes[d - 1]:
                ax.text(.5, .5, 'NOT_RUN / NOT_SCORED', ha='center', va='center', fontsize=7)
                ax.set_axis_off()
            axes[d - 1, 0].set_title(f'Pattern {d}: no sealed array', fontsize=9)
            continue
        epdir = root / row['endpoint_dir']
        with np.load(root / 'inputs' / 'private' / row['case_id'] / 'truth.npz', allow_pickle=False) as t:
            zt, tx, ty = t['z'], t['target_x_edges'], t['target_y_edges']
        with np.load(epdir / 'endpoint.npz', allow_pickle=False) as e:
            z = e['z']
        ct, chi = zt[12:61].reshape(7, 7), z[12:61].reshape(7, 7)
        lt, loss = zt[61:110].reshape(7, 7), z[61:110].reshape(7, 7)
        arrays = [ct, chi, chi - ct, lt, loss, loss - lt]
        bounds = [(0, 8), (0, 8), (-8, 8), (0, .25), (0, .25), (-.25, .25)]
        extent = [tx[0], tx[-1], ty[0], ty[-1]]
        for j, (arr, lim) in enumerate(zip(arrays, bounds)):
            ax = axes[d - 1, j]
            im = ax.imshow(arr, origin='lower', extent=extent, interpolation='nearest',
                           aspect='equal', vmin=lim[0], vmax=lim[1], cmap='coolwarm' if j in (2, 5) else 'viridis')
            ax.set_xticks([]); ax.set_yticks([])
            if d == 1: ax.set_title(['True contrast', 'Reconstructed contrast', 'Contrast error',
                                    'True loss', 'Reconstructed loss', 'Loss error'][j], fontsize=10)
            if j == 0:
                m = row['metrics']; ax.set_ylabel(f"Pattern {d}\nIoU {m['support_IoU']:.3f}", fontsize=9)
        nactual += 1
    # One fixed colorbar per column. Construct mappables even when row 1 is missing.
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    for j, lim in enumerate([(0, 8), (0, 8), (-8, 8), (0, .25), (0, .25), (-.25, .25)]):
        sm = ScalarMappable(norm=Normalize(*lim), cmap='coolwarm' if j in (2, 5) else 'viridis')
        fig.colorbar(sm, ax=axes[:, j], location='bottom', shrink=.85, pad=.01)
    scope_caption = '\nData-exposed B0 warm refinement; cumulative B0+GN costs; no independent blind claim' if any(r.get('not_cold_start') for r in digit_rows.values()) else ''
    fig.suptitle(f'{method}: {nactual}/9 sealed reconstructions scored from actual arrays\n'
                 'Shared contrast 0..8 / loss 0..0.25; signed errors ±8 / ±0.25; support threshold contrast ≥ 1' + scope_caption, fontsize=12)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix('.png'), dpi=160)
    fig.savefig(output.with_suffix('.pdf'))
    plt.close(fig)
    return dict(method=method, actual_reconstructed_rows=nactual, total_rows=9,
                png=str(output.with_suffix('.png')), pdf=str(output.with_suffix('.pdf')))


def plot_case(root, record):
    """Required post-seal material figures at fixed physical grid/scales."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path(root)
    require_method_freeze(root, record['case_id'])
    out = root / record['endpoint_dir']
    if not (out / 'ONLINE_SEALED.json').exists(): raise RuntimeError('No seal for private plots')
    with np.load(root / 'inputs' / 'private' / record['case_id'] / 'truth.npz', allow_pickle=False) as t:
        zt, tx, ty = t['z'], t['target_x_edges'], t['target_y_edges']
    with np.load(out / 'endpoint.npz', allow_pickle=False) as e: z = e['z']
    config = json.loads((root / 'inputs' / 'online' / record['case_id'] / 'config.json').read_text())
    wx = np.asarray(config.get('wall_x', [-1.5, -.5, .5, 1.5]))
    wy = np.asarray(config.get('wall_y', [.5, .875, 1.25]))
    title = f"{record['campaign']} / {record['method']} / {record['case_id']}"
    if record.get('scope'): title += '\n' + record['scope'] + '; B0 warm start' if record.get('not_cold_start') else '\n' + record['scope']
    files = []
    def render(name, arrays, edges, limits, labels, error=False):
        fig, axes = plt.subplots(1, len(arrays), figsize=(5 * len(arrays), 4.8), layout='constrained', squeeze=False)
        for ax, a, (x, y), lim, caption in zip(axes[0], arrays, edges, limits, labels):
            im = ax.pcolormesh(x, y, a, shading='flat', cmap='coolwarm' if error else 'viridis', vmin=lim[0], vmax=lim[1], rasterized=True)
            ax.set(xlabel='x / reference wavelength', ylabel='y / reference wavelength', title=caption)
            ax.set_aspect('equal'); fig.colorbar(im, ax=ax, shrink=.8)
        fig.suptitle(title + '\nActual cell arrays; fixed scales; no interpolation', fontsize=11)
        path = out / name; fig.savefig(path, dpi=170); plt.close(fig); files.append(str(path.relative_to(root)))
    target_edges = [(tx, ty)] * 2
    wall_edges = [(wx, wy)] * 2
    render('reconstructed_target.png', [z[12:61].reshape(7, 7), z[61:].reshape(7, 7)], target_edges,
           [(0, 8), (0, .25)], ['Reconstructed target contrast', 'Reconstructed target loss'])
    render('truth_target.png', [zt[12:61].reshape(7, 7), zt[61:].reshape(7, 7)], target_edges,
           [(0, 8), (0, .25)], ['True target contrast', 'True target loss'])
    render('reconstructed_wall.png', [z[:6].reshape(2, 3), z[6:12].reshape(2, 3)], wall_edges,
           [(1, 6), (0, .25)], ['Reconstructed wall epsilon', 'Reconstructed wall loss'])
    render('truth_wall.png', [zt[:6].reshape(2, 3), zt[6:12].reshape(2, 3)], wall_edges,
           [(1, 6), (0, .25)], ['True wall epsilon', 'True wall loss'])
    d = z - zt
    render('reconstruction_error.png', [d[12:61].reshape(7, 7), d[61:].reshape(7, 7), d[:6].reshape(2, 3), d[6:12].reshape(2, 3)],
           target_edges + wall_edges, [(-8, 8), (-.25, .25), (-5, 5), (-.25, .25)],
           ['Target contrast error', 'Target loss error', 'Wall epsilon error', 'Wall loss error'], error=True)
    return files


def plot_digit_triptych(root, records, method, campaign, output, quantity='contrast'):
    """Nine x three standalone panels; all scales fixed and cells un-interpolated."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    root, output = Path(root), Path(output)
    rows = {r['digit']: r for r in records if r.get('cohort') == 'digit' and r['method'] == method and r['campaign'] == campaign}
    fig, axes = plt.subplots(9, 3, figsize=(9, 23), layout='constrained')
    isloss = quantity == 'loss'
    material_slice = slice(61, 110) if isloss else slice(12, 61)
    vmax = .25 if isloss else 8.
    for digit in range(1, 10):
        row = rows.get(digit)
        if row is None:
            for ax in axes[digit - 1]:
                ax.text(.5, .5, 'NOT_RUN / NOT_SCORED', ha='center', va='center', fontsize=8)
                ax.set_axis_off()
            continue
        require_method_freeze(root, row['case_id'])
        with np.load(root / 'inputs' / 'private' / row['case_id'] / 'truth.npz', allow_pickle=False) as t:
            truth = t['z'][material_slice].reshape(7, 7); tx, ty = t['target_x_edges'], t['target_y_edges']
        with np.load(root / row['endpoint_dir'] / 'endpoint.npz', allow_pickle=False) as e:
            estimate = e['z'][material_slice].reshape(7, 7)
        for j, arr in enumerate((truth, estimate, estimate - truth)):
            ax = axes[digit - 1, j]
            ax.pcolormesh(tx, ty, arr, shading='flat', cmap='coolwarm' if j == 2 else 'viridis', vmin=-vmax if j == 2 else 0, vmax=vmax, rasterized=True)
            ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
            if j == 0: ax.set_ylabel(f'Pattern {digit}', fontsize=9)
    for j, caption in enumerate(('Truth', 'Reconstruction', 'Estimate minus truth')):
        axes[0, j].set_title(f'{caption} {quantity}', fontsize=10)
        fig.colorbar(ScalarMappable(norm=Normalize(-vmax if j == 2 else 0, vmax), cmap='coolwarm' if j == 2 else 'viridis'), ax=axes[:, j], location='bottom', shrink=.8)
    scope_caption = '\nData-exposed B0 warm refinement; cumulative B0+GN costs; no independent blind claim' if any(r.get('not_cold_start') for r in rows.values()) else ''
    fig.suptitle(f'{campaign} / {method}: {len(rows)}/9 actual sealed {quantity} reconstructions\n7×7 physical material cells; shared scales; no interpolation' + scope_caption, fontsize=10)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix('.png'), dpi=170); fig.savefig(output.with_suffix('.pdf')); plt.close(fig)
    return dict(campaign=campaign, method=method, quantity=quantity, actual_reconstructed_rows=len(rows), png=str(output.with_suffix('.png')), pdf=str(output.with_suffix('.pdf')))


def plot_quality_cost(root, records):
    """Final scatter plots from actual sealed endpoints and wrapper receipts only.

    Batch aggregates are audit comparators, never added to setup/warm partitions.
    Physical and CSI current-data objectives remain separately labeled samples.
    """
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path(root)
    actual = [r for r in records if r.get('cohort') == 'digit' and r.get('status') == 'POST_SEAL_ARRAY_SCORED']
    if not any('B4' in r.get('campaign', '') for r in actual):
        raise RuntimeError('Final quality/cost plots await actual scored B4 outputs')
    out = root / 'figures' / 'quality_cost'; out.mkdir(parents=True, exist_ok=True)
    rows, trace_rows, audit = [], [], []
    campaigns = sorted({r['campaign'] for r in actual})
    summaries = {}
    for campaign in campaigns:
        folder = root / 'results' / campaign
        setup_path = folder / 'PUBLIC_OPERATOR_SETUP_RECEIPT.json'
        batch_path = folder / 'BATCH_RECEIPT.json'
        setup = json.loads(setup_path.read_text()) if setup_path.exists() else {}
        batch = json.loads(batch_path.read_text()) if batch_path.exists() else {}
        group = sorted([r for r in actual if r['campaign'] == campaign], key=lambda r: r['case_id'])
        # Contract order controls the first cold case, never truth quality.
        contract_path = folder / 'FROZEN_BATCH_CONTRACT.json'
        contract = json.loads(contract_path.read_text()) if contract_path.exists() else {}
        order = contract.get('cases', [r['case_id'] for r in group])
        first = order[0] if order else None
        partition_cpu = sum(float(r.get('cost', {}).get('wrapper_process_cpu_seconds', 0)) for r in group)
        setup_cpu = setup.get('wrapper_process_cpu_seconds')
        summaries[campaign] = dict(actual_scored_cases=len(group), requested_digit_cases=9,
            actual_setup_cpu=setup_cpu, actual_case_wrapper_cpu_sum=partition_cpu,
            batch_aggregate_cpu=batch.get('wrapper_process_cpu_seconds'),
            batch_unpartitioned_cpu=batch.get('unpartitioned_process_cpu_seconds'),
            batch_aggregate_counted_in_totals=False,
            cost_scope=('Cumulative B0+GN display trajectories; shared upstream B0 receipt is not a new job or independent blind case' if any(r.get('cost_basis') == 'cumulative_B0_plus_GN_DATA_EXPOSED' for r in group) else ('Setup once plus warm partitions; aggregate is an alternative audit total' if setup else 'Independent standalone wrapper receipts')))
        audit.append(dict(campaign=campaign, **summaries[campaign]))
        for record in group:
            require_method_freeze(root, record['case_id'])
            epdir = root / record['endpoint_dir']
            if not (epdir / 'ONLINE_SEALED.json').exists(): raise RuntimeError('Unsealed quality/cost record')
            cost, metrics = record.get('cost', {}), record['metrics']
            warm = bool(cost.get('warm_case_only'))
            row = dict(campaign=campaign, method=record['method'], case_id=record['case_id'], digit=record['digit'],
                endpoint_dir=record['endpoint_dir'], status=record['online_status'],
                support_IoU=metrics['support_IoU'], complex_contrast_relative_L2=metrics['target_complex_contrast_relative_L2_f1'],
                contrast_centroid_error=metrics['contrast_centroid_error'], target_real_RMSE=metrics['target_real_contrast_RMSE'],
                target_loss_RMSE=metrics['target_loss_RMSE'], wall_epsilon_RMSE=metrics['wall_epsilon_RMSE'], wall_loss_RMSE=metrics['wall_loss_RMSE'],
                physical_full_objective_with_prior=record.get('objective'),
                full_objective_comparable=record.get('full_objective_comparable'),
                cost_basis=record.get('cost_basis', 'warm_case_partition' if warm else 'standalone_wrapper'), scope=record.get('scope'),
                actual_case_CPU=cost.get('wrapper_process_cpu_seconds'), actual_case_wall=cost.get('wall_seconds'),
                actual_case_FE=cost.get('full_fe'), actual_case_RHS=cost.get('rhs'),
                actual_shared_setup_CPU=setup_cpu, actual_shared_setup_FE=setup.get('full_fe'), actual_shared_setup_RHS=setup.get('rhs'),
                first_cold_CPU_setup_plus_first_warm=None,
                actual_setup_allocated_per_scored_case_CPU=None if setup_cpu is None else setup_cpu / len(group),
                heldout_status=metrics.get('prediction_heldout_status'),
                heldout_white_RMS=metrics.get('prediction_heldout_clean_white_RMS'),
                heldout_relative_L2=metrics.get('prediction_heldout_relative_L2'))
            if setup and record['case_id'] == first and row['actual_case_CPU'] is not None:
                row['first_cold_CPU_setup_plus_first_warm'] = setup_cpu + row['actual_case_CPU']
            with np.load(epdir / 'endpoint.npz', allow_pickle=False) as endpoint:
                row['physical_full_data_objective'] = float(endpoint['data_objective']) if 'data_objective' in endpoint else None
            rows.append(row)
            trace_path = epdir / 'TRACE.jsonl'
            if trace_path.exists():
                for index, line in enumerate(trace_path.read_text().splitlines()):
                    item = json.loads(line)
                    if item.get('kind') == 'CSI_ITERATION':
                        trace_rows.append(dict(campaign=campaign, case_id=record['case_id'], digit=record['digit'], trace_row=index,
                            start=item.get('start'), iteration=item.get('iteration'), kind='CSI_CURRENT_DATA',
                            objective=item.get('current_data_objective'), state_objective=item.get('state_objective')))
                    elif item.get('kind') == 'FULL_FEM_CHECKPOINT':
                        trace_rows.append(dict(campaign=campaign, case_id=record['case_id'], digit=record['digit'], trace_row=index,
                            start=item.get('start'), iteration=item.get('iteration'), kind='PHYSICAL_FULL_DATA',
                            objective=item.get('data_objective'), objective_with_prior=item.get('objective')))
    def write_csv(name, data):
        fields = sorted(set().union(*(r.keys() for r in data))) if data else ['status']
        with (out / name).open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(data)
    write_csv('ENDPOINT_COST_QUALITY_INPUTS.csv', rows)
    write_csv('CSI_VS_PHYSICAL_OBJECTIVE_INPUTS.csv', trace_rows)
    write_csv('COST_PARTITION_AUDIT_INPUTS.csv', audit)
    statuses = root / 'evaluation' / 'DIGITS_1_TO_9_METRICS.csv'
    if statuses.exists():
        with statuses.open(newline='') as stream:
            status_rows = [{k: r.get(k) for k in ('campaign', 'method', 'case_id', 'digit', 'status', 'online_status')} for r in csv.DictReader(stream)]
        write_csv('ACTUAL_AND_NOT_RUN_STATUS_INPUTS.csv', status_rows)
    status_caption = '; '.join(f'{c}: {summaries[c]["actual_scored_cases"]}/9 actual scored, {9-summaries[c]["actual_scored_cases"]} NOT_SCORED' for c in campaigns)
    colors = {c: plt.get_cmap('tab10')(i % 10) for i, c in enumerate(campaigns)}
    files = []
    def save(fig, name):
        fig.suptitle(name.replace('_', ' ') + '\n' + status_caption, fontsize=11)
        fig.savefig(out / f'{name}.png', dpi=170); fig.savefig(out / f'{name}.pdf'); plt.close(fig)
        files.append(str((out / f'{name}.png').relative_to(root)))
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), layout='constrained')
    for campaign in campaigns:
        rr = [r for r in rows if r['campaign'] == campaign and r['actual_case_CPU'] is not None]
        for ax, key, ylabel in zip(axes, ('complex_contrast_relative_L2', 'support_IoU', 'contrast_centroid_error'),
                                  ('Complex contrast relative L2', 'IoU at contrast ≥ 1', 'Centroid error / reference wavelength')):
            for basis in sorted({r['cost_basis'] for r in rr}):
                selected = [r for r in rr if r['cost_basis'] == basis and r[key] is not None]
                ax.scatter([r['actual_case_CPU'] for r in selected], [r[key] for r in selected], c=[colors[campaign]],
                           marker='s' if basis == 'warm_case_partition' else ('^' if basis == 'cumulative_B0_plus_GN_DATA_EXPOSED' else 'o'), label=campaign + ' / ' + basis)
                # Closely spaced warm runtimes can produce coincident labels; raw points stay exact.
                if basis != 'warm_case_partition':
                    for r in selected: ax.annotate(str(r['digit']), (r['actual_case_CPU'], r[key]), fontsize=7, xytext=(3, 3), textcoords='offset points')
            cold = [r for r in rr if r['first_cold_CPU_setup_plus_first_warm'] is not None and r[key] is not None]
            if cold: ax.scatter([r['first_cold_CPU_setup_plus_first_warm'] for r in cold], [r[key] for r in cold], marker='D', facecolors='none', edgecolors=[colors[campaign]], label=campaign + ' first cold setup+warm')
            ax.set(xlabel='Actual wrapper process CPU seconds', ylabel=ylabel); ax.grid(alpha=.2)
    axes[0].legend(fontsize=6); save(fig, 'ENDPOINT_QUALITY_VS_CPU')
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout='constrained')
    for campaign in campaigns:
        rr = [r for r in rows if r['campaign'] == campaign]
        for ax, x, y, xlabel, ylabel in ((axes[0], 'wall_epsilon_RMSE', 'complex_contrast_relative_L2', 'Wall epsilon RMSE', 'Target complex contrast relative L2'),
            (axes[1], 'wall_loss_RMSE', 'target_loss_RMSE', 'Wall loss RMSE', 'Target loss RMSE')):
            ax.scatter([r[x] for r in rr], [r[y] for r in rr], c=[colors[campaign]], label=campaign)
            for r in rr: ax.annotate(str(r['digit']), (r[x], r[y]), fontsize=7, xytext=(3, 3), textcoords='offset points')
            ax.set(xlabel=xlabel, ylabel=ylabel); ax.grid(alpha=.2)
    axes[0].legend(fontsize=7); save(fig, 'WALL_ERROR_VS_TARGET_ERROR')
    if trace_rows:
        fig, axes = plt.subplots(3, 3, figsize=(13, 11), layout='constrained')
        for digit, ax in enumerate(axes.ravel(), 1):
            selected = [r for r in trace_rows if r['digit'] == digit and r['objective'] is not None]
            for campaign in campaigns:
                for kind, marker in (('CSI_CURRENT_DATA', '.'), ('PHYSICAL_FULL_DATA', 'x')):
                    rr = [r for r in selected if r['campaign'] == campaign and r['kind'] == kind]
                    if rr: ax.scatter([r['trace_row'] for r in rr], [r['objective'] for r in rr], c=[colors[campaign]], marker=marker, s=12, label=campaign + ' ' + kind)
            if not selected: ax.text(.5, .5, 'NOT_RUN / NOT_SAVED', ha='center', va='center', transform=ax.transAxes)
            ax.set(title=f'Pattern {digit}', xlabel='Actual chronological trace row', ylabel='Data objective'); ax.grid(alpha=.2)
        handles, labels = next(((a.get_legend_handles_labels()) for a in axes.ravel() if a.get_legend_handles_labels()[0]), ([], []))
        if handles: fig.legend(handles, labels, loc='outside lower center', fontsize=6)
        save(fig, 'CSI_PSEUDODATA_VS_PHYSICAL_CHECKPOINTS')
    (out / 'COST_QUALITY_FIGURE_RECEIPT.json').write_text(json.dumps(dict(status='ACTUAL_ARRAY_SCATTER_PLOTS', files=files,
        rows=len(rows), traces=len(trace_rows), summaries=summaries, physics_calls=0,
        no_path_interpolation=True, no_truth_based_endpoint_or_method_selection=True,
        accounting='Case receipts used once; setup separated; batch aggregate audit never added to partitions'), indent=2) + '\n')
    return dict(files=files, actual_rows=len(rows), trace_rows=len(trace_rows), physics_calls=0)


def plot_heldout_cost(root, records):
    """Refresh only post-seal heldout columns/figure; fitting costs stay fixed."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path(root); out = root / 'figures' / 'quality_cost'; out.mkdir(parents=True, exist_ok=True)
    actual = [r for r in records if r.get('cohort') == 'digit' and r['metrics'].get('prediction_heldout_status') == 'SAVED_ACTUAL_ARRAY']
    rows = []
    for r in actual:
        metrics, cost, hcost = r['metrics'], r.get('cost', {}), r.get('heldout_cost', {})
        rows.append(dict(campaign=r['campaign'], method=r['method'], case_id=r['case_id'], digit=r['digit'], endpoint_dir=r['endpoint_dir'],
            heldout_status='SAVED_ACTUAL_ARRAY', heldout_white_RMS=metrics['prediction_heldout_clean_white_RMS'],
            heldout_relative_L2=metrics['prediction_heldout_relative_L2'],
            actual_fitting_wrapper_CPU=cost.get('wrapper_process_cpu_seconds'),
            fitting_cost_basis=r.get('cost_basis', 'warm_case_partition' if cost.get('warm_case_only') else 'standalone_wrapper'),
            heldout_prediction_FE=hcost.get('full_fe'), heldout_prediction_RHS=hcost.get('rhs'),
            heldout_prediction_partition_CPU=hcost.get('cpu_seconds'),
            heldout_prediction_wrapper_CPU=hcost.get('wrapper_process_cpu_seconds'),
            heldout_prediction_wall=hcost.get('wall_seconds')))
    fields = sorted(set().union(*(r.keys() for r in rows))) if rows else ['heldout_status']
    with (out / 'HELDOUT_COST_INPUTS.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(rows)
    # Refresh heldout columns in the existing endpoint CSV without redrawing old figures.
    endpoint_csv = out / 'ENDPOINT_COST_QUALITY_INPUTS.csv'
    if endpoint_csv.exists():
        with endpoint_csv.open(newline='') as stream: endpoints = list(csv.DictReader(stream))
        lookup = {r['endpoint_dir']: r for r in rows}
        for endpoint in endpoints:
            h = lookup.get(endpoint['endpoint_dir'])
            if h:
                for key in ('heldout_status', 'heldout_white_RMS', 'heldout_relative_L2', 'heldout_prediction_FE',
                            'heldout_prediction_RHS', 'heldout_prediction_partition_CPU', 'heldout_prediction_wrapper_CPU'):
                    endpoint[key] = h[key]
        fields = sorted(set().union(*(r.keys() for r in endpoints)))
        with endpoint_csv.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader(); writer.writerows(endpoints)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout='constrained')
    campaigns = sorted({r['campaign'] for r in rows})
    for i, campaign in enumerate(campaigns):
        selected = [r for r in rows if r['campaign'] == campaign and r['actual_fitting_wrapper_CPU'] is not None]
        for ax, key, ylabel in ((axes[0], 'heldout_relative_L2', 'Disjoint-frequency clean prediction relative L2'),
                                (axes[1], 'heldout_white_RMS', 'Disjoint-frequency clean prediction white RMS')):
            ax.scatter([r['actual_fitting_wrapper_CPU'] for r in selected], [r[key] for r in selected],
                       marker='s' if any(r['fitting_cost_basis']=='warm_case_partition' for r in selected) else 'o', label=campaign)
            if not any(r['fitting_cost_basis']=='warm_case_partition' for r in selected):
                for r in selected: ax.annotate(str(r['digit']), (r['actual_fitting_wrapper_CPU'], r[key]), fontsize=7, xytext=(3, 3), textcoords='offset points')
            ax.set(xlabel='Fitting CPU seconds (GN includes B0; B4 warm part)', ylabel=ylabel); ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    caption = '; '.join(f'{c}: {sum(r["campaign"]==c for r in rows)}/9 actual heldout arrays' for c in campaigns)
    fig.suptitle('Post-seal heldout f=0.925; no online use or endpoint reselection\n' + caption + '\nShared setup and heldout validation wrapper costs remain separate', fontsize=10)
    figure = out / 'HELDOUT_ERROR_VS_FITTING_CPU'
    fig.savefig(figure.with_suffix('.png'), dpi=170); fig.savefig(figure.with_suffix('.pdf')); plt.close(fig)
    validation_receipts = {}
    for shared_path in (root / 'results/heldout_predictions').glob('HELDOUT_PREDICTIONS_RECEIPT*.json'):
        payload = json.loads(shared_path.read_text()); receipt = payload.get('receipt')
        if receipt:
            validation_receipts[(receipt.get('job'), receipt.get('pid'))] = dict(receipt=receipt, recorded_prediction_rows=len(payload.get('rows', [])), source=str(shared_path.relative_to(root)))
    shared = list(validation_receipts.values())
    summary = dict(actual_heldout_rows=len(rows), figure=str(figure.with_suffix('.png').relative_to(root)),
        physics_calls=0, endpoint_reselection=False, fitting_costs_unchanged=True,
        shared_validation_wrapper_receipts=shared,
        accounting='Shared validation wrapper is authoritative; per-case prediction partitions are diagnostics and are not added again',
        per_campaign={c:dict(rows=sum(r['campaign']==c for r in rows),
            mean_white_RMS=float(np.mean([r['heldout_white_RMS'] for r in rows if r['campaign']==c])),
            mean_relative_L2=float(np.mean([r['heldout_relative_L2'] for r in rows if r['campaign']==c]))) for c in campaigns})
    (out / 'HELDOUT_REFRESH_RECEIPT.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


def plot_all_method_comparison(root, records):
    """Fixed truth|B0|B0+FullGN|B4 columns; all methods retained for every digit."""
    import csv
    import shutil
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    root = Path(root)
    campaigns = ('B0_FULL_FWI', 'FULL_GN_TV_SUPPLEMENT', 'B4_GSDA')
    indexed = {}
    for campaign in campaigns:
        rr = [r for r in records if r.get('cohort') == 'digit' and r['campaign'] == campaign]
        if len(rr) != 9 or len({r['digit'] for r in rr}) != 9:
            raise RuntimeError(f'Comparison requires exactly nine actual sealed rows for {campaign}')
        indexed[campaign] = {r['digit']: r for r in rr}
    arrays, inputs = {}, []
    for digit in range(1, 10):
        case = indexed[campaigns[0]][digit]['case_id']
        require_method_freeze(root, case)
        if any(indexed[c][digit]['case_id'] != case for c in campaigns): raise ValueError('Comparison case mismatch')
        with np.load(root / 'inputs/private' / case / 'truth.npz', allow_pickle=False) as t:
            truth, tx, ty = t['z'].copy(), t['target_x_edges'].copy(), t['target_y_edges'].copy()
        materials = [('Truth', truth)]
        for campaign in campaigns:
            record = indexed[campaign][digit]; epdir = root / record['endpoint_dir']
            if not (epdir / 'ONLINE_SEALED.json').exists(): raise RuntimeError('Comparison endpoint is not sealed')
            with np.load(epdir / 'endpoint.npz', allow_pickle=False) as endpoint:
                materials.append((campaign, endpoint['z'].copy()))
        arrays[digit] = (materials, tx, ty)
        for name, z in materials:
            contrast, loss = z[12:61].reshape(7, 7), z[61:].reshape(7, 7)
            for iy in range(7):
                for ix in range(7):
                    inputs.append(dict(digit=digit, case_id=case, method=name, row_ascending_y=iy, column=ix,
                        x_center=float((tx[ix]+tx[ix+1])/2), y_center=float((ty[iy]+ty[iy+1])/2),
                        contrast=float(contrast[iy,ix]), loss=float(loss[iy,ix])))
    input_csv = root / 'figures' / 'DIGITS_COMPARISON_ACTUAL_CELL_INPUTS.csv'
    with input_csv.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(inputs[0])); writer.writeheader(); writer.writerows(inputs)
    files = []
    for quantity, section, maximum in (('CONTRAST', slice(12,61), 8.), ('LOSS', slice(61,110), .25)):
        fig, axes = plt.subplots(9, 4, figsize=(12, 23), layout='constrained')
        for digit in range(1, 10):
            materials, tx, ty = arrays[digit]
            for j, (_, z) in enumerate(materials):
                ax = axes[digit-1,j]
                ax.pcolormesh(tx, ty, z[section].reshape(7, 7), shading='flat', cmap='viridis', vmin=0, vmax=maximum, rasterized=True)
                ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
                if j == 0: ax.set_ylabel(f'Pattern {digit}', fontsize=10)
        for j, title in enumerate(('Truth', 'B0 full FWI', 'B0 + FullGN\nData-exposed warm refinement', 'B4 GSDA')):
            axes[0,j].set_title(title, fontsize=10)
        fig.colorbar(ScalarMappable(norm=Normalize(0,maximum), cmap='viridis'), ax=axes.ravel().tolist(),
                     location='bottom', shrink=.8, label=quantity.lower() + f' (shared 0..{maximum:g})')
        fig.suptitle(f'All nine actual sealed reconstructions: {quantity.lower()}\nFixed columns; every method shown; 7×7 physical material cells without interpolation\nFullGN is data-exposed B0 warm refinement; its comparison costs include B0+GN', fontsize=11)
        output = root / f'DIGITS_1_TO_9_{quantity}'
        fig.savefig(output.with_suffix('.png'), dpi=170); fig.savefig(output.with_suffix('.pdf')); plt.close(fig)
        files.extend([output.with_suffix('.png').name, output.with_suffix('.pdf').name])
    for extension in ('png','pdf'):
        shutil.copyfile(root / f'DIGITS_1_TO_9_CONTRAST.{extension}', root / f'DIGITS_1_TO_9_RECONSTRUCTION.{extension}')
        files.append(f'DIGITS_1_TO_9_RECONSTRUCTION.{extension}')
    receipt = dict(status='27_ACTUAL_SEALED_METHOD_ENDPOINTS', distinct_digits=9, method_endpoint_rows=27,
        column_order=['Truth', 'B0 full FWI', 'B0+FullGN DATA_EXPOSED warm refinement', 'B4 GSDA'],
        files=files, actual_cell_csv=str(input_csv.relative_to(root)), cell_rows=len(inputs),
        contrast_scale=[0,8], loss_scale=[0,.25], physical_target_shape=[7,7],
        interpolation=False, display_does_not_add_resolution=True, truth_based_method_selection=False, physics_calls=0)
    (root / 'figures' / 'DIGITS_COMPARISON_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    return receipt


def plot_physics_evidence(root):
    """Compact existing E1/E2/E3 data figures; imports no forward solver."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(root); out=root/'figures/physics'; out.mkdir(parents=True,exist_ok=True)
    files=[]
    def table(name, rows):
        fields=sorted(set().union(*(r.keys() for r in rows))) if rows else ['status']
        with (out/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    def save(fig,name,title):
        fig.suptitle(title,fontsize=11)
        fig.savefig(out/f'{name}.png',dpi=170);fig.savefig(out/f'{name}.pdf');plt.close(fig)
        files.append(str((out/f'{name}.png').relative_to(root)))
    absorption=[]; born=[]
    for case in ('c010','c011','c012'):
        folder=root/'results/absorption_controls'/case
        summary=json.loads((folder/'CONTROL_RESULT.json').read_text())
        with np.load(folder/'full_state.npz',allow_pickle=False) as a:
            full=a['assembled_FEM_data'];native=a['data']
        data={'Native current':native}
        for name,filename in (('Absorbed','absorbed_state.npz'),('Reopened','reopened_state.npz')):
            with np.load(folder/filename,allow_pickle=False) as a:data[name]=a['data']
        with np.load(folder/'allwall_born_comparison.npz',allow_pickle=False) as a:
            data['Born']=a['wall_only_data']+a['Born_target_induced']
        for method,value in data.items():
            absorption.append(dict(case=case,method=method,relative_measurement_error=float(np.linalg.norm(value-full)/np.linalg.norm(full)),
                white_L2_error=float(np.linalg.norm((value-full)/.01)),scope='Offline truth-exposed mechanism control; not online imaging success'))
        born.append(dict(case=case,**summary['allwall_Born_comparison']))
    table('ABSORPTION_ACTUAL_INPUTS.csv',absorption);table('BORN_ACTUAL_INPUTS.csv',born)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    labels=('c010 ellipse','c011 nonconvex/air cell','c012 strong holed')
    for method,marker in (('Native current','o'),('Absorbed','s'),('Reopened','+'),('Born','^')):
        selected=[r for r in absorption if r['method']==method]
        axes[0].scatter(range(3),[r['relative_measurement_error'] for r in selected],marker=marker,label=method)
    axes[0].set(yscale='log',ylabel='Relative measurement discrepancy to assembled full FEM',xticks=range(3),xticklabels=labels)
    axes[0].legend(fontsize=7);axes[0].grid(alpha=.2)
    axes[1].bar(range(3),[r['whitened_RMS_independent_real_imag'] for r in born],color='#b85142')
    axes[1].set(ylabel='All-wall Born whitened RMS error',xticks=range(3),xticklabels=labels)
    for i,r in enumerate(born):axes[1].annotate(f"{r['whitened_RMS_independent_real_imag']:.3g}",(i,r['whitened_RMS_independent_real_imag']),xytext=(0,4),textcoords='offset points',ha='center')
    save(fig,'ABSORPTION_REOPEN_BORN','Actual offline full/absorbed/reopened/Born arrays; no target-recovery gate')
    spectral=json.loads((root/'results/spectral/SUMMARY.json').read_text())
    rows=spectral['rows'] if 'rows' in spectral else spectral['records']
    with np.load(root/'results/spectral/GD_MODES.npz',allow_pickle=False) as a:singular=a['singular_values']
    table('GD_MODE_INPUTS.csv',[dict(mode=i+1,singular_value=float(v)) for i,v in enumerate(singular)])
    sr=[]
    for r in rows:
        sr.append(dict(case=r['case'],kind='Kept projected48',strong_rank=None,
            white_error_to_fullFEM=r['projected_white_error_to_fullFEM'],schur_relative_to_projected=None,
            full_field_relative_error=r['projected_field_relative_error'], current_equation_relative_residual=r.get('current_equation_relative_residual')))
        for branch in r['schur_ablation']:
            sr.append(dict(case=r['case'],kind='Dropped feedback',strong_rank=branch['strong_rank'],
                white_error_to_fullFEM=branch['dropped_feedback_white_error_to_fullFEM'],
                schur_relative_to_projected=branch['schur_relative_to_full_projected']))
    table('SPECTRAL_ACTUAL_INPUTS.csv',sr)
    fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
    axes[0].scatter(np.arange(1,len(singular)+1),singular,s=15);axes[0].set(xlabel='Retained actual GD mode index',ylabel='Actual weighted GD singular value',yscale='log')
    positions=np.arange(len(sr));axes[1].scatter(positions,[r['white_error_to_fullFEM'] for r in sr],c=['#2b748e' if r['kind']=='Kept projected48' else '#b85142' for r in sr])
    captions=[r['case']+' '+('keep48' if r['strong_rank'] is None else 'drop'+str(r['strong_rank'])) for r in sr]
    axes[1].set(xticks=positions,xticklabels=captions,ylabel='Whitened L2 discrepancy to full FEM');axes[1].tick_params(axis='x',rotation=65,labelsize=7)
    field_rows=[r for r in sr if r['kind']=='Kept projected48']
    axes[2].bar(range(len(field_rows)),[r['full_field_relative_error'] for r in field_rows],color='#2b748e')
    axes[2].set(ylabel='48-mode total-field relative error to full FEM',xticks=range(len(field_rows)),xticklabels=[r['case'] for r in field_rows])
    for i,r in enumerate(field_rows):axes[2].annotate(f"{r['full_field_relative_error']:.3f}",(i,r['full_field_relative_error']),ha='center',xytext=(0,4),textcoords='offset points',fontsize=8)
    for ax in axes:ax.grid(alpha=.2)
    save(fig,'GD_MODES_AND_FEEDBACK','Actual GD modes and feedback ablations; full-FEM data and total-field discrepancies')
    mesh=json.loads((root/'results/independent_validation/SUMMARY.json').read_text())['rows']
    mr=[{k:r.get(k) for k in ('case','truth_white_mesh_distance','endpoint_white_mesh_distance','original_W_tolerance','within_original_tolerance','scope')} for r in mesh]
    table('MESH_VALIDATION_ACTUAL_INPUTS.csv',mr)
    fig,ax=plt.subplots(figsize=(8.5,4.8),layout='constrained')
    ax.bar(range(len(mr)),[r['truth_white_mesh_distance'] for r in mr],color=['#2b748e' if r['within_original_tolerance'] else '#b85142' for r in mr],label='Truth-array coarse/fine distance')
    for i,r in enumerate(mr):
        ax.annotate(f"{r['truth_white_mesh_distance']:.4f}",(i,r['truth_white_mesh_distance']),ha='center',xytext=(0,4),textcoords='offset points',fontsize=8)
        if r['endpoint_white_mesh_distance'] is not None:ax.scatter(i,r['endpoint_white_mesh_distance'],marker='x',s=60,color='black',label='c006 endpoint coarse/fine distance')
    threshold=mr[0]['original_W_tolerance'];ax.axhline(threshold,color='black',ls='--',lw=1,label=f'Original tolerance {threshold:g}')
    ax.set(xticks=range(len(mr)),xticklabels=[r['case'] for r in mr],ylabel='Whitened coarse/fine measurement L2 distance')
    ax.legend(fontsize=7);ax.grid(axis='y',alpha=.2)
    save(fig,'INDEPENDENT_MESH_CHECK','Actual further FEM refinement, including failures above 0.1; no continuum error bound')
    costs=[]
    for name in ('absorption_controls','spectral','independent_validation'):
        receipt=json.loads((root/'results'/name/'RECEIPT.json').read_text())
        costs.append(dict(job=name,**{k:receipt.get(k) for k in ('wrapper_process_cpu_seconds','wall_seconds','peak_rss_bytes','full_fe','rhs')}))
    table('PHYSICS_JOB_COST_INPUTS.csv',costs)
    fig,axes=plt.subplots(1,3,figsize=(12,4.5),layout='constrained')
    for ax,key,ylabel in zip(axes,('wrapper_process_cpu_seconds','full_fe','rhs'),('Wrapper process CPU seconds','Full FE reservations','Complex RHS reservations')):
        ax.bar(range(len(costs)),[r[key] for r in costs]);ax.set(xticks=range(len(costs)),xticklabels=[r['job'].replace('_','\n') for r in costs],ylabel=ylabel);ax.grid(axis='y',alpha=.2)
    save(fig,'PHYSICS_JOB_COSTS','Existing native job wrapper costs; diagnostics shown once, not added again to campaign accounting')
    result=dict(files=files,physics_calls=0,source_jobs=[r['job'] for r in costs],inputs='CSV files beside figures',scope='Existing measured data only; no interpolated curves or new scientific gates')
    (out/'PHYSICS_FIGURE_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def plot_native_controls(root, records, include_b0=False):
    """Actual paired native control matrix; cost partitions and actions preserved."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    root=Path(root); out=root/('figures/paired_controls_with_B0' if include_b0 else 'figures/native_controls');out.mkdir(parents=True,exist_ok=True)
    families=(('B0',) if include_b0 else ())+('B1','B1_RESET','B2','B3','B4','B5');cases=('c010','c011','c012')
    def campaign(family): return 'B0_PAIRED_CONTROLS' if family=='B0' else 'CONTROL_'+family
    lookup={(r['campaign'],r['case_id']):r for r in records if r.get('cohort')=='mechanism_control'}
    if any((campaign(family),case) not in lookup for family in families for case in cases):
        raise RuntimeError('Paired native figure requires all 18 actual scored control endpoints')
    costs=[];paired=[];cells=[];arrays={}
    for family in families:
        folder=root/'results'/campaign(family)
        if family=='B0':
            costs.append(dict(family=family,setup_CPU=None,setup_FE=None,setup_RHS=None,aggregate_CPU=None,aggregate_FE=None,aggregate_RHS=None,unpartitioned_CPU=None,aggregate_is_alternative_total=False,cost_scope='Three independent cold B0 child jobs; no shared setup'))
            continue
        setup=json.loads((folder/'PUBLIC_OPERATOR_SETUP_RECEIPT.json').read_text())
        aggregate=json.loads((folder/'BATCH_RECEIPT.json').read_text())
        costs.append(dict(family=family,setup_CPU=setup['wrapper_process_cpu_seconds'],setup_FE=setup['full_fe'],setup_RHS=setup['rhs'],
            aggregate_CPU=aggregate['wrapper_process_cpu_seconds'],aggregate_FE=aggregate['full_fe'],aggregate_RHS=aggregate['rhs'],
            unpartitioned_CPU=aggregate.get('unpartitioned_process_cpu_seconds'),aggregate_is_alternative_total=True))
    for case in cases:
        with np.load(root/'inputs/private'/case/'truth.npz',allow_pickle=False) as t:
            truth,tx,ty=t['z'].copy(),t['target_x_edges'].copy(),t['target_y_edges'].copy()
        materials=[('Truth',truth)]
        for family in families:
            r=lookup[campaign(family),case];epdir=root/r['endpoint_dir'];m=r['metrics'];cost=r['cost']
            seal=json.loads((epdir/'ONLINE_SEALED.json').read_text())
            with np.load(epdir/'endpoint.npz',allow_pickle=False) as endpoint:
                z=endpoint['z'].copy();data_objective=float(endpoint['data_objective']) if 'data_objective' in endpoint else None
            materials.append((family,z))
            paired.append(dict(case_id=case,family=family,method=r['method'],physical_full_objective_with_prior=r['objective'],physical_data_objective=data_objective,
                support_IoU=m['support_IoU'],complex_contrast_relative_L2=m['target_complex_contrast_relative_L2_f1'],centroid_error=m['contrast_centroid_error'],
                target_loss_RMSE=m['target_loss_RMSE'],wall_epsilon_RMSE=m['wall_epsilon_RMSE'],wall_loss_RMSE=m['wall_loss_RMSE'],
                case_fitting_CPU=cost['wrapper_process_cpu_seconds'],case_fitting_FE=cost['full_fe'],case_fitting_RHS=cost['rhs'],
                warm_case_CPU=cost['wrapper_process_cpu_seconds'] if cost.get('warm_case_only') else None,
                warm_case_FE=cost['full_fe'] if cost.get('warm_case_only') else None,warm_case_RHS=cost['rhs'] if cost.get('warm_case_only') else None,
                cost_basis=r.get('cost_basis'),whole_child_CPU=(r.get('end_to_end_cost') or {}).get('child_total_cpu_seconds'),
                whole_child_wall=(r.get('end_to_end_cost') or {}).get('external_wall_seconds'),whole_child_not_additive_to_Book=True,scope=r.get('scope'),
                case_fitting_wall=cost['wall_seconds'],warm_case_wall=cost['wall_seconds'] if cost.get('warm_case_only') else None,shared_setup_included_in_case=False,
                ABSORB_count=r['mechanism']['ABSORB_count'],REOPEN_count=r['mechanism']['REOPEN_count'],
                scheduled_absorptions=r['mechanism'].get('ABSORB_scheduled_count'),scheduled_reopens=r['mechanism'].get('REOPEN_scheduled_count'),
                automatic_reopens=r['mechanism'].get('REOPEN_automatic_conflict_count'),physical_reset_count=len(seal.get('physical_reset_events',[])),
                target_risk_status=r['mechanism']['target_risk_status'],heldout_status=m['prediction_heldout_status'],
                heldout_white_RMS=m.get('prediction_heldout_clean_white_RMS'),heldout_relative_L2=m.get('prediction_heldout_relative_L2'),
                guard_denied=len(r['guard']['denied']) if r.get('guard') else None,guard_source=r.get('guard_source')))
        arrays[case]=(materials,tx,ty)
        for method,z in materials:
            chi,loss=z[12:61].reshape(7,7),z[61:].reshape(7,7)
            for iy in range(7):
                for ix in range(7):cells.append(dict(case_id=case,method=method,row_ascending_y=iy,column=ix,contrast=float(chi[iy,ix]),loss=float(loss[iy,ix])))
    def table(name,rows):
        fields=sorted(set().union(*(r.keys() for r in rows)))
        with (out/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    table('PAIRED_FULL_OBJECTIVE_AND_COST.csv',paired);table('METHOD_SETUP_AND_AGGREGATE_COSTS.csv',costs);table('ACTUAL_CELL_INPUTS.csv',cells)
    files=[]
    for quantity,section,maximum in (('CONTRAST',slice(12,61),8.),('LOSS',slice(61,110),.25)):
        fig,axes=plt.subplots(3,len(families)+1,figsize=(18 if include_b0 else 16,8),layout='constrained')
        for i,case in enumerate(cases):
            materials,tx,ty=arrays[case]
            for j,(method,z) in enumerate(materials):
                ax=axes[i,j];ax.pcolormesh(tx,ty,z[section].reshape(7,7),shading='flat',cmap='viridis',vmin=0,vmax=maximum,rasterized=True)
                ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
                if i==0:ax.set_title(method.replace('_',' '),fontsize=11)
                if j==0:ax.set_ylabel({'c010':'c010 ellipse','c011':'c011 nonconvex\nwall air cell','c012':'c012 strong holed'}[case],fontsize=10)
        fig.colorbar(ScalarMappable(norm=Normalize(0,maximum),cmap='viridis'),ax=axes.ravel().tolist(),location='bottom',shrink=.8,label=f'{quantity.lower()} shared 0..{maximum:g}')
        fig.suptitle('Actual paired controls: Truth | '+('B0 | ' if include_b0 else '')+'B1 | B1 reset | B2 | B3 | B4 | B5\nSame physical worlds; 7×7 cells without interpolation; B0 is a frozen-method control extension',fontsize=11)
        path=out/f'NATIVE_CONTROLS_{quantity}';fig.savefig(path.with_suffix('.png'),dpi=170);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)
        files.append(str(path.with_suffix('.png').relative_to(root)))
    fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained')
    for i,case in enumerate(cases):
        rows=[next(r for r in paired if r['case_id']==case and r['family']==family) for family in families]
        for ax,key,label in ((axes[0],'physical_full_objective_with_prior','Full physical measured objective + frozen prior'),
                             (axes[1],'case_fitting_CPU','Fitting Book CPU seconds (B0 cold; CSI warm parts)'),
                             (axes[2],'complex_contrast_relative_L2','Target complex contrast relative L2')):
            ax.scatter(range(len(families)),[r[key] for r in rows],label=case)
            ax.set(xticks=range(len(families)),xticklabels=[f.replace('_',' ') for f in families],ylabel=label);ax.tick_params(axis='x',rotation=35);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);axes[0].set_yscale('log')
    fig.suptitle('Actual paired control endpoints; fixed measured-objective selection\nPublic setup charged once per method in separate CSV; aggregate is not added to partitions',fontsize=11)
    path=out/'NATIVE_CONTROLS_PAIRED_OBJECTIVE_COST';fig.savefig(path.with_suffix('.png'),dpi=170);fig.savefig(path.with_suffix('.pdf'));plt.close(fig);files.append(str(path.with_suffix('.png').relative_to(root)))
    result=dict(actual_native_control_endpoints=3*len(families),case_order=list(cases),method_order=list(families),files=files,physics_calls=0,
        cost_accounting='Warm-case rows plus setup once; aggregate alternative audit total, not added again',fixed_contrast_scale=[0,8],fixed_loss_scale=[0,.25],interpolation=False)
    (out/'NATIVE_CONTROL_FIGURE_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n');return result


def plot_native_heldout(root, records):
    """Refresh native-control heldout data only; retain already-viewed material PNGs."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(root);out=root/'figures/native_controls';out.mkdir(parents=True,exist_ok=True)
    families=('B1','B1_RESET','B2','B3','B4','B5');cases=('c010','c011','c012')
    selected=[r for r in records if r.get('campaign') in {'CONTROL_'+f for f in families}]
    rows=[]
    for r in selected:
        m=r['metrics'];hcost=r.get('heldout_cost',{})
        rows.append(dict(family=r['campaign'][8:],case_id=r['case_id'],heldout_status=m['prediction_heldout_status'],
            heldout_white_RMS=m.get('prediction_heldout_clean_white_RMS'),heldout_relative_L2=m.get('prediction_heldout_relative_L2'),
            heldout_FE=hcost.get('full_fe'),heldout_RHS=hcost.get('rhs'),heldout_partition_CPU=hcost.get('cpu_seconds'),
            warm_fitting_CPU=r['cost']['wrapper_process_cpu_seconds']))
    def write(path,values):
        fields=sorted(set().union(*(r.keys() for r in values)))
        with path.open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(values)
    write(out/'NATIVE_CONTROL_HELDOUT_INPUTS.csv',rows)
    paired_path=out/'PAIRED_FULL_OBJECTIVE_AND_COST.csv'
    if paired_path.exists():
        with paired_path.open(newline='') as stream:paired=list(csv.DictReader(stream))
        lookup={(r['family'],r['case_id']):r for r in rows}
        for row in paired:
            update=lookup.get((row['family'],row['case_id']))
            if update:row.update(update)
        write(paired_path,paired)
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    for case in cases:
        rr=[next((r for r in rows if r['family']==family and r['case_id']==case),None) for family in families]
        for ax,key,label in ((axes[0],'heldout_relative_L2','Heldout clean prediction relative L2'),(axes[1],'heldout_white_RMS','Heldout clean prediction whitened RMS')):
            good=[(i,r) for i,r in enumerate(rr) if r and r['heldout_status']=='SAVED_ACTUAL_ARRAY']
            ax.scatter([i for i,r in good],[r[key] for i,r in good],label=case)
            ax.set(xticks=range(6),xticklabels=[f.replace('_',' ') for f in families],ylabel=label);ax.tick_params(axis='x',rotation=30);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    actual=sum(r['heldout_status']=='SAVED_ACTUAL_ARRAY' for r in rows)
    fig.suptitle(f'Actual post-seal native-control prediction at f=0.925: {actual}/18 saved arrays\nPrediction validation only; no online use, endpoint reselection or new image resolution',fontsize=11)
    path=out/'NATIVE_CONTROL_HELDOUT';fig.savefig(path.with_suffix('.png'),dpi=170);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)
    summary=dict(actual_heldout_rows=actual,figure=str(path.with_suffix('.png').relative_to(root)),physics_calls=0,
        prior_material_figures_unchanged=True,endpoint_reselection=False,
        accounting='Fitting setup/parts retained separately; per-case validation CPU partitions are not added to the parent wrapper again',
        per_family={f:dict(rows=sum(r['family']==f and r['heldout_status']=='SAVED_ACTUAL_ARRAY' for r in rows),
            mean_relative_L2=float(np.mean([r['heldout_relative_L2'] for r in rows if r['family']==f and r['heldout_relative_L2'] is not None])),
            mean_white_RMS=float(np.mean([r['heldout_white_RMS'] for r in rows if r['family']==f and r['heldout_white_RMS'] is not None]))) for f in families})
    (out/'NATIVE_HELDOUT_REFRESH_RECEIPT.json').write_text(json.dumps(summary,indent=2)+'\n');return summary


def plot_tv_ablation(root, records):
    """Only the actual sealed same-data TV=.2/0 digit8 endpoints are compared."""
    import csv
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(root);case='c006';out=root/'figures/tv_ablation';out.mkdir(parents=True,exist_ok=True)
    campaigns=('B0_TV0_ABLATION','B0_FULL_FWI')
    chosen=[]
    for campaign in campaigns:
        rr=[r for r in records if r.get('campaign')==campaign and r['case_id']==case]
        if len(rr)!=1:raise RuntimeError('TV comparison requires both actual sealed c006 endpoints')
        chosen.append(rr[0])
    with np.load(root/'inputs/private'/case/'truth.npz',allow_pickle=False) as t:
        truth,tx,ty=t['z'].copy(),t['target_x_edges'].copy(),t['target_y_edges'].copy()
    arrays=[truth];rows=[]
    for r in chosen:
        epdir=root/r['endpoint_dir']
        if not (epdir/'ONLINE_SEALED.json').exists():raise RuntimeError('TV endpoint is unsealed')
        with np.load(epdir/'endpoint.npz',allow_pickle=False) as endpoint:
            arrays.append(endpoint['z'].copy());data_objective=float(endpoint['data_objective']) if 'data_objective' in endpoint else None
        m=r['metrics'];c=r['cost'];whole=r.get('end_to_end_cost') or {}
        rows.append(dict(campaign=r['campaign'],scope=r.get('scope'),case_id=case,digit=8,
            full_objective_with_different_TV_prior=r['objective'],common_physical_data_objective=data_objective,
            support_IoU=m['support_IoU'],complex_contrast_relative_L2=m['target_complex_contrast_relative_L2_f1'],
            real_RMSE=m['target_real_contrast_RMSE'],loss_RMSE=m['target_loss_RMSE'],centroid_error=m['contrast_centroid_error'],
            train_white_RMS=m['observation_white_RMS'],true_holes=m['true_topology']['holes_bg4'],predicted_holes=m['predicted_topology']['holes_bg4'],
            true_components8=m['true_topology']['components8'],predicted_components8=m['predicted_topology']['components8'],
            false_positive_cells=m['false_positive_cells'],offsupport_contrast_mass_fraction=m['offsupport_positive_contrast_mass_fraction'],
            true_holes_distinct=m['all_true_holes_distinct_enclosed'],heldout_status=m['prediction_heldout_status'],
            heldout_white_RMS=m.get('prediction_heldout_clean_white_RMS'),heldout_relative_L2=m.get('prediction_heldout_relative_L2'),
            full_fe=c['full_fe'],rhs=c['rhs'],book_CPU=c['wrapper_process_cpu_seconds'],whole_child_CPU=whole.get('child_total_cpu_seconds'),
            whole_child_not_additive_to_Book=True,selected_scaled_projected_gradient_inf=r.get('selected_scaled_projected_gradient_inf')))
    fields=sorted(set().union(*(r.keys() for r in rows)))
    with (out/'TV0_VS_TV02_ACTUAL_INPUTS.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    fig,axes=plt.subplots(2,3,figsize=(10,7),layout='constrained')
    for i,(section,maximum,label) in enumerate(((slice(12,61),8.,'contrast'),(slice(61,110),.25,'loss'))):
        for j,z in enumerate(arrays):
            ax=axes[i,j];im=ax.pcolormesh(tx,ty,z[section].reshape(7,7),shading='flat',cmap='viridis',vmin=0,vmax=maximum,rasterized=True)
            ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
            if i==0:
                caption=('Truth','Data-exposed TV=0 ablation','Frozen TV=.2 baseline')[j]
                if j:caption+=f"\nIoU {rows[j-1]['support_IoU']:.3f}; holes {rows[j-1]['true_holes']}→{rows[j-1]['predicted_holes']}"
                ax.set_title(caption,fontsize=10)
            if j==0:ax.set_ylabel(label)
        fig.colorbar(im,ax=axes[i,:],shrink=.8,label=f'{label}, shared 0..{maximum:g}')
    fig.suptitle('Actual same-data digit8 TV ablation; all110 material coordinates free\nFixed 7×7 cells/scales; no interpolation or truth-based selection\nDifferent TV priors make raw total objectives incomparable as an accuracy ranking',fontsize=10)
    path=out/'DIGIT8_TV0_VS_TV02';fig.savefig(path.with_suffix('.png'),dpi=170);fig.savefig(path.with_suffix('.pdf'));plt.close(fig)
    result=dict(actual_endpoints=2,figure=str(path.with_suffix('.png').relative_to(root)),physics_calls=0,
        material_selection_using_truth=False,raw_objectives_have_different_regularization=True,scope='DATA_EXPOSED_SAME_DATA_TV_ABLATION')
    (out/'TV_ABLATION_FIGURE_RECEIPT.json').write_text(json.dumps(result,indent=2)+'\n');return result

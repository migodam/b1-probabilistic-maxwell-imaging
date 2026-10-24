"""Separate full66 anchor-conditioned fit; never changes raw-y trainer semantics.

NPZ: x[N,66] labels, y[N,512] raw data, context_x[N,66] y-only anchor,
optional base_mean/base_chol, context_J[N,512,66] native whitened Jq.
--metric requires precomputed velocity_metric[N,66,66] and an explicit
truth-free TRAIN-frozen metric contract. No physics/metric construction here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import resource
import signal
import sys
import time

import numpy as np
import torch

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from src.p3r4 import anchor_density as anchor
    from src.p3r4.density import KINDS, make_model
    from src.p3r4.train_density import load_data, fit_normalization, normalized, write_json
else:
    from . import anchor_density as anchor
    from .density import KINDS, make_model
    from .train_density import load_data, fit_normalization, normalized, write_json


def load_anchor_data(path, j_summary='none', metric='none'):
    data = load_data(path)
    if 'offline_only' in Path(path).resolve().parts:
        raise ValueError('offline_only truth forbidden in training')
    with np.load(path, allow_pickle=False) as z:
        data['context_x'] = z['context_x']
        data['context_J'] = z['context_J'] if j_summary != 'none' else None
        data['velocity_metric'] = z['velocity_metric'] if metric != 'none' else None
    if data['parent'] is None:
        raise ValueError('anchor fits require explicit independent parent IDs')
    anchor.context_features(data['y'], data['context_x'], data['context_J'], j_summary)
    return data


def validate_metric(matrix, count):
    matrix = np.asarray(matrix)
    if matrix.shape != (count, 66, 66) or np.iscomplexobj(matrix) or not np.isfinite(matrix).all():
        raise ValueError('velocity_metric[N,66,66] finite real matrix required')
    if not np.allclose(matrix, matrix.transpose(0, 2, 1), atol=2e-6, rtol=0):
        raise ValueError('symmetric velocity-error metric required')
    lo, hi = np.inf, -np.inf
    for begin in range(0, count, 64):
        eig = torch.linalg.eigvalsh(torch.from_numpy(np.asarray(matrix[begin:begin+64], np.float64))).numpy()
        lo, hi = min(lo, float(eig.min())), max(hi, float(eig.max()))
    if lo < .5-1e-5 or hi > 2+1e-5:
        raise ValueError('velocity-error metric eigenvalues must stay within[0.5,2]')
    return matrix.astype(np.float32), dict(minimum_eigenvalue=lo, maximum_eigenvalue=hi)


def prepare(data, stats, args):
    features = anchor.context_features(data['y'], data['context_x'], data['context_J'], args.j_summary)
    mean, chol = anchor.base_parameters(data['context_x'], data['base_mean'], data['base_chol'], args.anchor_scale)
    center = mean if args.residual_center else np.zeros_like(mean)
    raw = dict(x=data['x']-center, y=features, parent=data['parent'], base_mean=mean-center, base_chol=chol)
    result = normalized(raw, stats, 'anchor' if args.kind == 'cfm' else 'standard', args.anchor_scale)
    result['velocity_metric'] = data['velocity_metric']
    return result


def cfm_error_loss(model, x1, y, generator, base_mean, base_chol, metric=None):
    x0 = model.base_sample(y, generator, base_mean, base_chol)
    t = torch.rand((len(y), 1), generator=generator, device='cpu', dtype=y.dtype).to(y.device)
    error = model.velocity((1-t)*x0+t*x1, t, y)-(x1-x0)
    if metric is None:
        return error.square().mean()
    # Weight the velocity ERROR only; never project or penalize velocity itself.
    return torch.einsum('bi,bij,bj->b', error, metric, error).mean()/66


def run(args):
    wall, cpu = time.monotonic(), time.process_time()
    torch.set_num_threads(1)
    torch.manual_seed(args.seed)
    generator = torch.Generator().manual_seed(args.seed+1103)
    rng = np.random.default_rng(args.seed)
    out = args.out
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('fresh separate anchor-fit destination required')
    out.mkdir(parents=True, exist_ok=True)
    train = load_anchor_data(args.train, args.j_summary, args.metric)
    val = load_anchor_data(args.validation, args.j_summary, args.metric)
    if set(train['parent']) & set(val['parent']):
        raise ValueError('TRAIN/validation parent overlap')
    if args.train_size:
        if not 2 <= args.train_size <= len(train['x']):
            raise ValueError('2<=train-size<=available required')
        train = {k: (v[:args.train_size] if v is not None else None) for k, v in train.items()}
    frame = None
    if args.frame:
        if 'locked' in str(args.frame).lower():
            raise ValueError('locked frame forbidden')
        if str(args.frame).endswith('.npy'):
            frame = np.load(args.frame, allow_pickle=False)
        else:
            with np.load(args.frame, allow_pickle=False) as z:
                frame = z['frame']
    train_features = anchor.context_features(train['y'], train['context_x'], train['context_J'], args.j_summary)
    mean, _ = anchor.base_parameters(train['context_x'], train['base_mean'], train['base_chol'], args.anchor_scale)
    residual = train['x']-mean if args.residual_center else train['x']
    stats = fit_normalization(residual, train_features, frame)
    metric_contract, metric_checks = None, None
    if args.metric != 'none':
        metric_contract = json.loads(args.metric_contract.read_text())
        if metric_contract.get('kind') != args.metric or metric_contract.get('coordinate') != 'network_standardized_full66_velocity':
            raise ValueError('metric contract kind/coordinate mismatch')
        if metric_contract.get('endpoint_label_dependent') is not False or metric_contract.get('TRAIN_only_tuning') is not True:
            raise ValueError('explicit endpoint-independent TRAIN-only metric tuning contract required')
        if metric_contract.get('bounds') != [.5, 2.]:
            raise ValueError('frozen metric bounds must be[0.5,2]')
        train['velocity_metric'], tcheck = validate_metric(train['velocity_metric'], len(train['x']))
        val['velocity_metric'], vcheck = validate_metric(val['velocity_metric'], len(val['x']))
        metric_checks = dict(train=tcheck, validation=vcheck)
    tr, va = prepare(train, stats, args), prepare(val, stats, args)
    context_dim = train_features.shape[1]
    config = dict(dimension=66, context_dim=context_dim, width=args.width)
    if args.kind == 'cfm':
        config.update(blocks=4, time_features=32, velocity_hidden=args.velocity_hidden)
    else:
        config.update(layers=args.layers, bins=args.bins, output_rank=args.output_rank,
                      weak_dim=args.weak_dim, bound=5., diag_hidden=args.diag_hidden)
    device = torch.device(args.device)
    if args.device == 'mps' and not torch.backends.mps.is_available():
        raise RuntimeError('MPS requested but unavailable')
    model = make_model(args.kind, **config).to(device)
    count = sum(p.numel() for p in model.parameters())
    if count > args.capacity_limit:
        raise ValueError('declared parameter cap exceeded')
    ac = dict(j_summary=args.j_summary, residual_center=args.residual_center, anchor_scale=args.anchor_scale,
              context_dim=context_dim, context='all raw512 + all Gaussianized context_x66'+(' +66 log1p column squared norms of Jx' if args.j_summary != 'none' else ''),
              J_source_coordinate='cached whitened native Jq; property columns multiply standard-normal pdf(context_x_eta)',
              translation='base_mean(y) if supplied, otherwise context_x(y); determinant1; no endpoint dependence',
              base_covariance='supplied fullrank raw Cholesky or anchor_scale*I; computational base, not qualified posterior covariance',
              physical_cost='all context states/Jacobians belong to pre-existing parent receipts; no context physics performed by this trainer',
              metric=args.metric, metric_contract=metric_contract, metric_eigenvalue_checks=metric_checks,
              metric_scope='network-standardized velocity ERROR; never velocity projection; supplied y-only matrices; provenance parent-owned')
    criterion = 'validation_latent_energy_score_raw_x66' if args.kind == 'cfm' else 'validation_nll_raw_x66'
    contract = dict(format=anchor.FORMAT, kind=args.kind, seed=args.seed, model_config=config,
                    train=str(args.train.resolve()), validation=str(args.validation.resolve()),
                    train_size=len(tr['x']), validation_size=len(va['x']), parameter_count=count,
                    parameters_per_training_parent=count/len(tr['x']), capacity_limit=args.capacity_limit,
                    full66_model=True, cfm_base='anchor' if args.kind == 'cfm' else 'standard',
                    anchor_contract=ac, learning_rate=args.lr, updates=args.updates, batch_size=args.batch_size,
                    eval_every=args.eval_every, validation_samples=args.validation_samples, validation_nfe=args.validation_nfe,
                    normalizer='TRAIN-only residual/frame and all condition columns; y512 never dropped',
                    frame=str(args.frame.resolve()) if args.frame else None,
                    split='unique independent parents, TRAIN/validation disjoint', criterion=criterion, update_zero_eligible=True,
                    reference_measure='raw Gaussianized dx66; translation determinant1 + normalizer/frame Jacobian',
                    input_provenance='context_x/base_mean/J/metric must be observation-derived and frozen; external cached evidence is parent-owned',
                    comparison_scope='new explicit anchor family; raw-y first-family checkpoints and semantics unchanged',
                    physics_calls=0, posterior_fidelity='NOT_ESTABLISHED', device=args.device)
    write_json(out/'RUN_CONTRACT.json', contract)
    np.savez_compressed(out/'NORMALIZATION.npz', **{k: v for k, v in stats.items() if isinstance(v, np.ndarray)})
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history, selected, selected_value, selected_state = [], 0, float('inf'), None
    update, training_loss, stop = 0, None, 'UPDATE_LIMIT'
    requested = [False]
    old = {}
    for sig in (signal.SIGINT, signal.SIGTERM):
        old[sig] = signal.signal(sig, lambda *_: requested.__setitem__(0, True))
    def batch(data, ids):
        return [torch.from_numpy(data[k][ids]).to(device) if data[k] is not None else None
                for k in ('x', 'y', 'base_mean', 'base_chol', 'velocity_metric')]
    def payload(state, chosen, value):
        return dict(format=anchor.FORMAT, kind=args.kind, seed=args.seed, state_dict=state,
                    model_config=config, normalization=stats, anchor_contract=ac, contract=contract,
                    selected_update=chosen, selection_metric=value, criterion=criterion, posterior_fidelity='NOT_ESTABLISHED')
    def evaluate():
        nonlocal selected, selected_value, selected_state
        model.eval()
        generator_val = torch.Generator().manual_seed(args.seed+7019)
        values = []
        fm = []
        for begin in range(0, len(va['x']), args.validation_batch_size):
            x, y, bmean, bchol, metric = batch(va, slice(begin, begin+args.validation_batch_size))
            with torch.no_grad():
                if args.kind != 'cfm':
                    value = -model.log_prob(x, y)+float(np.log(stats['x_scale']).sum())-stats['frame_logabsdet']
                else:
                    fm.append((float(cfm_error_loss(model, x, y, generator_val, bmean, bchol, metric).cpu()), len(x)))
                    sample = model.sample(y, args.validation_samples, generator_val, nfe=args.validation_nfe, base_mean=bmean, base_chol=bchol)
                    # Fixed orthogonal frame: raw Euclidean energy norm equals
                    # scaled frame-coordinate norm; translation cancels.
                    scale = torch.as_tensor(stats['x_scale'], dtype=x.dtype, device=device)
                    value = torch.linalg.vector_norm((sample-x[:, None])*scale, dim=-1).mean(1)
                    pair = x.new_zeros(len(x))
                    k = args.validation_samples
                    for i in range(k):
                        for j in range(i):
                            pair += torch.linalg.vector_norm((sample[:, i]-sample[:, j])*scale, dim=-1)
                    value = value-pair/(k*(k-1))
            values.append(value.cpu().numpy())
        values = np.concatenate(values).astype(np.float64)
        if not np.isfinite(values).all():
            raise FloatingPointError('nonfinite anchor validation score')
        value = float(values.mean())
        row = dict(update=update, train_minibatch_loss=training_loss, validation_metric=value,
                   validation_weighted_fm_mse=sum(v*n for v, n in fm)/len(va['x']) if fm else None,
                   learning_rate=args.lr, wall_seconds=time.monotonic()-wall)
        history.append(row)
        with (out/'loss_curve.jsonl').open('a') as file:
            file.write(json.dumps(row, allow_nan=False)+'\n')
        if value < selected_value:
            selected, selected_value = update, value
            selected_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            torch.save(payload(selected_state, selected, selected_value), out/'BEST.pt')
            np.savez_compressed(out/'validation_metric.npz', parent=va['parent'], metric=values)
        torch.save(dict(**payload({k:v.detach().cpu() for k,v in model.state_dict().items()}, update, value),
                        optimizer=optimizer.state_dict(), generator_state=generator.get_state(), numpy_rng_state=rng.bit_generator.state,
                        torch_rng_state=torch.get_rng_state()), out/'LAST.pt')
        print(json.dumps(dict(event='validation', **row, selected_update=selected)), flush=True)
    try:
        evaluate()
        for update in range(1, args.updates+1):
            if requested[0] or (args.max_wall and time.monotonic()-wall >= args.max_wall):
                update -= 1
                stop = 'SIGNAL' if requested[0] else 'WALL_LIMIT'
                break
            model.train()
            ids = rng.integers(0, len(tr['x']), min(args.batch_size, len(tr['x'])))
            x, y, bmean, bchol, metric = batch(tr, ids)
            optimizer.zero_grad(set_to_none=True)
            loss = cfm_error_loss(model, x, y, generator, bmean, bchol, metric) if args.kind == 'cfm' else -model.log_prob(x, y).mean()
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError('nonfinite anchor loss')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5., error_if_nonfinite=True)
            optimizer.step()
            training_loss = float(loss.detach().cpu())
            if update % args.eval_every == 0 or update == args.updates:
                evaluate()
        if history[-1]['update'] != update:
            evaluate()
    finally:
        for sig, handler in old.items():
            signal.signal(sig, handler)
    np.savez_compressed(out/'loss_curve.npz', update=np.array([r['update'] for r in history]), validation=np.array([r['validation_metric'] for r in history]))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7,4))
    ax.plot([r['update'] for r in history], [r['validation_metric'] for r in history], marker='o')
    ax.axvline(selected, ls='--', color='gray')
    ax.set(xlabel='Optimizer update (0 eligible)', ylabel=criterion, title=f'Anchor {args.kind}: train{len(tr["x"])}, parameters{count:,}')
    fig.tight_layout()
    fig.savefig(out/'LOSS_CURVE.png', dpi=150)
    plt.close(fig)
    result = dict(**contract, selected_update=selected, selected_initialization=selected==0,
                  selected_validation_metric=selected_value, completed_updates=update, stop_reason=stop,
                  wall_seconds=time.monotonic()-wall, cpu_seconds=time.process_time()-cpu,
                  peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024), history=history)
    write_json(out/'RESULT.json', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('train', 'validation', 'out'):
        parser.add_argument('--'+key, required=True, type=Path)
    parser.add_argument('--kind', choices=KINDS, required=True)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--updates', type=int, required=True)
    parser.add_argument('--device', choices=('cpu','mps'), default='cpu')
    parser.add_argument('--residual-center', action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument('--j-summary', choices=('none','log_column_norm'), default='none')
    parser.add_argument('--anchor-scale', type=float, default=1.)
    parser.add_argument('--metric', choices=('none','maxwell','decoder','random'), default='none')
    parser.add_argument('--metric-contract', type=Path)
    parser.add_argument('--frame', type=Path)
    parser.add_argument('--train-size', type=int)
    parser.add_argument('--width', type=int, default=128)
    parser.add_argument('--diag-hidden', type=int, default=655)
    parser.add_argument('--velocity-hidden', type=int, default=547, help='CFM feedforward hidden width; preserves width128 encoder/residual state and matches spline capacity')
    parser.add_argument('--layers', type=int, default=8)
    parser.add_argument('--bins', type=int, default=10)
    parser.add_argument('--output-rank', type=int, default=32)
    parser.add_argument('--weak-dim', type=int, default=16)
    parser.add_argument('--capacity-limit', type=int, default=2_000_000)
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--validation-batch-size', type=int, default=32)
    parser.add_argument('--validation-samples', type=int, default=4)
    parser.add_argument('--validation-nfe', type=int, default=16)
    parser.add_argument('--eval-every', type=int, default=100)
    parser.add_argument('--lr', type=float, default=2e-4)
    parser.add_argument('--weight-decay', type=float, default=1e-4)
    parser.add_argument('--max-wall','--max-seconds', type=float, default=0.)
    args = parser.parse_args(argv)
    if args.updates<0 or min(args.width,args.diag_hidden,args.velocity_hidden,args.layers,args.bins,args.output_rank,args.batch_size,args.validation_batch_size,args.eval_every,args.capacity_limit)<1:
        parser.error('nonnegative updates and positive dimensions/batches/cap required')
    if not np.isfinite([args.lr,args.weight_decay,args.max_wall,args.anchor_scale]).all() or args.lr<=0 or args.weight_decay<0 or args.max_wall<0 or args.anchor_scale<=0:
        parser.error('finite positive LR/anchor-scale; nonnegative decay/time required')
    if args.validation_samples<2 or args.validation_nfe<4 or args.validation_nfe%4:
        parser.error('validation-samples>=2 and RK4 NFE positive multiple4 required')
    if args.metric!='none' and (args.kind!='cfm' or args.metric_contract is None):
        parser.error('velocity-error metric requires CFM and --metric-contract')
    run(args)


if __name__ == '__main__':
    main()

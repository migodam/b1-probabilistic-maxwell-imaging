"""Bounded full66 density/CFM training on existing independent NPZ joint pairs.

Example: python -m src.p3r4.train_density --train TRAIN.npz --validation VAL.npz
    --out RUN --kind spline --seed 11 --updates 1000 --device cpu

NPZ contract: real raw Gaussianized x[N,66], real packed y[N,512]; optional
parent[N] IDs. CFM --base anchor additionally requires y-only base_mean[N,66],
and optionally base_chol[N,66,66]; otherwise anchor_scale*I is used. Normalizers
and optional fixed orthogonal --frame are fitted/selected using TRAIN only.
This entry point schedules no runs beyond its one explicitly requested fit.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time

import numpy as np
import torch

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from src.p3r4.density import KINDS, make_model
else:
    from .density import KINDS, make_model


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def load_data(path):
    if 'locked' in str(path).lower():
        raise ValueError('locked truth is forbidden in training/validation')
    with np.load(path, allow_pickle=False) as data:
        x, y = data['x'], data['y']
        parent_key = next((k for k in ('parent', 'parents', 'parent_id') if k in data), None)
        parent = data[parent_key].astype(str) if parent_key else None
        base_mean = data['base_mean'].copy() if 'base_mean' in data else None
        base_chol = data['base_chol'].copy() if 'base_chol' in data else None
    if x.ndim != 2 or x.shape[1] != 66 or y.shape != (len(x), 512) or len(x) < 2:
        raise ValueError('x[N,66], y[N,512] with N>=2 required')
    if np.iscomplexobj(x) or np.iscomplexobj(y) or not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError('finite real raw Gaussianized x and source-major y required')
    if parent is not None and (parent.shape != (len(x),) or len(set(parent)) != len(parent)):
        raise ValueError('unique independent parent IDs required')
    return dict(x=x, y=y, parent=parent, base_mean=base_mean, base_chol=base_chol)


def fit_normalization(x, y, frame=None, identity=False):
    if frame is None:
        frame = np.eye(66)
    frame = np.asarray(frame, dtype=np.float64)
    if frame.shape != (66, 66) or not np.isfinite(frame).all() or not np.allclose(frame.T@frame, np.eye(66), atol=1e-6, rtol=0):
        raise ValueError('frozen real orthogonal frame[66,66] required')
    coords = np.asarray(x, np.float64)@frame
    if identity:
        mean_x, scale_x = np.zeros(66), np.ones(66)
        mean_y, scale_y = np.zeros(512), np.ones(512)
    else:
        mean_x, scale_x = coords.mean(0), coords.std(0)
        mean_y, scale_y = np.asarray(y, np.float64).mean(0), np.asarray(y, np.float64).std(0)
        scale_x = np.where(scale_x < 1e-8, 1., scale_x)
        scale_y = np.where(scale_y < 1e-12, 1., scale_y)
    if not all(np.isfinite(v).all() for v in (mean_x, scale_x, mean_y, scale_y)):
        raise ValueError('nonfinite TRAIN normalization statistics')
    return dict(frame=frame, frame_inverse=np.linalg.inv(frame),
                frame_logabsdet=float(torch.linalg.slogdet(torch.from_numpy(frame)).logabsdet),
                x_mean=mean_x, x_scale=scale_x, y_mean=mean_y, y_scale=scale_y,
                fit='TRAIN only; orthogonal frame then per-coordinate centering/scaling',
                already_standardized=identity)


def normalized(data, stats, base='standard', anchor_scale=1.):
    x = ((data['x']@stats['frame']-stats['x_mean'])/stats['x_scale']).astype(np.float32)
    y = ((data['y']-stats['y_mean'])/stats['y_scale']).astype(np.float32)
    result = dict(x=x, y=y, parent=data['parent'], base_mean=None, base_chol=None)
    if base == 'anchor':
        mean = data['base_mean']
        if mean is None or mean.shape != data['x'].shape or not np.isfinite(mean).all():
            raise ValueError('--base anchor requires finite y-only base_mean[N,66] in both NPZ files')
        result['base_mean'] = ((mean@stats['frame']-stats['x_mean'])/stats['x_scale']).astype(np.float32)
        raw_chol = data['base_chol']
        if raw_chol is None:
            # An isotropic raw-coordinate base stays diagonal under any
            # orthogonal frame, followed by coordinate standardization.
            result['base_chol'] = np.broadcast_to(np.diag(anchor_scale/stats['x_scale']),
                                                   (len(x), 66, 66)).astype(np.float32).copy()
        else:
            if raw_chol.shape != (len(x), 66, 66) or not np.isfinite(raw_chol).all():
                raise ValueError('finite base_chol[N,66,66] required')
            if np.any(np.diagonal(raw_chol, axis1=1, axis2=2)<=0) or not np.allclose(raw_chol, np.tril(raw_chol)):
                raise ValueError('base_chol must be lower triangular with strictly positive diagonal')
            chol = np.empty_like(raw_chol, dtype=np.float32)
            transform = stats['frame'].T/stats['x_scale'][:, None]
            for begin in range(0, len(x), 128):
                a = transform@raw_chol[begin:begin+128]
                chol[begin:begin+128] = np.linalg.cholesky(a@a.transpose(0, 2, 1))
            result['base_chol'] = chol
    return result


def run(args):
    wall, cpu = time.monotonic(), time.process_time()
    torch.set_num_threads(args.threads)
    try:
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass
    torch.manual_seed(args.seed)
    rng = np.random.default_rng(args.seed)
    generator = torch.Generator(device='cpu').manual_seed(args.seed+1103)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'RESULT.json').exists() or (out/'BEST.pt').exists():
        raise FileExistsError('output contains prior results; choose a fresh --out')
    if Path(args.train).resolve() == Path(args.validation).resolve():
        raise ValueError('train and validation files must differ')
    train, val = load_data(args.train), load_data(args.validation)
    ids_checked = train['parent'] is not None and val['parent'] is not None
    if ids_checked and set(train['parent']) & set(val['parent']):
        raise ValueError('TRAIN/validation independent parents overlap')
    available = len(train['x'])
    if args.train_size:
        if args.train_size < 2 or args.train_size > available:
            raise ValueError('2<=train-size<=available required')
        # Parent ordering/subset is fixed in advance; no label-based selection.
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
    stats = fit_normalization(train['x'], train['y'], frame, args.already_standardized)
    tr, va = normalized(train, stats, args.base, args.anchor_scale), normalized(val, stats, args.base, args.anchor_scale)
    config = dict(dimension=66, context_dim=512, width=args.width, layers=args.layers,
                  bins=args.bins, output_rank=args.output_rank, weak_dim=args.weak_dim,
                  bound=args.bound, diag_hidden=args.diag_hidden)
    if args.kind == 'cfm':
        config = dict(dimension=66, context_dim=512, width=args.width, blocks=4, time_features=32,
                      velocity_hidden=args.velocity_hidden)
    dev = torch.device(args.device)
    if args.device == 'mps' and not torch.backends.mps.is_available():
        raise RuntimeError('MPS requested but unavailable; no silent device switch')
    model = make_model(args.kind, **config).to(dev)
    if args.kind == 'cfm':
        # The declared standard base is N(0,I) in RAW Gaussianized coordinates,
        # even though the velocity network works in training-standardized coords.
        with torch.no_grad():
            model.base_location.copy_(torch.as_tensor(-stats['x_mean']/stats['x_scale'], dtype=torch.float32, device=dev))
            model.base_scale.copy_(torch.as_tensor(1/stats['x_scale'], dtype=torch.float32, device=dev))
    parameter_count = sum(p.numel() for p in model.parameters())
    if parameter_count > args.capacity_limit:
        raise ValueError(f'{parameter_count} parameters exceed declared cap {args.capacity_limit}')
    criterion = 'validation_latent_energy_score_raw_x66' if args.kind == 'cfm' else 'validation_nll_raw_x66'
    contract = dict(format='p3r4_full66_density_v1', kind=args.kind, seed=args.seed,
                    train=str(Path(args.train).resolve()), validation=str(Path(args.validation).resolve()),
                    train_size=len(tr['x']), available_train_size=available, validation_size=len(va['x']),
                    parent_disjointness='CHECKED' if ids_checked else 'NOT_CHECKABLE_IDS_ABSENT',
                    model_config=config, parameter_count=parameter_count, capacity_limit=args.capacity_limit,
                    parameters_per_training_parent=parameter_count/len(tr['x']),
                    updates=args.updates, learning_rate=args.lr, learning_rate_check='positive finite fixed LR; no automatic LR search',
                    data_scale_check='TRAIN-only finite column scales; constant columns scale=1; raw y512 preserved',
                    capacity_check='actual parameter count versus cap; train-size and parameters/parent reported; no automatic enlargement',
                    batch_size=args.batch_size, eval_every=args.eval_every, validation_batch_size=args.validation_batch_size,
                    criterion=criterion, update_zero_eligible=True, device=args.device, threads=args.threads,
                    cfm_base=args.base, anchor_scale=args.anchor_scale,
                    base_provenance='external base_mean/base_chol must be y-only, frozen; content provenance is parent responsibility',
                    validation_samples=args.validation_samples, validation_nfe=args.validation_nfe,
                    reference_measure='raw Gaussianized dx66; full Jacobian of TRAIN normalizer included',
                    frame_path=str(Path(args.frame).resolve()) if args.frame else None,
                    frame_rule='fixed supplied TRAIN-only orthogonal frame; never depends on candidate x',
                    physics_calls=0, posterior_fidelity='NOT_ESTABLISHED',
                    cfm_density='full-trace numerical continuous-flow density only; not exact MH',
                    validation_fm_loss_is_not_density_quality=True)
    write_json(out/'RUN_CONTRACT.json', contract)
    np.savez_compressed(out/'NORMALIZATION.npz', **{k: v for k, v in stats.items() if isinstance(v, np.ndarray)})
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history, best_value, best_update, best_state = [], float('inf'), 0, None
    stop_reason, update, train_loss = 'UPDATE_LIMIT', 0, None
    stop_requested = [False]
    old_handlers = {}
    def on_signal(signum, frame):
        stop_requested[0] = True
    for sig in (signal.SIGTERM, signal.SIGINT):
        old_handlers[sig] = signal.signal(sig, on_signal)
    def tensors(data, ii):
        return tuple(torch.from_numpy(data[k][ii]).to(dev) if data[k] is not None else None
                     for k in ('x', 'y', 'base_mean', 'base_chol'))
    def evaluate():
        model.eval()
        losses, fmlosses = [], []
        vg = torch.Generator(device='cpu').manual_seed(args.seed+7019)
        # Fixed validation noise/times and sample indices at every checkpoint.
        for begin in range(0, len(va['x']), args.validation_batch_size):
            ii = slice(begin, begin+args.validation_batch_size)
            x, y, mean, chol = tensors(va, ii)
            with torch.no_grad():
                if args.kind != 'cfm':
                    loss = -model.log_prob(x, y)+float(np.log(stats['x_scale']).sum())-stats['frame_logabsdet']
                else:
                    fmlosses.append((float(model.loss(x, y, vg, mean, chol).cpu()), len(x)))
                    samples = model.sample(y, args.validation_samples, vg, nfe=args.validation_nfe,
                                           base_mean=mean, base_chol=chol)
                    k = args.validation_samples
                    raw_scale = torch.as_tensor(stats['x_scale'], dtype=x.dtype, device=dev)
                    first = torch.linalg.vector_norm((samples-x[:, None])*raw_scale, dim=-1).mean(1)
                    # U-statistic pair term: E||X-X'|| / 2, distinct draws only.
                    pair = x.new_zeros(len(x))
                    for i in range(k):
                        for j in range(i):
                            pair += torch.linalg.vector_norm((samples[:, i]-samples[:, j])*raw_scale, dim=-1)
                    loss = first-pair/(k*(k-1))
            losses.append(loss.cpu().numpy().astype(np.float64))
        rows = np.concatenate(losses)
        value = float(rows.mean())
        if not np.isfinite(rows).all():
            raise FloatingPointError('nonfinite validation metric')
        fm = sum(x*n for x, n in fmlosses)/len(va['x']) if fmlosses else None
        return value, rows, fm
    def checkpoint(state, selected_update, selected_metric):
        return dict(format='p3r4_full66_density_v1', kind=args.kind, model_config=config,
                    normalization=stats, state_dict=state, seed=args.seed,
                    selected_update=selected_update, selection_metric=selected_metric, criterion=criterion,
                    contract=contract, posterior_fidelity='NOT_ESTABLISHED')
    def validation_event():
        nonlocal best_value, best_update, best_state
        value, rows, fm = evaluate()
        row = dict(update=update, train_minibatch_loss=train_loss, validation_metric=value,
                   validation_fm_mse=fm, learning_rate=optimizer.param_groups[0]['lr'],
                   wall_seconds=time.monotonic()-wall, cpu_seconds=time.process_time()-cpu)
        row[criterion] = value
        history.append(row)
        with (out/'loss_curve.jsonl').open('a') as f:
            f.write(json.dumps(row, allow_nan=False)+'\n')
        if value < best_value:
            best_value, best_update = value, update
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            torch.save(checkpoint(best_state, best_update, best_value), out/'BEST.pt')
            np.savez_compressed(out/'validation_metric.npz', metric=rows,
                                parent=va['parent'] if va['parent'] is not None else np.arange(len(rows)))
        current = {k: v.detach().cpu() for k, v in model.state_dict().items()}
        torch.save(dict(**checkpoint(current, update, value), optimizer=optimizer.state_dict(),
                        torch_rng_state=torch.get_rng_state(), generator_state=generator.get_state(),
                        numpy_rng_state=rng.bit_generator.state), out/'LAST.pt')
        print(json.dumps(dict(event='validation', **row, selected_update=best_update)), flush=True)
    try:
        validation_event()  # update0 can win; never discard an initialization win.
        for update in range(1, args.updates+1):
            if stop_requested[0] or (args.max_seconds and time.monotonic()-wall >= args.max_seconds):
                update -= 1
                stop_reason = 'SIGNAL' if stop_requested[0] else 'WALL_LIMIT'
                break
            model.train()
            ii = rng.integers(0, len(tr['x']), size=min(args.batch_size, len(tr['x'])))
            x, y, mean, chol = tensors(tr, ii)
            optimizer.zero_grad(set_to_none=True)
            if args.kind == 'cfm':
                loss = model.loss(x, y, generator, mean, chol)
            else:
                loss = -model.log_prob(x, y).mean()
            if not bool(torch.isfinite(loss)):
                raise FloatingPointError('nonfinite training loss')
            loss.backward()
            grad_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), args.grad_clip, error_if_nonfinite=True)
            optimizer.step()
            train_loss = float(loss.detach().cpu())
            if update % args.eval_every == 0 or update == args.updates:
                validation_event()
        if history[-1]['update'] != update:
            validation_event()
    finally:
        for sig, handler in old_handlers.items():
            signal.signal(sig, handler)
    model.load_state_dict(best_state)
    result = dict(**contract, selected_update=best_update, selected_initialization=best_update==0,
                  selected_validation_metric=best_value, completed_updates=update, stop_reason=stop_reason,
                  wall_seconds=time.monotonic()-wall, cpu_seconds=time.process_time()-cpu,
                  peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
                  max_seconds=args.max_seconds, weight_decay=args.weight_decay,
                  history=history, checkpoint=str((out/'BEST.pt').resolve()))
    if args.kind != 'cfm':
        result['validation_nll'] = best_value
        result['validation_unit_prior_nll_raw_x66'] = float(.5*(val['x'].astype(np.float64)**2+np.log(2*np.pi)).sum(1).mean())
    np.savez_compressed(out/'loss_curve.npz', update=np.array([h['update'] for h in history]),
                        validation=np.array([h['validation_metric'] for h in history]),
                        training=np.array([np.nan if h['train_minibatch_loss'] is None else h['train_minibatch_loss'] for h in history]))
    # Diagnostic figure; numerical JSONL/NPZ is authoritative.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot([h['update'] for h in history], [h['validation_metric'] for h in history], marker='o', label=criterion)
    ax.axvline(best_update, color='gray', ls='--', label=f'selected update {best_update}')
    ax.set(xlabel='Optimizer updates (update0 eligible)', ylabel=criterion,
           title=f'{args.kind} | train {len(tr["x"])} | parameters {parameter_count:,}')
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out/'LOSS_CURVE.png', dpi=150)
    plt.close(fig)
    result.update(wall_seconds=time.monotonic()-wall, cpu_seconds=time.process_time()-cpu,
                  peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024))
    write_json(out/'RESULT.json', result)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--train', required=True, type=Path)
    parser.add_argument('--validation', required=True, type=Path)
    parser.add_argument('--out', '--output', required=True, type=Path)
    parser.add_argument('--kind', required=True, choices=KINDS)
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--updates', type=int, required=True)
    parser.add_argument('--device', choices=('cpu', 'mps'), default='cpu')
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--validation-batch-size', type=int, default=64)
    parser.add_argument('--eval-every', type=int, default=100)
    parser.add_argument('--validation-samples', type=int, default=4)
    parser.add_argument('--validation-nfe', type=int, default=16)
    parser.add_argument('--width', type=int, default=128)
    parser.add_argument('--layers', type=int, default=8)
    parser.add_argument('--bins', type=int, default=10)
    parser.add_argument('--output-rank', type=int, default=32)
    parser.add_argument('--weak-dim', type=int, default=16)
    parser.add_argument('--diag-hidden', type=int, default=655, help='two-layer diagonal head width; default matches spline parameter count')
    parser.add_argument('--velocity-hidden',type=int,default=None,help='optional CFM feedforward expansion; None preserves original model')
    parser.add_argument('--bound', type=float, default=5.)
    parser.add_argument('--capacity-limit', type=int, default=2_000_000)
    parser.add_argument('--train-size', type=int)
    parser.add_argument('--lr', '--learning-rate', type=float, default=2e-4)
    parser.add_argument('--weight-decay', type=float, default=1e-4)
    parser.add_argument('--grad-clip', type=float, default=5.)
    parser.add_argument('--max-seconds', '--max-wall', '--max-wall-seconds', type=float, default=0., help='0 delegates wall limit to parent scheduler')
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--already-standardized', action='store_true')
    parser.add_argument('--frame', type=Path, help='fixed TRAIN-only orthogonal frame66 NPZ frame or NPY')
    parser.add_argument('--base', choices=('standard', 'anchor'), default='standard')
    parser.add_argument('--anchor-scale', type=float, default=1.)
    args = parser.parse_args(argv)
    if args.updates < 0 or min(args.batch_size, args.validation_batch_size, args.eval_every, args.threads,
                              args.width, args.output_rank, args.diag_hidden, args.capacity_limit) < 1:
        parser.error('updates>=0 and positive batches/evaluation interval/threads/width/cap required')
    if not np.isfinite([args.lr, args.weight_decay, args.max_seconds, args.bound, args.grad_clip, args.anchor_scale]).all() or args.lr<=0 or args.weight_decay<0 or args.max_seconds<0 or args.bound<=0 or args.grad_clip<=0 or args.anchor_scale<=0:
        parser.error('finite positive lr/bound/clip/anchor-scale and nonnegative decay/time required')
    if args.layers<2 or args.bins<2 or args.validation_samples<2 or args.validation_nfe<4 or args.validation_nfe%4:
        parser.error('layers>=2,bins>=2,validation-samples>=2,RK4 nfe positive multiple of4 required')
    if args.kind != 'cfm' and args.base != 'standard':
        parser.error('anchor base applies only to CFM')
    run(args)


if __name__ == '__main__':
    main()

"""Supervised material-risk training, with train-only normalization and no physics."""
import copy
import json
import time
from pathlib import Path
import numpy as np
import torch
from torch import nn
from .risk_models import MaterialRegressor, RiskPredictor, observation_features


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load_data(path):
    if 'locked' in str(path).lower():
        raise ValueError('locked data are not accepted by the training interface')
    with np.load(path, allow_pickle=False) as data:
        parent_key = 'parent' if 'parent' in data else 'parents'
        obs, chi, parents = data['observation'], data['chi'], data[parent_key].astype(str)
    if chi.shape != (len(obs), 576) or len(parents) != len(obs):
        raise ValueError('material or parent dimensions do not match observations')
    if len(np.unique(parents)) != len(parents):
        raise ValueError('each split must contain independent unique parents')
    if not np.iscomplexobj(chi) or not np.isfinite(obs).all() or not np.isfinite(chi).all():
        raise ValueError('finite complex observations and materials required')
    return obs, chi, parents


def material_targets(chi, roi):
    c = chi[:, roi]
    return np.concatenate((c.real, c.imag), axis=1).astype('float32')


def risk_rows(prediction, truth):
    return np.mean(np.abs(prediction - truth) ** 2, axis=1).astype('float64')


def risk_summary(rows):
    return dict(complex_roi_mse=float(np.mean(rows)), mean_real_component_mse=float(np.mean(rows)/2),
                parent_sem=float(np.std(rows, ddof=1)/np.sqrt(len(rows))) if len(rows)>1 else 0., n_parent=len(rows))


def evaluate(model, x, y, batch_size, device):
    model.eval()
    total = 0.
    with torch.inference_mode():
        for start in range(0, len(x), batch_size):
            xb = torch.from_numpy(x[start:start+batch_size]).to(device)
            yb = torch.from_numpy(y[start:start+batch_size]).to(device)
            total += float(torch.sum((model(xb)-yb)**2).cpu())
    return total / y.size


def train_variant(kind, seed, train_obs, train_chi, train_parents, val_obs, val_chi, val_parents,
                  config, roi, output, epochs=100, batch_size=256, patience=12,
                  max_seconds=600, learning_rate=1e-3, weight_decay=1e-4, device='cpu'):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    torch.set_num_threads(1)
    if torch.get_num_interop_threads() != 1:
        torch.set_num_interop_threads(1)
    torch.manual_seed(seed)
    np.random.seed(seed)
    rng = np.random.default_rng(seed)
    start_wall, start_cpu = time.monotonic(), time.process_time()
    raw_x, raw_v = observation_features(train_obs), observation_features(val_obs)
    x_mean = raw_x.mean(0, dtype=np.float64).astype('float32')
    x_scale = raw_x.std(0, dtype=np.float64).astype('float32')
    x_scale[x_scale < 1e-12] = 1
    x, vx = (raw_x-x_mean)/x_scale, (raw_v-x_mean)/x_scale
    y, vy = material_targets(train_chi, roi), material_targets(val_chi, roi)
    y_mean = y.mean(0, dtype=np.float64).astype('float32')
    # A single scale preserves actual unweighted material MSE as the objective.
    y_scale = float(max(np.sqrt(np.mean((y-y_mean)**2)), 1e-8))
    payload = dict(format='p3r2_material_risk_v1', kind=kind, seed=int(seed),
                   x_mean=x_mean, x_scale=x_scale, y_mean=y_mean, y_scale=y_scale,
                   roi_indices=roi, host_chi=config['host_chi'], target_id=config.get('target_id'),
                   packing='source-major Re32 then Im32', training_parent=train_parents,
                   estimator='material conditional mean point estimator; posterior calibration NOT_ESTABLISHED')
    history, best_update, reason, model = [], 0, 'BASELINE', None
    # A rerun replaces the corresponding log instead of silently appending runs.
    (output/'epochs.jsonl').write_text('')
    if kind == 'ridge':
        # Ridge is regularized average MSE: (X'X/n + alpha I)^-1 X'Y/n.
        # Fit all alphas from the same train Gram matrix and select on validation.
        gram = x.T.astype('float64') @ x.astype('float64') / len(x)
        rhs = x.T.astype('float64') @ (y-y_mean).astype('float64') / len(x)
        values, vectors = np.linalg.eigh(gram)
        rotated = vectors.T @ rhs
        best_risk = float('inf')
        for alpha in (1e-4, 1e-3, 1e-2, 1e-1, 1., 10., 100.):
            coef = (vectors @ (rotated/(values[:,None]+alpha))).astype('float32')
            mse = float(np.mean((vx @ coef + y_mean - vy)**2))
            history.append(dict(alpha=alpha, validation_complex_roi_mse=2*mse))
            if mse < best_risk:
                best_risk, payload['coef'], payload['alpha'] = mse, coef.copy(), alpha
        reason = 'VALIDATION_SELECTED_RIDGE_ALPHA'
    elif kind != 'mean':
        dev = torch.device(device)
        model = MaterialRegressor(kind, len(y_mean)).to(dev)
        yn, vyn = (y-y_mean)/y_scale, (vy-y_mean)/y_scale
        optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
        scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, factor=0.5, patience=4, min_lr=1e-6)
        best_risk = evaluate(model, vx, vyn, batch_size, dev)
        best_state = copy.deepcopy(model.state_dict())
        history.append(dict(epoch=0, update=0, validation_complex_roi_mse=2*best_risk*y_scale**2,
                            train_complex_roi_mse=2*evaluate(model,x,yn,batch_size,dev)*y_scale**2,
                            learning_rate=learning_rate, wall_seconds=time.monotonic()-start_wall))
        with (output/'epochs.jsonl').open('a') as f:
            f.write(json.dumps(history[-1], allow_nan=False)+'\n')
        update, stale = 0, 0
        reason = 'EPOCH_LIMIT'
        for epoch in range(1, epochs+1):
            if time.monotonic()-start_wall >= max_seconds:
                reason = 'TIME_LIMIT'; break
            model.train()
            order = rng.permutation(len(x))
            batch_loss, seen = 0., 0
            interrupted = False
            for pos in range(0, len(order), batch_size):
                if time.monotonic()-start_wall >= max_seconds:
                    interrupted = True; break
                ii = order[pos:pos+batch_size]
                xb, yb = torch.from_numpy(x[ii]).to(dev), torch.from_numpy(yn[ii]).to(dev)
                optimizer.zero_grad(set_to_none=True)
                loss = nn.functional.mse_loss(model(xb), yb)
                if not torch.isfinite(loss):
                    raise RuntimeError('nonfinite training loss')
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), 5.)
                optimizer.step()
                update += 1; seen += len(ii); batch_loss += float(loss.detach().cpu())*len(ii)
            current = evaluate(model, vx, vyn, batch_size, dev)
            scheduler.step(current)
            row = dict(epoch=epoch, update=update, completed_epoch=not interrupted,
                       validation_complex_roi_mse=2*current*y_scale**2,
                       train_minibatch_complex_roi_mse=2*batch_loss/max(seen,1)*y_scale**2,
                       learning_rate=float(optimizer.param_groups[0]['lr']), wall_seconds=time.monotonic()-start_wall)
            history.append(row)
            with (output/'epochs.jsonl').open('a') as f:
                f.write(json.dumps(row, allow_nan=False)+'\n')
            if current < best_risk:
                best_risk, best_update, stale = current, update, 0
                best_state = copy.deepcopy(model.state_dict())
            else:
                stale += 1
            if interrupted:
                reason = 'TIME_LIMIT'; break
            if stale >= patience:
                reason = 'VALIDATION_PATIENCE'; break
        payload['state_dict'] = {k:v.detach().cpu() for k,v in best_state.items()}
    payload['best_update'] = best_update
    torch.save(payload, output/'BEST.pt')
    predictor = RiskPredictor(payload, device)
    val_prediction = predictor.predict(val_obs, batch_size)
    rows = risk_rows(val_prediction, val_chi[:,roi])
    train_rows = risk_rows(predictor.predict(train_obs, batch_size), train_chi[:,roi])
    np.savez_compressed(output/'validation_risks.npz', parent=val_parents, risk=rows,
                        prediction_roi=val_prediction, roi_indices=roi)
    summary = dict(kind=kind, seed=int(seed), train_size=len(train_obs), validation_size=len(val_obs),
                   selected_update=int(best_update), selected_initialization=bool(model is not None and best_update==0),
                   stop_reason=reason, validation=risk_summary(rows), train=risk_summary(train_rows),
                   parameter_count=sum(p.numel() for p in model.parameters()) if model else
                       int(y_mean.size + (payload['coef'].size if kind=='ridge' else 0)),
                   ridge_alpha=payload.get('alpha'), max_seconds=max_seconds, epochs_limit=epochs,
                   batch_size=batch_size, learning_rate=learning_rate, weight_decay=weight_decay,
                   training_device=str(device), torch_num_threads=1,
                   wall_seconds=float(time.monotonic()-start_wall), cpu_seconds=float(time.process_time()-start_cpu),
                   objective='actual complex ROI material MSE; training mean component MSE differs by factor 2 only',
                   normalization='observation train-only per feature; material train-only mean and single global scale',
                   selection='full validation each epoch and update zero; locked ground truth never read',
                   calibrated_posterior='NOT_ESTABLISHED', history=history)
    write_json(output/'RESULT.json', summary)
    return summary

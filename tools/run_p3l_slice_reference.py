"""Paid tensor quadrature on a declared, restricted 3-D conditional slice.

This executable is deliberately separate from FULL66 posterior inference. It
never changes TARGET.json, performs inner Jacobian queries, reads GT, or trains.
Merely importing this module performs no physics or campaign bookkeeping.
"""
from pathlib import Path
import gc
import json
import sys

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import log_ndtr, logsumexp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from p3l.core import OUT, Book, Target, write, save, plain
from p3l.conditional import make_frame

PARENT = 'DEV_000'
JOB = 'slice3_reference_DEV000'
ORDERS = ((9, 16), (17, 24), (25, 32))
MAX_ARRAY_BYTES = 400 * 1024**2


def log_normal_interval(lo, hi, mean):
    """Stable log Gaussian interval probability, including positive-tail boxes."""
    lower, upper = float(lo - mean), float(hi - mean)
    if not lower < upper:
        raise ValueError('nonempty prior interval required')
    if lower > 0:
        lower, upper = -upper, -lower
    a, b = float(log_ndtr(lower)), float(log_ndtr(upper))
    result = b + np.log(-np.expm1(a - b))
    if not np.isfinite(result):
        raise FloatingPointError('prior box mass not representable in log space')
    return float(result)


def freeze(path, protocol):
    """Exclusive registration: a prior attempt requires reviewed recovery."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(plain(protocol), stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        import os
        os.fsync(stream.fileno())


def grid(order_a, order_b, s):
    ax, aw = leggauss(order_a)
    bx, bw = leggauss(order_b)
    ax, aw, xi, xw = .5 * ax, .5 * aw, 6 * bx, 6 * bw
    aa1, aa2, xx = np.meshgrid(ax, ax, xi, indexing='ij')
    wa1, wa2, wx = np.meshgrid(aw, aw, xw, indexing='ij')
    return (np.column_stack((aa1.ravel(), aa2.ravel())), xx.ravel(),
            (wa1 * wa2 * wx * s).ravel(), ax, aw, xi, xw)


def evaluate_arrays(target, y, q0, Q3, s, a, xi, purpose):
    """Keep one Point at a time; Target persists each paid state immediately."""
    n = len(xi)
    # q, material, prediction and all scalar metadata. Do not retain LU/Points.
    required = n * (66 * 8 + 576 * 16 + 32 * 8 * 16 + 8 * 8 + 128)
    if required > MAX_ARRAY_BYTES:
        raise MemoryError('raw tensor arrays exceed 400 MiB policy')
    q = np.empty((n, 66))
    chi = np.empty((n, 576), complex)
    prediction = np.empty((n, 32, 8), complex)
    Phi = np.empty(n)
    logprior3 = np.empty(n)
    keys = np.empty(n, dtype='U96')
    mu = -Q3.T @ q0[:64]
    for i in range(n):
        x = np.r_[a[i], s * xi[i]]
        q[i] = q0
        q[i, :64] += Q3 @ x
        p = target.evaluate(q[i], purpose=purpose)
        chi[i] = p.state.chi
        prediction[i] = p.state.data
        Phi[i] = target.phi(p, y)
        logprior3[i] = -.5 * np.sum((x - mu)**2) - 1.5 * np.log(2 * np.pi)
        keys[i] = p.key
        del p
        # Receiver transfers are not used here; keep no obsolete state owners.
        if target.model is not None:
            target.model._transfers.clear()
    return dict(q=q, chi=chi, prediction=prediction, Phi=Phi,
                logprior3=logprior3, log_gamma=logprior3-Phi,
                a=a, xi=xi, b=s*xi, endpoint_keys=keys)


def summarize_order(raw, na, nb, s, aw, xw, log_box_mass):
    """All posterior moments use precisely the same tensor quadrature measure."""
    lg = raw['log_gamma']
    lw = lg + np.log(raw['quadrature_weights'])
    logZ = float(logsumexp(lw))
    weights = np.exp(lw-logZ)
    raw['normalized_joint_weights'] = weights
    mean = weights @ raw['chi']
    # Centered accumulation avoids cancellation near deterministic decoder cells.
    variance = np.zeros(mean.shape, float)
    for w, image in zip(weights, raw['chi']):
        variance += w * np.abs(image-mean)**2
    c = raw['chi'].reshape(na*na, nb, -1)
    g = lg.reshape(na*na, nb)
    lbw = np.log(s*xw)
    logZb = logsumexp(g+lbw[None, :], axis=1)
    conditional_weights = np.exp(g+lbw[None, :]-logZb[:, None])
    conditional_mean = np.einsum('ij,ijk->ik', conditional_weights, c)
    # A quadrature-node maximum is a discrete profile, never an optimized mode.
    mode_index = np.argmax(g, axis=1)
    row = np.arange(na*na)
    keys = raw['endpoint_keys'].reshape(na*na, nb)
    aweights = (aw[:, None]*aw[None, :]).ravel()
    weak_weights = np.exp(logZb+np.log(aweights)-logZ)
    result = dict(material_mean=mean, material_variance=variance,
                  conditional_b_material_mean=conditional_mean,
                  conditional_b_weights=conditional_weights,
                  weak_log_Zb_untruncated_prior3=logZb,
                  weak_log_Zb_box_prior=logZb-log_box_mass,
                  weak_marginal_density=np.exp(logZb-logZ),
                  weak_quadrature_weights=aweights,
                  weak_posterior_weights=weak_weights,
                  discrete_profile_b=raw['b'].reshape(na*na, nb)[row, mode_index],
                  discrete_profile_log_gamma=g[row, mode_index],
                  discrete_profile_material=c[row, mode_index],
                  discrete_profile_endpoint_keys=keys[row, mode_index],
                  discrete_profile_at_edge=(mode_index == 0) | (mode_index == nb-1))
    summary = dict(a_order=na, b_order=nb, states=len(lg),
                   logZ_untruncated_prior3_box_integral=logZ,
                   logZ_conditioned_prior3_box=logZ-log_box_mass,
                   log_prior_box_mass=log_box_mass,
                   maximum_normalized_quadrature_weight=float(weights.max()),
                   variance_min=float(variance.min()),
                   conditional_b_mode_status='DISCRETE_PROFILE_BEST_PAID_NODE_NOT_OPTIMIZED',
                   reference_status='FINITE_TENSOR_QUADRATURE_REFINEMENT_NOT_CERTIFIED',
                   credible_intervals='NOT_COMPUTED_NO_NODE_BIN_WIDENING',
                   Gaussian_volume='NOT_COMPUTED_OPTIONAL',
                   original_FULL_posterior=False)
    return result, summary


def main():
    config = OUT/'configs/B_SLICE3_REFERENCE.json'
    if config.exists():
        raise RuntimeError('slice3 reference already registered; refusing implicit replay')
    anchor_path = OUT/f'A/anchors/{PARENT}.npz'
    # Require the already-paid receiver cache. Never fill a missing anchor J
    # with new RHS while preparing the frozen experimental registration.
    with np.load(anchor_path, allow_pickle=False) as source:
        if 'receiver_T' not in source.files or not source['receiver_T'].size:
            raise RuntimeError('paid anchor receiver_T required; no new anchor Jacobian authorized')
    book = Book('B', JOB, max_wall=3600)
    status, error, rows = 'COMPLETE', None, []
    dest = OUT/'B/slice3'
    results_path = dest/'RESULTS.json'
    try:
        target = Target(book)
        y = target.observation(PARENT)
        anchor = target.load_cache(anchor_path, owner='P3L A paid DEV_000 anchor')
        # Cached receiver_T means this is only array contraction, no RHS.
        J = target.jacobian(anchor)[:, :64].copy()
        frame = make_frame(anchor.q, J, 2, kind='posterior_material_risk',
                           material_derivative=target.prior.derivative(anchor.q))
        restricted = J @ frame['Vs']
        H = restricted.T @ restricted
        _, vectors = np.linalg.eigh(.5*(H+H.T))
        v = frame['Vs'] @ vectors[:, -1]
        Q3 = np.column_stack((frame['Vw'], v))
        orth_error = float(np.linalg.norm(Q3.T@Q3-np.eye(3)))
        if orth_error > 1e-10:
            raise ValueError('Q3 is not orthonormal')
        q0 = anchor.q.copy()
        s = float(1/np.sqrt(1+np.sum((J@v)**2)))
        mu = -Q3.T@q0[:64]
        intervals = [(-.5, .5), (-.5, .5), (-6*s, 6*s)]
        log_box_mass = sum(log_normal_interval(lo, hi, m)
                           for (lo, hi), m in zip(intervals, mu))
        protocol = dict(job=JOB, stage='B', parent=PARENT, max_wall_seconds=3600,
            target_family='SLICE3_BOX', original_FULL_posterior=False,
            target_definition='Original FullMaxwell likelihood times original Gaussian3 conditional prior, conditioned on declared box; other 61 Gaussian coordinates and both native properties fixed at anchor',
            original_P1_target_unchanged=True, GT_used=False, training=False,
            frame='posterior_material_risk dimension 2; v dominant restricted Jz.T Jz eigenvector in Vs',
            anchor=str(anchor_path), q0=q0, Q3=Q3, prior3_mean=mu,
            prior3_covariance=np.eye(3), intervals=intervals, s=s,
            log_prior_box_mass=log_box_mass, orthogonality_error=orth_error,
            reference_measure='da1 da2 db; b=s*xi and quadrature includes Jacobian s',
            orders=ORDERS, all_new_states_paid_and_immediately_saved=True,
            inner_Jacobian_queries=0, memory_array_cap_MiB=400,
            frozen_b_section='b=0 once for each distinct weak a across all orders',
            mode='best paid b quadrature node: discrete profile, not optimized mode',
            confidence='numerical refinement only; no FULL posterior or reference qualification',
            credible_intervals='NOT_COMPUTED', optional_Gaussian_volume='NOT_COMPUTED')
        # No fresh full-state evaluations occurred above. Register every choice
        # before the first new likelihood query, including q0, frame and box.
        freeze(config, protocol)
        dest.mkdir(parents=True, exist_ok=True)
        save(dest/'FRAME.npz', q0=q0, Q3=Q3, Jz=J, Vw=frame['Vw'],
             v=v, s=s, prior3_mean=mu, prior3_covariance=np.eye(3),
             log_prior_box_mass=log_box_mass)
        del anchor, frame, J, restricted, H, vectors
        if target.model is not None:
            target.model._transfers.clear()
        gc.collect()
        write(results_path, dict(protocol=protocol, status='RUNNING', rows=rows))
        previous_mean, previous_logZ = None, None
        zero_cache = {}
        for na, nb in ORDERS:
            label = f'a{na}_b{nb}'
            a, xi, weights, ax, aw, bx, bw = grid(na, nb, s)
            raw = evaluate_arrays(target, y, q0, Q3, s, a, xi,
                                  purpose=f'p3l_slice3_box_{label}')
            raw['quadrature_weights'] = weights
            raw['log_gamma_box_prior'] = raw['log_gamma']-log_box_mass
            moments, row = summarize_order(raw, na, nb, s, aw, bw, log_box_mass)
            mean = moments['material_mean']
            row['refinement_from_previous_order'] = None if previous_mean is None else dict(
                material_mean_RMS_difference=float(np.linalg.norm(mean-previous_mean)/np.sqrt(mean.size)),
                logZ_difference=float(row['logZ_untruncated_prior3_box_integral']-previous_logZ))
            previous_mean = mean.copy()
            previous_logZ = row['logZ_untruncated_prior3_box_integral']
            save(dest/f'{label}_RAW.npz', **raw)
            save(dest/f'{label}_MOMENTS.npz', **moments, weak_nodes=ax,
                 xi_nodes=bx, xi_quadrature_weights=bw)
            del raw, moments, a, xi, weights
            gc.collect()
            # Deduplicate exactly shared weak nodes (not approximate neighbors).
            pairs = np.array([(x, z) for x in ax for z in ax])
            missing = [p for p in pairs if tuple(p) not in zero_cache]
            if missing:
                ma = np.asarray(missing)
                section = evaluate_arrays(target, y, q0, Q3, s, ma,
                                          np.zeros(len(ma)), purpose=f'p3l_slice3_fixed_b0_{label}')
                save(dest/f'{label}_B0_NEW.npz', **section)
                for i, p in enumerate(ma):
                    zero_cache[tuple(p)] = dict(chi=section['chi'][i].copy(),
                        Phi=float(section['Phi'][i]), log_gamma=float(section['log_gamma'][i]),
                        endpoint=str(section['endpoint_keys'][i]))
                del section, ma
            save(dest/f'{label}_B0.npz', a=pairs,
                 chi=np.asarray([zero_cache[tuple(p)]['chi'] for p in pairs]),
                 Phi=np.array([zero_cache[tuple(p)]['Phi'] for p in pairs]),
                 log_gamma=np.array([zero_cache[tuple(p)]['log_gamma'] for p in pairs]),
                 endpoint_keys=np.array([zero_cache[tuple(p)]['endpoint'] for p in pairs]))
            row['frozen_b0_new_states'] = len(missing)
            rows.append(row)
            write(results_path, dict(protocol=protocol, status='RUNNING', rows=rows))
            print(json.dumps(plain(row), allow_nan=False), flush=True)
            if target.model is not None:
                target.model._transfers.clear()
            gc.collect()
        write(results_path, dict(protocol=protocol, status='COMPLETE', rows=rows,
            conclusion='Reported refinement differences only; no certified quadrature accuracy or original FULL posterior reference'))
    except BaseException as exc:
        status, error = 'ERROR', repr(exc)
        write(dest/'INTERRUPTED.json', dict(status=status, error=error, completed_orders=rows,
            endpoints='Every completed paid state remains in the Book job endpoints directory',
            replay='No automatic replay; review saved endpoints and remaining existing authorization before explicit recovery registration'))
        raise
    finally:
        book.close(status, error)


if __name__ == '__main__':
    main()

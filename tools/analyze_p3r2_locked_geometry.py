"""Offline oracle geometry audit, only after all frozen online estimates close.

No forward, LU, adjoint, JVP or hidden online proposal is constructed here.
Existing receiver/field arrays give an already-paid material linear response.
"""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from p1_material.prior import MaterialPrior
from p3r2.core import OUT, write, save, pack
from infer_p3r2_locked import frozen_contract


def main():
    ap = argparse.ArgumentParser();ap.add_argument('--root', type=Path, default=OUT)
    a = ap.parse_args();out = a.root
    manifest, parents, _, _ = frozen_contract(out, out/'METHODS_FROZEN_FOR_LOCKED.json')
    for name in ('estimates', 'trajectory_estimates', 'evaluation', 'trajectory_evaluation'):
        d = json.loads((out/'L'/name/'COMPLETE.json').read_text())
        if d.get('status') != 'COMPLETE':raise RuntimeError('Frozen completed evaluation required')
    dest = out/'L/offline_geometry'
    if dest.exists():raise FileExistsError('No overwrite of offline geometry results')
    dest.mkdir(parents=True)
    cfg = json.loads((out/'configs/TARGET.json').read_text())
    points = json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
    prior = MaterialPrior(cfg, points);sigma = cfg['noise']['sigma']
    rows = [];singular = [];eigenvalues = []
    for parent in parents:
        with np.load(out/'L/contexts'/f'{parent}.npz') as z:
            q, chi, J, pred = (z[k].copy() for k in ('q0', 'chi0', 'J', 'prediction'))
        # Oracle information is opened only in this isolated offline evaluator.
        with np.load(out/'inputs/offline_only'/f'{parent}.npz') as z:
            qt, ct, pt = (z[k].copy() for k in ('q', 'chi', 'prediction'))
        D = prior.derivative(q)
        M = np.concatenate((D.real, D.imag), axis=0)/24
        U, s, vh = np.linalg.svd(M, full_matrices=False)
        rank = int(np.sum(s > max(s[0]*1e-8, 1e-12)))
        Q = U[:, :rank];T = vh[:rank].T/s[:rank]
        K = J@T
        H = K.T@K + T[:64].T@T[:64]
        val, V = np.linalg.eigh((H+H.T)/2)
        h = ct-chi;hm = np.r_[h.real, h.imag]/24
        coordinates = V.T@(Q.T@hm)
        total = float(hm@hm)
        tangent = float(np.sum(coordinates**2))
        fractions = {str(radius):float(np.sum(coordinates[val <= radius**-2]**2)/max(total, 1e-300))
            for radius in (.02, .04, .08)}
        actual = pack(pt-pred)*np.sqrt(2)/sigma
        latent_linear = J@(qt-q)
        with np.load(out/'L/contexts/cache'/f'{parent}.npz') as z:
            receiver, field = z['receiver_T'], z['field']
        C = np.einsum('qp,ps->qsp', receiver, field).reshape(256,576)
        material_linear = pack((C@h).reshape(32,8))*np.sqrt(2)/sigma
        decoder_error = material_linear-latent_linear
        maxwell_error = actual-material_linear
        with np.load(out/'L/estimates'/f'{parent}.npz') as z:
            baseline = float(np.mean(abs(z['probit']-ct)**2))
            learned = float(np.mean(abs(z['fused']-ct)**2))
        row = dict(parent=parent, oracle_geometry_diagnostic=True,
            material_RMS=float(np.sqrt(total)),material_tangent_rank=rank,
            tangent_energy_fraction=tangent/max(total,1e-300),
            material_secant_outside_tangent_fraction=max(0.,1-tangent/max(total,1e-300)),
            prior_aware_weak_energy_fraction=fractions,
            full_whitened_prediction_distance=float(np.linalg.norm(actual)),
            latent_linear_response_norm=float(np.linalg.norm(latent_linear)),
            material_linear_response_norm=float(np.linalg.norm(material_linear)),
            decoder_response_remainder_norm=float(np.linalg.norm(decoder_error)),
            Maxwell_response_remainder_norm=float(np.linalg.norm(maxwell_error)),
            response_remainder_cross_inner_product=float(decoder_error@maxwell_error),
            Gaussian_prior_energy_change=.5*float(qt[:64]@qt[:64]-q[:64]@q[:64]),
            actual_properties_on_boundary=bool(np.any(q[64:]<1e-6) or np.any(q[64:]>1-1e-6)),
            fused_material_risk_improvement=baseline-learned)
        rows.append(row);singular.append(s);eigenvalues.append(np.pad(val,(0,66-rank),constant_values=np.nan))
        if len(rows)%32==0:print(json.dumps(dict(completed=len(rows),new_RHS=0)),flush=True)
    correlation = {}
    for radius in (.02, .04, .08):
        x = [r['prior_aware_weak_energy_fraction'][str(radius)] for r in rows]
        y = [r['fused_material_risk_improvement'] for r in rows]
        stat = spearmanr(x,y)
        correlation[str(radius)] = dict(Spearman=float(stat.statistic),pvalue=float(stat.pvalue),
            scope='Exploratory offline association; no GT-directed online proposals or causal attribution')
    save(dest/'SPECTRA.npz',parent=np.asarray(parents),material_singular_values=np.stack(singular),
         prior_aware_task_eigenvalues=np.stack(eigenvalues))
    write(dest/'RESULTS.json',dict(status='COMPLETE',parents=parents,rows=rows,
        material_metric='[ReDg;ImDg]/24',Jacobian='native66 fully whitened original material decoder chain rule',
        weak_rule='eigenvalue<=1/radius^2 of K.TK+Tz.TTz in orthonormal material tangent',
        rank_rule='material singular value > max(smax*1e-8,1e-12)',
        probability_claim=False,exact_fiber_claim=False,GT_used_online=False,new_RHS=0,
        association=correlation,limitations='Finite secant projection at a single point is not a reachable path or posterior mass; properties near bounds have a tangent cone rather than a free linear tangent'))
    write(dest/'COMPLETE.json',dict(status='COMPLETE',parents=parents,new_RHS=0,offline_truth_only=True))


if __name__=='__main__':main()

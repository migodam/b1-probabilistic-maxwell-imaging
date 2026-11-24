"""Outward negative-space norm of already enclosed physical current arrays."""
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from flint import arb, ctx
from scipy.linalg import svd
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
NUMERICS = HERE.parent
STUDY = HERE.parents[2]
NATIVE = NUMERICS/'native'
ARRAY = 'pw360_32x32__pair2pct'
SAVED = 'u100_best_conditioned.npz'
PROJECTOR = 'u100_best_strict.json'
DISK = f'control_disk075__{ARRAY}'
PHYSICAL_IDS = tuple(f'optfresh_{i:02d}__{ARRAY}__u1' for i in (1,4,5))
sys.path.insert(0,str(STUDY))
from core.certificate import ball, frobenius_squared, gershgorin, safe_float, star


def elapsed(start):
    w,c = start
    return dict(wall_seconds=perf_counter()-w,cpu_seconds=process_time()-c)


def norm(a):
    return frobenius_squared(a).sqrt()


def endpoints(a, b):
    return [max(0.,safe_float(a,upper=False)),safe_float(b,upper=True)]


def residual_source(stem):
    receipt = HERE/'strict_endpoints'/f'{stem}.json'
    proof = json.loads(receipt.read_text())
    if stem == DISK:
        trial_path = receipt.with_suffix('.npz')
    else:
        trial_path = HERE/proof.get('source',f'search_runs/{stem}.npz')
    return receipt, trial_path, proof


def audit(stem, v, v_ball, leak, gmin, gmax):
    start = perf_counter(), process_time()
    receipt, trial_path, proof = residual_source(stem)
    with np.load(trial_path,allow_pickle=False) as f:
        material, current = f['material'], f['current_approx']
    if (material.shape!=(9216,) or current.shape!=(9216,32) or not np.isrealobj(material)
            or not np.isfinite(material).all() or material.min()<0 or material.max()>min(1.,proof['cap'])):
        raise ArithmeticError(f'Illegal material or geometry: {stem}')
    count = int(np.count_nonzero(material))
    physics = proof['forward']
    if count!=physics['active_cells'] or float(material.max())!=proof['peak']:
        raise ArithmeticError(f'Forward enclosure provenance disagrees with saved array: {stem}')
    current[material==0] = 0  # The existing forward proof used this exact approximation.
    error = (arb(physics['strict_inverse']['inverse_spectral_norm_upper'])
             * arb(physics['primal_defect_frobenius_upper'])).upper()
    jb = ball(current)
    approx = star(v_ball)*jb
    n = norm(jb)
    a = norm(approx)
    perturbation = leak*n.upper()
    flo = (a.lower()-perturbation)/gmax.sqrt()-error
    fhi = (a.upper()+perturbation)/gmin.sqrt()+error

    # A normalized coherent input probes peak energy. Its direction was selected
    # from the numerical product, but all evaluation and normalization use Arb.
    _, _, vh = svd(v.conj().T@current,full_matrices=False,check_finite=False)
    direction = ball(vh[0].conj()[:,None])
    dn = norm(direction)
    ca = norm(approx*direction)/dn
    cj = norm(jb*direction)/dn
    slo = (ca.lower()-leak*cj.upper())/gmax.sqrt()-error
    gram_max = gershgorin(approx*star(approx),upper=True)
    shi = min(fhi.upper(), ((gram_max.sqrt()+perturbation)/gmin.sqrt()+error).upper())
    if not flo>0 or not slo>0:
        raise ArithmeticError('Reference failed to prove a nonzero physical projection.')
    return dict(observation=proof['observation'], test_cap=1., material_peak=proof['peak'],
                active_cells=count, material_status='box_feasible_with_verified_regular_forward',
                noise_ball_status=proof['status'], observation_radius=proof['radius'],
                original_forward_residual_interval=[proof['distance_lower_for_fixed_candidate'],
                                                   proof['distance_upper_for_fixed_candidate']],
                arrays=str(trial_path.relative_to(HERE)), forward_receipt=str(receipt.relative_to(HERE)),
                full_current_frobenius_interval=endpoints(n.lower()-error,n.upper()+error),
                current_approximation_error_F_upper=safe_float(error,upper=True),
                ritz_projection_frobenius_interval=endpoints(a.lower(),a.upper()),
                negative_current_frobenius_interval=endpoints(flo,fhi),
                negative_current_frobenius_squared_interval=endpoints(flo.lower()**2,fhi.upper()**2),
                negative_current_spectral_norm_interval=endpoints(slo,shi),
                negative_current_spectral_norm_squared_interval=endpoints(slo.lower()**2,shi.upper()**2),
                **elapsed(start))


def main():
    total = perf_counter(), process_time()
    started = perf_counter(), process_time()
    ctx.prec,ctx.threads = 128,1
    replay = json.loads((NATIVE/PROJECTOR).read_text())
    leakage, eig = replay['projector'], replay['verified_state']
    with np.load(NATIVE/SAVED,allow_pickle=False) as f:
        v = f['negative_vec']
        values = f['negative_eig']
    if v.shape!=(9216,eig['exact_nu']) or eig['exact_nu']!=13 or np.any(values>=0):
        raise ArithmeticError('Unexpected saved negative basis.')
    v_ball = ball(v)
    gmin = arb(leakage['negative_projected_gram_min_lower'])
    gmax = gershgorin(star(v_ball)*v_ball,upper=True)
    ell = arb(leakage['ritz_positive_leak_upper'])
    if not gmin > 0:
        raise ArithmeticError('The saved projection bound does not imply full column rank.')
    setup = elapsed(started)
    rows = [audit(row,v,v_ball,ell,gmin,gmax) for row in (DISK,*PHYSICAL_IDS)]
    uniform_frob_floor = max(row['negative_current_frobenius_squared_interval'][0] for row in rows)
    uniform_source_matrix_floor = max(row['negative_current_spectral_norm_squared_interval'][0] for row in rows)
    result = dict(contract='96² 400 MHz 2 m full lossless real box [0,1], fixed 32 source and 32 receiver arrays',
                  current_units='Unweighted native discrete Euclidean norm over cells, with fixed total-fluence source normalization. Frobenius sums all 32 columns; spectral norm maximizes over unit coherent source inputs.',
                  projector_profile='scalar_grid_601', projection_center='uncentered', negative_rank=eig['exact_nu'],
                  arb_precision_bits=ctx.prec, negative_basis=f'../native/{SAVED}',
                  exact_projector_receipt=f'../native/{PROJECTOR}', ritz_positive_leak_frobenius_upper=safe_float(ell,upper=True),
                  negative_projected_gram_min_lower=safe_float(gmin,upper=False),
                  negative_projected_gram_max_upper=safe_float(gmax,upper=True),
                  necessary_uniform_full_box_frobenius_energy_cap_lower=uniform_frob_floor,
                  necessary_uniform_full_box_unit_coherent_energy_cap_lower=uniform_source_matrix_floor,
                  setup=setup, rows=rows, new_forward_solves=0, new_inverse_solves=0,
                  new_green_calls=0, optimizer_runs=0, **elapsed(total))
    (HERE/'hidden_current.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for row in rows:
        print(row['observation'], row['negative_current_frobenius_interval'],
              row['negative_current_spectral_norm_interval'],flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

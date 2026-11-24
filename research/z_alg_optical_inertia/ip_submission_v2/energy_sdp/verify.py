"""Arb proof for energy-capped four-region single-source witnesses on the 96² grid."""

import argparse
import csv
import json
import sys
from contextlib import ExitStack
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from probe import HERE, blocks, load
from local_moment import (ExactGeometry, acb_mat, arb, ball, certify_positive, ctx,
                          exact_local_weight_matrices, gershgorin,
                          h, hi, mid, norm_upper, star)

sys.path.insert(0, str(HERE.parent))
from energy_port_control import solve as port_only_bound

DEFAULT_GRID = (0, 10, 100, 250, 350, 400, 1000, 2000, 3000, 5000, 7000, 7500,
                8500, 10000, 15000, 26000, 1e5, 1e6, 1e8)


def eigen_upper(a):
    v = ball(eigh(h(mid(a)), check_finite=False)[1])
    gram_floor = gershgorin(star(v) * v, upper=False)
    if not gram_floor > 0:
        raise ArithmeticError('Unable to verify proposed source eigencolumns.')
    return max(arb(0), gershgorin(star(v) * a * v, upper=True)) / gram_floor


def verify_case(exact, x0, y0, z0, anchor_weight, fit_scale):
    wall, cpu = perf_counter(), process_time()
    nrays = z0.shape[1]
    dark = exact.project_to_data(z0, acb_mat(exact.s.nrows(), nrays))
    dark_sector, dark_optical = exact.forms(dark, exact.green(dark))
    rho = [[q[j,j].real for j in range(nrays)] for q in dark_sector]
    tau = [[q[j,j].real for j in range(nrays)] for q in dark_optical]
    y = arb(float(fit_scale))*ball(y0)
    x = exact.project_to_data(x0, y)
    sector, optical = exact.forms(x, exact.b + exact.green(x))
    fixed_weights = []
    if nrays > 4:
        identity = acb_mat(np.eye(y.ncols(), dtype=int).tolist())
        auxiliary = (arb(float(anchor_weight))/(nrays-4))*identity
        optical = [a + sum(tau[i][4:],arb(0))*auxiliary for i,a in enumerate(optical)]
        fixed_weights = [auxiliary]*(nrays-4)
    candidate_system = acb_mat(4,4,[tau[i][j] for i in range(4) for j in range(4)])
    weights = exact_local_weight_matrices(candidate_system, optical) + fixed_weights
    positive = [certify_positive(w) for w in weights]
    hull = []
    for i, q in enumerate(sector):
        for k, w in enumerate(weights):
            q += rho[i][k]*w
        hull.append(-certify_positive(-q))
    dark_vectors = star(dark) * dark
    covariance = acb_mat(y.ncols(), y.ncols())
    for k, w in enumerate(weights):
        covariance += dark_vectors[k, k].real * w
    mean_gram = star(x)*x
    full_moment = mean_gram + covariance
    mean_limit = eigen_upper(mean_gram)
    covariance_limit = eigen_upper(covariance)
    energy = eigen_upper(full_moment)
    energy_upper = hi(energy)
    evidence = {'weight_eigenvalue_floors': positive,
                     'local_optical_midpoints': [[float(q.mid()) for q in row] for row in tau],
                     'local_sector_midpoints': [[float(q.mid()) for q in row] for row in rho],
                     'interval_dark_receiver_residual_upper':norm_upper(exact.s*dark),
                     'fit_scale': fit_scale, 'fixed_covariance_coefficient_sum':anchor_weight,
                     'sector_eigenvalue_ceilings': hull,
                     'mean_current_energy_upper': hi(mean_limit),
                     'covariance_current_energy_upper': hi(covariance_limit),
                     'sum_of_separate_energy_uppers': hi(mean_limit+covariance_limit),
                     'uniform_witness_current_energy_upper':energy_upper,
                     'rank_upper':len(weights)+1,
                     'witness_second_moment_eigenvalues_proposal': np.linalg.eigvalsh(h(mid(full_moment))).tolist(),
                     'interval_mean_receiver_residual_upper':norm_upper(exact.s*x-y),
                     'columns_of_interval_green_evaluations': x.ncols()+dark.ncols(),
                     'receiver_null_projection_rhs_columns': x.ncols()+dark.ncols(),
                     'seconds_wall':perf_counter()-wall, 'seconds_cpu':process_time()-cpu}
    return energy_upper, evidence


def interpolation_upper(y_squared, cap, points):
    energy_eligible = [p for p in points if cap is None or p[0] <= cap]
    if any(p[1] == 1 for p in energy_eligible):
        return 0., 'VERIFIED_EXACT_ZERO'
    best = min((1-arb(scale))**2 for _, scale in energy_eligible)
    for e0, s0 in energy_eligible:
        for e1, s1 in points:
            if e1 <= cap:
                continue
            t = (arb(cap)-arb(e0))/(arb(e1)-arb(e0))
            candidate = (1-t)*(1-arb(s0))**2 + t*(1-arb(s1))**2
            best = min(best.upper(),candidate.upper())
    return hi((best*y_squared).sqrt()), 'BOUNDED_NOT_SOLVED'


def run(args):
    wall,cpu=perf_counter(),process_time()
    s,b,symbol,data,base=load()
    ctx.prec,ctx.threads = 160,1
    exact=ExactGeometry(s,b,symbol,blocks(2))
    cases = {}
    g_eigenvalues, g_eigenvectors = eigh(s@s.conj().T)
    ports = exact.gram, g_eigenvalues, g_eigenvectors
    with ExitStack() as ctxs:
        sources=[(file.name,ctxs.enter_context(np.load(file)))
                 for file in [args.proposals, *args.more]]
        for case,y in data.items():
            points = [(0.,0.)]
            certificates = {}
            for name,proposal in sources:
                z0 = proposal[case+'__rays'] if case+'__rays' in proposal else proposal['rays']
                anchor_weight = 0.
                if case+'__auxiliary_ray' in proposal:
                    z0 = np.column_stack((z0,proposal[case+'__auxiliary_ray']))
                    anchor_weight = float(proposal[case+'__auxiliary_weight'])
                scale = float(proposal[case+'__fit_scale']) if case+'__fit_scale' in proposal else 1.
                energy, details = verify_case(exact,proposal[case],y,z0,anchor_weight,scale)
                certificates[name] = details
                points.append((energy,scale))
                print(case,name,'uniform trace bound',energy,'weight floor',min(details['weight_eigenvalue_floors']),
                      'sector ceiling',max(details['sector_eigenvalue_ceilings']),flush=True)
            y_gram = star(ball(y))*ball(y)
            y_squared = hi(eigen_upper(y_gram))
            radius = next(row['inflated_radius'] for row in base['observations'] if row['id']==case)
            grid = []
            cap_grid = sorted(set(args.caps))
            if args.physical_caps:
                cap_grid = sorted(set(cap_grid + args.physical_caps[case.split('__')[0]]))
            for cap in cap_grid + [None]:
                upper, status = interpolation_upper(y_squared,cap,points)
                low = port_only_bound(ports,y,float(cap))['lower'] if cap is not None else 0.
                if cap == 0:
                    status = 'VERIFIED_PORT_NORM'
                if low > upper:
                    raise ArithmeticError('Inconsistent certified bounds.')
                grid.append({'energy_cap':cap, 'residual_lower':low, 'residual_upper':upper,
                             'status':status, 'residual_at_most_inflated_radius':upper <= radius,
                             'residual_strictly_exceeds_inflated_radius':low > radius})
            cases[case] = {'observed_port_energy_upper':y_squared,'inflated_radius':radius,
                           'proofs': certificates, 'energy_grid':grid}
    args.output.write_text(json.dumps({'partition':[2,2], 'material_upper':.75, 'interval_bits':ctx.prec,
                                       'energy_units':'unweighted discrete Euclidean squared current on the 96×96 grid',
                                       'objective_units':'receiver Euclidean residual norm',
                                       'coherent_max':'supremum over complex a with Euclidean norm 1',
                                       'exact_receiver_projection_solves':exact.receiver_solves,
                                       'exact_receiver_projection_rhs_columns':exact.receiver_rhs_columns,
                                       'full_material_forward_solves':0,
                                       'observations':cases,'seconds_wall':perf_counter()-wall,
                                       'seconds_cpu':process_time()-cpu}, indent=2, allow_nan=False)+'\n')
    if args.csv is not None:
        with args.csv.open('w',newline='') as output:
            names = ['observation',*next(iter(cases.values()))['energy_grid'][0].keys()]
            writer = csv.DictWriter(output, fieldnames=names)
            writer.writeheader()
            for case, details in cases.items():
                for row in details['energy_grid']:
                    writer.writerow({'observation':case, **row,
                                     'energy_cap':row['energy_cap'] if row['energy_cap'] is not None else 'inf'})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=type(HERE),default=HERE/'verified.json')
    parser.add_argument('--proposals',type=type(HERE),default=HERE/'alt025.npz')
    parser.add_argument('--more',type=type(HERE),nargs='*',default=[])
    parser.add_argument('--csv',type=type(HERE))
    parser.add_argument('--physical-caps',type=json.loads)
    parser.add_argument('--caps',type=float,nargs='+',default=DEFAULT_GRID)
    with threadpool_limits(1):
        run(parser.parse_args())

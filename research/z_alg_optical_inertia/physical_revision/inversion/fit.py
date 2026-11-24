"""Compare bounded full-wave inverse policies on frozen public observations."""

import argparse
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.sparse import coo_matrix
from threadpoolctl import threadpool_limits

from gate import HERE, STUDY, elapsed, clocks, original, write

SOURCE = STUDY.parent/'z_alg_xudong_cie/minpeak_native/search.py'
spec = spec_from_file_location('mature_full_wave_model', SOURCE)
native = module_from_spec(spec)
spec.loader.exec_module(native)
STAGES = ((16, 20), (32, 30), (48, 44))
RELEASED_CAP = 4.


class Grid:
    def __init__(self, side):
        self.side = side
        self.n = side**2
        z = np.linspace(0., side - 1., native.SIDE)
        k = np.minimum(z.astype(int), side - 2)
        v = z - k
        ids = np.arange(native.N)
        columns, values = [], []
        for dy, wy in enumerate((1-v, v)):
            for dx, wx in enumerate((1-v, v)):
                columns.append(((k[:, None]+dy)*side + k[None, :]+dx).ravel())
                values.append((wy[:, None]*wx[None, :]).ravel())
        self.matrix = coo_matrix((np.concatenate(values), (np.tile(ids, 4), np.concatenate(columns))),
                                 shape=(native.N, self.n)).tocsr()
        self.edge_n = 2*side*(side-1)

    def smoothness(self, x):
        shape = x.reshape(self.side, self.side)
        dx, dy = np.diff(shape, axis=1), np.diff(shape, axis=0)
        weight = .001/self.edge_n
        ridge = .000001/self.n
        cost = .5*weight*(np.sum(dx**2)+np.sum(dy**2))+.5*ridge*np.dot(x, x)
        gradient = ridge*shape
        gradient[:, 1:] += weight*dx
        gradient[:, :-1] -= weight*dx
        gradient[1:] += weight*dy
        gradient[:-1] -= weight*dy
        return cost, gradient.ravel()

    def from_previous(self, x):
        if x is None:
            return np.zeros(self.n)
        prev = int(np.sqrt(len(x)))
        points = np.linspace(0., 1., prev)
        broader = np.linspace(0., 1., self.side)
        coordinates = np.stack(np.meshgrid(broader, broader, indexing='ij'), axis=-1)
        return RegularGridInterpolator((points, points), x.reshape(prev, prev))(coordinates).ravel()


class StageLimit(Exception):
    pass


def search_stage(geometry, obs, initial, side, limit, cap):
    timer = clocks()
    model = native.Model(*geometry, obs.y)
    grid = Grid(side)
    previous = None
    trace = []
    best_loss, best_residual = None, None

    def objective(x):
        nonlocal previous, best_loss, best_residual
        if len(trace) >= limit:
            raise StageLimit
        m = np.clip(grid.matrix@x, 0., cap)
        try:
            data, gradient, discrepancy, state, ports = model.full(m, previous)
        except ArithmeticError as e:
            previous = None
            trace.append(dict(state_error=str(e)))
            if best_loss is None:
                raise
            d = x-best_loss['x']
            return 1.+.01*np.dot(d,d), .02*d
        reg, derivative = grid.smoothness(x)
        previous = state
        loss = float(data+reg)
        peak = float(m.max())
        trace.append(dict(residual=float(discrepancy), residual_over_radius=float(discrepancy/obs.padded_radius),
                          regularizer=float(reg), objective=loss, peak=peak))
        if best_loss is None or loss < best_loss['objective']:
            best_loss = dict(x=x.copy(), material=m, ports=ports, **trace[-1])
        if best_residual is None or discrepancy < best_residual['residual']:
            best_residual = dict(x=x.copy(), material=m, ports=ports, **trace[-1])
        return loss, grid.matrix.T@gradient+derivative

    x0 = np.clip(grid.from_previous(initial), 0., cap)
    try:
        result = minimize(objective, x0, method='L-BFGS-B', jac=True, bounds=[(0.,cap)]*grid.n,
                          options=dict(maxiter=limit, maxfun=limit, gtol=1e-11, ftol=2e-13,
                                       maxls=24, maxcor=12))
        stop = str(result.message)
    except StageLimit:
        stop = 'objective-call limit'
    if best_loss is None:
        raise ArithmeticError('No successful full-wave fit in stage.')
    row = dict(side=side, cap=cap, call_budget=limit, actual_calls=len(trace), counts=dict(model.count),
               selected_call=next(i for i,r in enumerate(trace)
                                  if r.get('objective')==best_loss['objective']),
               best_residual_call=next(i for i,r in enumerate(trace)
                                       if r.get('residual')==best_residual['residual']),
               trace=trace, stop=stop, **elapsed(timer))
    return best_loss, best_residual, row


def recheck(geometry, obs, selected):
    old, native.RTOL = native.RTOL, 1.e-10
    model = native.Model(*geometry, obs.y)
    timer = clocks()
    try:
        m = selected['material']
        value, ports, currents = model.full(m, gradient=False)
        rhs = m[:,None]*model.b
        defect = currents-m[:,None]*model.green(currents)-rhs
    finally:
        native.RTOL = old
    return ports, dict(norm_residual=float(value), residual_over_radius=float(value/obs.padded_radius),
                       rescan_port_difference=float(np.linalg.norm(ports-selected['ports'])),
                       state_equation_relative_defect=float(np.linalg.norm(defect)/max(1e-20,np.linalg.norm(rhs))),
                       state_equation_max_source_relative_defect=float(np.max(
                           native.norm_columns(defect)/np.maximum(1e-20,native.norm_columns(rhs)))),
                       counts=dict(model.count), **elapsed(timer))


def run(geometry, obs, gate, out, *, stages=STAGES, released_cap=RELEASED_CAP):
    timer = clocks()
    out.mkdir(parents=True, exist_ok=True)
    low_cap = gate['cap']
    stage_cache, paths = {}, {}
    low_stage = search_stage(geometry, obs, None, *stages[0], low_cap)
    stage_cache[('low',)] = low_stage
    ratio = low_stage[1]['residual_over_radius']
    discrepancy_release = ratio > 2.5
    policies = {'fixed_low': ('low','low','low'),
                'high_from_start': ('high','high','high'),
                'certificate': ('high','high','high') if gate['strict_exclusion'] else ('low','low','low'),
                'discrepancy': ('low','high','high') if discrepancy_release else ('low','low','low')}
    outputs = {}
    for name, sequence in policies.items():
        previous, path_stages = None, []
        for j, (which, (side, budget)) in enumerate(zip(sequence, stages, strict=True)):
            key = sequence[:j+1]
            if key == ('low',):
                item = low_stage
            else:
                if key not in stage_cache:
                    cap = low_cap if which=='low' else released_cap
                    item = search_stage(geometry, obs, previous, side, budget, cap)
                    stage_cache[key] = item
                    print(obs.id, key, item[1]['residual_over_radius'], item[2]['wall_seconds'], flush=True)
                item = stage_cache[key]
            loss, residual, row = item
            previous = loss['x']
            path_stages.append(row)
        if sequence not in paths:
            audit_ports, audit = recheck(geometry, obs, loss)
            if np.array_equal(loss['material'], residual['material']):
                min_ports, min_audit = audit_ports, audit
            else:
                min_ports, min_audit = recheck(geometry, obs, residual)
            paths[sequence] = (audit_ports, audit, min_ports, min_audit)
        audit_ports, audit, min_ports, min_audit = paths[sequence]
        outputs.update({name+'_material':loss['material'], name+'_ports':audit_ports,
                        name+'_minres_material':residual['material'], name+'_minres_ports':min_ports})
        profiles = tuple(policies)
        reports = outputs.setdefault('reports', {})
        reports[name] = dict(cap_schedule=[low_cap if term=='low' else released_cap for term in sequence],
                             shared_path=next(n for n in profiles if policies[n]==sequence),
                             stages=path_stages, audit=audit, audit_minres=min_audit,
                             peak=loss['peak'], minres_peak=residual['peak'],
                             charge_fit_wall=sum(r['wall_seconds'] for r in path_stages),
                             charge_fit_cpu=sum(r['cpu_seconds'] for r in path_stages),
                             charge_full_objective_calls=sum(r['actual_calls'] for r in path_stages))
    reports = outputs.pop('reports')
    arrays = out/(obs.id+'.npz')
    np.savez_compressed(arrays, Y=obs.y, **outputs)
    result = dict(id=obs.id, cap=low_cap, gate=gate, radius=obs.padded_radius,
                  discrepancy_stage1_minimum_ratio=ratio, discrepancy_releases=discrepancy_release,
                  distinct_prefixes=len(stage_cache), reports=reports, **elapsed(timer))
    write(arrays.with_suffix('.json'), result)
    print('DONE',obs.id,[(k,round(v['audit']['residual_over_radius'],4)) for k,v in reports.items()],flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append')
    parser.add_argument('--out', type=Path, default=HERE/'original_fits')
    parser.add_argument('--bootstrap', action='store_true', help='Explicit v2 data-only gate decisions for the three known eligible cases.')
    opt = parser.parse_args()
    geometry, original_obs = original()
    array = {r.id:r for r in original_obs}
    if opt.bootstrap:
        best = {}
        from gate import inherited_distances
        for oid, cap, d, src in inherited_distances(original_obs):
            if d > best.get((oid,cap),{}).get('strict_distance_lower',0.):
                r = array[oid].padded_radius
                best[(oid,cap)] = dict(id=oid, cap=cap, strict_distance_lower=d,
                                       strict_ratio=d/r, strict_exclusion=d>r, source='v2/'+src,
                                       requested_window_impossible=d>2.5*r)
        candidates = list(best.values())
    else:
        candidates = json.loads((HERE/'original_gate.json').read_text())['decisions']
    picked = {}
    for row in candidates:
        if row['strict_exclusion'] and not row['requested_window_impossible']:
            oid = row['id']
            if (not opt.case or oid in opt.case) and row['cap']>picked.get(oid,{}).get('cap',0.):
                picked[oid] = row
    for name, candidate in picked.items():
        if (opt.out/(name+'.npz')).exists() and (opt.out/(name+'.json')).exists():
            continue
        print('START', name, candidate['cap'], candidate['strict_ratio'], flush=True)
        run(geometry, array[name], candidate, opt.out)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

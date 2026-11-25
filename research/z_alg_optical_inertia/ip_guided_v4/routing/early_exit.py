#!/usr/bin/env python3
"""Compare one strict current-iterate neighborhood gate with plain discrepancy and restart."""
import sys
from argparse import ArgumentParser

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from routing_run import HERE, clock, elapsed, other_suite, module, read, write
from cie import branch


def verifier():
    tools = HERE.parent/'boxes'
    box_run = module('routing_box_run',tools/'run.py')
    previous = sys.modules['run']
    try:
        sys.modules['run'] = box_run
        audit = module('routing_box_verify',tools/'verify.py')
    finally:
        sys.modules['run'] = previous
    return box_run,audit


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--halfwidth',type=float,default=.04)
    argv = parser.parse_args()
    ctx.prec,ctx.threads = 128,1
    box_api, ver = verifier()
    out = HERE/('cie_early_exit' if argv.halfwidth==.04 else 'cie_early_exit_narrow')
    out.mkdir(exist_ok=True)
    for geo,obs in other_suite('mild_mid'):
        if obs.object_id not in ('guided_mildmid_03','guided_mildmid_04'):
            continue
        archive = HERE/'online'/(obs.id+'.npz')
        with np.load(archive,allow_pickle=False) as a:
            born = a['born_only_material']
        born_cost = read(archive.with_suffix('.json'))['born_linear_fit']
        zero = np.zeros(len(geo.b))
        r_stop, m_stop = branch(geo,obs,'first_four_steps',zero,8,4,8.,mf=3)
        r_full,m_full = branch(geo,obs,'continue_to_sixteen',zero,8,16,8.,mf=3)
        t = clock()
        x = m_stop['4_material']
        lower,upper = np.maximum(0.,x-argv.halfwidth),np.minimum(20.,x+argv.halfwidth)
        box_geo = box_api.Geometry(geo.s,geo.b,geo.side,G=geo.g)
        physical = box_api.Observation(obs.id,box_geo,obs.y,obs.radius,obs.source)
        grid = ver.ProofGrid(box_geo)
        mu = 10*float(np.linalg.norm(geo.s,ord=2))**2
        arrays = {}
        try:
            strict,arrays = ver.prove(physical,grid,np.arange(len(x)),lower,upper,mu)
        except (ArithmeticError,ValueError,ZeroDivisionError,np.linalg.LinAlgError) as exc:
            strict = dict(failure=str(exc),strictly_excludes=False)
        strict.update(total_new_cost=elapsed(t))
        r_restart,m_restart = branch(geo,obs,'born_restart_twelve_after_four',born,8,12,8.,mf=3)
        point = r_stop['audits'][-1]['native']['residual']
        r_simple = 'restart' if point>obs.radius else 'continue'
        r_l4 = 'restart' if strict['strictly_excludes'] else 'continue'
        def charged(named,include_certificate):
            costs = dict()
            for k in ('wall_seconds','cpu_seconds'):
                base = r_stop[k]+r_restart[k]+born_cost[k] if named=='restart' else r_full[k]
                costs[k] = base+(strict['total_new_cost'][k] if include_certificate else 0.)
            return costs
        final = m_full[str(r_full['step'])+'_material']
        outside = (final<lower)|(final>upper)
        row=dict(id=obs.id,object_id=obs.object_id,radius=obs.radius,
                 material_box_halfwidth=argv.halfwidth,
                 material_box_restriction='144 all active in the symmetric box clipped to [0,20]',
                 arms=dict(first_checkpoint=r_stop,continue_to_16=r_full,born_restart=r_restart),
                 born_linear_fit_cost=born_cost,
                 strict_box= strict,one_point_relative_discrepancy=point/obs.radius,
                 unrestricted_continuation_linf_from_checkpoint=float(np.max(abs(x-final))),
                 unrestricted_continuation_cells_outside_box=int(outside.sum()),
                 simple_discrepancy=dict(selected=r_simple,charged=charged(r_simple,False)),
                 l4_rule=dict(selected=r_l4,charged=charged(r_l4,True)),
                 deterministic_original_checkpoint_replayed=True,
                 scope='A bounded local neighborhood probe. Future unrestricted CIE iterates may leave this neighborhood.')
        all_arrays = {label+'_'+k:v for label,raw in [('prefix',m_stop),('continue',m_full),('restart',m_restart),
                                                       ('strict',arrays)] for k,v in raw.items()}
        np.savez_compressed(out/(obs.id+'.npz'),**all_arrays)
        write(out/(obs.id+'.json'),row)
        print(obs.id,'L4',strict.get('distance_over_radius_lower'),r_l4,'simple',r_simple,
              'neumann',strict.get('center',{}).get('neumann_residual_over_radius_lower'),flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

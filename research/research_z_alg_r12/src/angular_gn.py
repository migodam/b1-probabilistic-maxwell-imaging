"""One fixed coherent exposure followed by a matched physical GN reconstruction."""
from argparse import ArgumentParser
import json
import os
from pathlib import Path
import sys
from time import perf_counter, process_time

for name in ('OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from angular_gap import RINGS, cost, layouts, norm
from weak_r4 import DEFAULT_ARCHIVE, HERE, PARENT_IDS, R4, archived_material_score, read_npz


def fit(Model, reg, ops, m, target, scale, view_noise, max_steps):
    t0, c0 = perf_counter(), process_time()
    model = Model(ops)
    state = model.forward(m)
    damping = 3e-9
    trace = []
    original_source_count = 12

    def score(z):
        residual = z['y']-target
        prior = reg.evaluate(z['m'])[0]
        old = norm(residual[:,:original_source_count])/(view_noise*np.sqrt(original_source_count))
        new = norm(residual[:,original_source_count:])/view_noise if target.shape[1]>original_source_count else None
        return dict(phi=.5*norm(residual/scale)**2 + prior, old_noise_multiple=old, new_noise_multiple=new)

    val = score(state)
    trace.append(val)
    for _ in range(max_steps):
        J = model.jacobian(state)/scale
        r = Model.pack(state['y']-target)/scale
        _, g_prior, w_prior = reg.evaluate(state['m'])
        g = J.T@r + g_prior
        matrix = J.T@J + (reg.D.T @ reg.D.multiply(w_prior[:,None])).toarray()
        accepted = False
        for _ in range(4):
            p = cho_solve(cho_factor(matrix+damping*np.eye(m.size), check_finite=False), -g, check_finite=False)
            for a in (1.,.5,.25,.125):
                trial = model.forward(np.clip(m+a*p,.02,20.))
                trial_value = score(trial)
                if trial_value['phi'] < val['phi']:
                    state, val, m, accepted = trial, trial_value, trial['m'], True
                    damping = max(1e-12, .35*damping)
                    break
            if accepted:
                break
            damping *= 10
        if not accepted:
            break
        trace.append(val)
    return m, dict(accepted_steps=len(trace)-1, trace=trace, **val, cost=cost(t0,c0,**model.counts))


def run(args):
    sys.path.insert(0,str(args.archive.resolve()))
    from research.research_z_alg_r4_diversity.src.physics import Model, Regularizer
    history = args.archive/R4
    ops = read_npz(history/'data/operators.npz')
    weak_summary = json.loads((HERE/'results/weak_r4.json').read_text())
    by_case = {r['parent_id']:r for r in weak_summary['cases']}
    acquisition_manifest = json.loads((history/'data/manifests/evaluation.json').read_text())
    records = {r['obs_id']:r for r in acquisition_manifest['records']}
    reg = Regularizer(ops['grid_shape'], smooth=1e-5,tv=1e-5,epsilon=.01)
    geometry = layouts(ops['centers'],float(ops['k']))
    outputs = {}
    result = dict(config=dict(cases=args.cases, one_exposure=True, emitters=24,
                               common_start='archived all-scales discrepancy GN anchor',
                               extra_warm_start='best fit of two saved feasible candidate points to new exposure only',
                               fixed_iterations=args.steps, damping0=3e-9, prior=dict(smooth=1e-5,tv=1e-5,epsilon=.01),
                               scaling='same ||old Y|| for all data blocks and arms, fixed absolute complex noise per exposure',
                               truth_usage='Archived GT loaded only after fixed source selection, used for synthetic new observation and scoring. Same perturbation vector replayed in every exposure within each case.'), cases={})
    with np.load(args.feeds,allow_pickle=False) as feeds, np.load(args.materials,allow_pickle=False) as materials:
        for case in args.cases:
            # These feeds are frozen stored numerical choices from angular_gap.
            fields, arms = [], []
            for R in RINGS:
                for method in ('single','coherent_l2_2','coherent_l2_4'):
                    arm = f'{R:.2f}_24_{method}'
                    fields.append(geometry[R,24]['E0']@feeds[case+'_'+arm])
                    arms.append(arm)
            fields = np.column_stack(fields)
            row = records[by_case[case]['obs_id']]
            truth = read_npz(history/'data'/row['truth_path'])
            old, anchor = (materials[case+'_'+s] for s in ('old_Y','anchor'))
            candidates = [materials[case+'_'+s] for s in ('plus','minus')]
            scale = norm(old)
            per_view = .001*scale/np.sqrt(old.shape[1])
            rng = np.random.default_rng([2026111911,int(case.split('p')[-1])])
            noise = rng.normal(size=old.shape[0])+1j*rng.normal(size=old.shape[0])
            noise *= per_view/norm(noise)
            t0, c0 = perf_counter(), process_time()
            simulation = Model(dict(ops,E0=fields))
            data = simulation.forward(truth['m'])['y']+noise[:,None]
            predictions = np.stack([simulation.forward(candidate)['y'] for candidate in candidates])
            case_result = dict(simulation_and_candidate_selection_cost=cost(t0,c0,**simulation.counts), arms={})
            old_m, old_stats = fit(Model, reg, ops, anchor, old, scale, per_view, args.steps)
            outputs[case+'_old_only'] = old_m
            case_result['old_only'] = dict(**old_stats, **archived_material_score(old_m,truth),
                                           initial=archived_material_score(anchor,truth))
            for i, name in enumerate(arms):
                Y = np.column_stack((old, data[:,i]))
                arm_ops = dict(ops,E0=np.column_stack((ops['E0'],fields[:,i])))
                discrepancies = np.linalg.norm(predictions[:,:,i]-data[:,i],axis=1)/per_view
                winner = int(np.argmin(discrepancies))
                info = dict(discrepancies=discrepancies.tolist(), selected=('plus','minus')[winner],
                            winner_only=archived_material_score(candidates[winner],truth))
                for label, m0 in [('anchor_gn',anchor),('winner_gn',candidates[winner])]:
                    m, stats = fit(Model,reg,arm_ops,m0,Y,scale,per_view,args.steps)
                    outputs[f'{case}_{name}_{label}'] = m
                    info[label] = dict(**stats,**archived_material_score(m,truth))
                case_result['arms'][name] = info
            result['cases'][case] = case_result
            print(case,'old-only and far/near coherent (L2≤4) common-init NRMSE',
                  round(case_result['old_only']['nrmse'],4),
                  *(round(case_result['arms'][f'{R:.2f}_24_coherent_l2_4']['anchor_gn']['nrmse'],4) for R in RINGS), flush=True)
    args.output.with_suffix('.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(args.output.with_suffix('.npz'),**outputs)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('--archive',type=Path,default=Path(os.environ.get('ZALG_WEAK_ARCHIVE',DEFAULT_ARCHIVE)))
    p.add_argument('--feeds',type=Path,default=HERE/'results/angular_gap.npz')
    p.add_argument('--materials',type=Path,default=HERE/'results/weak_r4.npz')
    p.add_argument('--output',type=Path,default=HERE/'results/angular_gn')
    p.add_argument('--cases',nargs='+',default=list(PARENT_IDS))
    p.add_argument('--steps',type=int,default=10)
    run(p.parse_args())

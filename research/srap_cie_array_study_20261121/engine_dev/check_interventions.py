"""Development-only branch experiments using public training fields, never truth/heldouts."""
import argparse
import json
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import engine
from engine import legacy, np

CASE_NAMES = ('la180_austria_rot90_eps3p5','la120_austria_eps3p5',
              'la180_coupled_triplet_eps4p5','fresnel_twin_8GHz_half180_tx')
ARMS = ('dual','reset','full_ap','srap','graph')
SAMPLE = (0,1,4,8)
HISTORICAL = HERE.parents[2] / 'limited_aperture_routed_dr_20261120/.runs/public'
if not HISTORICAL.is_dir():
    HISTORICAL = HERE.parent.parent / 'limited_aperture_routed_dr_20261120/.runs/public'


def read_public(name, directory):
    with np.load(directory / f'{name}.npz', allow_pickle=False) as f:
        p = {k: f[k] for k in ('E','S','Y','centers','radius','k','delta','case_json','rx_ids','valid_mask') if k in f}
    case = json.loads(str(p['case_json']))
    if name in CASE_NAMES[:2] or name == 'la180_austria_rot0_eps3p5':
        case['public_schedule'] = 'legacy_low'
    return p,case


def model_for(p, case):
    return legacy.PaperModel(p,case,legacy.Clock(600.))


def dense_comparison(directory):
    name = 'la180_austria_rot0_eps3p5'
    p, case = read_public(name,directory)
    model = model_for(p,case)
    state = {}
    def collect(event,step,stage,a):
        if event == 'before' and not state:
            state.update({k:v.copy() for k,v in a.items()})
    model.checkpoint_hook = collect
    for _,step,_,_ in engine.run_cie_arm(model,p['Y'],case,'dual',np.zeros(model.n),p['delta']):
        if state:
            break
    weak,strong = engine.route_basis(model)
    dense_model = model_for(p,case)
    q0 = model.clock.wall()
    sparse = engine.BallProjector(model,p['Y'],p['delta'],np.column_stack((strong,weak)))
    sparse_setup = model.clock.wall()-q0
    t = dense_model.clock.wall()
    dense = legacy.make_projector(dense_model,p['Y'],p['delta'])
    dense_setup = dense_model.clock.wall()-t
    gx,gq = state['J']/6.,state['Q']
    cx,cq,_ = dense.material(gx,gq)
    projections = []
    for label,x,q in [('first_cycle',cx,cq),('air',np.zeros_like(cx),model.E)]:
        t0 = time.perf_counter()
        xd,qd = dense.affine_ball(x,q)
        td = time.perf_counter()-t0
        t0 = time.perf_counter()
        xs,qs = sparse.affine_ball(x,q)
        ts = time.perf_counter()-t0
        t_s = 6*model.ports.S
        # Stationarity on the propagation graph, with the shared-port active multiplier.
        mu = sparse.last_A_info['multiplier']
        probe = sparse.H(xs) - x - sparse.BH(q-model.E) + mu*t_s.conj().T @ (t_s @ xs - p['Y'])
        normalizer = np.linalg.norm(x+sparse.BH(q-model.E)) + np.linalg.norm(mu*t_s.conj().T@p['Y'])
        rxs = weak@(weak.conj().T @ (xs-gx)); rxd = weak@(weak.conj().T @ (xd-gx))
        rqs = weak@(weak.conj().T @ (qs-gq)); rqd = weak@(weak.conj().T @ (qd-gq))
        sj_s,_,m_s = sparse.material(gx+rxs,gq+rqs)
        sj_d,_,m_d = dense.material(gx+rxd,gq+rqd)
        jd,gr = engine.dpcs(model,p['Y'],15,False)
        dd = legacy.FourierCurrents(model.shape,5,gr,model.clock.counts)
        ss = engine.FourierCurrents(model.shape,5,gr,model.clock.counts)
        ad = dd.project_coefficients(6*sj_d-jd)
        aa = ss.project_coefficients(6*sj_s-jd)
        projections.append(dict(input=label,cycle=step,
            x_relative=engine.ratio(xs-xd,xd),q_relative=engine.ratio(qs-qd,qd),
            projected_material_relative=engine.ratio(m_s-m_d,m_d),
            fourier_coefficient_relative=engine.ratio(aa-ad,ad),
            reinjected_J_relative=engine.ratio(ss.apply(aa)-dd.apply(ad),jd+dd.apply(ad)),
            data_ball_relative=float(np.linalg.norm(t_s @ xs-p['Y'])/p['delta']),
            normal_equation_relative=float(np.linalg.norm(probe)/max(normalizer,1e-300)),
            multiplier_relative=float(abs(mu-dense.last_A_info['multiplier'])/max(1.,abs(mu))),
            dense_s=td,cg_s=ts))
    out = dict(case=name, dense_setup_s=dense_setup,fft_cg_setup_s=sparse_setup,
               cg_max_true_H_residual=float(model.clock.counts['H_max_true_relative_residual']),
               cg_rhs_iterations=model.clock.counts['H_CG_rhs_iterations'],
               scipy_dense_cholesky=dense_model.clock.counts['H_cholesky'],projections=projections)
    (HERE/'dense_comparison.json').write_text(json.dumps(out,indent=2))
    return out


def interventions(name,directory):
    p, case = read_public(name,directory)
    prefix = None
    summary = []
    folder = HERE / 'checkpoints'
    folder.mkdir(exist_ok=True)
    for method in ARMS:
        model = model_for(p,case)
        capture = {}
        first = None
        start = time.perf_counter()
        def hook(event,count,stage,arrays):
            nonlocal first
            if event == 'before' and first is None:
                first=count
            if count == first:
                for k,v in arrays.items():
                    capture[f'{event}_{k}'] = v.copy()
        model.checkpoint_hook = hook
        history = []
        for m,step,stage,force in engine.run_cie_arm(model,p['Y'],case,method,np.zeros(model.n),p['delta']):
            if first is None:
                continue
            if step-first in SAMPLE:
                phi,res,pred = model.objective(m,p['Y'])
                capture[f't{step-first}_m'] = m.copy()
                capture[f't{step-first}_pred'] = pred.copy()
                history.append(dict(steps_after=step-first,stage=stage,phi=phi,
                                    physical_train_residual=res))
            if step-first == max(SAMPLE):
                break
        before = {k:v for k,v in capture.items() if k.startswith('before_')}
        same = True if prefix is None else all(np.array_equal(v,prefix[k]) for k,v in before.items())
        if prefix is None:
            prefix = before
        if not same:
            raise ArithmeticError(f'Branch checkpoint differs for {name}/{method}')
        if method == 'dual':
            np.savez_compressed(folder/f'{name}__shared.npz',**before)
        np.savez_compressed(folder/f'{name}__{method}.npz',**{k:v for k,v in capture.items()
                                                               if not k.startswith('before_')})
        summary.append(dict(case=name,method=method,shared_checkpoint_bitwise=same,step=first,
                            experiment_wall_s=time.perf_counter()-start,counts=dict(model.clock.counts),
                            trajectory=history,trace=model.trace_route))
        print(name,method,first,history[-1],flush=True)
    (HERE/f'{name}.json').write_text(json.dumps(summary,indent=2))
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--public',type=Path,default=HISTORICAL)
    parser.add_argument('--cases',nargs='*',default=CASE_NAMES)
    parser.add_argument('--skip-dense',action='store_true')
    args=parser.parse_args()
    if not args.skip_dense:
        print(json.dumps(dense_comparison(args.public),indent=2),flush=True)
    for name in args.cases:
        interventions(name,args.public)

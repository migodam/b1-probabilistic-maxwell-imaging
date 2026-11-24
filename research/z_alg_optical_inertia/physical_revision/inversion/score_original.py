"""Score sealed full-wave maps on material and offset synthetic ports."""

import argparse
import json

import numpy as np
from scipy.special import hankel1, jv
from threadpoolctl import threadpool_limits

from fit import native
from gate import HERE, STUDY, clocks, elapsed, original, write


def probes(held):
    width,side,k,receiver_radius = 2.,96,2*np.pi*400e6/3e8,3.75
    x = (np.arange(side)+.5)*(width/side)-width/2.
    yy,xx = np.meshgrid(x,x,indexing='ij')
    points = np.c_[xx.ravel(),yy.ravel()]

    def angular_ports(n):
        theta = np.deg2rad((np.arange(n)+.5*held)*360/n)
        return np.c_[np.cos(theta),np.sin(theta)]

    tx = angular_ports(10)
    rx = receiver_radius*angular_ports(20)
    radius = width/side/np.sqrt(np.pi)
    d = np.linalg.norm(rx[:,None]-points[None],axis=2)
    s = .5j*np.pi*k*radius*jv(1,k*radius)*hankel1(0,k*d)
    return s,np.exp(1j*k*(points@tx.T))


def predict(geometry, shape):
    model = native.Model(*geometry,np.ones((20,10)))
    _, ports, currents = model.full(shape, gradient=False)
    rhs = shape[:,None]*model.b
    defect = currents-shape[:,None]*model.green(currents)-rhs
    return ports, dict(state_equation_relative_defect=float(np.linalg.norm(defect)/np.linalg.norm(rhs)),
                       counts=dict(model.count))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=type(HERE),default=HERE/'original_fits')
    args = parser.parse_args()
    native.RTOL = 1.e-10
    start = clocks()
    geometry, observations = original()
    s,b = probes(False)
    check = dict(train_S_max_difference=float(np.max(np.abs(s-geometry[1]))),
                 train_B_max_difference=float(np.max(np.abs(b-geometry[2]))))
    if max(check.values())>1.e-12:
        raise ArithmeticError('Reconstructed original training geometry does not match.')

    held = (geometry[0],*probes(True))
    cache = {}
    scored = {}
    for obs in observations:
        path = args.input/(obs.id+'.npz')
        if not path.exists():
            continue
        t = clocks()
        with np.load(path,allow_pickle=False) as predicted:
            if not np.array_equal(predicted['Y'],obs.y):
                raise ArithmeticError('Prediction does not match public observation.')
            reference_name = obs.object_id+'.npz'
            if reference_name not in cache:
                f = STUDY/'synthetic/evaluation'/reference_name
                with np.load(f,allow_pickle=False) as reference:
                    truth = reference['m_fine']
                ref_train, info_train = predict(geometry,truth)
                ref_held, info_held = predict(held,truth)
                cache[reference_name] = truth,ref_train,ref_held,info_train,info_held
            truth,ref_train,ref_held,info_train,info_held = cache[reference_name]
            receipt = json.loads(path.with_suffix('.json').read_text())
            rows = {}
            inside = truth>0
            reused = {}
            for policy,value in receipt['reports'].items():
                if policy != value['shared_path']:
                    rows[policy] = rows[value['shared_path']]
                    continue
                material = predicted[policy+'_material']
                pred_train = predicted[policy+'_ports']
                pred_held, info = predict(held,material)
                delta = material-truth
                rows[policy] = dict(
                    peak=float(np.max(material)), full_relative_l2=float(np.linalg.norm(delta)/np.linalg.norm(truth)),
                    inside_relative_l2=float(np.linalg.norm(delta[inside])/np.linalg.norm(truth[inside])),
                    background_rms=float(np.sqrt(np.mean(delta[~inside]**2))),
                    train_clean_port_relative_l2=float(np.linalg.norm(pred_train-ref_train)/np.linalg.norm(ref_train)),
                    born_at_truth_relative_l2=float(np.linalg.norm(geometry[1]@(truth[:,None]*geometry[2])-ref_train)/
                                                     np.linalg.norm(ref_train)),
                    held_clean_relative_l2=float(np.linalg.norm(pred_held-ref_held)/np.linalg.norm(ref_held)),
                    held_audit=info)
                reused[policy] = (pred_train,pred_held)
            if receipt['discrepancy_stage1_minimum_ratio']<=2.5:
                raise ArithmeticError('Original exact discrepancy requires an independent replay.')
            rows['exact_discrepancy'] = rows['discrepancy']
            scored[obs.id] = dict(
                target_peak=float(truth.max()),target_inside_cells=int(np.count_nonzero(inside)), policies=rows,
                original_exact_r_control='Same observed first-stage threshold crossing, same cap schedule and identical saved path',
                train_truth_audit=info_train,held_truth_audit=info_held,**elapsed(t))
            np.savez_compressed(args.input/(obs.id+'_offline_ports.npz'), truth_train=ref_train,
                             truth_held=ref_held,
                             **{name+'_train':ports[0] for name,ports in reused.items()},
                             **{name+'_held':ports[1] for name,ports in reused.items()})
        print(obs.id, 'scored',flush=True)
    write(HERE/'original_offline_scores.json',dict(scope='96² tight numerical full solver with offset artificial hold-outs',
                                                   training_geometry_match=check,results=scored,**elapsed(start)))


if __name__ == '__main__':
    with threadpool_limits(1):
        main()

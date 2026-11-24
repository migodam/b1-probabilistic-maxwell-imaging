"""Numerical all-pixel frame dual probe using low positive-state witnesses."""
from time import perf_counter,process_time

import numpy as np
from scipy.linalg import eigh, qr
from scipy.optimize import linprog
from threadpoolctl import threadpool_limits

from outward import clipped_best
from replay import HERE, configure, herm, save_json, times


def main():
    start=perf_counter(),process_time()
    name="u100_best"
    geo,native,_,_=configure(name)
    solved=np.load(HERE/f"{name}_conditioned.npz",allow_pickle=False)
    eig=solved["negative_eig"]
    neg=solved["negative_vec"]
    n=neg/np.sqrt(-eig)
    gn=solved["negative_gn"]
    w=solved["frame_weight_upper_numeric"]
    nu=len(eig)
    native.modes(nu+24)
    basis=native.eig_basis[:,nu:]
    basis-=neg@(neg.conj().T@basis)
    basis,_=qr(basis,mode="economic")
    pos,rot=eigh(herm(basis.conj().T@native.numeric_action(basis)))
    if pos.min()<=0:
        raise ArithmeticError("Positive-frame proposal has an unresolved positive mode.")
    basis=basis@rot
    trial=basis/np.sqrt(pos)
    trial_g=geo.numerical_green(trial)
    # For any test vector the independently stronger endpoint can be
    # assigned at each pixel. Each resulting weighted density is <= 1.
    d=np.maximum(abs(trial)**2,abs(trial-native.w.cap*trial_g)**2).T
    sums=[]
    for mode,rot in enumerate([np.eye(nu),solved["homogeneous_subspace"]]):
        for k in (1,2,4,6,9,nu):
            a,b=n@rot[:,:k],gn@rot[:,:k]
            x=clipped_best(a,b,native.w.cap)
            obj=np.sum(abs(a-x[:,None]*b)**2,axis=1)
            t=perf_counter(), process_time()
            proposal=linprog(-obj,A_ub=d,b_ub=np.ones(len(d)),
                             bounds=np.column_stack((np.zeros(len(w)),w)),method="highs")
            sums.append(dict(dimension=k,basis=mode,status=proposal.message,
                            trace_upper_numeric=-float(proposal.fun) if proposal.success else None,**times(*t)))
            print(sums[-1],flush=True)
    outfile=HERE/"u100_best_all_pixel_probe"
    save_json(outfile.with_suffix(".json"),dict(rows=sums,positive_witness_count=len(pos),
        positive_ritz=pos.tolist(),operations=dict(geo.cost),**times(*start)))
    np.savez_compressed(outfile.with_suffix(".npz"), positive_trial=trial,
                        positive_trial_green=trial_g, frame_densities=d)


if __name__=="__main__":
    with threadpool_limits(1):
        main()

"""Search for a common two-frequency opponent using weighted worst coherent waveforms."""
from pathlib import Path
from time import perf_counter, process_time
import argparse
import json

import numpy as np
from scipy.optimize import minimize
from threadpoolctl import threadpool_limits
from probe import Circuit, HERE, R, at_load, material_task, summarize_pair


def run(side, t, sectors, steps, start=None, order=14):
    clock = perf_counter(), process_time()
    ix, shells = material_task(side)
    circuits = [Circuit(side, f, ix) for f in (400, 460)]
    xy = circuits[0].disk.centers[ix] - [.125, .0625]
    angle = np.arctan2(xy[:, 1], xy[:, 0]) % (2*np.pi)
    sector = np.minimum((sectors*angle/(2*np.pi)).astype(int), sectors-1)
    group = np.where(shells, 1 + (shells-1)*sectors + sector, 0)
    exterior = [group==g for g in np.unique(group[group>0])]
    num = len(exterior)
    source = start or f"common_{side}_{t}_shared.npz"
    with np.load(HERE/source) as p:
        r=float(p["r"])
        v0=np.array([p[key][mask].mean() for key in ("plus","minus") for mask in exterior])
    loads = np.exp(2j*np.pi*np.arange(4)/4)
    prefix=f"joint_wave_{side}_{t}_s{sectors}"
    trace = []
    def fields(v):
        plus=np.zeros(len(ix))
        plus[group==0]=t
        minus=np.zeros(len(ix))
        for j,mask in enumerate(exterior):
            plus[mask]=v[j]
            minus[mask]=v[num+j]
        return plus,minus
    def objective(v):
        value=0.
        grad=np.zeros_like(v)
        peak=[]
        pair=fields(v)
        for circuit in circuits:
            states=[]
            for x in pair:
                s, left, right=circuit.solve(x, jac=True)
                derivatives=np.array([-2j*R*(left[:,m]@right[m]) for m in exterior])
                states.append((s, derivatives))
            largest=0.
            for load in loads:
                delta=(at_load(states[0][0],load)-at_load(states[1][0],load))/(2*r)
                u,sv,vh=np.linalg.svd(delta,full_matrices=False)
                largest=max(largest,float(sv[0]))
                value+=sum(sv**order)
                mode=(u*(sv**(order-1)))@vh
                for j,(s,deriv) in enumerate(states):
                    k=load/(1-load*s[24,24])
                    l=np.eye(16,25,dtype=complex)
                    b=l.T.copy()
                    l[:,24]=k*s[:16,24]
                    b[24,:]=k*s[24,:16]
                    dd=(l@deriv@b)/(2*r)
                    grad[j*num:(j+1)*num]+=(1-2*j)*order*np.einsum("ij,kij->k",mode.conj(),dd).real
            peak.append(largest)
        pnorm=(value/(len(loads)*2))**(1/order)
        trace.append([float(pnorm),*peak])
        return pnorm,grad*(pnorm/(order*value))
    out=minimize(objective,v0,method="L-BFGS-B",jac=True,bounds=[(0.,4.5)]*len(v0),
                 options=dict(maxiter=steps, ftol=1e-13, gtol=1e-8, maxcor=18))
    plus,minus=fields(out.x)
    data=dict(ix=ix, group=group, v=out.x, plus=plus, minus=minus, r=r)
    np.savez_compressed(HERE/(prefix+".npz"),**data)
    bands={}
    for circuit in circuits:
        bands[str(circuit.band)]=summarize_pair(circuit.solve(plus),circuit.solve(minus),r)
    receipt=dict(side=side,t=t,start=source,sectors=sectors,order=order,iterations=int(out.nit),
                 evaluations=int(out.nfev),stop=str(out.message),history=trace,diagnostics=bands,
                 calls=[w.count for w in circuits],wall_s=perf_counter()-clock[0],cpu_s=process_time()-clock[1])
    (HERE/(prefix+".json")).write_text(json.dumps(receipt,indent=2)+"\n")
    print(prefix,[b["all_loads_grid_sup_transfer_over_2r"] for b in bands.values()],flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("side",type=int)
    parser.add_argument("t",type=float)
    parser.add_argument("sectors",type=int)
    parser.add_argument("--steps",type=int,default=100)
    parser.add_argument("--start")
    parser.add_argument("--order",type=int,default=14)
    args=parser.parse_args()
    with threadpool_limits(1):
        run(**vars(args))

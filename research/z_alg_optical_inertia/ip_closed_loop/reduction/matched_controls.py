"""Verify actual Neumann divergence and the scalar positive-state obstruction."""
import numpy as np
from flint import arb,ctx

from benchmark import (HERE, ProofGrid, clock, elapsed, endpoints, freeze_config, load_model,
                       neumann_probe, prepare_candidates, threadpool_limits, write)
from interval_backend import diagonal, row_scaled, sub
from verify import dark_obstruction, spectral_moment


def run():
    freeze = freeze_config()
    matched = []
    for side,idx in ((32,17),(96,87)):
        started = clock()
        spec = freeze['scales'][str(side)]
        geom,ids,(g,s,b),y,r,base,sites,name = load_model(side,spec['candidate_base_upper'])
        w = np.full(len(ids),spec['weight'])
        candidates = prepare_candidates(side,ids,sites,base,w,1000,spec['seed'])
        lo,up = endpoints(candidates[idx],np.zeros(len(base)),base)
        _,(eg,es,_) = ProofGrid(geom).restricted(ids)
        md,rd = (lo+up)/2,(up-lo)/2
        md_exact = [(arb(float(l))+arb(float(u)))/2 for l,u in zip(lo,up,strict=True)]
        rd_exact = [(arb(float(u))-arb(float(l)))/2 for l,u in zip(lo,up,strict=True)]
        mg = md[:,None]*g
        ident = np.eye(len(g))
        a,rg = ident-mg,rd[:,None]*g
        i_acb = diagonal([1]*len(g))
        emg = row_scaled(eg,md_exact)
        exact_a,exact_rg = i_acb-emg,row_scaled(eg,rd_exact)
        neumann = neumann_probe(g,s,b,y,lo,up,r,eg)
        # The upper-endpoint perturbation is legal on the same original material box.
        spectral = spectral_moment(a,rg,exact_a,exact_rg)
        # Uniform sector weights have a dark direction that no residual weight changes.
        u,nz = lo*up,np.flatnonzero(lo*up)
        zz = g[nz]
        c = ident - mg - mg.conj().T + zz.conj().T@(u[nz,None]*zz)
        # Only these two nonzero products enter the expanded sector matrix.
        gr = sub(eg,nz,range(len(g)))
        products = [arb(float(lo[i]))*arb(float(up[i])) for i in nz]
        ec = i_acb-emg-emg.conjugate().transpose()+gr.conjugate().transpose()*row_scaled(gr,products)
        dark = dark_obstruction(c,s,ec,es)
        out = dict(side=side,observation=name,candidate_index=idx,neumann=neumann,
                   actual_upper_endpoints_legal_perturbation=spectral,uniform_sector_psd_obstruction=dark,
                   scopes={'neumann':'2-norm center contraction plus actual all-upper-endpoints spectral radius',
                           'positive_state':'every positive scalar times this uniform physical sector weight, same receiver matrix',
                           'exclusions':'matched candidate orientations chosen among the independent replay sample'},
                   cost=elapsed(started))
        matched.append(out)
        write(HERE/'matched_controls.json',matched)
        print(side, spectral, dark,flush=True)


if __name__ == '__main__':
    ctx.prec,ctx.threads = 128,1
    with threadpool_limits(limits=1):
        run()

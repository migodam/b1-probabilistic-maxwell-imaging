"""Task-restricted four-block numerical identity on actual TM worlds.

This is a solver qualification fixture, not a wall ambiguity witness or ROM.
It uses the existing two seeds, fixed frequencies and declared finite targets.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.linalg import lu_factor,lu_solve,svdvals
from .accounting import Book
from .vie import VIE
from .geometry import antennas,material_grid,probes

def main(root):
    root=Path(root);out=root/'results/four_blocks';out.mkdir(parents=True,exist_ok=True)
    book=Book(root/'results');tx,rx=antennas();rows=[];arrays={}
    for seed in [2026102201,2026102202]:
        for f in [.8,1.2]:
            p=VIE(64,f,tx,rx,book,source_cell_average=True)
            w=material_grid(p,seed);base=p.solve(w)
            eps3=material_grid(p,seed,probe=probes()[8]);ports=np.flatnonzero(eps3-w)
            active=base.active
            tic=time.perf_counter()
            rhs=(w[active]-1)[:,None]*p.green_block(active,ports)
            port_response,_=p._solve(base.factor,rhs,label='fourblock_port_rhs')
            A=p.S[:,ports]+p.S[:,active]@port_response
            B=base.field[ports]
            C=p.green_block(ports,ports)+p.green_block(ports,active)@port_response
            setup=time.perf_counter()-tic
            for probe_index in [8,9]:
                target=material_grid(p,seed,probe=probes()[probe_index]);v=(target-w)[ports]
                core=np.eye(len(ports))-C*v[None,:]
                book.reserve('validation','fourblock_target_core_LU',factorizations=1,freq=f,active_cells=len(ports))
                tic=time.perf_counter();factor=lu_factor(core,check_finite=False)
                book.reserve('validation','fourblock_target_core_rhs',rhs=4,freq=f)
                E=lu_solve(factor,B,check_finite=False);core_time=time.perf_counter()-tic
                reconstructed=base.data+(A@(v[:,None]*E)).T
                born=base.data+(A@(v[:,None]*B)).T
                direct=p.solve(target)
                induced=direct.data-base.data
                key=f's{seed}_f{f}_p{probe_index}'
                row=dict(seed=seed,frequency=f,probe_index=probe_index,port_rank=len(ports),
                         full_field_equivalence_relative=float(np.linalg.norm(reconstructed-direct.data)/np.linalg.norm(direct.data)),
                         equivalence_whitened=float(np.linalg.norm(reconstructed-direct.data)/.01),
                         born_induced_relative_error=float(np.linalg.norm(born-direct.data)/np.linalg.norm(induced)),
                         core_sigma_min=float(svdvals(core)[-1]),core_condition=float(np.linalg.cond(core)),
                         port_setup_seconds=setup,core_seconds=core_time,direct_state_seconds=direct.timings['state_seconds'],
                         port_rhs=len(ports),direct_rhs=4,small_core_rhs=4,
                         status='DISCRETE_IDENTITY_ONLY_UNQUALIFIED_CONTINUUM_GRID',
                         full_query_signature=False)
                rows.append(row);print(json.dumps(row),flush=True)
                for label,value in [('full',direct.data),('cal',base.data),('born',born),('fourblock',reconstructed),('A',A),('B',B),('C',C)]:arrays[f'{key}_{label}']=value
    (out/'RESULTS.json').write_text(json.dumps(rows,indent=2))
    np.savez_compressed(out/'ARRAYS.npz',**arrays)

if __name__=='__main__':main(Path(__file__).resolve().parents[2])

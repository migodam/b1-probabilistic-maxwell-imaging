"""Actual FEM wall/target Schur identity, full rescattering and Born control.

Uses the already-audited coarse base geometry solely for this discrete identity.
The L2-projected target-field coordinates make the analytic incident forcing
exact in the discrete equation. No nodal-incident interpolation is substituted.
"""
from pathlib import Path
import json,gc,time
import numpy as np
from scipy.linalg import cho_factor,cho_solve,lu_factor,lu_solve
from skfem import asm
from wslam.accounting import BudgetBook,json_default
from wslam.forward import Evaluator
from wslam.world import World,polygon_inside
from wslam.fem_backend import _contrast_mass,_forcing
from wslam.physics import paired_feedback_residual,target_feedback_terms
root=Path(__file__).resolve().parents[1];out=root/'results/e0/four_blocks';out.mkdir(parents=True,exist_ok=False)
book=BudgetBook();solver=Evaluator(book,'E0E1',h=.22,n=128,grid_refinement=1)
world=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
wall=solver.evaluate(world,include_target=False);p,s0=wall.backend,wall.state
mask=polygon_inside(s0.material_points,world.pose,world.radius,solver.n)
mass=asm(_contrast_mass,p.material_basis,contrast=(p.k**2*mask).reshape(p._material_shape))[p._free_dofs][:,p._free_dofs].tocsc()
idx=np.unique(mass.nonzero()[0]);m=len(idx)
M=mass[idx][:,idx].toarray().real
if not m:raise RuntimeError('Empty target port')
chol=cho_factor(M,lower=True,check_finite=True)
forcing=np.column_stack([asm(_forcing,p.material_basis,load=(p.k**2*mask*s0.field[:,s]).reshape(p._material_shape))[p._free_dofs][idx] for s in range(len(p.tx))])
B=cho_solve(chol,forcing)
book.reserve('E0E1','wall_green_target_ports',rhs=int(m),factorizations=0,freq=1.,backend='existing_fem_factor',role='four_block_setup')
G=np.empty((len(p._free_dofs),m),complex);maxres=0.;start=time.perf_counter()
for startcol in range(0,m,24):
    last=min(m,startcol+24);rhs=np.zeros((len(p._free_dofs),last-startcol),complex,order='F')
    rhs[idx[startcol:last],np.arange(last-startcol)]=1
    g=s0.factor.solve(rhs)
    residual=float(np.linalg.norm(s0._matrix@g-rhs,axis=0).max())
    if residual>p.residual_tolerance:raise RuntimeError('Green port residual failed')
    G[:,startcol:last]=g;maxres=max(maxres,residual)
    p._event('backend_solve',label='wall_green_target_ports',rhs=last-startcol,residual=residual,transpose='N')
C=G[idx];A=np.asarray(p._receiver_probe[:,p._free_dofs]@G)
np.savez_compressed(out/'BLOCKS.npz',A=A,B=B,C=C,mass=M,wall_data=wall.data,free_indices=idx)
rows=[]
for contrast in (.5,3.,5.):
    target=World(world.support,world.epsinf,world.ell,world.pose,contrast=contrast)
    full=solver.evaluate(target)
    chi=contrast+1j*target.target_loss
    V=chi*M;core=np.eye(m)-C@V
    corefactor=lu_factor(core);q=lu_solve(corefactor,B)
    delta=A@V@q
    predicted=wall.data+delta.T
    born=wall.data+(A@V@B).T
    exact_delta=full.data-wall.data
    feedback=paired_feedback_residual(p,s0,full.state)
    row=dict(target_contrast=contrast,target_ports=m,setup_rhs=m,direct_sources=len(p.tx),
        relative_full_identity=float(np.linalg.norm(predicted-full.data)/np.linalg.norm(full.data)),
        target_induced_relative_identity=float(np.linalg.norm(delta.T-exact_delta)/np.linalg.norm(exact_delta)),
        Born_target_induced_relative_error=float(np.linalg.norm(born-full.data)/np.linalg.norm(exact_delta)),
        feedback_equation=feedback,core_solve_relative_residual=float(np.linalg.norm(core@q-B)/np.linalg.norm(B)),
        minimum_mass_eigenvalue=float(np.linalg.eigvalsh(M)[0]),
        setup_wall_seconds=time.perf_counter()-start,
        source_field_coordinate='B=M_target^{-1} integral_target k2 phi E_wall; exact L2 projected forcing, not interpolated nodal incident',
        scope='same discrete FEM exact four-block identity; no continuum certificate or speedup')
    rows.append(row);(out/'RESULTS.json').write_text(json.dumps(rows,indent=2,default=json_default))
    np.savez_compressed(out/f'contrast_{contrast}.npz',full=full.data,predicted=predicted,born=born,wall=wall.data,core=core,solution=q)
    print(json.dumps(row,default=json_default),flush=True)
del G;gc.collect()

#!/usr/bin/env python3
"""Paid E0 metrics only. Run through tools/run_job.py; never self-launches jobs."""
from pathlib import Path
from dataclasses import replace
import gc
import json
import time
import numpy as np
from scipy.special import hankel1, jv, jvp, h1vp
from wslam.accounting import BudgetBook, json_default
from wslam.geometry import enumerate_supports
from wslam.world import World
from wslam.forward import Evaluator, volume_extraction
from wslam.fem_backend import PMLConfig
from wslam.physics import pack

ROOT = Path(__file__).resolve().parents[1]


def line_mie(tx, rx, k, center, radius, ep):
    """Graf addition, actual unit exterior point sources; scattering only.

    Independently evaluated Bessel series copied from the historical mathematical
    expression, without importing any historical job/accounting/solver module.
    """
    m = np.arange(-32, 33)
    n = np.sqrt(complex(ep))
    z = k*radius
    a = (n*jvp(m,n*z)*jv(m,z)-jv(m,n*z)*jvp(m,z))/(jv(m,n*z)*h1vp(m,z)-n*jvp(m,n*z)*hankel1(m,z))
    rs, rr = tx-np.asarray(center), rx-np.asarray(center)
    incoming = .25j*hankel1(m[:,None],k*np.linalg.norm(rs,axis=1)[None,:])*np.exp(-1j*m[:,None]*np.arctan2(rs[:,1],rs[:,0])[None,:])
    outgoing = hankel1(m[:,None],k*np.linalg.norm(rr,axis=1)[None,:])*np.exp(1j*m[:,None]*np.arctan2(rr[:,1],rr[:,0])[None,:])
    return (outgoing.T@(a[:,None]*incoming)).T


def white(a):
    return float(np.linalg.norm(pack(a)))


def main():
    if not (ROOT/'EXPERIMENT_CONTRACT.json').exists():
        raise RuntimeError('Freeze contract first')
    out = ROOT/'results/e0/qualification'
    out.mkdir(parents=True,exist_ok=False)
    book = BudgetBook()
    rows = []
    started = time.perf_counter()
    def record(row):
        rows.append(row)
        (out/'RESULTS.json').write_text(json.dumps(dict(
            status='METRICS_ONLY_PARENT_SCIENTIFIC_JUDGMENT_REQUIRED',
            rows=rows,counts=book.totals('E0E1'),wall_seconds=time.perf_counter()-started,
            scope='fixed E0 configurations; no automatic continuum or campaign PASS'),
            indent=2,default=json_default,allow_nan=False))
        print(json.dumps(row,default=json_default,allow_nan=False),flush=True)
    def observation(label, ev, fields=False):
        arrays=dict(data=ev.data,volume_data=volume_extraction(ev),
                    tx=ev.backend.tx,rx=ev.backend.rx)
        if fields:
            arrays.update(points=ev.state.material_points,weights=ev.state.material_weights,
                          epsilon=ev.state.epsilon,field=ev.state.field,
                          current=ev.state.current,scattered_dofs=ev.state.scattered_dofs,
                          mesh_p=ev.backend.mesh.p,mesh_t=ev.backend.mesh.t)
        np.savez_compressed(out/(label+'.npz'),**arrays)
    cases=[]
    for i,(mask,prior) in enumerate(enumerate_supports(2)):
        cases.append((f'support_{i:02d}',World(mask,3.75,.125,(.2,2.)),1.))
    for name,e,loss in [('weak',1.5,0.),('strong',6.,.25)]:
        for f in (.85,1.15):
            for j,pose in enumerate(((-1.,1.6),(1.,2.6))):
                cases.append((f'{name}_f{f}_corner{j}',World(np.ones((2,2),bool),e,loss,pose),f))
    ring=np.ones((3,3),bool);ring[1,1]=False
    hetero=np.where(np.indices((3,3)).sum(axis=0)%2,2.5,4.5)
    cases.append(('ring_3x3',World(ring,hetero,.1,(.2,2.)),1.))
    for index,(label,world,freq) in enumerate(cases):
        coarse_solver=Evaluator(book,'E0E1',h=.22,n=128,grid_refinement=1)
        fine_solver=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2)
        coarse=coarse_solver.evaluate(world,freq)
        fine=fine_solver.evaluate(world,freq)
        observation(label+'_primary',coarse,fields=(label=='support_00'))
        observation(label+'_independent',fine,fields=(label=='support_00'))
        point_delta=white(coarse.data-fine.data)
        record(dict(case=label,frequency=freq,world=world.metadata(),
            primary=dict(h=.22,n=128,grid_refinement=1,ndofs=coarse.backend.ndofs,
                         residual=coarse.state.residual,point_volume_white=white(volume_extraction(coarse)-coarse.data)),
            independent=dict(h=.17,n=256,grid_refinement=2,ndofs=fine.backend.ndofs,
                             residual=fine.state.residual,point_volume_white=white(volume_extraction(fine)-fine.data)),
            point_grid_white=point_delta,point_grid_exceeds_point1=point_delta>.1))
        coarse_solver.clear();fine_solver.clear();del coarse,fine;gc.collect()
    # Independent empty-wall cylinder: exact circle Mie versus polygon FEM,
    # retaining the geometrical polygon approximation as an explicit error.
    empty=World(np.zeros((2,2),bool),3.75,.125,(.2,2.))
    fine_solver=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2)
    ev=fine_solver.evaluate(empty,1.)
    reference=line_mie(ev.backend.tx,ev.backend.rx,ev.backend.k,empty.pose,empty.radius,
                       1+empty.contrast+1j*empty.target_loss)
    observation('empty_wall_line_mie',ev,fields=True)
    np.savez_compressed(out/'MIE_REFERENCE.npz',scattered=reference,direct=ev.state.direct)
    record(dict(case='independent_unit_line_Mie',point_scattering_white=white(ev.state.scattered-reference),
                volume_scattering_white=white(volume_extraction(ev)-ev.state.direct-reference),
                circle_polygon_n=256,series_orders=[-32,32],residual=ev.state.residual,
                scope='analytic circular cylinder compared to fitted polygon FEM'))
    fine_solver.clear();del ev;gc.collect()
    # Frozen PML comparison, both alternative strength and thickness, no tuning.
    world=World(np.ones((2,2),bool),6.,.25,(-1.,1.6))
    primary=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2)
    ev=primary.evaluate(world,1.15);base_data=ev.data.copy();observation('pml_baseline',ev)
    primary.clear();del ev;gc.collect()
    for label,pm in [('strength6',PMLConfig(strength=6)),('thickness1p25',PMLConfig(thickness=1.25))]:
        solver=Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2,pml=pm)
        ev=solver.evaluate(world,1.15);observation('pml_'+label,ev)
        delta=white(ev.data-base_data)
        record(dict(case='pml_'+label,frequency=1.15,point_difference_white=delta,
                    exceeds_point1=delta>.1,residual=ev.state.residual,pml=pm.as_dict()))
        solver.clear();del ev;gc.collect()
    # Full physical seven-column Jacobian, with distinct left/right wall eps.
    base=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
    solver=Evaluator(book,'E0E1',h=.22,n=128,grid_refinement=1)
    ev=solver.evaluate(base,1.,jacobian=True,kind='left_right')
    np.savez_compressed(out/'JACOBIAN.npz',data=ev.data,jacobian=ev.jacobian,names=np.asarray(ev.names))
    observation('derivative_base',ev,fields=True)
    def shifted(name,step):
        if name in ('eps_left','eps_right'):
            e=np.array(base.epsinf,copy=True);e[:,0 if name=='eps_left' else 1]+=step
            return replace(base,epsinf=e)
        if name in ('pose_x','pose_y'):
            pose=base.pose.copy();pose[0 if name=='pose_x' else 1]+=step
            return replace(base,pose=pose)
        return replace(base,**{dict(target_contrast='contrast',target_loss='target_loss',target_radius='radius')[name]:getattr(base,dict(target_contrast='contrast',target_loss='target_loss',target_radius='radius')[name])+step})
    for col,name in enumerate(ev.names):
        steps=([.002,.001,.0005] if name.startswith('pose_') else [.0005] if name=='target_radius'
               else [.002,.001] if name.startswith('eps_') else [.002])
        for step in steps:
            fdsolver=Evaluator(book,'E0E1',h=.22,n=128,grid_refinement=1)
            plus=fdsolver.evaluate(shifted(name,step));minus=fdsolver.evaluate(shifted(name,-step))
            fd=(plus.data-minus.data).ravel()/(2*step);j=ev.jacobian[:,col]
            np.savez_compressed(out/f'FD_{name}_{step}.npz',plus=plus.data,minus=minus.data,
                                derivative=fd,jacobian=j)
            record(dict(case='central_FD',parameter=name,step=step,
                relative_error=float(np.linalg.norm(fd-j)/max(np.linalg.norm(fd),np.finfo(float).tiny)),
                derivative_difference_white=white(fd-j),
                maximum_residual=max(plus.state.residual,minus.state.residual)))
            fdsolver.clear();del plus,minus;gc.collect()
    p,s=ev.backend,ev.state
    rng=np.random.default_rng(2026102501)
    coeff=base.material_directions(s.material_points,p.freq,128,'left_right')
    weights=np.array([.13,-.21,.17,.003])
    direction=coeff@weights
    cot=rng.normal(size=s.data.shape)+1j*rng.normal(size=s.data.shape)
    product=p.jvp(s,lambda points:direction)
    adjoint=p.vjp(s,cot)
    lhs=float(np.real(np.vdot(cot,product)))
    rhs=float(np.real(np.vdot(adjoint,direction)))
    record(dict(case='real_VJP_dot',lhs=lhs,rhs=rhs,
                relative_defect=abs(lhs-rhs)/max(abs(lhs),abs(rhs),np.finfo(float).tiny),
                direction='left/right eps and target contrast/loss coefficient combination'))
    # Direct raw observations, plotted without interpretation or smoothing.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,figsize=(8,5),sharex=True)
    for source,data in enumerate(ev.data):
        axes[0].plot(p.rx[:,0],data.real,'o-',label=f'TX{source}')
        axes[1].plot(p.rx[:,0],data.imag,'o-')
    axes[0].set_ylabel('Re total E');axes[1].set_ylabel('Im total E')
    axes[1].set_xlabel('Receiver x');axes[0].legend();fig.tight_layout()
    fig.savefig(out/'RAW_POINT_OBSERVATIONS.png',dpi=160);plt.close(fig)
    solver.clear()
    record(dict(case='driver_complete',overall_gate='NOT_ASSIGNED_PARENT_REVIEW',
                source_plot='RAW_POINT_OBSERVATIONS.png',
                scope='metrics and fixed configurations only; no alternate configuration chosen'))


if __name__=='__main__':
    main()

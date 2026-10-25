#!/usr/bin/env python3
"""Fixed E1 local mechanism experiment; run only after parent E0 release.

No automatic gate, no A0 observation selection, no algorithmic damping as prior.
"""
from pathlib import Path
from dataclasses import replace
import gc,json,time
import numpy as np
from scipy.linalg import block_diag,eigh
from scipy.optimize import least_squares
from wslam.accounting import BudgetBook,json_default
from wslam.world import World
from wslam.forward import Evaluator
from wslam.physics import acquisition,pack,packj
from wslam.information import (target_structure,generalized_retention,schur,
    dynamic_information,svd_audit)
ROOT=Path(__file__).resolve().parents[1]


def white(a):return float(np.linalg.norm(pack(a)))


def main():
    gate=json.loads((ROOT/'results/e0/gate.json').read_text())
    if gate['decision']!='PASS':raise RuntimeError('E0 not qualified; E1 is not authorized')
    # Parent owns E0 release and launches this paid driver explicitly.
    out=ROOT/'results/e1/dynamic';out.mkdir(parents=True,exist_ok=False)
    book=BudgetBook();tx,rx=acquisition(outer=True);rows=[];started=time.perf_counter()
    base=World(np.ones((2,2),bool),np.array([[3.,4.],[3.,4.]]),.08,(.2,2.))
    def solver(**kw):return Evaluator(book,'E0E1',h=.17,n=256,grid_refinement=2,**kw)
    def record(row):
        rows.append(row)
        (out/'RESULTS.json').write_text(json.dumps(dict(status='METRICS_ONLY_PARENT_REVIEW',rows=rows,
            counts=book.totals('E0E1'),wall_seconds=time.perf_counter()-started,
            overall_gate='NOT_ASSIGNED',scope='known local base E1; no broad-prior ambiguity or A0 selection'),
            indent=2,default=json_default,allow_nan=False))
        print(json.dumps(row,default=json_default,allow_nan=False),flush=True)
    def save_audit(label,result,extra=None):
        arrays={}
        def flatten(value,prefix=''):
            if isinstance(value,dict):
                for k,v in value.items():flatten(v,prefix+k+'__')
            elif isinstance(value,np.ndarray):arrays[prefix[:-2]]=value
        flatten(result)
        if extra:arrays.update(extra)
        np.savez_compressed(out/(label+'.npz'),**arrays)
    raw={}
    for name,world,f in [('pose0_f1',base,1.),('pose0_f085',base,.85),
                         ('pose0_f115',base,1.15),('pose1_f1',replace(base,pose=np.array([.3,2.05])),1.)]:
        sol=solver(tx=tx,rx=rx);ev=sol.evaluate(world,f,jacobian=True,kind='left_right')
        raw[name]=dict(data=ev.data.copy(),J=ev.jacobian.reshape(len(tx),len(rx),7).copy())
        np.savez_compressed(out/(name+'.npz'),data=ev.data,J=ev.jacobian,names=np.asarray(ev.names),tx=tx,rx=rx)
        record(dict(case='physical_linearization',label=name,world=world.metadata(),frequency=f,
                    residual=ev.state.residual,ndofs=ev.backend.ndofs,jacobian_names=ev.names))
        sol.clear();del ev;gc.collect()
    def block(name,source,receivers=None):
        r=np.arange(len(rx)) if receivers is None else np.asarray(receivers)
        item=raw[name];j=item['J'][source,r];f=item['data'][source,r]
        return dict(A=packj(j[:,:2]),B=packj(j[:,2:4]),C=packj(j[:,4:7]),
                    gamma=packj(np.column_stack((f,1j*f))),F=f,J=j)
    def diagnostic(label,blocks,motion=False):
        aa=[b['A'] for b in blocks];bb=[b['B'] for b in blocks]
        a=np.vstack(aa);shared_b=np.vstack(bb);ind_b=block_diag(*bb)
        independent=target_structure(aa,bb,structure='independent')
        shared=target_structure(aa,bb,structure='shared')
        ri=generalized_retention(a,ind_b);rs=generalized_retention(a,shared_b)
        result=dict(independent=independent,shared=shared,retention_independent=ri,retention_shared=rs)
        extra=dict(A=a,B_shared=shared_b,B_independent=ind_b)
        if len(blocks)==2:
            d=np.column_stack((-np.eye(2),np.eye(2)));q=.03**2*np.eye(2)
            dyn=target_structure(aa,bb,structure='motion',d=d,q=q)
            result['motion']=dyn;extra.update(D=d,Q=q)
        save_audit(label,result,extra)
        record(dict(case='information_structure',label=label,
            scope='one complex RX port mechanism only' if 'minimal' in label else 'full16RX local information',
            K0=rs['K0'],Keff_shared=shared['information'],Keff_independent=independent['information'],
            rho_shared=rs['eigenvalues'],rho_independent=ri['eigenvalues'],
            angles_shared=rs['principal_angles'],angles_independent=ri['principal_angles'],
            shared_absolute_eigenvalues=np.linalg.eigvalsh(shared['information']),
            independent_absolute_eigenvalues=np.linalg.eigvalsh(independent['information']),
            motion_formula=result.get('motion',{}).get('formula'),
            motion_interpretation=('actual known mean displacement(.1,.05)' if motion else
              'hypothetical zero-mean nuisance-equality precision; no actual trajectory claim') if len(blocks)==2 else None,
            motion_information=result.get('motion',{}).get('information'),
            C_shared=shared['C'],R_norm_shared=float(np.linalg.norm(shared['R'])),
            DC=result.get('motion',{}).get('DC'),cross_covariance=None,
            covariance_status='No wall prior; covariance only in separately labeled proper Gaussian diagnostic'))
        return result,a,shared_b,ind_b
    left=block('pose0_f1',0)
    diagnostic('full_single_left',[left])
    for name,extra in [('base_right',block('pose0_f1',1)),('outer_left',block('pose0_f1',2)),
                       ('frequency115',block('pose0_f115',0))]:
        diagnostic('full_joint_'+name,[left,extra])
    diagnostic('full_actual_motion',[left,block('pose1_f1',0)],motion=True)
    ml,mr=block('pose0_f1',0,[8]),block('pose0_f1',1,[8])
    diagnostic('minimal_single_left',[ml])
    minimal,a,bshared,bind=diagnostic('minimal_joint_twoTX',[ml,mr])
    # Anchored translation nuisance affects antenna coordinates only. The
    # material and target world remain fixed for this paid +/- .002 finite FD.
    nav=[]
    for axis in range(2):
        delta=np.zeros(2);delta[axis]=.002
        values=[]
        for sign in (1.,-1.):
            sol=solver(tx=tx+sign*delta,rx=rx+sign*delta)
            ev=sol.evaluate(base,1.);values.append(ev.data.copy())
            np.savez_compressed(out/f'NAV_axis{axis}_sign{sign:+}.npz',data=ev.data,tx=ev.backend.tx,rx=ev.backend.rx)
            sol.clear();del ev;gc.collect()
        nav.append((values[0]-values[1])/.004)
    nav=np.stack(nav,axis=-1)
    np.savez_compressed(out/'NAV_JACOBIAN.npz',J=nav,step=.002,
                        names=np.asarray(['anchored_nav_translation_x','anchored_nav_translation_y']))
    full=[block('pose0_f1',s) for s in (0,1)]
    af=np.vstack([b['A'] for b in full]);pose=np.vstack([b['B'] for b in full]);c=np.vstack([b['C'] for b in full]);gamma=np.vstack([b['gamma'] for b in full])
    navreal=np.vstack([packj(nav[s]) for s in (0,1)])
    for label,nuisance,precision in [
        ('pose_only',pose,np.zeros((2,2))),
        ('pose_intrinsic_shared',np.column_stack((pose,c)),np.zeros((5,5))),
        ('pose_intrinsic_gamma_anchorednav',np.column_stack((pose,c,gamma,navreal)),
         np.diag([0.]*5+[1/.005**2]*4))]:
        result=schur(af,nuisance,precision);ret=generalized_retention(af,nuisance,prior_precision=precision)
        save_audit('confusion_'+label,dict(schur=result,retention=ret),dict(A=af,B=nuisance,statistical_nuisance_precision=precision))
        # Proper local Gaussian diagnostic uses an explicitly separate weak
        # wall/pose/intrinsic prior; no GN/LM damping is reinterpreted as prior.
        joined=np.column_stack((af,nuisance));weak_sd=100.
        prior=np.diag(np.r_[np.full(2,weak_sd**-2),np.full(nuisance.shape[1],weak_sd**-2)])
        prior[2:,2:]+=precision
        h=joined.T@joined+prior;cov=np.linalg.solve(h,np.eye(len(h)))
        np.savez_compressed(out/('LOCAL_GAUSSIAN_DIAGNOSTIC_'+label+'.npz'),precision=h,covariance=cov,
                            cross_covariance=cov[:2,2:],declared_prior_precision=prior)
        record(dict(case='confusion_model',label=label,K0=ret['K0'],Keff=result['information'],
            data_only_Keff=result['data_only_information'],rho=ret['eigenvalues'],
            angles=ret['principal_angles'],statistical_nuisance_precision=precision,
            gaussian_scope='LOCAL_GAUSSIAN_DIAGNOSTIC; weak sd100 physical coordinates, gamma/nav sd.005; no posterior certificate',
            gaussian_cross_covariance=cov[:2,2:],gaussian_wall_covariance=cov[:2,:2],
            algorithmic_damping_as_prior=False))
    # Finite compensation direction is fixed by local generalized information
    # gain, then max-absolute normalized; no endpoint result chooses direction.
    k0=a.T@a;gain=minimal['shared']['information']-minimal['independent']['information']
    values,vectors=eigh((gain+gain.T)/2,(k0+k0.T)/2)
    direction=vectors[:,-1];direction/=np.max(np.abs(direction))
    if direction[np.argmax(abs(direction))]<0:direction=-direction
    record(dict(case='profile_direction_frozen',wall_direction=direction,
                generalized_gain_eigenvalue=float(values[-1]),steps=[-.03,-.01,-.002,.002,.01,.03],
                source='minimal twoTX shared-minus-independent information',scope='known E1 local base only'))
    # Each profile uses fixed wall change, continuous bounded poses, and actual
    # noiseless baseline data. Independent poses are separate per transmitter.
    for mode in ('single_left','shared_twoTX','independent_twoTX'):
        source_ids=[0] if mode=='single_left' else [0,1]
        mini_tx=tx[source_ids];mini_rx=rx[[8]]
        target=np.stack([raw['pose0_f1']['data'][i,[8]] for i in source_ids])
        for step in (-.03,-.01,-.002,.002,.01,.03):
            eps=np.asarray(base.epsinf).copy();eps[:,0]+=step*direction[0];eps[:,1]+=step*direction[1]
            changed=replace(base,epsinf=eps)
            independent=mode=='independent_twoTX'
            cache={};solvers=[solver(tx=mini_tx[[i]],rx=mini_rx) for i in range(2)] if independent else [solver(tx=mini_tx,rx=mini_rx)]
            def evaluate_poses(x):
                key=np.asarray(x).tobytes()
                if key in cache:return cache[key]
                data=[];js=[]
                for i,sl in enumerate(solvers):
                    pose=x[2*i:2*i+2] if independent else x
                    ev=sl.evaluate(replace(changed,pose=np.asarray(pose)),1.,jacobian=True,kind='left_right')
                    data.append(ev.data.copy());js.append(ev.jacobian[:,2:4].copy())
                d=np.concatenate(data,axis=0) if independent else data[0]
                j=block_diag(*[packj(z) for z in js]) if independent else packj(js[0])
                # Independent residual packing must follow the same per-source
                # real/imag blocks as block_diag Jacobian.
                r=np.concatenate([pack(data[i]-target[[i]]) for i in range(len(data))]) if independent else pack(d-target)
                cache.clear();cache[key]=(r,j,d);return cache[key]
            initial=np.tile(base.pose,2 if independent else 1)
            lo=np.tile([-1.,1.6],2 if independent else 1);hi=np.tile([1.,2.6],2 if independent else 1)
            try:
                fit=least_squares(lambda x:evaluate_poses(x)[0],initial,
                    jac=lambda x:evaluate_poses(x)[1],bounds=(lo,hi),max_nfev=8)
                residual,jac,data=evaluate_poses(fit.x)
                ablocks=[ml['A']] if mode=='single_left' else [ml['A'],mr['A']]
                bblocks=[ml['B']] if mode=='single_left' else [ml['B'],mr['B']]
                aprofile=np.vstack(ablocks);bprofile=block_diag(*bblocks) if independent else np.vstack(bblocks)
                delta_pose=fit.x-initial
                predicted=aprofile@(step*direction)+bprofile@delta_pose
                actual_block=np.concatenate([pack(data[[i]]-target[[i]]) for i in range(len(source_ids))])
                result=dict(case='finite_compensation_profile',mode=mode,wall_step=step,
                    endpoint_pose=fit.x,wall_epsinf=eps,success=bool(fit.success),status=int(fit.status),
                    message=str(fit.message),nfev=fit.nfev,njev=fit.njev,
                    endpoint_white=float(np.linalg.norm(residual)),local_predicted_white=float(np.linalg.norm(predicted)),
                    actual_minus_local_white=float(np.linalg.norm(actual_block-predicted)),
                    bounds=[lo,hi],finite_nonconvergence_retained=not bool(fit.success),
                    scope='single complex RX8 mechanism; not full-aperture ambiguity')
                label=f'PROFILE_{mode}_{step:+.3f}'
                np.savez_compressed(out/(label+'.npz'),data=data,target=target,pose=fit.x,wall_epsinf=eps,
                                    residual=residual,source_block_residual=actual_block,local_prediction=predicted,jacobian=jac)
                if abs(step)==.01 and mode in ('single_left','shared_twoTX'):
                    fine=Evaluator(book,'E0E1',h=.17,n=384,grid_refinement=3,tx=mini_tx,rx=mini_rx)
                    fe=fine.evaluate(replace(changed,pose=fit.x),1.)
                    result['independent_grid_white']=white(fe.data-data)
                    result['independent_grid_exceeds_point1']=result['independent_grid_white']>.1
                    np.savez_compressed(out/(label+'_INDEPENDENT.npz'),data=fe.data,pose=fit.x,wall_epsinf=eps)
                    fine.clear();del fe;gc.collect()
                record(result)
            except Exception as error:
                record(dict(case='finite_compensation_profile_failure',mode=mode,wall_step=step,
                            error=type(error).__name__+': '+str(error),retained=True))
                # Budget exhaustion must stop, never be converted to a skipped
                # profile that permits unbudgeted subsequent physical actions.
                from wslam.accounting import BudgetExceeded
                if isinstance(error,BudgetExceeded):raise
            finally:
                for sl in solvers:sl.clear()
                gc.collect()
    record(dict(case='driver_complete',overall_gate='NOT_ASSIGNED_PARENT_REVIEW',
                algorithm_performance_claim=False,true_posterior_certificate=False))


if __name__=='__main__':main()

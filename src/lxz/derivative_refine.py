"""Final five-state derivative qualification; no new search or target change."""
import json,sys,time
import numpy as np
from .local2d import Book,MaterialPrior,prior,OUT,pack,dump

def run():
    dest=OUT/'results/E2_derivative_refinement';dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE.json').exists():raise RuntimeError('final audit already completed')
    book=Book();start=time.perf_counter();assets=OUT/'delegated/e0_2d';name='250033_weak'
    meta=json.loads((assets/f'PUBLIC_TRAIN_CALIBRATION_{name}.json').read_text())
    z=np.load(assets/f'PUBLIC_TRAIN_CALIBRATION_{name}.npz');sigma=float(z['sigma'])
    chart=MaterialPrior(meta['generator'],'weak');u=chart.coordinates();x,row=chart.decode(u)
    sys.path.insert(0,str(prior.A30/'src'));from a30.physics import FrozenMaxwellAdapter
    model=FrozenMaxwellAdapter.from_background_npz(prior.INPUTS/'background/weak.npz',book,purpose_prefix='LXZ_DERIVATIVE_FINAL')
    book.context=dict(anchor=name,method='derivative_refinement',role='base_rebuild_charged')
    base=model.full_state(x,purpose='derivative_audit_base');T=model.receiver_transfer(base)
    J=np.einsum('qp,ps->qsp',T,base.field).reshape(256,576)
    D=chart.analytic_derivative(u);Dfd=chart.derivative(u);A=J@D
    DM=np.r_[D.real,D.imag]/24;_,s,V=np.linalg.svd(DM,full_matrices=False);rank=int(np.sum(s>s[0]*1e-7))
    basis=V[:rank].T/s[:rank];_,sv,W=np.linalg.svd(np.r_[A.real,A.imag]@basis,full_matrices=False)
    v=basis@W[-1];v/=np.sqrt(np.mean(abs(D@v)**2));h=min(1e-3,.02/np.linalg.norm(v));rows=[]
    arrays=dict(chi=x,prediction=base.data,analytic_Dg=D,finite_Dg=Dfd,J_material=J,weak_direction=v,
        material_singular_values=s,prior_u=u,response_per_material_singular_values=sv)
    for k,step in enumerate([h,h/10]):
        fields=[]
        for sign in [-1,1]:
            xp,r=chart.decode(u+sign*step*v);book.context=dict(anchor=name,method='derivative_refinement',step=step,sign=sign)
            ss=model.full_state(xp,purpose='actual_centered_difference');fields.append(ss.data)
            key=f'h{k}_{sign}';arrays[key+'_chi']=xp;arrays[key+'_prediction']=ss.data
            np.savez_compressed(dest/'RAW_ARRAYS.npz',**arrays)
        fd=(fields[1]-fields[0])/(2*step);jv=(A@v).reshape(32,8)
        rows.append(dict(step=step,parameter_step_norm=float(step*np.linalg.norm(v)),relative_error=float(np.linalg.norm(fd-jv)/np.linalg.norm(jv)),
            analytic_Jv_norm=float(np.linalg.norm(jv)),whitened_Jv_norm=float(np.linalg.norm(jv)*np.sqrt(2)/sigma)))
    result=dict(anchor=name,scope='derivative qualification only; does not retroactively rerun method trajectories',
        old_finite_difference_Dg_relative_error=float(np.linalg.norm(D-Dfd)/np.linalg.norm(D)),
        old_FD_error_in_new_weak_direction=float(np.linalg.norm(J@(Dfd-D)@v)/np.linalg.norm(J@D@v)),
        exact_peak_normalization_derivative=True,rank=rank,material_singular_values=s,audits=rows,
        total_budget=book.totals(),wall_seconds=time.perf_counter()-start,search_states_added=0,paid_states_added=5,paid_RHS_added=72)
    dump(dest/'COMPLETE.json',result);book.append(dict(event='JOB_END',controlled_wall_seconds=time.perf_counter()-start,counts=dict(book.counts)))
    print(json.dumps({k:v for k,v in result.items() if k not in ['material_singular_values','audits']}),flush=True)

if __name__=='__main__':run()

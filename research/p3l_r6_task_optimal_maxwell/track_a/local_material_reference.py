"""Noise-adapted, volume-correct material-coordinate quadrature.

Different integration path from the global source panels. Retains exact prior
fibre and BOTH overlapping source labels. This is Q1/Q2 refinement evidence,
not a global interval certificate. Gaussian whitening is a mature baseline.
"""
import sys,time,json
from pathlib import Path
import numpy as np
from numpy.polynomial.hermite import hermgauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from shared.physics import OUT,Book,Target,write,save,SIGMA
from track_a.stratified import fibre,PH,RLO,RHI,ILO,IHI,quotient_representative

def prior_tasks(V):
    R,I=V;f,phase=fibre(R,I,96)
    pure=PH/((RHI-RLO)*(IHI-ILO)) if RLO<R<RHI and ILO<I<IHI else 0.
    return f+pure,(np.array([0.,f,pure]),phase+pure)

def main():
    write(OUT/'track_a/LOCAL_REFINEMENT_PROTOCOL.json',dict(reason='global panels resolve host seam but miss narrow regular likelihood; compare ordinary GN-whitened quadrature, not another chart',before_queries=time.time(),cases=['REGULAR','PURE_OVERLAP'],orders=[15,21],new_RHS_limit_remaining=13448,all_hyperparameters_fixed=True,teacher_scope='3D fixed submodel material pushforward only'))
    b=Book('A','local_material_GN_reference',role='LOW_DIMENSIONAL_MODEL_REFERENCE',max_wall=1200)
    status='COMPLETE';error=None
    try:
        target=Target(b);obs=np.load(OUT/'track_a/OBSERVATIONS.npz');records=[]
        for name in ['REGULAR','PURE_OVERLAP']:
            index=list(obs['names']).index(name);y=obs['observation'][index];sigma=float(obs['sigma'][index])
            target.sigma=sigma
            # Choose initial material from previous paid response nodes solely
            # by this observation's likelihood, never by generating latent.
            best_phi=np.inf;V=None
            for path in (OUT/'track_a/raw').glob('S3*.npz'):
                with np.load(path) as z:
                    ph=np.sum(np.abs(z['prediction']-y)**2,axis=(1,2))/sigma**2;j=int(np.argmin(ph))
                    if ph[j]<best_phi:V=np.array([z['V'][j].real,z['V'][j].imag]);best_phi=float(ph[j])
            trace=[];current=None
            for iteration in range(20):
                current=target.evaluate(quotient_representative(*V),'local_material_mode')
                Q=float(target.prior.decode(current.q)[1]['fraction'][target.prior.mask][0])
                J=target.jacobian(current);JV=np.column_stack((J[:,64]/Q,J[:,65]/(.15*Q)))
                r=target.white(current.state.data-y);step=-np.linalg.solve(JV.T@JV,JV.T@r)
                trace.append(dict(iteration=iteration,V=V,phi=float(.5*r@r),step=step))
                if np.linalg.norm(step)<1e-11:break
                if np.linalg.norm(step)>.01:step*=.01/np.linalg.norm(step)
                accepted=False
                for alpha in [1.,.5,.25,.125]:
                    proposal=V+alpha*step
                    try:q=quotient_representative(*proposal)
                    except ValueError:continue
                    candidate=target.evaluate(q,'local_material_mode_line_search')
                    if target.phi(candidate,y)<target.phi(current,y):V=proposal;accepted=True;break
                if not accepted:break
            # Recompute at final point (all attempts, including repeated points,
            # are explicitly paid; no phantom cache discount).
            current=target.evaluate(quotient_representative(*V),'local_material_whitening')
            Q=float(target.prior.decode(current.q)[1]['fraction'][target.prior.mask][0]);J=target.jacobian(current)
            JV=np.column_stack((J[:,64]/Q,J[:,65]/(.15*Q)));H=JV.T@JV;C=np.linalg.inv(H);chol=np.linalg.cholesky(C)
            center_phi=target.phi(current,y);logdet=float(np.linalg.slogdet(C)[1]);base_before=b.totals()['RHS']
            for order in [15,21]:
                gh,w=hermgauss(order);samples=[];weight=[];labels=[];phase=[];phis=[]
                for i,a in enumerate(gh):
                    for j,c in enumerate(gh):
                        z=np.sqrt(2)*np.array([a,c]);v=V+chol@z
                        try:q=quotient_representative(*v)
                        except ValueError:continue
                        density,(bran,event)=prior_tasks(v)
                        p=target.evaluate(q,'material_prior_fibre_GH');phi=target.phi(p,y)
                        logq=-np.log(2*np.pi)-.5*logdet-.5*(z@z)
                        ww=w[i]*w[j]/np.pi*np.exp(center_phi-phi-logq)*density
                        samples.append(v);weight.append(ww);labels.append(bran/max(density,1e-300));phase.append(event/max(density,1e-300));phis.append(phi)
                samples=np.asarray(samples);weight=np.asarray(weight);prob=weight/weight.sum();mu=prob@samples
                sd=float(np.sqrt(np.sum(prob*np.sum((samples-mu)**2,axis=1))))
                rows=dict(observation=name,order=order,mean_increment=mu,material_ROI_SD=sd,
                          log_evidence=float(np.log(weight.sum())-center_phi),source_layer_mass=prob@np.asarray(labels),
                          latent_phase_gt_half=float(prob@phase),material_real_increment_gt_002=float(prob@(samples[:,0]>.02)),
                          normalizer_quadrature='prior-fibre-correct Gaussian importance quadrature; not empirical-node TV',
                          material_Gaussian_mode=V,Gaussian_mode_mean_gap_SD=float(np.linalg.norm(V-mu)/max(sd,1e-30)),
                          total_additional_RHS_so_far=b.totals()['RHS']-base_before)
                records.append(rows);save(OUT/f'track_a/raw/{name}_material_GH{order}.npz',V=samples,weights=weight,phi=np.asarray(phis),source_probability=np.asarray(labels),mode=V,C=C)
                write(OUT/'track_a/LOCAL_REFINEMENT_PARTIAL.json',dict(rows=records));print(json.dumps(dict(case=name,order=order,mean=mu.tolist(),RHS=b.totals()['stages']['A']['RHS'])),flush=True)
            write(OUT/f'track_a/{name}_MODE_TRACE.json',dict(trace=trace,H=H,mode=V,GT_read=False))
        write(OUT/'track_a/LOCAL_REFINEMENT_RESULTS.json',dict(status='COMPLETE_Q1_Q2_LOCAL_REFINEMENT',rows=records,global_tail_certificate=False))
    except BaseException as exc:status='FAILED';error=repr(exc);write(OUT/'track_a/LOCAL_REFINEMENT_FAILURE.json',dict(error=error));raise
    finally:b.close(status,error)

if __name__=='__main__':main()

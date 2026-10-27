"""Original-coordinate and gauge-integrated P1 submodels, Full TM24 only.

No response-chart inversion. The 3D uniform-field submodel is integrated both
in its original three coordinates and in 2D material amplitude plus an EXACT
one-dimensional prior fibre. All quadrature diagnostics are Q1/Q2, not interval
certificates. Finite node measures are never assigned TV against a density.
"""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.special import ndtr,ndtri,betaincinv
from scipy.stats import norm
from numpy.polynomial.legendre import leggauss
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from shared.physics import OUT,Book,Target,SIGMA,write,save,pack
from p1_material.prior import saturation

LO,HI=-4.,4.
Z0=ndtr(HI)-ndtr(LO);A=1/1.3
PH=(ndtr(-A)-ndtr(LO))/Z0
ULO,UHI=ndtr(LO),ndtr(HI)
RLO,RHI=.8+ULO,.8+UHI
ILO,IHI=.03+.15*ULO,.03+.15*UHI

def gl(n,a,b):
    t,w=leggauss(n);return a+(b-a)*(t+1)/2,w*(b-a)/2

def inverse_fraction(q):
    return (2*betaincinv(3,3,np.asarray(q))-1)/1.3

def fibre(R,I,n=48):
    """Transition material density and phase>1/2 numerator at V=R+iI.

    dV Jacobian is q^2. Integrating z directly avoids the apparent q=1
    inverse-derivative singularity. Both properties are truncated-Gaussianized
    uniforms; their exact widths, not untruncated intervals, appear below.
    """
    lo=max(R/RHI,I/IHI);hi=min(R/RLO,I/ILO,1.)
    if not 0<lo<hi:return 0.,0.
    zl,zh=inverse_fraction(lo),inverse_fraction(hi)
    z,w=gl(n,float(zl),float(zh));q=saturation(1.3*z)
    density=float(np.sum(w*norm.pdf(z)/Z0/q**2))/((RHI-RLO)*(IHI-ILO))
    if zh<=0:phase=0.
    elif zl>=0:phase=density
    else:
        z,w=gl(n,0.,float(zh));q=saturation(1.3*z)
        phase=float(np.sum(w*norm.pdf(z)/Z0/q**2))/((RHI-RLO)*(IHI-ILO))
    return density,phase

def quotient_representative(R,I):
    lo=max(R/RHI,I/IHI);hi=min(R/RLO,I/ILO,1.)
    if not 0<lo<=hi:raise ValueError('outside uniform material image')
    q=.5*(lo+hi);z=float(inverse_fraction(q))
    native=np.r_[np.zeros(64),R/q-.8,(I/q-.03)/.15];native[0]=z
    return native

class UniformPhysics:
    def __init__(self,target):self.target=target;self.cache={};self.calls=0
    def __call__(self,V,native=None):
        # Cache only exact requested floating material amplitudes, no rounding.
        V=complex(V);key=np.complex128(V).tobytes()
        if key not in self.cache:
            q=(np.r_[-2.,np.zeros(63),.5,.5] if V==0 else quotient_representative(V.real,V.imag)) if native is None else native
            point=self.target.evaluate(q,'stratified_uniform_material')
            expected=self.target.prior.host+V*self.target.prior.mask
            error=float(np.max(np.abs(point.state.chi-expected)))
            if error>2e-13:raise RuntimeError('representative is not requested material')
            self.cache[key]=(point.state.data.copy(),error);self.calls+=1
        return self.cache[key][0]

def source_nodes(dimension,nz,nu,seam_scale=SIGMA**(1/3),scaled=True):
    """Disjoint host/pure/transition prior integration. Host is one atom."""
    v=[0j];w=[PH];layer=[0];phase=[0.]
    ur,wr=gl(nu,ULO,UHI);wr=wr/(UHI-ULO)
    ui,wi=(gl(nu,ULO,UHI) if dimension==3 else (np.array([.5]),np.array([1.])))
    if dimension==3:wi=wi/(UHI-ULO)
    # The full-Gaussian phase coordinate was explicitly truncated, not silently
    # dropped. Pure offset gauge has exactly PH probability.
    for r,a in zip(ur,wr):
        for im,b in zip(ui,wi):v.append(complex(.8+r,.03+.15*im));w.append(PH*a*b);layer.append(2);phase.append(1.)
    cuts=sorted(set([0.,min(2*seam_scale,2*A),min(5*seam_scale,2*A),A,2*A])) if scaled else [0.,A,2*A]
    for left,right in zip(cuts[:-1],cuts[1:]):
        h,wh=gl(nz,left,right);z=-A+h;weights=wh*norm.pdf(z)/Z0
        for zz,c in zip(z,weights):
            Q=float(saturation(1.3*zz))
            for r,a in zip(ur,wr):
                for im,b in zip(ui,wi):
                    v.append(Q*complex(.8+r,.03+.15*im));w.append(c*a*b);layer.append(1);phase.append(float(Q>.5))
    return dict(V=np.asarray(v),prior_weights=np.asarray(w),layer=np.asarray(layer),phase_numerator=np.asarray(w)*np.asarray(phase))

def quotient_nodes(nr,nt,nf=48,seam_scale=SIGMA**(1/3)):
    """Cubic radial coordinate cancels the q^(-5/3) area-density singularity."""
    v=[0j];w=[PH];layer=[0];phase=[0.]
    # Piece boundaries include all changes in the rectangle's angular fibre.
    angles=sorted(set([ILO/RHI,ILO/RLO,IHI/RHI,IHI/RLO]))
    tmax=RHI**(1/3)
    for th0,th1 in zip(angles[:-1],angles[1:]):
        theta,wa=gl(nt,th0,th1)
        for th,b in zip(theta,wa):
            # At these radii the fibre touches q=1 or switches a property
            # boundary. Explicit panels avoid integrating hidden step/kinks.
            rmax=min(RHI,IHI/th)
            radial=[0.,rmax,RLO,ILO/th,(2*seam_scale)**3,(5*seam_scale)**3]
            cuts=sorted(set([float(r**(1/3)) for r in radial if 0<=r<=rmax]))
            for left,right in zip(cuts[:-1],cuts[1:]):
                t,wt=gl(nr,left,right)
                for u,au in zip(t,wt):
                    R=float(u**3)
                    I=th*R;density,event=fibre(R,I,nf)
                    if density>0:
                        jac=3*u**5*au*b
                        v.append(complex(R,I));w.append(jac*density);layer.append(1);phase.append(jac*event)
    # Exact pure offset probability, ordinary smooth property integral. It is
    # added to transition density on overlapping material images, not discarded.
    rr,wr=gl(nr,RLO,RHI);ii,wi=gl(nt,ILO,IHI)
    for R,a in zip(rr,wr):
        for I,b in zip(ii,wi):
            p=PH*a*b/((RHI-RLO)*(IHI-ILO))
            v.append(complex(R,I));w.append(p);layer.append(2);phase.append(p)
    return dict(V=np.asarray(v),prior_weights=np.asarray(w),layer=np.asarray(layer),phase_numerator=np.asarray(phase))

def evaluate_nodes(physics,nodes):
    y=np.stack([physics(v) for v in nodes['V']]);return dict(nodes,prediction=y)

def summarize(nodes,observation,sigma,reference_phi):
    pred=nodes['prediction'];difference=pred-observation
    phi=np.sum(np.abs(difference)**2,axis=(1,2))/sigma**2
    # Reference shift includes exact host; no underflow is interpreted as proof.
    L=np.exp(np.clip(reference_phi-phi,-745,700));a=nodes['prior_weights']*L;Z=float(a.sum());prob=a/Z
    V=nodes['V'];mean=np.sum(prob*V);var=float(np.sum(prob*np.abs(V-mean)**2))
    branches=[float(prob[nodes['layer']==k].sum()) for k in (0,1,2)]
    phase=float(np.sum(nodes['phase_numerator']*L)/Z)
    event=float(np.sum(prob*(V.real>.02)))
    fw=np.stack([pack(x)*np.sqrt(2) for x in pred]);fm=prob@fw
    centeredV=np.column_stack((V.real-mean.real,V.imag-mean.imag))
    cov=(centeredV.T*prob)@(fw-fm)
    between=np.zeros_like(cov)
    for k,mass in enumerate(branches):
        if mass>1e-290:
            keep=nodes['layer']==k;weights=prob[keep]/mass
            vm=np.sum(weights*V[keep]);yk=weights@fw[keep]
            between+=mass*np.outer([vm.real-mean.real,vm.imag-mean.imag],yk-fm)
    return dict(log_evidence=float(np.log(Z)-reference_phi),scaled_evidence=Z,mean_increment=[float(mean.real),float(mean.imag)],
                material_ROI_SD=float(np.sqrt(var)),source_layer_mass=branches,latent_phase_gt_half=phase,
                material_real_increment_gt_002=event,prior_mass_quadrature=float(nodes['prior_weights'].sum()),
                within_task_response_covariance_norm=float(np.linalg.norm(cov-between)),
                between_task_response_covariance_norm=float(np.linalg.norm(between)),
                total_task_response_covariance_norm=float(np.linalg.norm(cov)),
                posterior_mean_observation_derivative_norm=float(np.linalg.norm(cov)/sigma**2),
                covariance_qualification='quadrature conditional covariance of low-dimensional model only',
                underflow_nodes=int(np.sum(reference_phi-phi<-700)),node_count=len(V))

def compensation(target,physics):
    q=np.r_[np.zeros(64),.4,.45];q[0]=1.
    base=target.evaluate(q,'pure_seam_compensation_base');J=target.jacobian(base)
    delta=complex(1.2,.0975);c=-2.74625*(delta.real*J[:,64]+delta.imag/.15*J[:,65])
    rows=[]
    for dims in ([64],[64,65]):
        Aop=J[:,dims];Q,_=np.linalg.qr(Aop);perp=c-Q@(Q.T@c)
        rows.append(dict(strong_columns=dims,c_norm=float(np.linalg.norm(c)),normal_norm=float(np.linalg.norm(perp)),
                         normal_ratio=float(np.linalg.norm(perp)/np.linalg.norm(c)),singular_values=np.linalg.svd(Aop,compute_uv=False)))
    trials=[]
    for h in [.08,.04,.02,.01,.005]:
        native=q.copy();native[0]=A-h;fraction=float(saturation(1.3*native[0]));un=target.evaluate(native,'uncompensated_upper_seam')
        adjusted=native.copy();adjusted[64]=delta.real/fraction-.8;adjusted[65]=(delta.imag/fraction-.03)/.15
        comp=target.evaluate(adjusted,'exact_material_fibre_compensation')
        trials.append(dict(h=h,fraction=fraction,uncompensated_white_norm=float(np.linalg.norm(target.white(un.state.data-base.state.data))),
                           compensated_white_norm=float(np.linalg.norm(target.white(comp.state.data-base.state.data))),
                           compensated_material_RMS=float(np.linalg.norm(comp.state.chi-base.state.chi)/24),
                           original_prior_penalty_change=float(.5*(native[:64]@native[:64]-q[:64]@q[:64])),
                           legal=True,interpretation='same-material gauge, not a second material world'))
    write(OUT/'track_a/COMPENSATION.json',dict(operators=rows,finite_trials=trials,
          consequence='two adjustable properties remove the whole cubic material coefficient; exact QD gauge, no claim of a new exact material fibre'))

def main():
    cfg=dict(created_before_queries=time.time(),prior='active standard normals truncated[-4,4], remainder zero',
             models={'S1':[0,64],'S3':[0,64,65]},observations='controlled host and regular worlds, not an SBC population',
             source_levels_S1=[[12,12],[24,16],[32,24]],source_levels_S3=[[10,8],[16,9]],
             quotient_levels=[[12,6,48],[20,8,64]],noise_scales=[1.,.5,.25],
             max_RHS=160000,reference_status='refinement must qualify per task; no global interval certificate')
    write(OUT/'track_a/FROZEN_PROTOCOL.json',cfg)
    b=Book('A','stratified_uniform_submodels',role='LOW_DIMENSIONAL_MODEL_REFERENCE',max_wall=3600)
    status='COMPLETE';error=None;tick=time.perf_counter()
    try:
        target=Target(b);physics=UniformPhysics(target);host=physics(0j)
        rng=np.random.default_rng(202610270611);innovation=(rng.normal(size=(32,8))+1j*rng.normal(size=(32,8)))/np.sqrt(2)
        observations=[]
        for scale in cfg['noise_scales']:
            sigma=SIGMA*scale;y=host+sigma*innovation;observations.append((f'HOST_sigma{scale}',y,sigma))
        # These are mechanism observations with fixed public worlds, not an
        # online anchor or unbiased prior-predictive population assessment.
        for name,V in [('REGULAR',complex(.12,.015)),('PURE_OVERLAP',complex(1.2,.0975))]:
            mean=physics(V);observations.append((name,mean+SIGMA*innovation,SIGMA))
        save(OUT/'track_a/OBSERVATIONS.npz',observation=np.stack([a[1] for a in observations]),sigma=np.array([a[2] for a in observations]),names=np.array([a[0] for a in observations]))
        rows=[]
        schedules=[]
        for nz,nu in cfg['source_levels_S1']:schedules.append((f'S1_scaled_{nz}_{nu}',lambda nz=nz,nu=nu:source_nodes(2,nz,nu)))
        schedules.append(('S1_unscaled_32_24',lambda:source_nodes(2,32,24,scaled=False)))
        for nz,nu in cfg['source_levels_S3']:schedules.append((f'S3_source_{nz}_{nu}',lambda nz=nz,nu=nu:source_nodes(3,nz,nu)))
        for nr,nt,nf in cfg['quotient_levels']:schedules.append((f'S3_quotient_{nr}_{nt}',lambda nr=nr,nt=nt,nf=nf:quotient_nodes(nr,nt,nf)))
        for name,make in schedules:
            start=time.perf_counter();before=Book.totals()['RHS'];nodes=evaluate_nodes(physics,make())
            save(OUT/'track_a/raw'/f'{name}.npz',**nodes)
            for oname,y,sigma in observations:
                ref=float(np.sum(np.abs(host-y)**2)/sigma**2)
                # For nonhost observations a finite non-underflow reference
                # shift is the minimum over available exact Full responses.
                ref=min(ref,float(np.min(np.sum(np.abs(nodes['prediction']-y)**2,axis=(1,2))/sigma**2)))
                summary=summarize(nodes,y,sigma,ref)
                rows.append(dict(method=name,observation=oname,sigma=sigma,**summary,
                                 additional_shared_RHS=Book.totals()['RHS']-before,method_wall_seconds=time.perf_counter()-start))
            write(OUT/'track_a/PARTIAL_RESULTS.json',dict(rows=rows,used_RHS=Book.totals()['stages']['A']['RHS']))
            print(json.dumps(dict(event='A_QUADRATURE',method=name,nodes=len(nodes['V']),total_RHS=Book.totals()['RHS'],stage_RHS=Book.totals()['stages']['A']['RHS'])),flush=True)
        compensation(target,physics)
        write(OUT/'track_a/RESULTS.json',dict(status='COMPLETED_WITH_PER_CASE_REFERENCE_QUALIFICATION_PENDING',rows=rows,
                    new_RHS=Book.totals()['stages']['A']['RHS'],wall_seconds=time.perf_counter()-tick,
                    source_model_scope='2D/3D fixed submodels, never Full66 marginal',new_4D_model='NOT_RUN_UNTIL_2D_3D_REFERENCE_GATES_PASS'))
    except BaseException as exc:status='FAILED';error=repr(exc);write(OUT/'track_a/FAILURE.json',dict(error=error));raise
    finally:b.close(status,error)

if __name__=='__main__':main()

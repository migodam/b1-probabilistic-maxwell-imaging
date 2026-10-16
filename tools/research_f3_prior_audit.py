"""Independent quotient-prior audit: decoder and cached C/D algebra only.
Reads online generator JSON, cached F1 p/C/D, never truth, LU, field or RHS.
"""
from pathlib import Path
import copy, json, sys, time, os
import numpy as np
from scipy.integrate import quad
from scipy.stats import kstest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from prior_constrained import prior
A=np.log(9/4); LO=np.r_[np.full(4,-.13),np.full(2,.05),.4,.4]; SP=np.r_[np.full(4,.26),.045,.045,.5,.5]
STEPS=(1e-3,5e-4,2.5e-4)

def density(r):
    r=np.asarray(r); z=np.abs(r)
    return np.where(z<=A,np.maximum(0,(.81*np.exp(-z)-.16*np.exp(z))/.5),0.)
def cdf(r):
    r=np.asarray(r); z=np.minimum(np.abs(r),A)
    v=(.81*(-np.expm1(-z))-.16*np.expm1(z))/.5
    return np.clip(.5+np.sign(r)*v,0,1)
def q_of_p(p):
    w=LO[6:]+SP[6:]*p[6:]
    return np.r_[p[:6],(np.log(w[0]/w[1])+A)/(2*A)]
def physical(q):
    r=2*A*q[6]-A
    return np.r_[LO[:6]+SP[:6]*q[:6],.6*np.exp(r/2),.6*np.exp(-r/2)]
def row_from_a(meta,a):
    r=copy.deepcopy(meta['generator']);d=r['geometry']
    d['gaussian_centers'][:2]=a[:4].reshape(2,2).tolist(); d['gaussian_widths'][:2]=a[4:6].tolist();d['gaussian_weights'][:2]=a[6:].tolist()
    return r

def decode(meta,a):return prior.decode(row_from_a(meta,a),background_id=meta['background_id'])
def derivative(meta,q):
    a=physical(q);d=meta['generator']['geometry'];p=prior.geometry_points()-np.array(d['center']);c,s=np.cos(d['rotation']),np.sin(d['rotation'])
    xx=(c*p[:,0]+s*p[:,1])/d['scale'];yy=(-s*p[:,0]+c*p[:,1])/(d['scale']*d['axis_ratio'])
    vals=np.zeros(len(p)); dv=np.zeros((len(p),7))
    for i in range(2):
        dx=xx-a[2*i];dy=yy-a[2*i+1];rr=dx*dx+dy*dy;w=a[6+i];width=a[4+i];g=w*np.exp(-rr/(2*width**2));vals+=g
        dv[:,2*i]=g*dx/width**2*SP[2*i];dv[:,2*i+1]=g*dy/width**2*SP[2*i+1];dv[:,4+i]=g*rr/width**3*SP[4+i];dv[:,6]+=(1 if i==0 else -1)*A*g
    im=np.argmax(vals);peak=vals[im]
    return 1.2*prior.V*(dv/peak-vals[:,None]*dv[im]/peak**2)
def mass(D):return np.r_[D.real,D.imag]/24
def rel(a,b):return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-300))
def svrank(D):
    u,s,v=np.linalg.svd(D,full_matrices=False);n=int(sum(s>s[0]*1e-7));return u[:,:n],s,v[:n].T/s[:n],n

def main():
    start=time.perf_counter();rng=np.random.default_rng(20261017);out=ROOT/'research/research_f3_quotient_exploration';out.mkdir(parents=True,exist_ok=True)
    weights=rng.uniform(.4,.9,(400000,2));rs=np.log(weights[:,0]/weights[:,1]);ks=kstest(rs,cdf);integ=quad(lambda x:float(density(x)),-A,A,points=[0],epsabs=1e-13)
    eps=np.array([1e-2,1e-3,1e-4,1e-5]);checks={'scope':'ZERO_PHYSICS_DECODER_AND_CACHED_C_D_ONLY','seed':20261017,'monte_carlo_n':len(rs),'quadrature_mass':integ[0],'quadrature_error_estimate':integ[1],'cdf_endpoints':cdf(np.array([-A,0,A])).tolist(),'cdf_derivative_max_error':float(max(abs((cdf(x+1e-6)-cdf(x-1e-6))/2e-6-density(x)) for x in np.linspace(-A+1e-3,A-1e-3,1000))),'monte_carlo_KS':{'statistic':float(ks.statistic),'pvalue':float(ks.pvalue),'DKW_99pct_bound':float(np.sqrt(np.log(200)/2/len(rs)))},'boundary':{'epsilon':eps.tolist(),'density_over_distance':(density(A-eps)/eps).tolist(),'limit':1.44,'tail_mass_over_distance_squared':((1-cdf(A-eps))/eps**2).tolist(),'tail_limit':.72},'normalized_t_density_jacobian':float(2*A),'t_density_integral':quad(lambda t:float(2*A*density(2*A*t-A)),0,1,points=[.5])[0],'canonical_endpoint_weights':[physical(np.r_[np.full(6,.5),t])[6:].tolist() for t in (0,1)],'cases':{},'reads':[],'new_full_states':0,'new_LU_factorizations':0,'new_LU_solves':0,'new_complex_RHS':0,'new_Mie_calls':0,'GT_arrays_read':False}
    legacy=Path(os.environ.get('B1_F3_LEGACY_ROOT',str(ROOT.parent/'b1-research-f/research/research_f_intrinsic_transport')))
    old=legacy/'results/F1'
    for bg in ('weak','ring'):
        case=f'p250073_{bg}_a1_snr20';mp=legacy/'inputs/online'/f'{case}.json';meta=json.loads(mp.read_text());checks['reads'].append(str(mp));rr=[];endpoint={}
        for j in (0,11):
            sp=old/case/f'state_{j:02}.npz';jp=sp.with_name(sp.stem+'_JAC.npz')
            with np.load(sp) as z:p=z['p']
            with np.load(jp) as z:C=z['C'];D=z['D']
            checks['reads'] += [f'{sp}:p',f'{jp}:C,D'];q=q_of_p(p);Dq=derivative(meta,q);L=np.zeros((8,7));L[:6,:6]=np.eye(6);w=LO[6:]+SP[6:]*p[6:];L[6,6]=A*w[0]/.5;L[7,6]=-A*w[1]/.5
            go=np.r_[np.zeros(6),w/.5];Q,sv,E,n=svrank(mass(D));Qq,sq,Eq,nq=svrank(mass(Dq));B=C@D@E;Bq=C@Dq@Eq
            # Singular values of complex response alone are not the real-material response. Pack real and imaginary before SVD.
            sb=np.linalg.svd(np.r_[B.real,B.imag],compute_uv=False);sbq=np.linalg.svd(np.r_[Bq.real,Bq.imag],compute_uv=False)
            fd=[];direction=rng.normal(size=7);direction/=np.linalg.norm(direction)
            # Shape-only at endpoints avoids claiming two-sided support at a ratio boundary.
            direction[6]=0;direction/=np.linalg.norm(direction)
            for h in STEPS:
                fd.append({'h':h,'relative_error':rel((decode(meta,physical(q+h*direction))-decode(meta,physical(q-h*direction)))/(2*h),Dq@direction)})
            endpoint[str(j)]={'q':q.tolist(),'image_equivalence_relative_error':rel(decode(meta,physical(q)),decode(meta,LO+SP*p)),'analytic_chain_relative_error':rel(Dq,D@L),'old_material_rank':n,'quotient_material_rank':nq,'old_material_singular_values':sv.tolist(),'quotient_material_singular_values':sq.tolist(),'old_retained_material_condition':float(sv[0]/sv[6]),'quotient_material_condition':float(sq[0]/sq[-1]),'material_gauge_relative_defect':float(np.linalg.norm(D@go)/np.linalg.norm(D)),'cached_response_gauge_relative_defect':float(np.linalg.norm(C@D@go)/np.linalg.norm(C@D)),'tangent_subspace_relative_defect':rel(Q@(Q.T@Qq),Qq),'material_metric_response_sv_relative_error':rel(sbq,sb),'FD_shape_only':fd}
        for k in range(256):
            p=rng.uniform(0,1,8);q=q_of_p(p);rr.append(rel(decode(meta,physical(q)),decode(meta,LO+SP*p)))
        # Full seven-coordinate derivative at a strict interior, fixed independently.
        q=np.r_[rng.uniform(.2,.8,6),.5];Dq=derivative(meta,q);dire=rng.normal(size=7);dire/=np.linalg.norm(dire)
        fullfd=[{'h':h,'relative_error':rel((decode(meta,physical(q+h*dire))-decode(meta,physical(q-h*dire)))/(2*h),Dq@dire)} for h in STEPS]
        # One-sided second-order boundary derivative respecting support.
        bfd=[]
        for t,sgn in [(0,1),(1,-1)]:
            q=np.r_[[.3,.4,.7,.6,.35,.65],t];f0=decode(meta,physical(q));e=np.eye(7)[6]*sgn;dq=derivative(meta,q)@e
            bfd.append({'t':t,'steps':[{'h':h,'relative_error':rel((-3*f0+4*decode(meta,physical(q+h*e))-decode(meta,physical(q+2*h*e)))/(2*h),dq)} for h in STEPS]})
        same_component_errors=[rel(decode(meta,physical(np.r_[np.full(6,.5),t])),decode(meta,physical(np.r_[np.full(6,.5),.5]))) for t in (0.,.25,.75,1.)]
        checks['cases'][bg]={'coincident_components_ratio_image_errors':same_component_errors,'cached_endpoints':endpoint,'random_support_n':len(rr),'random_support_max_image_relative_error':max(rr),'FD_full_seven_interior':fullfd,'FD_ratio_boundary_one_sided':bfd}
    checks['backend_modules_loaded']=[m for m in sys.modules if any(x in m for x in ('a30.physics','finite_feedback','lxz.local2d','research_f.common'))]
    checks['elapsed_seconds']=time.perf_counter()-start
    checks['derivative_and_cached_algebra_pass']=bool(all(all(e['analytic_chain_relative_error']<1e-12 and e['tangent_subspace_relative_defect']<1e-11 and e['material_metric_response_sv_relative_error']<1e-10 and e['old_material_rank']==7 and e['quotient_material_rank']==7 and max(s['relative_error'] for s in e['FD_shape_only'])<1e-4 for e in v['cached_endpoints'].values()) and max(s['relative_error'] for s in v['FD_full_seven_interior'])<1e-4 and max(s['relative_error'] for b in v['FD_ratio_boundary_one_sided'] for s in b['steps'])<1e-4 for v in checks['cases'].values()))
    checks['numeric_consistency_pass']=bool(checks['derivative_and_cached_algebra_pass'] and abs(integ[0]-1)<1e-12 and ks.statistic<checks['monte_carlo_KS']['DKW_99pct_bound'] and all(v['random_support_max_image_relative_error']<1e-12 for v in checks['cases'].values()) and not checks['backend_modules_loaded'])
    (out/'PRIOR_AUDIT.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n')
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4));ax.hist(rs,bins=100,density=True,alpha=.35,label='400,000 original weight draws');x=np.linspace(-A,A,501);ax.plot(x,density(x),label='Exact quotient density');ax.axhline(1/(2*A),color='gray',ls='--',label='Uniform ratio law (different prior)');ax.set(xlabel=r'$\rho=\log(w_0/w_1)$',ylabel='Density');ax.legend();fig.tight_layout();fig.savefig(out/'quotient_prior.png',dpi=180);plt.close(fig)
    print(json.dumps({'numeric_consistency_pass':checks['numeric_consistency_pass'],'elapsed_seconds':checks['elapsed_seconds'],'KS':checks['monte_carlo_KS']}))
if __name__=='__main__':main()

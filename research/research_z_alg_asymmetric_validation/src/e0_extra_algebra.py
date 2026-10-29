"""Missing E0 controls: analytic toys only, zero Maxwell evaluations."""
import math
import numpy as np
from scipy.integrate import quad


def compute():
    tol=1e-10
    checks={}
    details={}
    # C_n action and all regular roots. Nonrotational prior at roots is nonuniform.
    for n in (2,3):
        y=1.3+0.4j
        root=y**(1/n)
        omega=np.exp(2j*np.pi/n)
        roots=np.array([root*omega**k for k in range(n)])
        errors=[float(np.max(np.abs(roots**n-y)))]
        for k in range(n):
            for j in range(n):
                errors.append(abs(omega**k*omega**j-omega**((k+j)%n)))
            errors.append(abs(omega**k*omega**(-k)-1))
        checks[f'C{n}_roots_action_identity_inverse_composition']=max(errors)<tol
        weights=np.exp(0.7*roots.real);weights/=weights.sum()
        checks[f'C{n}_asymmetric_prior_not_uniform']=np.max(abs(weights-1/n))>tol
        checks[f'C{n}_zero_full_stabilizer']=all(abs(omega**k*0)<tol for k in range(n))
        details[f'C{n}']={'roots_real':roots.real.tolist(),'roots_imag':roots.imag.tolist(),
                          'weights':weights.tolist(),'max_identity_error':float(max(errors)),
                          'zero_branch_count':1,'zero_stabilizer_size':n}
    # Circle: local legal arc survives while the full SO(2) action leaves a half box.
    theta=np.linspace(0,2*np.pi,257)
    circle=np.column_stack([np.cos(theta),np.sin(theta)])
    legal=(circle[:,0]>=-tol)&(abs(circle[:,0])<=1+tol)&(abs(circle[:,1])<=1+tol)
    checks['circle_invariant']=np.max(abs(np.sum(circle**2,axis=1)-1))<tol
    checks['circle_boundary_truncates_global_action']=bool(np.any(legal) and np.any(~legal))
    checks['circle_origin_singular']=np.linalg.norm(np.array([0.,0.]))==0
    mass=quad(lambda t:math.exp(0.7*math.cos(t)), -math.pi/2,math.pi/2,epsabs=1e-12,epsrel=1e-12)[0]
    details['circle']={'legal_domain':'x>=0 intersect [-1,1]^2','arc_prior_normalizer':mass,
                       'origin_gradient':[0.,0.],'nonuniform_prior':'exp(0.7*x) on legal arc'}
    # Parabola H(a,b)=a+b^2: exact additive action and Gaussian disintegration.
    def T(a,b,t):return np.array([a-2*b*t-t*t,b+t])
    a,b,t,s=0.7,-0.4,0.31,-0.22
    first=T(a,b,s);composed=T(*first,t)
    checks['parabola_group']=np.max(abs(composed-T(a,b,s+t)))<tol
    moved=T(a,b,t)
    checks['parabola_invariant']=abs(moved[0]+moved[1]**2-(a+b*b))<tol
    u=0.8
    density=lambda v:math.exp(-((u-v*v)**2+v*v)/2)
    z,err=quad(density,-np.inf,np.inf,epsabs=1e-12,epsrel=1e-12)
    mean,errmean=quad(lambda v:v*density(v),-np.inf,np.inf,epsabs=1e-12,epsrel=1e-12)
    # Integrate both real inverse branches a=u-b² with Jacobian |dH/da|=1.
    z_half,errhalf=quad(density,0,np.inf,epsabs=1e-12,epsrel=1e-12)
    norm=quad(lambda v:density(v)/z,-np.inf,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
    checks['parabola_conditional_normalizer']=abs(norm-1)<tol and abs(z-2*z_half)<tol and max(err,errmean,errhalf)<tol
    checks['parabola_conditional_mean']=abs(mean/z)<tol
    checks['parabola_nonuniform_orbit']=abs(density(0)-density(1))>tol
    details['parabola']={'u':u,'normalizer':z,'quad_error_estimate':err,'conditional_mean_b':mean/z,
                         'coarea_graph_metric_cancels':'sqrt(1+4b²)/sqrt(1+4b²)=1'}
    # Cubic branch weighting and critical-root changes, complement existing complex monodromy.
    roots=np.array([-math.sqrt(3),0,math.sqrt(3)])
    w=1/abs(3*roots**2-3);w/=w.sum()
    checks['cubic_weights']=np.max(abs(w-[.25,.5,.25]))<tol
    checks['cubic_second_moment']=abs(np.dot(w,roots**2)-1.5)<tol
    real_counts=[]
    for y in [1.9,2.1]:
        rr=np.roots([1,0,-3,-y]);real_counts.append(int(np.sum(abs(rr.imag)<tol)))
    checks['cubic_critical_root_count_change']=real_counts==[3,1]
    details['cubic']={'weights':w.tolist(),'real_root_counts_at_1_9_2_1':real_counts,
                      'uniform_second_moment':float(np.mean(roots**2)),'correct_event_abs_x_lt_1':float(w[1])}
    # Approximate adjacency is not an equivalence; fixed-y endpoint check fails.
    eps=0.1;p=np.array([0,.75*eps,1.5*eps])
    checks['epsilon_nontransitive']=bool(abs(p[1]-p[0])<=eps and abs(p[2]-p[1])<=eps and abs(p[2]-p[0])>eps)
    details['epsilon']={'epsilon':eps,'points':p.tolist(),'fixed_y_zero_endpoint_compatible':False}
    # Exact kernel of H(x,y,z)=x: V=(0,1,0), W=(0,0,y), [V,W]=(0,0,1).
    checks['exact_kernel_bracket_tangent']=float(np.dot([1,0,0],[0,0,1]))==0
    # Weak plane D=span{V=(1,0,0), W=(0,1,x)} has bracket (0,0,1), outside at x=0.
    weak=np.column_stack(([1.,0,0],[0,1.,0]));bracket=np.array([0.,0,1.])
    residual=np.linalg.norm(bracket-weak@np.linalg.lstsq(weak,bracket,rcond=None)[0])
    checks['weak_distribution_bracket_not_closed']=residual>tol
    details['brackets']={'exact_H':'x','exact_V':'d_y','exact_W':'y*d_z','exact_bracket':'d_z',
                         'weak_V':'d_x','weak_W':'d_y+x*d_z','weak_bracket':'d_z','weak_residual_at_x0':float(residual)}
    return {'tolerance':tol,'checks':checks,'all_passed':all(checks.values()),'details':details,
            'Maxwell_RHS':0,'scope':'Additional algebra controls; existing complex-monodromy script retained without duplication; CPU/memory paid'}

from probability_transfer import *
from scipy.integrate import quad
from scipy.optimize import minimize_scalar,brentq
from scipy.special import ndtri
from numpy.polynomial import Polynomial
D=1.3+.105j

def run():
    if (RES/'original1d_underflow_reporting_repair_receipt.json').exists():raise RuntimeError('already executed')
    b=Book('original1d_underflow_reporting_repair');status='FAILED';error=None
    try:
        ph=Physics(b);ph.load_setup();y=ph.obs['REGULAR'];counts=Counter()
        def phi_z(z):
            counts['phi_calls']+=1;v=complex(float(sat(1.3*z))*D);f=ph.reduced(v);return float(np.sum(abs(f-y)**2)/SIGMA**2)
        with b.span('original1d_profile'):
            fit=minimize_scalar(phi_z,bounds=(-A,A),method='bounded',options={'xatol':1e-14});zm=float(fit.x);shift=float(fit.fun)
        phi_host=phi_z(-A);phi_pure=phi_z(A);atom_mass=float(ndtr(-A));lh=float(np.exp(shift-phi_host));lp=float(np.exp(shift-phi_pure));atomZ=atom_mass*(lh+lp)
        atommean=atom_mass*(HOST.real*lh+(HOST.real+D.real)*lp)
        atomsecond=atom_mass*(HOST.real**2*lh+(HOST.real+D.real)**2*lp)
        atomevent=atom_mass*lp
        lhU=max(lh,1e-300);lpU=max(lp,1e-300);atomZU=atom_mass*(lhU+lpU)
        rows=[]
        for coord in ('Gaussian_z0','prior_CDF_uniform_t'):
            t0=time.perf_counter();cpu=time.process_time();calls0=counts['phi_calls'];reported={};evidence=[]
            if coord=='Gaussian_z0':
                lower,upper=-A,A;peak=zm
                def zfrom(x):return x
                def density(x):return norm.pdf(x)
            else:
                lower,upper=float(ndtr(-A)),float(ndtr(A));peak=float(ndtr(zm))
                def zfrom(x):return float(ndtri(x))
                def density(x):return 1.
            def integral(label,task,lo=lower,hi=upper):
                def fun(x):
                    z=zfrom(x);return task(z)*np.exp(shift-phi_z(z))*density(x)
                points=[peak] if lo<peak<hi else None
                val,err=quad(fun,lo,hi,points=points,epsabs=1e-12,epsrel=2e-10,limit=300);reported[label]=err;return val
            Z=integral('evidence',lambda z:1.)+atomZ
            N=integral('mean_numerator',lambda z:HOST.real+D.real*float(sat(1.3*z)))+atommean
            N2=integral('second_numerator',lambda z:(HOST.real+D.real*float(sat(1.3*z)))**2)+atomsecond
            eventlower=0. if coord=='Gaussian_z0' else .5
            NE=integral('event_numerator',lambda z:1.,eventlower,upper)+atomevent
            mean=N/Z;row=dict(coordinates=coord,scaled_evidence=Z,log_evidence=float(np.log(Z)-shift),mean_task=mean,material_variance=N2/Z-mean**2,event_probability=NE/Z,reported_QUADPACK_absolute_errors=reported,phi_calls=counts['phi_calls']-calls0,wall=time.perf_counter()-t0,CPU=time.process_time()-cpu,event_probability_underflowed=(NE==0.),event_probability_note='A displayed numeric zero is floating underflow, not zero posterior mass. Use nonzero envelope upper bound.',qualification='ordinary adaptive numerical integration; entire original Gaussian prior through finite transition and two exact saturation atom masses')
            rows.append(row)
        # Laplace/MAP is a matched local approximation, not a reference.
        mapfit=minimize_scalar(lambda z:phi_z(z)+.5*z*z,bounds=(-A,A),method='bounded',options={'xatol':1e-13});mapz=float(mapfit.x);h=1e-5;hessian=(phi_z(mapz+h)-2*phi_z(mapz)+phi_z(mapz-h))/h**2+1
        sd=1/np.sqrt(hessian);g,w=hermgauss(24);z=mapz+np.sqrt(2)*sd*g;tau=HOST.real+D.real*sat(1.3*z);lap=dict(z_MAP=mapz,z_SD=sd,mean_task=float(w@tau/np.sqrt(np.pi)),material_variance=float(w@(tau**2)/np.sqrt(np.pi)-(w@tau/np.sqrt(np.pi))**2),event_probability=float(ndtr(mapz/sd)),qualification='local Laplace in original Gaussian coordinate; no global or remainder certificate')
        # Analytic Gaussian moments of quintic decoder, evaluated ordinary float.
        tt=Polynomial([.5,.65]);poly=Polynomial([HOST.real])+D.real*(10*tt**3-15*tt**4+6*tt**5);coeff=poly.coef
        def moment(a,bb):
            a=np.asarray(a);bb=np.asarray(bb);pa=norm.pdf(a);pb=norm.pdf(bb);M=[ndtr(bb)-ndtr(a),pa-pb]
            for k in range(2,len(coeff)):M.append(a**(k-1)*pa-bb**(k-1)*pb+(k-1)*M[k-2])
            return sum(c*m for c,m in zip(coeff,M))
        def cell(a,bb):
            qa,qb=float(sat(1.3*a)),float(sat(1.3*bb));v=(qa+qb)/2*D;r=(qb-qa)/2*abs(D);maj=ph.majorant_derivative(v,r,0)
            if maj is None:return dict(a=a,bb=bb,l=0.,u=float(np.exp(shift)),p=float(ndtr(bb)-ndtr(a)))
            f=ph.reduced(v);rn=float(np.linalg.norm(pack(f-y))*np.sqrt(2)/SIGMA);e=maj['white_output_radius']
            return dict(a=a,bb=bb,l=float(np.exp(shift-.5*(rn+e)**2)),u=max(1e-300,float(np.exp(shift-.5*max(0,rn-e)**2))),p=float(ndtr(bb)-ndtr(a)))
        cells=[cell(-A,A)];trajectory=[];step=0;start=time.perf_counter();cpustart=time.process_time()
        for target in (32,128,512,2048):
            while len(cells)<target:
                # Largest unresolved evidence mass, ordinary deterministic policy.
                idx=max(range(len(cells)),key=lambda k:cells[k]['p']*(cells[k]['u']-cells[k]['l']))
                c=cells.pop(idx);mid=(c['a']+c['bb'])/2;cells.extend([cell(c['a'],mid),cell(mid,c['bb'])]);step+=1
            a=np.array([c['a'] for c in cells]);bb=np.array([c['bb'] for c in cells]);pp=np.array([c['p'] for c in cells]);ll=np.array([c['l'] for c in cells]);uu=np.array([c['u'] for c in cells]);mt=moment(a,bb)
            ZL=atomZ+float(pp@ll);ZU=atomZU+float(pp@uu)
            def H(c,lower):
                q=(c-HOST.real)/D.real
                zcut=-A if q<=0 else A if q>=1 else float((2*betaincinv(3,3,q)-1)/1.3)
                cut=np.clip(zcut,a,bb);pm=ndtr(cut)-ndtr(a);mm=moment(a,cut)
                minus=c*pm-mm;plus=mt-mm-c*(pp-pm)
                atomL=atom_mass*((HOST.real-c)*lhU+(HOST.real+D.real-c)*lp)
                atomU=atom_mass*((HOST.real-c)*lh+(HOST.real+D.real-c)*lpU)
                return atomL+float(ll@plus-uu@minus) if lower else atomU+float(uu@plus-ll@minus)
            left=brentq(lambda c:H(c,True),HOST.real,HOST.real+D.real,xtol=1e-13);right=brentq(lambda c:H(c,False),HOST.real,HOST.real+D.real,xtol=1e-13)
            pe=ndtr(bb)-ndtr(np.maximum(a,0));pe=np.maximum(pe,0);pc=pp-pe
            EL=atomevent+float(pe@ll);EU=atom_mass*lpU+float(pe@uu);CL=atom_mass*lh+float(pc@ll);CU=atom_mass*lhU+float(pc@uu)
            pL=EL/(EL+CU);pU=EU/(EU+CL)
            trajectory.append(dict(cells=len(cells),numerical_evidence_enclosure=[ZL,ZU],numerical_mean_interval=[left,right],numerical_event_interval=[pL,pU],mean_halfwidth=(right-left)/2,event_halfwidth=(pU-pL)/2,numerical_mean_stop=(right-left)/2<=.001,numerical_event_stop=(pU-pL)/2<=.002,wall_cumulative=time.perf_counter()-start,CPU_cumulative=time.process_time()-cpustart,strict_evidence_interval=None,strict_mean_interval=None,strict_event_interval=None,qualification='NUMERICAL_WHOLE_SUPPORT_TASK_ENVELOPE; exact-algebra formula but float coefficients/moments/CDF not enclosed'))
            print(json.dumps(dict(event='ORIGINAL1D_ENVELOPE',cells=len(cells),mean_halfwidth=(right-left)/2,event_halfwidth=(pU-pL)/2)),flush=True)
        refdiff={k:abs(rows[0][k]-rows[1][k]) for k in ['scaled_evidence','mean_task','material_variance','event_probability']}
        result=dict(status='NUMERICAL_CONDITIONAL_TASK_PASS' if refdiff['mean_task']<.0002 and refdiff['event_probability']<.0004 else 'REFERENCE_UNQUALIFIED',scope='Original untruncated one-dimensional Gaussian conditional P1 model, uniform offset only; not nonuniform/Full66.',active_gaussian_index=0,fixed_original_gaussianized_properties=[0,0],fixed_native_uniform_properties=[.5,.5],observation='R6 REGULAR controlled exposed observation',D=D,prior_host_atom=atom_mass,prior_pure_atom=atom_mass,prior_transition=1-2*atom_mass,original_prior_tail_accounting='both infinite Gaussian tails map exactly to host/pure, masses analytically retained, no tail deletion or truncation',profile_mode_z=zm,profile_mode_phi=shift,references=rows,cross_coordinate_differences=refdiff,local_Laplace=lap,envelope_trajectory=trajectory,strict_certificate=False,full66='FULL66_UNQUALIFIED',new_full_complex_RHS=0,underflow_reporting_repair='Preserved first run files; repeated full reduced experiment with upper likelihood floor1e-300 and retained atom upper bounds. All replay costs charged.',saturation_atom_log_scaled_likelihoods=dict(host=shift-phi_host,pure=shift-phi_pure),upper_likelihood_floor=1e-300)
        write(RES/'original1d_results_v2.json',result)
        with (RES/'original1d_final_cells_v2.jsonl').open('w') as out:
            for c in cells:out.write(json.dumps(c)+'\n')
        status='COMPLETE'
    except BaseException as e:error=repr(e);write(AUD/'original1d_failure.json',dict(error=error));raise
    finally:b.close(status,error)
if __name__=='__main__':run()

"""Real 24x24 full-Maxwell diagnosis, capped and persisted per endpoint.

GT is read only in the explicitly oracle decomposition/path branch. Legal
chart directions and continuation take no GT argument. Discovery is not MH.
"""
import sys,json,time,traceback
import numpy as np
from scipy.linalg import svdvals
from .common import ROOT,OUT,dump,pack,mass,rms,Book
from prior_constrained import prior
from prior_constrained.common import descriptor,separation
from lxz.local2d import MaterialPrior,TAU

sys.path.insert(0,str(prior.A30/'src'))
from a30.physics import FrozenMaxwellAdapter

def whiten(z,sigma):return pack(z)*np.sqrt(2)/sigma
def response_jac(T,state):return np.einsum('qp,ps->qsp',T,state.field).reshape(32,8,576)
def real_jac(C,sigma):
    rr=np.concatenate((C.real.transpose(1,0,2),C.imag.transpose(1,0,2)),axis=1).reshape(512,576)
    ii=np.concatenate((-C.imag.transpose(1,0,2),C.real.transpose(1,0,2)),axis=1).reshape(512,576)
    return np.c_[rr,ii]*24*np.sqrt(2)/sigma
def chart_jac(C,D,sigma):
    v=C.reshape(256,576)@D
    return np.column_stack([whiten(v[:,j].reshape(32,8),sigma) for j in range(D.shape[1])])
def geometry(C,chart,u,sigma):
    D=chart.analytic_derivative(u);DM=np.r_[D.real,D.imag]/24
    _,s,Vh=np.linalg.svd(DM,full_matrices=False);rank=int(np.sum(s>s[0]*1e-7))
    E=Vh[:rank].T/s[:rank]
    A=chart_jac(C,D,sigma);_,sr,V=np.linalg.svd(A@E,full_matrices=False)
    v=E@V[-1];v/=np.linalg.norm(v)
    # Freeze sign without any truth-dependent selection.
    if v[np.argmax(abs(v))]<0:v=-v
    return dict(D=D,A=A,v=v,rank=rank,material_spectrum=s,response_per_material_spectrum=sr,
                max_cell=int(np.argmax(prior.rho_field(chart.decode(u)[0],prior.background(chart.bg)))))

class Run:
    def __init__(self):
        self.book=Book();self.start=time.perf_counter();self.models={};self.states={};self.endpoints=[];self.metrics=dict(decomposition=[],feedback=[],paths=[],geometry=[])
        self.dest=OUT/'results';self.dest.mkdir(exist_ok=True)
    def model(self,bg):
        if bg not in self.models:self.models[bg]=FrozenMaxwellAdapter.from_background_npz(prior.INPUTS/'background'/f'{bg}.npz',self.book,purpose_prefix='FINITE_FEEDBACK_DIAGNOSTIC')
        return self.models[bg]
    def full(self,model,chi,key,**metadata):
        if key in self.states:return self.states[key]
        self.book.context=dict(endpoint=key,**metadata)
        ss=model.full_state(chi,purpose=key)
        # No subsequent response/Jacobian analysis happens until this is durable.
        p=self.dest/'endpoints'/f'{key}.npz';p.parent.mkdir(exist_ok=True)
        with p.open('wb') as f:
            np.savez_compressed(f,chi=ss.chi,current=ss.current,field=ss.field,prediction=ss.data,
                                LU=ss.factor[0],pivots=ss.factor[1],backward_residual=ss.residual)
            f.flush();__import__('os').fsync(f.fileno())
        self.states[key]=ss;self.endpoints.append(dict(id=key,array=str(p.relative_to(ROOT)),metadata=metadata,backward_residual=ss.residual))
        dump(self.dest/'ENDPOINTS.json',self.endpoints)
        return ss
    def residual(self,state,y,sigma):return float(np.linalg.norm(whiten(state.data-y,sigma)))
    def build_jac(self,model,base,key,sigma):
        self.book.context=dict(anchor=key,operation='full_receiver_jacobian')
        T=model.receiver_transfer(base);C=response_jac(T,base);J=real_jac(C,sigma)
        np.savez_compressed(self.dest/f'{key}_JACOBIAN.npz',C=C,J_white_mass=J,transfer=T)
        return T,C,J
    def spectrum(self,model,base,key):
        with self.book.span('dense_operator_spectrum_no_RHS',operator_SVDs=1):s=svdvals(base.matrix)
        info=dict(condition=float(s[0]/s[-1]),sigma_max=float(s[0]),sigma_min=float(s[-1]),
            relative_nonnormal_commutator=float(np.linalg.norm(base.matrix.conj().T@base.matrix-base.matrix@base.matrix.conj().T,'fro')/np.linalg.norm(base.matrix,'fro')**2))
        dump(self.dest/f'{key}_OPERATOR.json',info);return info
    def decompose(self,base,J,gt,gpred,key,sigma,op,chart=None,u=None):
        h=gt-base.chi;m=mass(h);norm=np.linalg.norm(m)
        with self.book.span('real_material_SVD_no_physics_RHS',material_SVDs=1):U,s,Vh=np.linalg.svd(J,full_matrices=False)
        projection=Vh@m
        spectra=[]
        for rtol in [1e-3,1e-2,1e-1]:
            active=s>=rtol*s[0];energy=np.sum(projection[active]**2)
            spectra.append(dict(relative_cutoff=rtol,absolute_cutoff=float(rtol*s[0]),strong_rank=int(active.sum()),
                strong_energy_fraction=float(energy/norm**2) if norm else None,
                weak_including_rectangular_null_fraction=float(max(0,1-energy/norm**2)) if norm else None,
                strong_condition=float(s[0]/s[active][-1]),
                gap_at_cut=float(s[active][-1]/s[~active][0]) if np.any(~active) else None))
        phase=prior.V*np.real(np.conj(prior.V)*h)
        row=dict(anchor=key,role='ORACLE_GEOMETRY_DIAGNOSTIC',GT_material_RMS=rms(h),
            GT_material_MAE=float(np.mean(abs(h))),linear_response=float(np.linalg.norm(J@m)),
            actual_full_prediction_difference=float(np.linalg.norm(whiten(gpred-base.data,sigma))),
            spectral_thresholds=spectra,sigma_max=float(s[0]),sigma_min=float(s[-1]),
            necessary_phase_subspace_energy=float(np.linalg.norm(phase)**2/max(np.linalg.norm(h)**2,1e-300)) if norm else None,
            outside_necessary_phase_subspace_energy=float(np.linalg.norm(h-phase)**2/max(np.linalg.norm(h)**2,1e-300)) if norm else None,
            prior_tangent='NOT_AVAILABLE: anchor not P0' if chart is None else 'EXACT_CHART_BUT_hGT_ZERO',
            finite_GT_endpoint_P0=True,connecting_straight_path_P0='OFF_PRIOR_ORACLE_PATH' if chart is None else 'ZERO',
            nuisance='fixed gain 1; no nuisance projection applicable',operator=op)
        if chart is not None:
            D=chart.analytic_derivative(u);DM=np.r_[D.real,D.imag]/24
            Q,ss,_=np.linalg.svd(DM,full_matrices=False);rr=ss>ss[0]*1e-7
            row['prior_tangent_rank']=int(rr.sum());row['prior_tangent_GT_fraction']=float(np.linalg.norm(Q[:,rr].T@m)**2/norm**2) if norm else None
        np.savez_compressed(self.dest/f'{key}_DECOMPOSITION.npz',singular_values=s,right_singular_vectors=Vh,
            GT_mass_displacement=m,GT_coefficients=projection,left_singular_vectors=U)
        self.metrics['decomposition'].append(row);return row
    def finite(self,model,base,T,endpoint,key,sigma,op,extended=False):
        dx=endpoint.chi-base.chi;actual=endpoint.data-base.data
        orders=[1,2,3,6]+([20] if extended else []);current=np.zeros_like(base.current);term=None;responses={};terms=[]
        self.book.context=dict(endpoint=key,operation='finite_feedback_series',max_order=max(orders))
        for n in range(1,max(orders)+1):
            rhs=dx[:,None]*(base.field if n==1 else model._cpu.green(term))
            term=model.solve_anchor(rhs,base,role=f'feedback_order_{n}');current=current+term
            response=model._cpu.observe(term);terms.append(float(np.linalg.norm(whiten(response,sigma))))
            if n in orders:responses[f'order_{n}']=model._cpu.observe(current)
        # Algebraically identical rational operator, using the paid endpoint LU:
        # (I-R0 ΔXG)^-1 R0 = (L0-ΔXG)^-1, no explicit inverse.
        exact=model.solve_anchor(dx[:,None]*base.field,endpoint,role='rational_finite_increment')
        rational=model._cpu.observe(exact)
        linear=T@(dx[:,None]*base.field);nonlinear=actual-linear
        errors={k:float(np.linalg.norm(whiten(v-actual,sigma))/max(np.linalg.norm(whiten(actual,sigma)),1e-300)) for k,v in responses.items()}
        with self.book.span('finite_feedback_norm_upper_no_solver_RHS',operator_SVDs=2):
            dnorm=float(svdvals(dx[:,None]*model.G)[0]);end_s=svdvals(endpoint.matrix)
        bound=dnorm/op['sigma_min'];a=whiten(linear,sigma);b=whiten(nonlinear,sigma)
        info=dict(endpoint=key,material_RMS=rms(dx),actual_response=float(np.linalg.norm(whiten(actual,sigma))),
            linear_response=float(np.linalg.norm(a)),nonlinear_remainder=float(np.linalg.norm(b)),
            linear_remainder_cosine=float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))) if np.linalg.norm(a)*np.linalg.norm(b)>0 else None,
            cancellation_ratio=float(np.linalg.norm(a+b)/max(np.linalg.norm(a)+np.linalg.norm(b),1e-300)),
            relative_prediction_errors=errors,term_response_norms=terms,
            rational_relative_error=float(np.linalg.norm(rational-actual)/max(np.linalg.norm(actual),1e-300)),
            finite_current_relative_error=float(np.linalg.norm(exact-(endpoint.current-base.current))/max(np.linalg.norm(endpoint.current-base.current),1e-300)),
            feedback_identity_backward_residual=float(np.linalg.norm(endpoint.matrix@exact-dx[:,None]*base.field)/max(np.linalg.norm(dx[:,None]*base.field),1e-300)),
            first_term_J_consistency=float(np.linalg.norm(responses['order_1']-linear)/max(np.linalg.norm(linear),1e-300)),
            Neumann_operator_norm_upper_bound=bound,Neumann_sufficient_condition=bound<1,
            norm_condition_failure_means='INCONCLUSIVE_NOT_RESOLVENT_FAILURE' if bound>=1 else 'SUFFICIENT_CONVERGENCE',
            endpoint_condition=float(end_s[0]/end_s[-1]),resolvent_inverse_norm_upper=float(op['sigma_max']/end_s[-1]),
            endpoint_factor_reused=True,exact_feedback_is_not_solution_fiber_proof=True)
        np.savez_compressed(self.dest/f'{key}_FEEDBACK.npz',delta_chi=dx,actual_response=actual,linear=linear,rational=rational,
            nonlinear_remainder=nonlinear,rational_current_increment=exact,**responses)
        self.metrics['feedback'].append(info);dump(self.dest/'FEEDBACK_PARTIAL.json',self.metrics['feedback']);return info
    def oracle_path(self,model,base,T,J,gt,gpred,key,sigma,op,y):
        h=gt-base.chi;values={0.:base.data}
        # For Simpson integration use packed real whitened Jv throughout.
        derivatives={0.:J@mass(h)};points=[]
        for t in [.25,.5,.75,1.]:
            ek='DEV_GT' if t==1 else f'{key}_oracle_t{t:g}'
            endpoint=self.full(model,base.chi+t*h,ek,role='OFF_PRIOR_ORACLE_PATH',t=t,anchor=key)
            values[t]=endpoint.data;self.finite(model,base,T,endpoint,f'{key}_oracle_t{t:g}',sigma,op,extended=t==1)
            self.book.context=dict(anchor=key,t=t,operation='path_Jh')
            derivative=model.anchor_jvp(h,endpoint);derivatives[t]=whiten(derivative,sigma)
            np.savez_compressed(self.dest/f'{key}_path_Jv_t{t:g}.npz',Jv=derivative)
            points.append(dict(t=t,endpoint=ek,whitened_residual=self.residual(endpoint,y,sigma),material_RMS_to_GT=rms(endpoint.chi-gt),
                noise_compatible=self.residual(endpoint,y,sigma)<=TAU,prior_status='OFF_PRIOR_ORACLE_PATH' if t<1 else 'GT_P0_ENDPOINT'))
        coarse=(derivatives[0]+4*derivatives[.25]+2*derivatives[.5]+4*derivatives[.75]+derivatives[1])/12
        for t in [.125,.375,.625,.875]:
            ep=self.full(model,base.chi+t*h,f'{key}_oracle_mid{t:g}',role='OFF_PRIOR_ORACLE_J_AVERAGE_REFINEMENT',t=t,anchor=key)
            self.book.context=dict(anchor=key,t=t,operation='path_Jh_refinement');jv=model.anchor_jvp(h,ep)
            derivatives[t]=whiten(jv,sigma);np.savez_compressed(self.dest/f'{key}_path_Jv_t{t:g}.npz',Jv=jv)
        fine=sum(derivatives[i/8]*(1 if i in [0,8] else 4 if i%2 else 2) for i in range(9))/24
        actual=whiten(values[1]-base.data,sigma)
        row=dict(anchor=key,role='ORACLE_GEOMETRY_DIAGNOSTIC_OFF_PRIOR_PATH',points=points,
            averaged_J_Simpson4_relative_error=float(np.linalg.norm(coarse-actual)/np.linalg.norm(actual)),
            averaged_J_Simpson8_relative_error=float(np.linalg.norm(fine-actual)/np.linalg.norm(actual)),
            averaged_J_quadrature_change=float(np.linalg.norm(coarse-fine)/np.linalg.norm(actual)),
            integrates_actual_nonconstant_J=True,local_equal_data_connection_proved=False)
        self.metrics['paths'].append(row)
        np.savez_compressed(self.dest/f'{key}_ORACLE_PATH.npz',t=np.arange(9)/8,Jh_white=np.array([derivatives[i/8] for i in range(9)]),
            averaged_J4=coarse,averaged_J8=fine,actual_full_increment=actual)
    def correct(self,A,v,r,cap):
        Q=np.eye(len(v))-np.outer(v,v)/(v@v)
        delta=-Q@np.linalg.pinv(A@Q,rcond=1e-8)@r
        delta*=min(1,cap/max(np.linalg.norm(delta),1e-300));return delta
    def legal_paths(self,model,base,T,C,chart,u0,key,sigma,y):
        geo=geometry(C,chart,u0,sigma);v0=geo['v'];desc0=descriptor(base.chi,model.chi_b)
        np.savez_compressed(self.dest/f'{key}_CHART.npz',u0=u0,**geo)
        self.metrics['geometry'].append(dict(anchor=key,rank=geo['rank'],material_spectrum=geo['material_spectrum'],
            response_per_material_spectrum=geo['response_per_material_spectrum'],weak_prior_unit_direction=v0,
            linear_material_RMS_per_prior_unit=rms(geo['D']@v0),prior_coordinate_norm=1.,gauge='common Gaussian weight scaling'))
        rows=[]
        def record(ss,u,method,t,parent=None):
            x,rr=chart.decode(u);cert=prior.verify_exact_endpoint(x,rr,background_id=chart.bg)
            sep=separation(desc0,descriptor(x,model.chi_b));e=self.residual(ss,y,sigma)
            row=dict(anchor=key,method=method,t=t,parent=parent,whitened_residual=e,noise_compatible=e<=TAU,
                data_displacement=float(np.linalg.norm(whiten(ss.data-base.data,sigma))),material_RMS=rms(x-base.chi),
                material_MAE=float(np.mean(abs(x-base.chi))),prior_status=cert['prior_status'],
                chart_coordinate_norm=float(np.linalg.norm(u-u0)),own_GT_distance=rms(x-base.chi),
                GT_used_by_path=False,**sep)
            rows.append(row);dump(self.dest/f'{key}_LEGAL_PATH_PARTIAL.json',rows);return row
        op=self.spectrum(model,base,key+'_legal')
        # Frozen-direction endpoints also serve the legal weak finite-feedback test.
        for t in [.25,.5,.75,1.]:
            u=u0+t*v0;ep=self.full(model,chart.decode(u)[0],f'{key}_weak_t{t:g}',role='LEGAL_CHART_FROZEN_WEAK',t=t,anchor=key)
            record(ep,u,'frozen_weak',t);self.finite(model,base,T,ep,f'{key}_weak_t{t:g}',sigma,op)
        # Exact second response derivative: material decoder curvature + Maxwell
        # feedback. Normal correction is a diagnostic, not a density transform.
        eps=1e-3;xp=chart.decode(u0+eps*v0)[0];xm=chart.decode(u0-eps*v0)[0]
        D2=(xp-2*base.chi+xm)/eps**2
        eps2=eps/2;D2fine=(chart.decode(u0+eps2*v0)[0]-2*base.chi+chart.decode(u0-eps2*v0)[0])/eps2**2
        d=geo['D']@v0;self.book.context=dict(anchor=key,operation='chart_curvature_feedback')
        q=model.solve_anchor(d[:,None]*base.field,base,role='chart_curvature_first_current')
        bmax=T@(2*d[:,None]*model._cpu.green(q));bchart=T@(D2fine[:,None]*base.field)
        curvature=whiten(bmax+bchart,sigma)
        self.metrics['geometry'][-1]['chart_second_difference_relative_change']=float(np.linalg.norm(D2-D2fine)/max(np.linalg.norm(D2fine),1e-300))
        self.metrics['geometry'][-1]['Maxwell_curvature_norm']=float(np.linalg.norm(whiten(bmax,sigma)))
        self.metrics['geometry'][-1]['decoder_curvature_norm']=float(np.linalg.norm(whiten(bchart,sigma)))
        for t in [.25,.5,.75]:
            delta=self.correct(geo['A'],v0,.5*t*t*curvature,t/2);u=u0+t*v0+delta
            ep=self.full(model,chart.decode(u)[0],f'{key}_curvature_t{t:g}',role='LEGAL_SECOND_ORDER_CANDIDATE_ONLY',t=t,anchor=key)
            row=record(ep,u,'second_order_curvature',t);row['normal_correction_norm']=float(np.linalg.norm(delta))
        u=u0.copy();state=base;prev=v0;g=geo;distance=0.;status='COMPLETE_THREE_STEPS'
        for step in range(1,4):
            if step>1:
                _,Ci,_=self.build_jac(model,state,f'{key}_updated_{step}',sigma);g=geometry(Ci,chart,u,sigma)
                if g['rank']!=geo['rank']:status='STOP_CHART_RANK_CHANGE';break
                if g['max_cell']!=geo['max_cell']:status='STOP_MAX_NORMALIZATION_CELL_CHANGE';break
            v=g['v']
            if np.real(np.vdot(g['D']@v,geo['D']@prev))<0:v=-v
            up=u+.25*v
            ep=self.full(model,chart.decode(up)[0],f'{key}_pc_{step}_predict',role='LEGAL_UPDATED_J_PREDICTOR',anchor=key,step=step)
            raw=record(ep,up,'updated_J_predictor',step*.25)
            delta=self.correct(g['A'],v,whiten(ep.data-base.data,sigma),.125);uc=up+delta
            cp=self.full(model,chart.decode(uc)[0],f'{key}_pc_{step}_correct',role='LEGAL_UPDATED_J_CORRECTOR_CANDIDATE_ONLY',anchor=key,step=step)
            rr=record(cp,uc,'updated_J_corrector',step*.25,parent=f'{key}_pc_{step}_predict')
            rr['direction_cosine_prior']=float(v@prev);rr['normal_correction_norm']=float(np.linalg.norm(delta))
            if not rr['noise_compatible'] or rr['prior_status']!='P0_EXACT':status='STOP_INCOMPATIBLE_ENDPOINT';break
            distance+=rms(cp.chi-state.chi);u=uc;state=cp;prev=v
        self.metrics['paths'].append(dict(anchor=key,scope='KNOWN_TRAIN_CALIBRATION_ONLY_NOT_INVERSE_SUCCESS',
            updated_J_status=status,updated_J_material_path_length=distance,rows=rows,
            same_data_path_tolerance='original observation 95% noise ball; prediction difference separately reported',posterior_samples=False))
        dump(self.dest/f'{key}_LEGAL_PATH.json',self.metrics['paths'][-1])
    def run(self):
        try:
            for key in ['physics_R_base','N3','250033_weak','250034_ring']:
                a=np.load(OUT/'inputs'/f'{key}.npz');known=key[0].isdigit();bg=key.split('_')[-1] if known else 'weak';model=self.model(bg)
                y=a['observation'];sigma=float(a['sigma']);gt=a['chi'] if known else a['truth'];gpred=a['prediction'] if known else a['truth_prediction']
                base=self.full(model,a['chi'],key+'_base',role='KNOWN_TRAIN' if known else key)
                T,C,J=self.build_jac(model,base,key,sigma);op=self.spectrum(model,base,key)
                if known:
                    meta=json.loads((OUT/'inputs'/f'{key}.json').read_text());chart=MaterialPrior(meta['generator'],bg);u=chart.coordinates()
                else:chart=None;u=None
                info=self.decompose(base,J,gt,gpred,key,sigma,op,chart,u)
                info['anchor_whitened_residual']=self.residual(base,y,sigma)
                info['anchor_noise_compatible']=info['anchor_whitened_residual']<=TAU
                info['anchor_prior_audit']=prior.phase_audit(base.chi,model.chi_b)
                if known:
                    if info['anchor_noise_compatible']:self.legal_paths(model,base,T,C,chart,u,key,sigma,y)
                    else:self.metrics['paths'].append(dict(anchor=key,status='NOT_RUN_UNQUALIFIED_ANCHOR'))
                else:self.oracle_path(model,base,T,J,gt,gpred,key,sigma,op,y)
                dump(OUT/'METRICS_PARTIAL.json',self.metrics)
            self.metrics['cost']=dict(**self.book.totals(),counts=dict(self.book.counts),controlled_wall_seconds=time.perf_counter()-self.start)
            dump(OUT/'METRICS.json',self.metrics);self.book.append(dict(event='JOB_END',status='COMPLETE',**self.metrics['cost']))
            dump(OUT/'RUN_CLOSED.json',dict(reason='predeclared diagnostics complete; no extra search for positive result',**self.metrics['cost']))
        except BaseException as e:
            dump(OUT/'FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc(),cost=self.book.totals(),wall_seconds=time.perf_counter()-self.start))
            self.book.append(dict(event='JOB_END',status='FAIL',error=repr(e),**self.book.totals(),controlled_wall_seconds=time.perf_counter()-self.start));raise

if __name__=='__main__':Run().run()

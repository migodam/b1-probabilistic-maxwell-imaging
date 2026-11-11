"""Public measured/benchmark input adapters. Never opens material truth.

Raw Fresnel files are local official observations, never regenerated or replaced
by simulated data. Defaults locate existing sibling R8 data without copying it.
"""
from pathlib import Path
import numpy as np
from research.research_max_struct_r8.src.fresnel_adapter import prepare
from research.research_max_struct_r8.src.paper_physics import acquisition,make_model


def existing_r8_root():
    repo=Path(__file__).resolve().parents[3]
    choices=[repo/'research/research_max_struct_r8',
             repo.parent/'r8_worktree/research/research_max_struct_r8']
    for root in choices:
        if (root/'data_local/fresnel/dielTM_dec8f.exp').is_file():return root
    raise FileNotFoundError('Official local R8 Fresnel raw data absent; supply raw_path. Never synthesize a replacement.')


def prepare_fresnel(case='dielTM_dec8f', *, raw_path=None, frequency=2.,side=24,
                     width=.16,delta_relative=.15):
    if case not in ('dielTM_dec8f','twodielTM_8f','FoamTwinDielTM'):
        raise ValueError('unsupported actual measurement package')
    path=Path(raw_path) if raw_path is not None else existing_r8_root()/'data_local/fresnel'/f'{case}.exp'
    if not path.is_file():raise FileNotFoundError(path)
    d=prepare(path,frequency=frequency,side=side,width=width,
              illumination='multipole',matrix_free=False)
    m=d['model']
    ports=list(np.moveaxis(m.Ss,2,0))
    metadata=dict(case=case,source_path=str(path),frequency_GHz=frequency,
        side=side,width_m=width,components=1,raw_original_redistributed=False,
        incident='R8 order2 outgoing multipoles fitted only to empty-target measurements',
        calibration=d['calibration'],training_sources=d['train_sources'].tolist(),
        train_mask=d['train_mask'],withheld_mask=d['withheld_mask'],
        training_channels=int(d['Y_train'].size),
        withheld_channels=int(d['withheld_mask'].sum()),
        source_dependent_ports=True,noise_covariance_known=False,
        calibrated_receive_antenna_response=False,
        delta_relative=float(delta_relative),delta_status='PREDECLARED_MODEL_DISCREPANCY_DIAGNOSTIC_NOT_NOISE_CERTIFICATE',
        material_truth_available=False,heldout_used_for_selection=False)
    return dict(model=m,G=m.G,E=m.E,S=ports,Y=d['Y_train'],
        delta=float(delta_relative*np.linalg.norm(d['Y_train'])),components=1,
        lower=0.,upper=20.,metadata=metadata,context=d)


def prepare_austria(public_path=None):
    """Load existing independent-grid classic observation, never its material GT."""
    path=Path(public_path) if public_path else existing_r8_root()/'RESULTS/paper_benchmarks/austria_eps4_400MHz/PUBLIC_DATA.npz'
    with np.load(path,allow_pickle=False) as z:
        side=int(z['inverse_side']);frequency=float(z['frequency']);Y=z['Y'].copy()
    public=acquisition(side,frequency);m=make_model(public)
    return dict(model=m,G=m.G,E=m.E,S=m.S,Y=Y,components=1,lower=0.,upper=20.,
        metadata=dict(source_path=str(path),inverse_side=side,frequency_Hz=frequency,
            benchmark='R8 Austria adapted physical acquisition',forward_side=150,
            truth_access='NONE',raw_observation='existing noisy PUBLIC_DATA',
            material_convention='real chi=epsilon_r-1',noise='historical 10% fixed-Frobenius complex noise',
            qualification='equal-area disk pulse VIE, crossgrid observation; not literal author reproduction'))


from research.research_max_struct_r8.src.fresnel_adapter import SourceReceiverModel


class SourceSpecificPublicModel(SourceReceiverModel):
    """Portable measured package with already selected source illuminations.

    Required public keys: G,E0,S,grid_shape,source_receiver_indices. S is the
    original common receiver matrix; each source uses its actual index list.
    Y and raw measurements are never loaded by this constructor.
    """
    components=1
    def __init__(self,public):
        if 'S_ports' in public:
            ports=np.asarray(public['S_ports'],complex)
            p=dict(public,S=ports[0])
            indices=np.broadcast_to(np.arange(ports.shape[1]),(ports.shape[0],ports.shape[1]))
            super().__init__(p,np.arange(public['E0'].shape[1]),indices)
            self.Ss=np.moveaxis(ports,0,2).copy()
        else:
            super().__init__(public,np.arange(public['E0'].shape[1]),public['source_receiver_indices'])
        self.projector_S=self.projection_ports=list(np.moveaxis(self.Ss,2,0))
        self.metadata=dict(components=1,source_dependent_ports=True,
            objective='0.5 Frobenius relative training residual plus inherited R8 prior',
            noise_covariance_known=False,receive_antenna_calibration_certified=False)


    def data_objective(self,residual,Y):
        return .5*float(np.vdot(residual,residual).real/np.vdot(Y,Y).real)

    def data_terms(self,m,Y):
        return self.data_objective(self.state(m)['Y']-Y,Y),float(self.reg.evaluate(m)[0])


def run_fresnel_auxiliary(model,Y,m0,method,check,audit_callback,audit_every=50,
                          *,cost_callback=None):
    """Unchanged R8 source-specific CIE loop with per-call paid audit hooks.

    A private function-globals copy wraps only the inherited material readout,
    SVD and Counter constructors. Material updates/PR state are unchanged.
    Every audit_every completed material updates adds one original physical
    audit; all original initial/stage-end audits remain. A finally receipt is
    retained as model.r13_auxiliary_cost even when caller BudgetEnd interrupts.
    """
    import time
    from types import FunctionType
    from collections import Counter
    from research.research_max_struct_r8.src import fresnel_solver as inherited
    arm={'cie':'cie_single','dual_cie':'cie_dual','global_dual':'global_dual'}.get(method)
    if arm is None:raise ValueError('No qualified measured source-dependent G19/A+ transplant supplied')
    if int(audit_every)!=audit_every or audit_every<1:raise ValueError('audit_every must be a positive integer')
    Y=np.asarray(Y,complex)
    if not hasattr(model,'Ss') or model.Ss.shape!=(Y.shape[0],model.E.shape[0],model.E.shape[1]):
        raise ValueError('Actual source-specific S_s with matching receive/source dimensions required')
    if Y.shape[1]!=model.E.shape[1] or not np.isfinite(model.Ss).all():
        raise ValueError('Invalid actual source-specific ports/observations')
    started,wall=time.process_time(),time.perf_counter();before=dict(model.cost)
    counts=Counter(material_updates_completed=0,periodic_physical_audits=0,
        original_initial_and_stage_physical_audits=0,source_SVD_completed=0)
    captured=[];result=None;status='RUNNING';error=None
    def physical_audit(m,metadata):
        counts['physical_audit_attempts']+=1
        check();value=model.physical_score(m,Y);state=model.state(m)
        counter=('periodic_physical_audits' if metadata['audit_kind']=='periodic_material_update'
            else 'original_initial_and_stage_physical_audits')
        counts[counter]+=1
        audit_callback(m.copy(),state,value,dict(measured_native_method=arm,
            material_updates=counts['material_updates_completed'],**metadata))
        check();return value
    class Proxy:
        def __getattr__(self,key):return getattr(model,key)
        def physical_score(self,m,Y):
            return physical_audit(m,dict(paid_stage_candidate=True,audit_kind='inherited_initial_or_stage_end'))
    def regularizer(m):
        check();result=model.reg.evaluate(m)[:2];check();return result
    def observed_readout(*args,**kwargs):
        check();answer=inherited.run.__globals__['_baseline_readout'](*args,**kwargs)
        counts['material_updates_completed']+=1
        if counts['material_updates_completed']%audit_every==0:
            physical_audit(answer[0],dict(paid_periodic_candidate=True,audit_kind='periodic_material_update',
                audit_every=int(audit_every)))
        check();return answer
    def observed_svd(*args,**kwargs):
        check();answer=inherited.run.__globals__['svd'](*args,**kwargs)
        counts['source_SVD_completed']+=1;check();return answer
    def captured_counter(*args,**kwargs):
        counter=Counter(*args,**kwargs);captured.append(counter);return counter
    scope=dict(inherited.run.__globals__)
    scope.update(_baseline_readout=observed_readout,svd=observed_svd,Counter=captured_counter)
    run=FunctionType(inherited.run.__code__,scope,inherited.run.__name__,
        inherited.run.__defaults__,inherited.run.__closure__)
    run.__kwdefaults__=inherited.run.__kwdefaults__
    config=dict(host=.1,material_box=[0,20],dpcs_rtol=.001,coefficient_tol=0,
        miniter=50,material_inner_maxiter=60,material_inner_ftol=1e-12,material_inner_gtol=1e-8,
        regularizer=regularizer,stages=[dict(mf=mf,beta_material=beta,
            beta_current=(1. if arm in ('cie_dual','global_dual') else beta),maxiter=it)
            for mf,beta,it in [(3,1.,50),(6,.2,50),(10,.2,50),(12,.2,150)]])
    try:
        check();result=run(Proxy(),Y,m0,arm,config);status='COMPLETE'
        return result
    except BaseException as exc:
        status='INTERRUPTED_OR_FAILED';error=type(exc).__name__;raise
    finally:
        receipt=dict(status=status,error_type=error,method=arm,audit_every=int(audit_every),
            operation_counts=dict(captured[0]) if captured else{},observer_counts=dict(counts),
            physical_model_delta={k:v-before.get(k,0) for k,v in model.cost.items()},
            CPU_seconds=time.process_time()-started,wall_seconds=time.perf_counter()-wall,
            setup_complete=counts['source_SVD_completed']==model.E.shape[1],
            setup_SVD_scope='one original receiver SVD per source; interrupted setup is retained',
            PR_counter_scope='legacy PR_CG_steps increments after readout; interruption inside paid periodic audit may leave one paid update counted only in material_updates_completed',
            physics='actual supplied source-specific ports and analytic receive adjoint',
            implementation='private per-call copy of unchanged R8 fresnel_solver.run; no shared globals or PR restart')
        model.r13_auxiliary_cost=receipt
        if result is not None:result['r13_audit_receipt']=receipt
        if cost_callback is not None:cost_callback(receipt)


def classic_auxiliary_config(public,method,model,Y,out=None):
    """Exact inherited Austria stage/rank contract for the common runner.

    Uses public acquisition metadata only. Caller adds budget observers,
    corrector and the explicitly new G19 material readout when appropriate.
    """
    stages=[]
    for mf,beta,it,tol in zip(public['baseline_mf'],public['baseline_beta'],
            public['baseline_maxiter'],public['baseline_coefficient_tol']):
        stages.append(dict(mf=int(mf),beta_material=float(beta),
            beta_current=float(beta) if method=='cie' else 8.,maxiter=int(it),
            coefficient_tol=float(tol)))
    arm={'cie':'cie_single','dual_cie':'cie_dual','global_dual':'global_dual'}.get(method,'a_plus')
    # The original Austria constructor deliberately bypasses R6's fixed
    # host.1/[.02,20]/mf3,6,10 contract. Preserve that real paper adapter.
    config=dict(arm=arm,stages=stages,material_box=[0.,20.],
        method='cie',sourcewise_cg=method in ('cie','dual_cie'),
        physics_family='equal_area_disk',dpcs_rank=15,coefficient_tol=0.,miniter=2,
        fft_geometry=dict(coords=public['coords'],radius=float(public['radius']),k=float(public['k'])),
        regularizer=lambda m:model.reg.evaluate(m)[:2],aplus_regularizer=model.reg,
        physical_sigma=float(np.linalg.norm(Y)),material_inner_maxiter=60,
        cache_material_value_gradient=True,
        initialize_physical_current=method.startswith('g19'))
    if out is not None:
        config['receipt_path']=str(Path(out)/'AUXILIARY_EXECUTION.json')
    return config

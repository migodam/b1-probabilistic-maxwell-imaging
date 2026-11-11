"""Current-family controls over the *unchanged* R8 FFT/PR alternating loop.

No global monkeypatch, GT, remote backend, new current optimizer, or physical
corrector lives here. Full modes use every native Fourier coefficient from
stage one. ``full_current`` omits only the DPCS projection on those columns.
This is an NFFT-space transplant into G19, not an original NFFT-SOM solver.
"""
from collections import Counter
from pathlib import Path
from types import FunctionType
import time
import numpy as np
from scipy.fft import fftn
from scipy.linalg import svdvals
from research.research_max_struct_r8.src import cie_fft as inherited
from research.research_max_struct_r6.src.fft_tsom_ops import NativeFourierAPCS
from research.research_max_struct_r11.src.material_readouts import prepare_full, full_material_readout

MODES = ('restricted', 'full_apcs', 'full_current')


def _physical_audit_due(config, iteration, audit_every):
    due = config.get('physical_audit_due')
    return bool(due()) if due is not None else bool(audit_every and iteration % audit_every == 0)


class _Proxy:
    def __init__(self, target, **overrides):
        self._target = target
        self.__dict__.update(overrides)

    def __getattr__(self, name):
        return getattr(self._target, name)


def _local_function(function, **overrides):
    """Bind frozen function code to a private globals dictionary per invocation.

    This provides a basis factory seam without editing/reimplementing the PR
    body or changing module/process globals. Imported callees stay original.
    """
    scope = dict(function.__globals__)
    scope.update(overrides)
    result = FunctionType(function.__code__, scope, function.__name__,
                          function.__defaults__, function.__closure__)
    result.__kwdefaults__ = function.__kwdefaults__
    return result


def run_full_current(model, Y, m0, config, mode='restricted', *,
                     budget_check=None, audit_callback=None, audit_every=50,
                     g19=True):
    """Run inherited G19 with one changed admissible current family.

    Required config is the R12/R8 ``arm_config('a_plus', ...)`` dictionary
    including fft_geometry and receipt_path (or a prepared matrix_denominator).
    All stages, beta values, prior, weights, bounds and PR settings are retained;
    only ``mf`` is set to the full mask for the two full modes. The caller owns
    process CPU caps and any subsequent original-data corrector.

    ``budget_check()`` may raise any exception, including BaseException; it is
    checked before/after setup and physical audits and inside the original loop.
    ``audit_callback(m,state,value,metadata)`` receives already-paid physical
    candidates; metadata records phase, current iteration and physical cost.
    No corrector or extra forward runs inside callbacks. An over-budget forward
    is charged but is excluded when the post-forward budget check raises.

    ``g19=False`` reuses the same audit scheduler for inherited CIE/Dual/Global
    or an explicitly supplied A+ readout. It never substitutes G19 material
    readout into CIE. For a strict algorithm control retain the original config.
    On budget exceptions, callbacks already received all eligible incumbents;
    final receipt is still delivered through config['cost_callback'].

    Optional config['physical_audit_due'] is a no-argument callable replacing
    the interval for periodic audits only. It is checked after a completed
    material update; inherited initial/stage-end audits remain unchanged.
    """
    if mode not in MODES:
        raise ValueError('Unknown current family: '+str(mode))
    if int(audit_every) != audit_every or audit_every < 0:
        raise ValueError('audit_every must be a nonnegative integer')
    cfg = dict(config)
    if cfg.get('physical_audit_due') is not None and not callable(cfg['physical_audit_due']):
        raise TypeError('physical_audit_due must be a no-argument callable')
    physics_family = cfg.get('physics_family', 'original_disk_pulse')
    backend = inherited
    physics_overrides = {}
    if physics_family == 'equal_area_disk':
        from research.research_max_struct_r8.src import paper_cie_fft
        from research.research_max_struct_r8.src.paper_physics import EqualAreaGreen, PaperOneHopAbsorption
        if g19 or cfg.get('initialize_physical_current', False):
            # The historical paper loop has no physical-current warm-start
            # block. Retain R8's explicit G19 warm-start loop and swap only
            # compatible physical actions; its global stop rule stays active.
            physics_overrides = dict(DiskPulseToeplitzGreen=EqualAreaGreen,
                                     OneHopAbsorption=PaperOneHopAbsorption)
        else:
            # Mature cold CIE/Dual retains historical per-stage stopping.
            backend = paper_cie_fft
    elif physics_family != 'original_disk_pulse':
        raise ValueError('Unsupported explicit auxiliary physics family: '+str(physics_family))
    cfg['stages'] = [dict(stage) for stage in config['stages']]
    if cfg.get('current_anchor') is not None:
        raise ValueError('Current-anchor translation is a separate historical control')
    if g19:
        if cfg.get('initialize_physical_current', True) is not True:
            raise ValueError('R13 G19 current controls retain the physical-current warm start')
        cfg['initialize_physical_current'] = True
    arm = cfg.get('arm', 'cie_dual')
    if g19 and arm != 'a_plus':
        raise ValueError('G19 requires original a_plus/global PR config; use g19=False for inherited CIE')
    n, ns = model.E.shape
    shape = tuple(map(int, model.shape))
    if mode != 'restricted':
        if len(shape) != 2 or shape[0] != shape[1] or shape[0] % 2:
            raise ValueError('Full native mask currently qualified only for even square 2D grids')
        for stage in cfg['stages']:
            stage['mf'] = shape[0]//2
    wall, cpu = time.perf_counter(), time.process_time()
    model_before = dict(model.cost)
    diagnostics = Counter()
    projection = None
    spaces, audits = [], []
    progress = dict(global_iteration=0, stage=-1, iteration=0)
    last_basis = None
    caller_cost = cfg.get('cost_callback')
    caller_iteration = cfg.get('iteration_callback')
    original_check = cfg.get('budget_check')

    def check():
        if original_check is not None:
            original_check()
        if budget_check is not None and budget_check is not original_check:
            budget_check()

    cfg['budget_check'] = check
    prepared_here = g19 and 'matrix_denominator' not in cfg
    try:
        check()
        if g19:
            if prepared_here:
                if not cfg.get('receipt_path'):
                    raise ValueError('Unprepared G19 requires receipt_path for denominator setup artifacts')
                Path(cfg['receipt_path']).parent.mkdir(parents=True, exist_ok=True)
                cfg = prepare_full(model, cfg)
            cfg['material_readout'] = full_material_readout
        setup_cpu = time.process_time()-cpu
        check()
    except BaseException as exc:
        failed_setup = dict(status='SETUP_FAILED' if isinstance(exc, Exception) else 'SETUP_INTERRUPTED_BY_CALLER',
            error=repr(exc), current_control=mode, truth_access='NONE',
            whole_adapter_CPU_seconds=time.process_time()-cpu,
            whole_adapter_wall_seconds=time.perf_counter()-wall,
            operation_counts={}, physical_model_delta={}, eligible_physical_audits=0)
        if cfg.get('receipt_path'):
            import json
            path = Path(cfg['receipt_path'])
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(failed_setup, indent=2)+'\n')
        if caller_cost is not None:
            caller_cost(failed_setup)
        raise
    reg = cfg.get('regularizer', lambda x: model.reg.evaluate(x)[:2])

    def audit(m, phase):
        check()
        state = model.state(m)
        residual = state['Y']-Y
        value = .5*float(np.vdot(residual, residual).real/np.vdot(Y, Y).real)+float(reg(m)[0])
        check()
        meta = dict(phase=phase, **progress, adapter_CPU_seconds=time.process_time()-cpu,
                    relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)),
                    physical_cost=dict(model.cost), objective=value)
        row = dict(m=np.asarray(m).copy(), objective=value, metadata=meta)
        audits.append(row)
        if audit_callback is not None:
            audit_callback(m.copy(), state, value, meta)
        check()
        return state

    proxy_model = _Proxy(model, state=lambda m: audit(m, 'initial' if not audits else 'stage'))

    def iteration(row):
        progress.update(global_iteration=progress['global_iteration']+1,
                        stage=int(row['stage']), iteration=int(row['iteration']))
        if _physical_audit_due(cfg, progress['global_iteration'], audit_every):
            audit(row['m'], 'periodic')
            diagnostics['periodic_physical_audits'] += 1
        if caller_iteration is not None:
            caller_iteration(row)
        check()

    cfg['iteration_callback'] = iteration

    def basis_factory(grid_shape, mf, Vd, cost):
        nonlocal last_basis
        check()
        actual_vd = np.zeros((n, 0), complex) if mode == 'full_current' else Vd
        basis = NativeFourierAPCS(grid_shape, mf, actual_vd, cost)
        # Actual independent complex spatial rank, rather than Fourier column
        # count. For P F_I: rank=k-r+rank(F_complement^* Vd).
        if mode == 'full_current':
            rank = n
        elif basis.k == n:
            rank = n-Vd.shape[1]
        else:
            spectrum = fftn(Vd.reshape(shape+(Vd.shape[1],)), axes=(0, 1), norm='ortho').reshape(n, -1)
            selected = np.ravel_multi_index(basis.frequencies.T, shape)
            complement = np.ones(n, bool)
            complement[selected] = False
            singular = svdvals(spectrum[complement], check_finite=False)
            rank = basis.k-Vd.shape[1]+int(np.count_nonzero(singular > 1e-12))
            diagnostics['current_rank_diagnostic_FFT_2D_transforms'] += Vd.shape[1]
            diagnostics['current_rank_diagnostic_small_SVD'] += 1
        info = dict(stage=len(spaces), inherited_mf=int(config['stages'][len(spaces)]['mf']),
                    actual_mf=int(mf), coefficient_columns=basis.k,
                    independent_complex_spatial_rank=int(rank),
                    actual_real_source_current_DOF=int(2*rank*ns), material_real_DOF=n,
                    receiver_DPCS_rank=int(Vd.shape[1]), DPCS_locked=mode != 'full_current',
                    rank_rule='rank(P F_I)=k-r+rank(F_complement* Vd), tol 1e-12; full masks exact n-r')
        spaces.append(info)
        last_basis = (basis, Vd)
        check()
        return basis

    def observed_lstsq(A, b, *args, **kwargs):
        nonlocal projection
        check()
        result = np.linalg.lstsq(A, b, *args, **kwargs)
        if projection is None:
            state = getattr(model, '_current', None)
            if state is None or not np.array_equal(state['m'], m0):
                raise RuntimeError('Warm projection must reuse the initial paid physical state')
            jp = state['j']
            remainder = A@result[0]-b
            projected = jp+remainder
            Q = model.E+model.G@projected
            diagnostics['initial_projection_diagnostic_dense_G_RHS'] += ns
            port_change = model.S@remainder
            diagnostics['initial_projection_diagnostic_dense_S_RHS'] += ns
            _, Vd = last_basis
            diagnostics['initial_projection_diagnostic_DPCS_RHS'] += ns
            projection = dict(relative_current_loss=float(np.linalg.norm(remainder)/max(np.linalg.norm(jp), 1e-300)),
                constitutive_defect=float(np.linalg.norm(projected-m0[:, None]*Q)/max(np.linalg.norm(projected), 1e-300)),
                relative_port_change=float(np.linalg.norm(port_change)/max(np.linalg.norm(Y), 1e-300)),
                DPCS_projection_loss=float(np.linalg.norm(Vd@(Vd.conj().T@remainder))/max(np.linalg.norm(jp), 1e-300)),
                LS_rank=int(result[2]), singular_values=result[3].tolist(), truth_access='NONE',
                semantics='Original initial physical current projected into this affine current family; LS unchanged')
            spaces[0]['warm_LS_observed_rank'] = int(result[2])
        check()
        return result

    numpy_proxy = _Proxy(np, linalg=_Proxy(np.linalg, lstsq=observed_lstsq))
    loop = _local_function(backend._run_cie, NativeFourierAPCS=basis_factory, np=numpy_proxy,
                           **physics_overrides)
    runner = _local_function(backend.run_cie_fft, _run_cie=loop)

    def receipt_callback(receipt):
        # Original wrapper's physical delta remains exact; setup diagnostics are
        # separately charged and whole-adapter timing includes denominator setup.
        if receipt['status'] == 'RUNNING':
            receipt['status'] = 'INTERRUPTED_BY_CALLER'
        receipt.update(current_control=mode, physics_family=physics_family,
            inherited_loop_entrypoint=backend.__name__+'._run_cie',
            auxiliary_Green='equal-area disk pulse FFT at unchanged supplied radius' if physics_family == 'equal_area_disk' else receipt['auxiliary_Green'],
            effective_stage_coefficient_tolerances=[float(s.get('coefficient_tol',cfg.get('coefficient_tol',1e-3)))
                if backend is not inherited else float(cfg.get('coefficient_tol',1e-3)) for s in cfg['stages']],
            stopping_policy='historical paper per-stage thresholds' if backend is not inherited else 'R8 global threshold; supplied per-stage threshold fields are inactive',
            algorithm='inherited R8 FFT/PR; G19 readout' if g19 else 'inherited '+arm,
            qualification='NFFT full-Fourier SPACE transplant; not original NFFT-SOM reproduction',
            Fourier='full native Fourier, no DPCS projection' if mode == 'full_current' else 'native Fourier with original DPCS projection',
            current_spaces=spaces, initial_current_projection=projection,
            adapter_diagnostic_cost=dict(diagnostics), denominator_setup_CPU_seconds=setup_cpu,
            denominator_setup_cost=dict(G_SVD=int(prepared_here), S_SVD=int(prepared_here),
                projected_basis_SVD=int(prepared_here),
                Green_basis_and_Galerkin_setup=int(prepared_here),
                dense_Green_complex_RHS=int(cfg['matrix_denominator'].rank) if prepared_here else 0),
            whole_adapter_CPU_seconds=time.process_time()-cpu,
            whole_adapter_wall_seconds=time.perf_counter()-wall,
            audit_every_material_updates=int(audit_every), eligible_physical_audits=len(audits),
            periodic_physical_audit_policy='caller physical_audit_due after completed material updates' if cfg.get('physical_audit_due') is not None else 'fixed material-update interval',
            audit_every_is_fallback_only=cfg.get('physical_audit_due') is not None,
            audit_schedule='inherited initial/stage ends plus periodic; no PR restart',
            initialize_physical_current=bool(cfg.get('initialize_physical_current', False)),
            best_audit=min((a['objective'] for a in audits), default=None),
            whole_adapter_physics_delta={k:model.cost.get(k, 0)-model_before.get(k, 0) for k in set(model.cost)|set(model_before)},
            truth_access='NONE')
        # Outer frozen wrapper writes its receipt before invoking this callback.
        # Write augmented fields to the same explicit new-run destination.
        if cfg.get('receipt_path'):
            import json
            Path(cfg['receipt_path']).write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
        if caller_cost is not None:
            caller_cost(receipt)

    cfg['cost_callback'] = receipt_callback
    result = runner(proxy_model, Y, m0, cfg)
    result['physical_audit_candidates'] = audits
    result['current_spaces'] = spaces
    result['initial_current_projection'] = projection
    if audits:
        best = min(audits, key=lambda row: row['objective'])
        result['m'] = best['m'].copy()
        result['selection'] = dict(policy='lowest original physical Phi among all eligible paid audits',
                                   objective=best['objective'], metadata=best['metadata'])
    return result


def run_audited_cie(model, Y, m0, config, **kwargs):
    """Same paid physical-audit scheduler, preserving inherited material update."""
    return run_full_current(model, Y, m0, config, g19=False, **kwargs)

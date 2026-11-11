"""Qualified NFFT Eq(3)–(6)/joint-L-BFGS adaptation, not author reproduction.

Complete unitary Fourier correction J=JD+F alpha is represented by the spatial
correction Z=F alpha. This orthogonal change of coordinates preserves the global
Euclidean joint metric; no current/material scaling or preconditioning is added.
The 2018 conference paper does not specify its CG/material/initialization rules
fully. Real material, common box/prior and joint L-BFGS-B are explicit adaptations.
"""
from collections import Counter
from pathlib import Path
import json
import time
import numpy as np
from scipy.linalg import svd
from scipy.optimize import minimize


def _physical_audit_due(config, iteration, audit_every):
    due = config.get('physical_audit_due')
    return bool(due()) if due is not None else bool(iteration % audit_every == 0)


def joint_value_gradient(G, E, S, Y, JD, Z, m, wd, ws, regularizer):
    """Analytic real/imag-current and shared-real-material LSIE derivatives."""
    J = JD+Z
    Q = E+G@J
    D = S@J-Y
    e = J-m[:, None]*Q
    ew = e*ws[None, :]
    data = float(np.sum(abs(D)**2*wd[None, :]))
    state = float(np.sum(abs(e)**2*ws[None, :]))
    prior, gp = regularizer(m)[:2]
    current_gradient = 2*(S.conj().T@(D*wd[None, :])+ew-G.conj().T@(m[:, None]*ew))
    material_gradient = -2*np.real(np.sum(ew.conj()*Q, axis=1))+np.asarray(gp, float)
    return data+state+float(prior), current_gradient, material_gradient, dict(data=data, state=state, prior=float(prior))


def run_nfft(model, Y, m0, config, check, audit_callback, audit_every=50):
    """Run the declared LSIE equations control on public scalar G/E/S/Y.

    config: material_box, dpcs_rank OR dpcs_rtol, regularizer(m)->(R,g),
    maxiter (default 3000), ftol (1e-12), gtol (1e-8), maxls (20), optional
    receipt_path and cost_callback(receipt). No half factors: joint objective
    is sum_s ||SJ_s-Y_s||²/||Y_s||² + ||J_s-mQ_s||²/||JD_s||² + R(m).
    It is not the physical candidate-selection objective, which stays the
    supplied model's original physical Phi. No GT or physical adjoint is used.

    check() owns the absolute process budget and may raise BaseException.
    audit_callback(m,paid_state,physical_Phi,metadata) receives eligible paid
    initial, periodic accepted-update, and terminal materials. No PR restart
    or outer corrector is introduced. Finally receipts survive interruption.

    Optional config['physical_audit_due'] is a no-argument callable replacing
    the interval for periodic audits only. It runs after an accepted joint
    update; native initial/terminal audits remain unchanged.
    """
    if int(audit_every) != audit_every or audit_every < 1:
        raise ValueError('Positive integer audit interval required')
    wall, cpu = time.perf_counter(), time.process_time()
    before = dict(model.cost)
    cost = Counter()
    status, error, result, opt = 'RUNNING', None, None, None
    audits, joint_history = [], []
    receiver_info = None
    last_audited = None
    cfg = dict(config)
    try:
        if cfg.get('physical_audit_due') is not None and not callable(cfg['physical_audit_due']):
            raise TypeError('physical_audit_due must be a no-argument callable')
        check()
        G, E, S = [np.asarray(x, complex) for x in (model.G, model.E, model.S)]
        Y, m0 = np.asarray(Y, complex), np.asarray(m0, float)
        n, ns = E.shape
        if G.shape != (n, n) or S.ndim != 2 or S.shape[1] != n or Y.shape != (S.shape[0], ns):
            raise ValueError('Expected common scalar G/E/S/Y acquisition')
        if m0.shape != (n,) or hasattr(model, 'Ss') or getattr(model, 'components', 1) != 1:
            raise ValueError('This qualified control supports real scalar material and common receiver ports only')
        lo, hi = [np.broadcast_to(np.asarray(v, float), (n,)).copy() for v in cfg.get('material_box', (.02, 20.))]
        if np.any(lo > hi) or np.any(m0 < lo) or np.any(m0 > hi) or not np.isfinite(m0).all():
            raise ValueError('Invalid material box/initializer')
        reg = cfg.get('regularizer', lambda m: model.reg.evaluate(m)[:2])
        yd = np.sum(abs(Y)**2, axis=0)
        if np.any(yd <= 0) or not np.isfinite(yd).all():
            raise ValueError('Nonzero finite observations required by source normalization')
        u, singular, vh = svd(S, full_matrices=False, check_finite=False)
        cost['receiver_SVD'] += 1
        check()
        if ('dpcs_rank' in cfg) == ('dpcs_rtol' in cfg):
            raise ValueError('Specify exactly one frozen public DPCS rank rule')
        rank = int(cfg['dpcs_rank']) if 'dpcs_rank' in cfg else int(np.count_nonzero(singular > float(cfg['dpcs_rtol'])*singular[0]))
        if not 1 <= rank <= len(singular) or singular[rank-1] <= np.finfo(float).eps*max(S.shape)*singular[0]:
            raise ValueError('A supported nonzero DPCS rank is required for the Eq6 determined-current normalization')
        Vd = vh[:rank].conj().T
        JD = Vd@((u[:, :rank].conj().T@Y)/singular[:rank, None])
        jd_energy = np.sum(abs(JD)**2, axis=0)
        if np.any(jd_energy <= 0) or not np.isfinite(jd_energy).all():
            raise ValueError('Eq6 state normalization requires nonzero finite determined current per source; no epsilon rescue')
        wd, ws = 1/yd, 1/jd_energy
        receiver_info = dict(rank=rank, source_determined_current_energy=jd_energy.tolist(),
            data_weights=wd.tolist(), state_weights=ws.tolist(),
            rank_source='caller-frozen public rank rule; no truth selection',
            Fourier='complete unitary Fourier space represented by its spatial correction',
            complex_current_DOF=n*ns, real_current_DOF=2*n*ns, material_real_DOF=n)
        size = n*ns

        def pack(Z, m):
            return np.r_[Z.real.ravel(), Z.imag.ravel(), m]

        def unpack(x):
            return (x[:size]+1j*x[size:2*size]).reshape(n, ns), x[2*size:]

        def evaluate(x):
            check()
            Z, m = unpack(x)
            value, current_gradient, material_gradient, parts = joint_value_gradient(G, E, S, Y, JD, Z, m, wd, ws, reg)
            gradient = pack(current_gradient, material_gradient)
            cost['joint_objective_gradient_evaluations'] += 1
            cost['dense_G_forward_calls'] += 1
            cost['dense_G_forward_complex_RHS'] += ns
            cost['dense_G_adjoint_calls'] += 1
            cost['dense_G_adjoint_complex_RHS'] += ns
            cost['receiver_forward_calls'] += 1
            cost['receiver_forward_complex_RHS'] += ns
            cost['receiver_adjoint_calls'] += 1
            cost['receiver_adjoint_complex_RHS'] += ns
            cost['joint_material_prior_evaluations'] += 1
            if not np.isfinite(value) or not np.isfinite(gradient).all():
                raise FloatingPointError('Nonfinite qualified NFFT equations objective/gradient')
            joint_history.append(dict(evaluation=cost['joint_objective_gradient_evaluations'],
                CPU_seconds=time.process_time()-cpu, objective=value, **parts,
                gradient_inf=float(np.linalg.norm(gradient, np.inf)),
                eligibility_scope='joint trial evaluation only; never a physical selection candidate'))
            check()
            return value, gradient

        def audit(m, reason):
            nonlocal last_audited
            check()
            paid = model.state(m)
            cost['physical_audit_calls'] += 1
            residual = paid['Y']-Y
            data = model.data_objective(residual, Y) if hasattr(model, 'data_objective') else .5*float(np.vdot(residual, residual).real/np.vdot(Y, Y).real)
            physical_value = float(data)+float(reg(m)[0])
            cost['physical_selection_prior_evaluations'] += 1
            try:
                check()
            except BaseException:
                cost['overrun_physical_audits_excluded'] += 1
                raise
            metadata = dict(method='NFFT_EQUATIONS_JOINT_LBFGSB_ADAPTATION',
                qualification='NOT_ORIGINAL_NFFT_SOM', reason=reason,
                accepted_joint_updates=cost['accepted_joint_updates'],
                objective=physical_value, relative_residual=float(np.linalg.norm(residual)/np.linalg.norm(Y)),
                adapter_CPU_seconds=time.process_time()-cpu, physical_cost=dict(model.cost))
            audits.append(dict(m=m.copy(), objective=physical_value, metadata=metadata))
            last_audited = m.copy()
            if audit_callback is not None:
                audit_callback(m.copy(), paid, physical_value, metadata)
            check()

        def accepted(x):
            cost['accepted_joint_updates'] += 1
            if _physical_audit_due(cfg, cost['accepted_joint_updates'], audit_every):
                _, m = unpack(x)
                audit(m, 'periodic_accepted_joint_update')
            check()

        audit(m0, 'common_initial_physical_material')
        x0 = pack(np.zeros_like(JD), m0)
        bounds = [(None, None)]*(2*size)+list(zip(lo, hi))
        options = dict(maxiter=int(cfg.get('maxiter', 3000)), ftol=float(cfg.get('ftol', 1e-12)),
            gtol=float(cfg.get('gtol', 1e-8)), maxls=int(cfg.get('maxls', 20)))
        if 'maxfun' in cfg:
            options['maxfun'] = int(cfg['maxfun'])
        opt = minimize(evaluate, x0, jac=True, method='L-BFGS-B', bounds=bounds,
                       callback=accepted, options=options)
        final_Z, final_m = unpack(opt.x)
        if last_audited is None or not np.array_equal(final_m, last_audited):
            audit(final_m, 'optimizer_terminal_material')
        best = min(audits, key=lambda row: row['objective'])
        result = dict(status='OPTIMIZER_TERMINATED', m=best['m'].copy(),
            physical_audit_candidates=audits, joint_history=joint_history,
            final_joint_m=final_m.copy(), final_joint_J=JD+final_Z,
            optimizer=dict(message=str(opt.message), success=bool(opt.success),
                accepted_updates=int(opt.nit), evaluations=int(opt.nfev),
                final_joint_objective=float(opt.fun), options=options),
            selection='lowest original physical Phi among eligible paid audits only',
            qualification='NOT_ORIGINAL_NFFT_SOM')
        status = 'COMPLETE'
        return result
    except BaseException as exc:
        status = 'FAILED' if isinstance(exc, Exception) else 'INTERRUPTED_BY_CALLER'
        error = dict(type=type(exc).__name__, message=str(exc))
        raise
    finally:
        receipt = dict(status=status, error=error, truth_access='NONE',
            method='NFFT_EQUATIONS_JOINT_LBFGSB_ADAPTATION', qualification='NOT_ORIGINAL_NFFT_SOM',
            primary_source='https://usncursi.org/archive/aps-ursi/2018/papers/1927.pdf Eq3–6',
            objective='sum_s wd_s ||S J_s-Y_s||² + ws_s ||J_s-m(E_s+GJ_s)||² + R(m); no halves',
            normalization='wd=1/||Y_s||²; ws=1/||JD_s||², not incident-field energy',
            explicit_adaptations=['real shared material', 'caller common material box/prior',
                'joint L-BFGS-B replaces unspecified original CG details',
                'data denominator interpreted as observed scattered-data energy'],
            initialization='spatial correction zero, hence J=JD; material exactly supplied m0',
            coordinate_metric='unscaled global Euclidean; unitary full Fourier equivalent to spatial correction',
            receiver_setup=receiver_info, operation_counts=dict(cost),
            physical_model_delta={k:model.cost.get(k, 0)-before.get(k, 0) for k in set(before)|set(model.cost)},
            eligible_physical_audits=len(audits), best_physical_Phi=min((a['objective'] for a in audits), default=None),
            periodic_physical_audit_policy='caller physical_audit_due after accepted joint updates' if cfg.get('physical_audit_due') is not None else 'fixed accepted-update interval',
            audit_every_fallback=int(audit_every),
            whole_adapter_CPU_seconds=time.process_time()-cpu,
            whole_adapter_wall_seconds=time.perf_counter()-wall,
            timer_scope='adapter setup/current/material/paid physical audits; caller process timer additionally includes imports, final serialization and exit')
        model.r13_nfft_cost = receipt
        if result is not None:
            result['execution_receipt'] = receipt
        if cfg.get('receipt_path'):
            path = Path(cfg['receipt_path']);path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
        if cfg.get('cost_callback') is not None:
            cfg['cost_callback'](receipt)

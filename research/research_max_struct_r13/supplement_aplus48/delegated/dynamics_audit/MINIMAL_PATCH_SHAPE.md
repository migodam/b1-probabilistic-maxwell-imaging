# Minimal recommended patch shape — diagnostic only, not applied

This is an implementation sketch for the parent, not executable production code. It does not authorize editing frozen R8/R7 or C2 files. Primary A+48 gets **none** of these changes.

## 1. Existing-hook-only first option

A wrapper supplied as config.material_readout captures actual fixed-current inputs and the original return, calling the original once:

```python
def observed_readout(blocks, J, Q, m, beta, Os, Ls, rho, bounds,
                     reg, config, cost):
    # Copying is diagnostic work. Do not call a physics/objective routine here.
    token = trace.begin_readout(J=J, Q=Q, m_before=m, beta=beta,
                               Os=Os, Ls=Ls, rho=rho, bounds=bounds)
    new_m, inner = original_readout(
        blocks, J, Q, m, beta, Os, Ls, rho, bounds, reg, config, cost)
    trace.end_readout(token, m_after=new_m, native_inner=inner)
    return new_m, inner                 # original objects, original arithmetic
```

Trace methods own independent copies and never mutate inputs. If original_readout raises, no completed end record is invented. Existing iteration_callback can additionally save m/J/A; stage_callback saves m/J/Q/A. This option lacks PR arrays and native history.

## 2. Complete outer-update hook in new diagnostic loop copy

Inside an otherwise exact copy of R8 _run_cie, after the existing history.append at source:306–319, before its existing iteration_callback:

```python
# Only new diagnostic read/copy operations; no new contractions or physics.
observer = config.get("diagnostic_observer")
if observer is not None:
    observer.outer_update(
        stage=stage_index, iteration=it + 1,
        global_iteration=iterations + 1,
        m_before=m_before_readout, m_after=m,
        J=J, Q=Q, A=A, gradient=grad, direction=direction,
        alpha=alpha, PR_plus=pr, restarted=restarted,
        native_history=history[-1], operation_counts=cost,
        phase="material_complete_before_native_callback_and_step_counter")
# Existing callback and all later statements remain in exactly the old order.
```

The observer copies **existing** scalar values and arrays; it must not compute new state_value, regularizer, projection or gradient calls. Native history already contains their frozen computed values.

At setup and stage entry add corresponding copy-only records for actual JD/Vd/metrics/bounds/rho/frequencies/betas/A/J/Q. Save old_grad/old_direction only when a checkpoint is taken before their assignment; at this post-update hook they equal grad/direction.

Use private per-invocation binding (as in full_current._local_function), with a backend proxy pointing to this new diagnostic loop. Leave frozen shared module globals and runner source untouched. Do not use inspect/trace/sys.settrace on every line or a new orchestration framework.

## 3. DR copy hook in a new diagnostic runner copy

Preserve the original material and reflected-affine call and addition expressions. Keep pre-step arrays alive only when tracing (the original X/Q updates rebind):

```python
pre_x, pre_q = x, q
cx, cq, cm = P.material(x, q)           # existing call, do not repeat
ax, aq = P.affine_ball(2*cx-x, 2*cq-q) # existing call, do not repeat
x, q = x+ax-cx, q+aq-cq                # existing expression/order
if observer is not None:
    observer.dr_update(
        index=k, pre_X=pre_x, pre_Q=pre_q,
        cX=cx, cQ=cq, material=cm,
        reflected_A_X=ax, reflected_A_Q=aq,
        post_X=x, post_Q=q,
        A_info=P.last_A_info, projector_cost=P.costs)
# Existing periodic physical audit stays after this point, unchanged.
```

No “last shadow” is recomputed after the loop. Copy last_A_info before any later A projection. Distinguish cm(k) and post_iterator(k+1) in output names.

## 4. Diagnostic owner finally and receipts

```python
try:
    run_diagnostic_copy(...)
finally:
    trace.flush_completed_records()
    trace.write_receipt(
        role="POST_SCORE_OFFLINE_DIAGNOSTIC",
        primary_selection_changed=False,
        completed_records_only=True,
        stop_exception=original_stop,
        # measured whole CPU/wall, nested trace CPU, bytes and missing coverage
    )
```

Finally only flushes observations; it never continues a step or creates/audits a candidate. All failures and partial files are retained. If flushing fails, report failed instrumentation and paid cost, without discarding the original algorithm exception.

Optional inner-color/root hooks should be added only for a fixed parent question; emit the arrays already computed immediately before the existing helper return. Do not recompute roots or enable inherited_gtol to obtain diagnostics.

No patch in this document has been applied, syntax-tested as code, or executed. Parent owns the concrete implementation and review.


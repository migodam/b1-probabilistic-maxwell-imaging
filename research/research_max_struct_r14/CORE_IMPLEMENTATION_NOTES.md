# R14 stateful finite-PR core

Implemented only `src/dualmetric_core.py`; historical R8/R7 and physical-model modules are unchanged. No experiments or truth scoring are run by import or construction. The parent owns mathematical interpretation, physical audits, candidate selection, and verification.

## Stable API

- `StatefulDualMetric(model, Y, m0, config, pair='01')`.
- Public live state: `m`, `A`/`coefficients`, `J`, `Q`, `cost`/`counts` (`Counter`), `iterations`/`global_iteration`, `stage`, `iteration`, `old_grad`/`old_gradient`, `old_direction`, `stage_pending`, `done`, `status`, `history`, `stage_states`.
- `step()` returns exactly one completed current+material update. Returned row includes all inherited scalar history fields and copied `m`, `m_before`, `J`, `Q`, `coefficients`, current `gradient`/`old_grad`, `direction`, `previous_gradient`, `previous_direction`, `alpha`, `PR`, `stage_end`, `stage_stop`, and `operation_counts`. Only scalar native rows are retained in `history`.
- `state_dict()`, `checkpoint(path)`, `from_checkpoint(model,Y,config,path,pair=None,reset_memory=False)`, `from_state_dict(...)`, `clone(...)`, `switch_pair(pair,reset_memory=False)`.
- Optional paid, cached-state diagnostics: `value_current(m,J,Q,include_prior=False)`, `value_material(m,J,Q)`, and full `material_value_gradient(m,J,Q)`. These do not solve Maxwell or run extra FFT propagation. Explicit `J` requires explicit cached `Q`.

## Inherited operations and pair semantics

The core retains native `NativeFourierAPCS` and `DiskPulseToeplitzGreen`, SVD-derived DPCS, fixed source metrics and Cholesky factors, global/sourcewise real PR+ axis, original non-descent restart rule, exact fixed-material quadratic line search, native squared-coefficient stopping, nested-mask zero-padding transport, and per-stage PR reset. Pair `01` is the original cold-DPCS A+-Arc current/material pairing when supplied the original R7 final flags.

`beta0` is stage `beta_current`, falling back exactly as inherited to `beta_material`/`beta`; `beta1` is stage `beta_material`/`beta`. Both scalar blocks use `beta0`; both onehop blocks use `beta1`. P1 forward calls the actual `OneHopAbsorption.apply`; its adjoint is `conj(p)*(z+Goff* (conj(a)*z))` with inherited block actions. Material P1 calls the original R7 `mobius_material_readout`, including the supplied Arc/MM flags and full original material dependence. Scalar material calls the inherited `_baseline_readout` with `beta0` and the original sourcewise/full-metric whitening branch.

`onehop_self_only=True` zeros only the newly owned `blocks.goff` copy after original extraction. It retains the self denominator and all other state/flags. This is the parent's single self-versus-offdiagonal control.

## Checkpoint and stage boundaries

Checkpoints contain material, coefficients, J/Q, old PR gradient/direction, old frequency coefficient mapping, complete stage/iteration/global position, stage-end pending flag, complete native history and stage states, DPCS/receiver basis and singular values, JD/GJD/SJD, source metrics/factors, frozen numerical configuration flags, and cached finite-Green spectrum and onehop blocks. Restore rebuilds existing action objects from those cached arrays; no SVD, FFT setup, propagation, warm projection, or physical solve is repeated. Numerical flags in the checkpoint take precedence over corresponding new configuration fields; runtime callbacks, regularizer, and `aplus_regularizer` objects must be supplied again. Caller process timing includes serialization and restoration work.

Stage-ending updates finish transport bookkeeping before callbacks, so a callback budget exception leaves a complete restartable state. The next `step()` enters the next mask and resets PR. `stage` still identifies the just-completed stage in a stage-end checkpoint. A complete last stage sets `done=True`; further `step()` calls raise `StopIteration`. `max_iterations` is a resumable caller budget: exhaustion sets `status='ITERATION_BUDGET'` and raises `StopIteration` without falsifying a completed stage event.

`budget_check` runs before/after setup, before/after stage entry, and at step entry/before work/after completion. It may raise `BaseException`. Original material readouts remain atomic as in R8. Optional `iteration_callback` receives a shallow row dictionary whose state arrays are already detached from live state; `deepcopy_iteration_callback=True` also detaches every nested item. Stage callbacks receive detached stage state. The core does no automatic initial, periodic, or stage-end physical audit. Cold DPCS is default; only explicit `initialize_physical_current=True` permits the historical initial physical solve if a paid matching state/current has not been supplied.

The checkpoint is a real inspectable compressed NPZ, loaded with `allow_pickle=False`. Principal arrays use explicit names (`m`, `A`, `J`, `Q`, `old_grad`, `old_direction`, `Vd`, `JD`, `Od`, `Os`); nested arrays also have readable names. A plain Unicode `__state_tree_json__` member records safe typed list/tuple/dict structure, including tuple dictionary keys, complex scalars and array references. No Python pickle or object arrays are accepted. Supplied observations, incident field, receiver operator, grid, native mask and coefficient shape are checked on restore. The caller must supply the same physical Green model and original prior. All four modes extract the same original onehop block cache at setup, permitting an objective-only same-state switch; this common setup cost is recorded even for scalar-only `00`.

## Verification boundary

Syntax compilation passed. Parent verification should compare `01` with original R8+R7 over full trajectories, test checkpoint continuation inside stages and exactly at mask boundaries, check fixed-material P1 adjoint and gradients, and verify clone isolation. No claim of empirical parity or successful physical recovery is made by this implementation note.

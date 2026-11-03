# Independent passive-adapter derivative review

Read `src/passive_continuation.py` and `experiments/passive_path.py`.
At fixed gamma, c=m+i gamma and dc/dm=1 for real material parameter m.
From `(I-diag(c)G)j=diag(c)E`, a real variation gives
`delta j=L_c^-1 diag(Q)delta m`. The inherited source-adjoint expression
`Re sum_s conj(a_s)Q_s`, with `L_c* a=S* residual / ||Y||²`, therefore remains
the correct real parameter gradient. The unchanged prior acts on real m.
No missing derivative of gamma exists because each stage holds gamma fixed.

The adapter assigns the new factor/state into the parent physics object and
retains an explicit `state['lu']`. The m-only state cache is valid in the
executed contract because a new model instance is constructed for each gamma.
Mutating `model.gamma` after a state is cached is unsupported and would need
explicit cache invalidation. This is a scope caveat, not an observed run fault.

Within each passive stage the corrector selects by that passive model's data
objective plus the same real-m prior. Only the final gamma0 comparison uses
the original physical objective. Calling every intermediate passive evaluation
an original-physics validation would be inaccurate; the final original
selection contract is explicit in the path runner.

The single requested 4x4 gamma1 directional check is recorded in
`audit/PASSIVE_GRADIENT_WIRING.json`: analytic0.001419016821450081,
central0.0014190168179384342, absolute difference3.51e-12. Step1e-5 and one
direction were used without tuning. Total cost4physicalLU/30physicalRHS
(24source+6adjoint), including the disclosed physical target root. This is
an adapter wiring check, not a new identity suite or proof of imaging benefit.

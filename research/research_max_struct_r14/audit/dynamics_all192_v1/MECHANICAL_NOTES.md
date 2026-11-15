# Mechanical dynamics notes — no GT

Analyzed 192 registered trajectories, 519199 recorded native updates, and 7 already-paid intervention checkpoints.

F_before uses the previous native material prior, including across stage boundaries; the first uniform prior is0. F_after is evaluated in the current block metric after readout. Stage changes can change the potential definition, so stagewise work sums do not establish one global Lyapunov function.

TRAJECTORIES/NATIVE_STAGE_SUMMARY retain raw maximum positive total/current/material changes and counts relative to a128eps-scaled arithmetic reference. That reference is neither a validated floating error bound nor a scientific quality gate. Iterations within each trajectory are dependent.

INTERVENTION_SNAPSHOTS compare originalPhi, data/prior, auxiliary current-metric potential and constitutive defect after0/1/8/32 native updates. Every comparison restores the same saved full state; reset branches additionally clearPR memory. Fixed update counts do not equal CPU.

Finite JVP/additivity tables copy existing paid central differences, actual feasible epsilons and branch/active-set flags. Small finite defects do not prove differentiability, contraction, a basin, or material recovery. A false recorded nonsmooth flag only means available flags did not change; unrecorded branches remain unknown.

The only available angle is the paid latest material displacement versus negative originalPhi gradient at the previous material. A raw constitutive-defect vector angle is unavailable and stays null.

Final counters may exceed last recorded native counters because budget interruption can occur inside a step or after its completion. Missing last rows are not added as synthetic iterations. Full receipt costs and per-event deltas are nested views, not additive ledgers.

No truth, family/regime label or recovery threshold was read/applied. No Maxwell/adjoint/optimizer/NN call was made. Scientific interpretation remains with parent. Re-run with a new --out after remaining36 downloads; never overwrite this pilot artifact.

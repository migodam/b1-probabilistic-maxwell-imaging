# R3 implementation/reuse audit and competing explanations

2026-11-03. Read-only audit of published R1/R2 and mature local source, followed
by the parent's explicitly authorized minimal development diagnostics/control
plan. All new writes are confined to `research/research_max_struct_r3/`.
R2 GT may be used only as published **OFFLINE_DIAGNOSTIC** reference. It cannot
become an unseen confirmation target or an online candidate selector.

## Actual model and APIs

The applicable backend is scalar-TM complex128 nonlinear disk-pulse VIE, not
the R1 vector DDA. Reuse
`research.research_fbg_r3_experiments.src.physics.MaxwellModel(public)`, with
`G:(576,576)`, `E:(576,12)`, `S:(48,576)`, real shared contrast `m:(576,)`.
It imports `research_z_alg_r4_diversity.src.physics.Model`, whose retained
factorization solves `L=I-diag(m)G`, `L j=diag(m)E` for all12 sources together.
`state(m)` returns j/Q/Y/L/lu, preserving the factor for `data_value_gradient`
and receiver-transpose Jacobian. Exact cached-material equality avoids a new
factorization; changing m invalidates it.

`model.data_value_gradient(m,Y)` computes original normalized data plus the
unchanged prior with one source adjoint (12 RHS) after the matching physical
forward. `model.reg.evaluate(m)` returns value/gradient/edge weights directly;
`model.regularizer(m)` additionally constructs an unused dense Hessian for a
plain value/gradient-only inner call. Reusing the former through a tuple
wrapper can be an implementation optimization only if prior/gradient remain
identical; this audit does not silently change the mature baseline.

Mature current entry is
`research.research_fbg_r3_experiments.src.cie_tsom_baseline.run_cie(model,Y,m0,cfg)`.
It builds receiver SVD DPCS and freezes
`JD=VD diag(sigma_D)^-1 UD* Y`, then uses
`J=JD+(I-VD VD*)F_mf A` with native projected Fourier coefficient coordinates.
APCS is generally not receiver nullspace; its data term is retained.
`fourier_apcs(shape,mf,VD,cost)` returns the native projected columns and actual
rank, not an arbitrary orthonormal replacement. Replacing its coordinates
without transporting PR geometry would change the baseline algorithm.

`_readout(J,Q,m,beta,weights,rho,bounds,reg,has_reg,cfg,cost)` fits one shared
real modified contrast per cell and applies the same physical-m prior by
bounded L-BFGS-B. Both current and material beta are explicitly transported
through m. Historical dual-beta is a two-metric alternation, not common-objective
monotone block descent. Fixed weights are `wd_s=1/||Y_s||²` and
`ws_s=1/||E_s||²`; the material prior is smooth1e-4 plus smoothedTV1e-4 on
differences divided by19.5, epsilon=.01, box[.02,20].

Strong physical corrector is
`full_maxwell_corrector.corrector(model,Y,m0,{'method':'lbfgs','max_nfev':1000,'box':(.02,20.)})`.
Here max_nfev names maxiter1000, maxfun4000, ftol1e-12, gtol1e-8, maxls30.
It records every physical incumbent and returns the smallest original objective,
not the last iterate or a GT-ranked reconstruction.

R2 `src/cie_global.py::run_cie` reuses the already-existing global PR flag,
adds full fixed source metrics and global Frobenius stopping. R2 wiring
diagnostic confirms covariance, not superior inversion. R2
`src/profile_current.py::run_profile` is known finite-penalty variable
projection; its **uniform total-energy** metrics differ from native per-source
metrics. It cannot by itself isolate hard DPCS or PR error. Exact native-metric
current LS generally needs12 differently weighted spatial solves; one shared
Cholesky is justified only for proportional source metrics.

## Source and raw-evidence review

Read R1 `src/dynamics_audit.py`, `src/tensor_families.py`, vector backend audit,
source/provenance and D1/D3 reports; R2 native/global/profile/direct runner,
capacity control, endpoint diagnostics and all relevant curved-bar raw
auxiliary/corrector histories. `R2_HISTORY_READBACK.json` extracts the decisive
rows without new Maxwell operations.

For the curved bar (`obs_b20d5a55677d4eca`), native source-normalized auxiliary
data term stays near3.61e-6 while current-state term falls1.757→.03317 and
material-state term at final beta.2 is .009666. The selected physical seed
before corrector has relative residual .41574; the corrector ends at .100735
and MAXITER. Global starts its corrector at .15007 and ends .107831 with a
relative-function-reduction stop. Its projected-gradient max2.73e-6 is above
the gtol1e-8 criterion. Native's final projected-gradient max2.91e-7 also
exceeds it. These are not certified local minima/basins.

The mf12 posthoc control reduces final current-state term to .02647 and
physical residual to .097693, while published material NRMSE worsens from
.510695 to .558744. Full APCS rank559 plus fixedDPCS17 spans576 linear spatial
directions, but the17 measured coefficients remain fixed. Thus failure of
mf12 does not falsify a hard-DPCS explanation. It does falsify the unsupported
claim that merely exposing all APCS directions necessarily restores the image.

At the native/global corrected endpoints, best mf10 fixed-DPCS L2 current
projection errors are .2576/.3020; mf12 leaves .04590/.02637. Freeing DPCS at
mf12 reduces ordinary L2 error to roundoff. The corresponding projected ports
can have residual about7.7e-4 while the actual physical material has residual
about.1. These are affine-current approximation diagnostics, not physical
recovery and not H(curl). Saved stage2 J transported to a corrected material
does not represent the original stage2 state defect.

Atlas source snapshots inspected: R1
`research/delegated/theory/ATLAS_READ_EVIDENCE.md`, sections06 and10, including
F-OPT-003/004/009 and F-BI-020/021/022/024. F-BI-021 already explicitly proposes
strong-mode corrections `[VD,Q]` rather than noisy DPCS hard locking. F-OPT-003
and F-BI-022 already cover finite-penalty current elimination; F-BI-024 gives
physical defect correction through the full resolvent. None is a new R3
identity. Physical source paths are recorded in R1 `audit/D_INPUT_PROVENANCE.json`.

Additional implementation inspected for reuse: FBG
`joint_current_model.py` (correct off-manifold derivative includes state-weight
derivative), `current_normal.py`, `profile_schur.py`, `controlled_release.py`,
`finite_response.py`; R4 `src/physics.py` and `src/solver.py`; Z-ALG asymmetric
`src/physics.py`, `representation_audit.py`, `endpoint_rank_review.py`.
Z-ALG's six-real-cell/four-real-observation families are not this576-material
imaging backend, and exact fixed-host or4/6-cell response kernels require their
stated fixed-outside assumptions. Do not replace this full-grid problem with
those constraints or claim current dark directions are legal material fibers.

## Two competing hypotheses and smallest discriminating diagnostic

**H1: noisy strong-mode hard locking is the operative incompatibility.**
The physical current at the true material need not satisfy the noisy fixed
DPCS coefficients. Full APCS cannot repair that mismatch. At fixed published
GT/native/global m, use the parent's exact minimum-state-energy calculation
under `D J=d_obs`, with `D=VD*`, `H=P_beta L` and
`A=D H^-1`. Compute `D L^-1` with17 transpose RHS from one retained physical
LU, then scale its columns by `P_beta^-1`. The minimum state energy for each
source is `b*(AA*)^-1 b`, where `b=d_obs-D j_phys`. Compare this floor against
the **same-anchor** weighted defect of stored auxiliary currents. At GT,
separate actual noise contribution from physical material mismatch. This
distinguishes an operative floor from merely observing a nonzero theoretical
noise conflict. Nearzero floor relative to actual defect weakens H1 as the
dominant quantitative explanation.

**H2: finite-penalty/resolvent shielding dominates even after freeing DPCS.**
The same-metric exact full-current profile can undervalue a poor physical
material, because relaxed currents absorb physical inconsistency. On the
same three fixed materials, compute exact full-current profile with **native
wd/ws and beta_current1**, original Phi and physical correction amplification.
No nonlinear material optimization is needed. A substantial discrepancy or
ranking reversal between exact relaxed profile and physical Phi supports this
mechanism; small discrepancy with high inner optimality instead directs the
next investigation to the nonlinear material path/termination. The existing
R2 uniform-total-energy profile is not sufficient evidence either way.

These hypotheses can coexist. The diagnostic should quantify floor, excess
defect and effective residual attenuation separately, then let parent theory
review choose the next pivot. Neither outcome licenses declaring a true basin
from capped iterations, nor repeating broad parameter searches.

## Parent-authorized first controls (not executed by this audit)

One ordinary strong-mode-release run on the disclosed R2 curved bar: keep
native stages3/6 unchanged; at final mf10 use `[VD,Q]`, with17 strong correction
coefficients initialized exactly zero and carried Fourier coefficients intact.
Keep the data residual, shared material, prior, beta schedule,500×3 steps and
same1000-iteration corrector. This is the already-known F-BI-021 control, not
a novel solver claim. Record rank, actual coefficients, current/material DOF,
matrix actions, inner costs and physical validation separately.

One clean-data native counterfactual reuses the stored published cleanY and
the identical schedule/init/prior/corrector. Both runs are
**OFFLINE_DIAGNOSTIC / POSTHOC_DEVELOPMENT**, not fresh confirmation. No online
GT selection is allowed. Their exact dispositions, costs and failures must
survive regardless of outcome. Parent owns the next scientific decision.

## First authorized execution readback

`results/diagnostics/FIXED_MATERIAL_H1_H2.json` records the GT hard-DPCS
minimum weighted state term **.0289064**, with coefficient mismatch.481806
entirely caused by actual noise (nonnoise mismatch4.86e-13). This is not
negligible beside the original native final auxiliary state term.03317, but
those quantities concern different materials and do not establish dominance.
At published GT/native/global corrected materials, exact same-native-metric
free-current profile+prior is1.14257e-5/.0105956/.00570939. GT is substantially
preferable, so these comparisons show no GT-versus-bad-endpoint ranking
reversal. Relaxed/physical native data-energy ratios are.843/.171/.0798:
shielding is stronger at the bad endpoints, without proving that the profile
globally prefers them. Fixed-material work costs3physicalLU/87physicalRHS,
36currentCholesky/36currentRHS,3smallGramCholesky/36GramRHS; wall.228610s.

Both first controls are sealed under `results/development/` and `predictions/`.
The strong-release wrapper reuses the unchanged mature source in memory,
releases final-stage17DPCS coordinates with zero seed, and restores its hook
in `finally`. Noisy release ends at relative physical residual.100662 versus
R2native.100735; clean native ends.107822 relative to stored cleanY. Both
performed1500 alternating steps plus a MAXITER corrector. Projected-gradient
maxima4.49e-7/1.26e-7 exceed1e-8. Release lowers current-state auxiliary term
to.01030 but raises auxiliary data term to6.82e-4; clean native auxiliary
data reaches1.19e-10 while its physical material remains poor. Thus noisy
hard locking is not necessary for the observed failure, and releasing strong
coordinates only at the final stage is not sufficient. This motivates the
parent's next material-path/termination pivot; it is not a universal rejection
of DPCS release or proof of a true basin. Control costs2064physicalLU/
49392physicalRHS,22.6430601wall/22.195747CPU seconds. Independent material
scores are delegated after seal and do not select or trigger repeats here.

## Second parent-authorized mechanism controls

After the parent observed a passive auxiliary-medium path reduce the original
residual to.0510, `experiments/known_controls.py` executed exactly two familiar
controls from the identical R2native corrected seed, with no new noise, grid,
beta search or GT selection. Known GlobalDual received up to1000iterations
then the same original physical corrector up to1000iterations. Both stopped
by relative function reduction; final residual.10073486244 and projected
gradient2.91e-7 are essentially unchanged. Paid cost104physicalLU/2484physical
RHS,100small source Gram factors,wall.785387s. Source Gram RHS are distinct
from physical LU RHS; each solve uses48receiver columns in a12×12 Gram.

The original physical objective and gradient were then both multiplied by
fixed1e6, with gtol correspondingly multiplied to.01, ftol1e-12 unchanged,
and1750maximum iterations. This keeps the physical minimizers and equivalent
gradient tolerance, while avoiding the unscaled absolute function-reduction
floor near objective.005. It ran745iterations/748evaluations and stopped by
PGTOL; original projected-gradient maximum9.915e-9 now meets1e-8. Yet relative
data residual remains.1007336867. Paid cost748physicalLU/17964physicalRHS,
wall5.023931s. No additional corrector or inverse run followed this control.

Consequently the earlier absolute-ftol stop was a real numerical termination
limitation, but fixing it does not remove this poor data fit. First-order
stationarity is stronger evidence than the prior capped endpoint, but it is
not a Hessian/local-minimum/global-basin certificate. Neither this simple
restart nor the tested known GlobalDual continuation explains the passive
path's separate improvement. The parent owns its interpretation and any new
mechanism-driven decision. Stored prediction scoring remains independent.

## Formula-driven finite coordinate/transport controls

The parent next authorized exactly two controls from the scaled-original
stationary endpoint, implemented in `experiments/finite_coordinate_paths.py`.
For one cell the exact Woodbury data response is `Delta F=d/(1-cd) o_i q_i`.
Each scan enumerates legal bounds,0 and real roots of the specified quadratic
data derivative, evaluates the **unchanged full prior** on every nonduplicate
candidate, and physically forward-validates the best predicted improvement.
This is a finite candidate rule; it does not find all stationary points of
the prior-augmented one-dimensional objective. Legal-interval denominator
minimum and the candidate's roots/pole distance are recorded.

Two-cell swaps exchange the values at each of1104fixed4-neighbor edges.
The exact2×2 Schur update uses all12source RHS and the same full prior.
Identity swaps need no small solve and are counted separately. Any proposed
improvement must strictly lower actual originalPhi after a full forward;
the initial incumbent is retained. Each accepted move would rebuild the
factor/resolvent, capped at10accepts, before the same physical corrector1000.
Full `L^-1` is explicitly paid as576RHS per anchor, not described as free
low-dimensional physics. Small swap solves are separate auxiliary work.

One-cell control evaluated2436candidates across2anchors. Its sole accepted
change was cell419, d=-2.26486e-4, with originalPhi decrease only2.47e-13;
predicted versus actual objective differed1.04e-17. The next scan found no
enumerated improvement. Residual stayed.100733686709. Neighbor control tested
all1104edges:909identities and195nontrivial2×2solves. The best nontrivial
objective change was **+5.41e-9**; no swap was accepted. Its residual stayed
unchanged. Correctors made no substantive progress.

Thus these **declared** single-cell-data-root and nearest-neighbor exchange
proposals do not escape this stationary development endpoint. They do not
exclude prior-induced interior coordinate optima missed by the finite rule,
larger/nonlocal transport, uphill intermediate paths or other initial states.
No global coordinate-optimality, true-basin or novel Woodbury claim follows.
Costs: one-cell2physicalLU/1200physicalRHS including1152resolventRHS;
swaps1physicalLU/612physicalRHS including576resolventRHS, plus195small solves/
2340smallRHS. Walls.181578/.048836s. All traces, candidate tables and selected
predictions are preserved. A subsequent system-Python readback failed to import
NumPy before science; it was repeated only as saved-array inspection in the
existing environment, with0LU/RHS and failure details in
`audit/FINITE_SCAN_READBACK.json`.

## One structural extension: all unordered swaps

The parent authorized replacing the nearest-neighbor dictionary by **all
165600unordered distinct-cell pairs**, with no distance/gradient/GT filter.
`experiments/all_pair_swaps.py` uses the same seed/objective/full prior and
same accepted-move/corrector budgets. A batch size8192 partitions memory only.
Per anchor it forms `A=O* r Q*`, `BO=O*O`, `BQ=Q Q*`; each nonidentity pair's
2×2 transfer K gives data norm change by the exact rank-two Gram contractions.
The full historical prior increment is computed on incident edges, using the
actual Regularizer.D coefficients and smooth/TV/epsilon values. A shared edge
gets the combined endpoint increment once. Thus neither a local surrogate
prior nor a different spatial regularizer is introduced.

The first complete scan had119881identity and45719nonidentity exchanges.
All were nonimproving; best nonidentity pair(1,81) had objective change
**+4.53247e-10**. The arm therefore stopped greedy search after one scan,
accepted0moves, and its common physical corrector did not alter the endpoint.
Relative residual.100733686717 and projected-gradient maximum9.915e-9 remain.
This is finite floating enumeration of the declared transposition neighborhood,
not global material/permutation optimality or a no-continuous-path theorem.

Costs:1physicalLU/612fullRHS including576resolventRHS;45719small2×2solves/
91438smallRHS in21batches; one cross Gram plus two norm Grams. Wall.118616s.
Process lifetime peak140.8MB; separately listed batch-array peak1.89MB and
three shared score Grams15.93MB are not a complete allocator profiler.
The initial script launch had a syntax typo before imports/numerics, recorded
with0LU/RHS and unavailable pre-instrumented CPU/wall in
`audit/ALL_PAIR_STARTUP_FAILURE.json`; it was fixed before this sole scan.

`audit/ALL_PAIR_SCORE_READBACK.json` compares only already-paid saved scores:
new Gram/incident-prior scores agree with previous direct receiver-response/
full-prior nearest-neighbor scores on all1104overlapping edges to1.39e-17,
at exactly the same stored anchor. No new physical or candidate solve was
added for this changed-wiring inspection. Prediction/receipt/candidate array
are sealed under `results/development/all_pair_swaps/` for independent scoring.

## Paired-ray and dressed-interaction refinements

The parent identified an omission of endpoint swaps: a bounded paired material
path can have an improving interior data stationary point, and equal values
can jointly increase despite giving an identity exchange. The next single
dictionary therefore used both directions(1,1)/(1,-1) on1104fixed nearest
pairs. `experiments/paired_material_rays.py` derives the exact response
numerator `t diag(1,s)-s t² adj(C_pair)` over quadratic determinant p(t),
builds quartic N/D with rank-two Grams, maps each legal tinterval to[-1,1],
and constructs the structurally degree≤6 numerator N'D−ND'. It enumerates
numerical real roots/bounds/0, then evaluates the same full material prior.
No tgrid, root-completeness certificate or exact totalPhi-coordinate optimum
is claimed. Numerical roots and their backward residuals are saved.

Nearest-pair run:2208direction cases,985zero-width intervals,1223degree6root
problems,2470finite candidates. The best nonzero proposal was(0,1),sign+,
t3.81845e-5, with **positive** objective difference6.34e-14. Its actual full
forward agreed with finite prediction to2.60e-18, so no finite move or new
corrector optimization was accepted. Cost3physicalLU/624physicalRHS including
576resolventRHS; wall.649186s, processpeak167.1MB. Maximum numerical-root
backward residual6.09e-14 does not prove root completeness.

The parent then specified one physics-selected graph rather than all-pair
ray enumeration: each cell chooses one partner maximizing the squared mixed
second response `Oi Cij Qj+Oj Cji Qi`. Its exact Gram norm is evaluated once
over the fixed matrix; unordered union gives at most576pairs without a distance
or gradient threshold. `experiments/dressed_neighbor_rays.py` reuses the
paired-ray engine on this graph. A fresh retained LU and576RHS resolvent are
paid; no historical C/O cache reuse is asserted.

Selected graph has567pairs,229nearest and338nonlocal; median grid Manhattan
distance3,max15; physical median.08944lambda,max.52612lambda. The1134directions
contain621nonzero-width intervals and1250finite candidates. Pair(223,246),
sign+,t=-3.937e-5 strictly reduced actualPhi by only **1.0573e-13**, with
finite/fullforward difference1.30e-17. The same corrector immediately meets
PGTOL, and physical relative residual remains.100733686715. This is numerical
first-order tolerance-scale adjustment, not a substantive escape. Cost2LU/
624fullRHS including576resolvent, wall.409758s, processpeak156.5MB.

Both predictions, graph distances, polynomial coefficients/roots, complete
finite candidate tables and failure/cost receipts are sealed. The parent
requested analysis after this graph, so no further dictionary, inverse run
or optimizer tuning followed. The negative outcomes distinguish these finite
pair structures from the separate passive-path behavior without converting
them into a global material/barrier theorem.

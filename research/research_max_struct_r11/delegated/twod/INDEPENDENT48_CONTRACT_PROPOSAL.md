# 48-object independent confirmation: root review requested, not generated

Three concrete JSON contracts are in configs/twod_independent48_{generator,solver,scoring}.json. Seed202611071107 differs from R10 seed2026110710. Each family has independent high and mixed draws at χ8/12/18:8×2×3=48; no low-only cases. Geometry/noise use independent SeedSequence children; there is no matching-geometry shortcut across regimes or contrast levels and no forward-based rejection.

Physics: inherited24×24 disk-pulse scalar-TM G/E/S,576 independent real contrast unknowns,12 exterior sources/48 receivers, same-grid finite truth, material host.1, box[.02,20], unchanged inherited physical prior. Independent circular complex noise has E|N|²=σ² with σ=.01||Yclean||F/√(48×12). Public files only operators,Y,σ,noise convention/grid/bounds; family/regime/amplitude/seed/GT remain generator/scoring-only. This phase is not cross-grid evidence.

## Geometry ranges (unchanged R10 draw_material, independently redrawn)

Each draw first has translation x/y uniform[-.08,.08], angle uniform[0,2π], axis scales uniform[.9,1.1] in fractional ROI halfwidth coordinates. Existing topologies, not eight variants of curved bars:

|Family|Parameters/ranges|Topology|
|---|---|---|
|annulus|outer[.56,.69], inner[.27,.36]|closed shell with hole|
|slotted_shell|outer[.60,.70],inner[.34,.40],slot halfwidth[.12,.20]|open shell|
|multi_island|separation[.30,.39],radius[.19,.25]|three disconnected bodies|
|curved_bars|radius[.46,.56],halfwidth[.13,.17]|two curved bodies|
|comb|shift[-.035,.035],halfwidth[.10,.13]|four teeth and spine|
|coated_inclusion|major[.56,.68],minor[.40,.49],fixed inner ellipse .27/.21 offset(.09,-.04)|nested two-material inclusion in mixed regime|
|perforated_star|phase[-.2,.2],radius[.51,.60],five-lobe modulation.22,holes radii.16/.12|perforated five-lobe boundary|
|bridge|body radius[.26,.32],separation[.37,.43],neck halfwidth[.10,.14]|two bodies with thin connecting neck|

High cases set the entire support to χ8/12/18. Mixed sets support to2.5, and high region is u+.23v>0 (coated inclusion uses its specified inner ellipse). Only legality checks: support nonempty, mixed both phases nonempty, finite box material. Fail without reroll for any illegal draw; no GT quality gate or shape selection.

## All six arms, one mature schedule

1. CIE-single: inherited sourcewise PR, beta-current=beta-material.
2. Dual-CIE: inherited sourcewise PR, beta-current=1.
3. Global-Dual: inherited global PR, scalar mature material readout.
4. fullG19: global PR + R10 full moving-denominator Green material readout.
5. freezeA: same global PR + R11 lagged_material_readout/prepare_lagged.
6. contractedT: same global PR + R11 contracted_lagged_material_readout/prepare_contracted_lagged.

All use host initialization.1, public noise-ball99 DPCS, original source metrics, FFT-APCS stages mf3/6/10 at500 iterations each and mf12 at1500, miniter50/coefficient_tol0. Material beta1/.2/.2/.2; current beta all1 except CIE-single follows material beta. Candidate inner allowance8 per R10/R11 readout, fixed before confirmation. Physical corrector original L-BFGS-B maxiter3000/maxfun12000, ftol1e-12/gtol.01/scaling1e6/maxls30. No wall cutoff, no GT iterate selection, no budget padding.

Each child has480 whole-child CPU seconds, independent checkpoints60/120/240/480, complete imports/setup/aux/corrector/failures/tail charged. Atomic overrun is charged and excluded from earlier-cap prediction. Natural stopping retains actual spent CPU. Report actual GetProcessTimes, original physicalLU/RHS frontier, auxiliary counts separately, and each object's all-method minimum achieved LU endpoint. Same allowedCPU is not equal actual expenditure or equal LU.

Primary material thresholds are inherited NRMSE≤.15, strict supportχ>1 IoU≥.7, and legalbox; continuous scores remain primary numerical evidence. NRMSE normalizer=max(GT)-.1. Region errors include χ8 explicitly. All288 terminal outcomes before any new GT scoring, all failures/missing checkpoints retained, separate root-authorized offline scorer. No confirmation GT used to choose an arm or parameters.

Only after root review: set frozen statuses, generate48 once, then immediately execute CIE/Dual/Global baseline manifest on authorized P8/P10. Candidate six-arm configuration freezes together; no waiting for exposure/confirmation GT to choose parameters.

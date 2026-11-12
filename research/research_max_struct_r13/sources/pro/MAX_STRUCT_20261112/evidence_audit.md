# R12 paired evidence: independent source and saved-state audit

Date: 2026-11-12, Asia/Singapore. Frozen repository commit: `ec3a37da6f3ee90beb61e753067d4f639f4118f5`.

## 1. What was actually read and computed

Read the original R12 task, R12 continuity instructions, the complete paired report, `run_paired.py`, the exact imported `run_cie_fft` implementation, its inherited R6 `arm_config`, `NativeFourierAPCS`, the field-current projector and runner, and all `offline_v2` CSVs. Read R8 `NEURAL_WARM_STARTS.md` and `WARM_DEVELOPMENT_ANALYSIS.md` to check whether the proposed explanations had already been discovered or tested.

The new computation uses the eleven actual paired starts, their frozen G/E/S and observations, saved physical initial/raw/primary arrays, saved receiver singular vectors and DPCS currents, 44 run costs, and offline truth only for retrospective outcome analysis. It additionally retrieves ten full saved candidate pools from the public release archive and scores their already-paid physical material candidates offline. No optimizer was rerun. These analyses require **zero new Maxwell forward solves**. The corresponding code and numerical tables are in `audit_saved_evidence.py`, `plot_evidence_audit.py`, `analyze_candidate_pool.py`, and `evidence_analysis/`.

This is a retrospective mechanism analysis of exposed development cases. It is not a new independent confirmation set, a new algorithm trial, or a validated online classifier.

## 2. The most consequential implementation detail

The supplement invokes

`solver.arm_config('a_plus', common, Y, model.E, [dict(mf=12,maxiter=1500)])`.

The last argument is a **continuation stage**, not a replacement of the stage list. R6 `arm_config` prepends the frozen stages

\[
(m_f,\beta_m,\beta_J,N)=(3,1,1,500),\ (6,.2,1,500),\ (10,.2,1,500),
\]

then appends \((12,.2,1,1500)\). `initialize_physical_current` acts only at the first stage, so the initial legal least-squares projection uses **36 Fourier columns**, not the full 576-column Fourier grid. This is also recorded independently in all eleven `G19_EXECUTION.json` files as `initial_projection_columns=36` and `initial_projection_rank=36`.

However, **all eleven runs later reached the full \(m_f=12\) stage**. A and C completed 3,000 auxiliary steps. The B runs completed respectively 1,954, 1,983, 1,921, and 2,275 steps. The first three stages account for 1,500 steps, so each B run had already spent substantial work in the full Fourier stage before its proposal budget ended.

On a 24×24 grid, the four-corner mask with \(m_f=12\) includes every Fourier mode. Its projected basis spans the entire orthogonal complement of the fixed DPCS subspace. Therefore:

- Initial Fourier truncation can substantially change the trajectory.
- It is **not** an explanation that the final G19 auxiliary family permanently lacked those Fourier directions.
- The final affine family still fixes the retained receiver coordinates to \(J_D\). This is a different restriction and persists at full Fourier rank.

Any proposed mechanism needs to explain why later expansion does or does not repair the trajectory, rather than treating the initial 36-column defect as a permanent rank limitation.

## 3. Exact decomposition, and what the eleven starts show

Let \(J_p\) be the saved physical current at the shared starting material. Let \(P_D=V_DV_D^*\) be the retained receiver-current projector, and let \(P_U\) project onto the actual projected \(m_f=3\) Fourier columns. Their ranges are orthogonal. The legal initialization is

\[
J_{\mathrm{init}}=J_D+P_UJ_p.
\]

Consequently

\[
\Delta J=J_{\mathrm{init}}-J_p
=\underbrace{J_D-P_DJ_p}_{d_D}
-\underbrace{(I-P_D-P_U)J_p}_{d_F},
\qquad \langle d_D,d_F\rangle_F=0.
\]

This decomposition was already derived in R8; it is **not a new theorem**. Its application to the eleven latest paired starts is new analysis.

| Start | Relative total current loss | Fourier share of lost current energy | Fixed-DPCS share | Fixed-DPCS port change / ‖Y‖ | Fourier port change / ‖Y‖ |
|---|---:|---:|---:|---:|---:|
| A ring 0001 | .818064 | 99.9736% | .0264% | .136590 | .003918 |
| A shell 0004 | .955891 | 94.0707% | 5.9293% | .589725 | .000976 |
| A bridge 0022 | .845417 | 99.6432% | .3568% | .168206 | .002534 |
| B 0001 | .805530 | 99.1060% | .8940% | .182944 | .002729 |
| B 0003 | .798085 | 98.9573% | 1.0427% | .217339 | .005673 |
| B 0010 | .908098 | 98.7670% | 1.2330% | .315194 | .003101 |
| B 0022 | .980478 | 90.7035% | 9.2965% | .814907 | .000909 |
| C 0001 | .437273 | 81.7386% | 18.2614% | .892003 | .002642 |
| C 0004 | .497748 | 97.4271% | 2.5729% | .561383 | .009115 |
| C 0022 | .716089 | 62.7607% | 37.2393% | .842722 | .000116 |
| C 0027 | .494110 | 94.0768% | 5.9232% | .786479 | .017181 |

The median Fourier share is **97.4271%**. Nevertheless, the fixed DPCS correction supplies at least **99.9178% of the induced received-field change energy** in every case. A ring is particularly instructive: fixed DPCS accounts for only 0.0264% of lost current energy but causes a port shift 34.9 times larger than the Fourier component. C0022 gives a ratio exceeding 7,200. Because the receiver-SVD components are orthogonal also after applying S, this is not a cancellation artifact. The DPCS change largely corrects an initial measurement mismatch, so this observation alone does not establish that it is harmful.

This directly demonstrates why **current energy, receiver visibility, and material recovery cannot be substituted for each other**. It does not demonstrate that the large Fourier component is harmless: propagation and material readout can amplify and mix it.

At the starting material \(m\), the induced constitutive defect is exactly

\[
e=(I-D_mG)\Delta J.
\]

The two orthogonal current components need not remain orthogonal after multiplication by \(I-D_mG\). Even a material-oriented diagnostic should act on this transformed defect and the common-source material tangent \(D_hQ\), rather than ranking directions only by \(\|\Delta J\|\) or \(\|G\Delta J\|\). A rowwise local shared-real material pullback has the form

\[
h_i^{\mathrm{LS}}=\frac{\Re\langle Q_i,e_i\rangle}{\|Q_i\|^2},
\]

for an interior physical starting state and fixed-current least squares. This is a descriptive constitutive tangent, **not the actual nonlocal G19 material step**, and not a finite-distance success certificate.

### The scalar projection defect is not an established success predictor

For the eleven cases, the retrospective Spearman correlation between initial current loss and final G19−DR NRMSE is .3818. Restricting to the seven actual bad endpoints in A+B gives **−.1071**. The sample is small and exposed; the associated p-values are descriptive and no significance or forecasting claim is justified.

The concrete counterexamples matter more than fitting a threshold:

- B0022 has the largest loss, .9805, and DR rescues it.
- B0010 has loss .9081 and remains a failure for every compared method.
- A0001 has loss .8181 yet G19 reaches a more accurate final material than DR.
- C0022 has loss .7161 and both G19 and DR finish accurately, with G19 slightly better under primary selection.

R8 had already invalidated the simpler claim that reducing the initial constitutive defect guarantees better recovery: all eight development defects decreased under the physical-complement anchor, but several materials worsened. Its 42-object independent extension also failed to turn this interface change into a stable replacement for NN→direct.

## 4. Hard DPCS equality and actual good physical states

For every saved DR physical candidate we computed

\[
\frac{\|J_D-P_DJ_{\rm physical}\|_F}{\|J_{\rm physical}\|_F},\qquad
\frac{\|S(J_D-P_DJ_{\rm physical})\|_F}{\delta}.
\]

Here \(\delta=\sigma\sqrt{\chi^2_{2N_y}(.99)/2}\), exactly matching `run_paired.py`. The NPZ `noise_radius` field is normalized, so using it directly as an absolute field radius would produce a misleading ratio; the supplied audit code uses the correct absolute convention.

Ten of eleven DR primary candidates lie inside the declared data ball. Their physical-current discrepancy with the hard DPCS affine family is **0.2894%–1.1131% of current norm**, while its received-field component occupies **0.3498–0.6028 data-ball radii**. These are nonzero, measurable discrepancies, but modest in total current norm. The failed B0010 instead has 10.10% discrepancy and about 30.06 radii in the retained coordinates.

The correct statement is that exact hard retained-data equality excludes those ten candidates' **exact physical triples** from the exact lifted feasible intersection. It does **not** prove that G19 cannot output those materials: G19 permits nonzero constitutive residual during its search, returns a material candidate, and applies a separate physical corrector. Nor does a nonzero DPCS mismatch prove that relaxing it will improve the search. This distinction is necessary to avoid turning an algebraic incompatibility into an unsupported global algorithmic impossibility claim.

## 5. Corrector audit: data feasibility does not resolve material selection

The 44 actually executed method/case transitions contain eight increases in NRMSE between raw and primary. Twenty B bypass objects must not be counted four times as additional experiments. DR has two harmful transitions out of eleven, both in the four NN development objects:

| DR start | Raw NRMSE | Primary NRMSE | Raw residual / δ | Primary residual / δ | Material step NRMSE | Change in predicted data / δ |
|---|---:|---:|---:|---:|---:|---:|
| C0004 | .074037 | .078824 | .922251 | .916391 | .090955 | .078946 |
| C0022 | .055806 | .067208 | .960022 | .955784 | .055116 | .449022 |

C0004 moves a considerable distance in material space while barely changing the predicted data relative to the declared noise radius. C0022 is stronger evidence against a simple objective explanation: its data term decreases by \(4.46\times10^{-9}\), its prior decreases by \(1.236\times10^{-6}\), and the material becomes less accurate. Both candidates remain data-feasible.

Therefore a constrained corrector can prevent leaving a prescribed ball, but cannot promise preservation of GT accuracy **within** the ball. The already registered secondary rule—minimum prior among data-feasible paid candidates—also selects the same harmful C0022 corrected candidate. A proposal to replace \(D+R\) with `min R subject to discrepancy` has useful estimator semantics, but this dataset already contains a counterexample to interpreting it as a material-accuracy guarantee.

Conversely, abandoning all correction also loses genuine gains. The three A cases improve markedly on average, as do three effective B rescues. The scientifically relevant separation is between generating diverse physical candidates, identifying data-distinguishable differences, and applying an explicit material-selection assumption. A universally safe GT-preserving corrector cannot be inferred from the current objective or a noise-ball safeguard.

### Incumbent semantics that affect interpretation

`raw` is the minimum original physical \(\Phi\) among the initial state and all physically audited proposal candidates. It is not necessarily the last DR/AP shadow or last G19 material readout. `corrected` and `primary` are global incumbents across all eligible candidates, rather than necessarily the terminal corrector iterate. Thus nonincreasing \(\Phi\) is built into selection; it is not independent evidence that the selected material improved.

The projection runner audits every 50 iterations. Any claim about unobserved intermediate material quality would require further saved states or a rerun; the current scalar projection history does not contain every iterate's material.

### New candidate-pool test: did AP generate a good material that Phi discarded?

Ten full `CANDIDATES.npz` arrays were extracted from the original public paired-results ZIP using normal HTTP Range requests. The extraction transferred 23,538,520 bytes, validated every selected entry's CRC and uncompressed length, and did not fetch the entire 176,584,550-byte archive. The retrieval receipt is in `paired_selected_candidates/retrieval_receipt.json`. Each pool was aligned with its original `HISTORY.json`; only eligible initial/proposal candidates were included before the corrector.

The oracle here is **only a retrospective diagnostic**, not a new output rule: it asks whether an accurate material ever appeared among the already-audited proposals.

| Exposed case | Best generated AP proposal NRMSE | Best generated DR proposal NRMSE | Published DR raw NRMSE |
|---|---:|---:|---:|
| A shell 0004 | .514830 | .173264 | .173264 |
| B 0001 | .423836 | .117463 | .117501 |
| B 0003 | .424947 | .167737 | .167742 |
| B 0022 | .295946 | .063050 | .063057 |
| C 0001 | .273914 | .107503 | .107522 |

Each method/case has 62 saved pre-corrector physical candidates, including the duplicated initial audit and 60 regular proposal audits. Even the oracle AP proposal is much worse than the DR proposal in all five cases. The best AP proposals remain far outside the ball; their residuals are respectively 58.94, 23.08, 32.29, 77.69, and 83.51 radii. DR's published raw selection differs from the proposal oracle by less than .00004 NRMSE in these cases.

Thus **the observed AP-versus-DR advantage in these five paired examples is not explained by Phi rejecting an equally good AP proposal that was already in the saved candidate pool**. The difference lies earlier in candidate generation along the matched full-state iterations. This sharpens the evidence for a dynamics effect. It still does not prove that no unsaved intermediate AP iterate was useful, that no longer AP run can succeed, or that DR is globally superior.

### Another fairness factor: physical candidate audit opportunities

G19 and DR use the same original physical Phi for final selection, but their schedules expose very different numbers of candidates to that selector. Directly counted from `HISTORY.json`:

- G19 has five proposal-phase physical rows in A/C: the initial material plus four stage-boundary audits.
- G19 has only four proposal-phase physical rows in B: initial plus the completed mf3, mf6 and mf10 stages.
- AP and DR have 61 proposal-phase physical rows: duplicate initial plus 60 regular audits.

The B G19 runs did enter mf12, but the proposal cap interrupted that stage before an eligible end-of-stage physical audit. Their published raw candidates therefore came from the initial material or the three earlier stage endpoints. This is a **finite-budget candidate-auditing confound**, distinct from both the eventual Fourier span and the DR reflection mechanism. Costs are correctly charged, so it does not invalidate the reported algorithm-versus-algorithm outcome. It does prevent attributing the entire difference to one dynamical or representation ingredient. The cleanest next G19-versus-DR mechanism test should match paid physical audit opportunities, or explicitly compare schedules as a separate factor. AP versus DR remains better isolated because their audit schedules are identical.

## 6. Costs and what is actually paired

| Cohort | G19 added CPU | DR added CPU | AP added CPU | Direct added CPU |
|---|---:|---:|---:|---:|
| A, 3 starts | 518.2344 s | 426.8594 s | 234.1094 s | 237.2031 s |
| B, 4 triggered rescues | 355.0000 s | 336.8750 s | 239.3281 s | 62.5000 s |
| C, 4 NN starts | 790.3594 s | 544.3281 s | 550.3750 s | 374.2969 s |

These are the completed external whole-child costs from `offline_v2/CHILD_COST_LEDGER.csv`. Common upper limits are not equal actual expenditure. The B DR strategy retains all twenty nontriggered objects and adds 336.875 s to the original 4,290.1875 s G19 cost, an increase of 7.8522%. Its NRMSE change .143891→.106261 includes every object. Those twenty retained cases are a bypass policy, not proof that running DR on them would be harmless.

The ten good-data DR endpoints should not be equated with ten successes under every material metric. In the original B24, all twenty untriggered objects have NRMSE below .15, but eleven fail the joint label because IoU is below .7. Their average NRMSE is .096506. That is evidence of remaining structural error at low residual, not proof of exact nonidentifiability for each object.

For the C cases, primary means are G19 .086768 and DR .086094; registered secondary deployment means are G19 .080104 and DR .082901. The ranking changes with the predefined material estimator. This is another reason to separate search competence from selection quality.

## 7. Most useful implications for the next experiment

1. **Use a factorial intervention, not total current loss as the mechanism.** Keep dynamics and readout fixed, separately change early Fourier continuation, the fixed DPCS equality, and the preservation of the starting physical state. Compare equal actual costs and retain the original schedule. The R8 physical-complement intervention is mandatory prior evidence, not a new proposed method.
2. **Test persistent restrictions on material-relevant tangents.** A complete final Fourier basis removes its own truncation obstruction, but the receiver lock persists. Measure its effect after physical propagation and shared-material pullback. A large norm or a tiny data-space angle alone does not establish that the relevant branch is reachable.
3. **Compare AP and DR in the same full state and selection protocol.** Current evidence gives a credible dynamics difference because AP and DR use the same projectors, metric, initial physical fields, and 3,000 steps. The newly scored candidate pools show that five AP failures are not rescued by an offline oracle over their observed proposal materials. This still does not isolate DR versus G19, because their representations, objectives, and physical candidate-audit opportunities differ.
4. **Treat in-ball ambiguity separately.** Adding a discrepancy constraint cannot cure the C0022 counterexample. Any claimed selection improvement needs an independent prior, withheld physical information, acquisition change, or a stated uncertainty/solution-family output.
5. **Do not call this new confirmation.** The mechanism figure and correlations use eleven exposed starts. Freeze any proposed predictor before new cases. Include failures B0010, G19-favored ring A0001, and both harmful-corrector examples in the next falsification set.

## 8. Source map

All paths are relative to the frozen repository:

- `research/research_max_struct_r12/supplement_paired/src/run_paired.py`: exact initialization, budgets, candidate selection, correction, saved-state semantics.
- `research/research_max_struct_r8/src/benchmark_solver.py`: inherited arm-config dispatch.
- `research/research_max_struct_r6/src/benchmark_solver.py`: first three stages and appended continuation.
- `research/research_max_struct_r8/src/cie_fft.py`: receiver lock, warm-current least-squares projection, stage progression.
- `research/research_max_struct_r6/src/fft_tsom_ops.py`: exact four-corner Fourier mask.
- `research/research_max_struct_r12/supplement_paired/runtime/field_current.py`: identical AP/DR projectors and distinct update dynamics.
- `research/research_max_struct_r12/supplement_paired/RESULTS/offline_v2/`: scores, corrector decompositions, receiver-current projection defects, actual completed process costs.
- `research/research_max_struct_r12/supplement_paired/inputs/` and `RESULTS/*/{initial,raw,primary}.npz`: actual fixed operators and physical states used in new analysis.
- `research/research_max_struct_r8/NEURAL_WARM_STARTS.md`, `WARM_DEVELOPMENT_ANALYSIS.md`: earlier exact split, fixed physical-complement anchor, 8-object development counterexamples and 42-object confirmation.

The source and implementation findings are independent of whether a new noise-compatible geometry ultimately succeeds. They constrain the claims that any such proposal can honestly make.

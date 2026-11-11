# Reduced true3D exploratory material inversion

Three public development observations: matched mixed coated cube, matched connected unequal inclusions, and the same coated cube generated on an8³ grid and inverted on6³. All inverse grids have216 real chi unknowns and648 complex vector-current components, with60 common dipole sources and180 electric receiver channels. Bounds chi[0,4] (epsilon[1,5]), uniform chi=.25 and1% raw iid circular-complex noise are fixed. This is a reduced dyadic CM+RR DDA adaptation, not a literal2019 reproduction.

Baseline methods directly call inherited R8 VectorPR3D.optimize_pr: CIE_sourcewise, Dual_sourcewise, CIE_global and Aplus_global. Their120 auxiliary cycles use inherited mf1→2 and beta6→1 (Dual current beta3 late), inner12. AP/DR use scale6, 3000 projections/audit50/320CPU cap with one raw iid99% noise ball and no whitening. Physical candidate selection and all common correctors use the same inherited source-normalized physical Phi and same quadratic prior, with corrector maxiter120. The projection feasibility ball and candidate-scoring objective have explicitly different data normalization.

All18 predictions were sealed before scorer GT read. Selection never uses GT. Crossgrid truth for scoring is the same analytic geometry averaged onto the inverse cells. No independent fresh confirmation or high-epsilon8 claim.

| Case | Method | Final NRMSE | Raw NRMSE | Physical residual | Phi | Paid CPU s |
|---|---|---:|---:|---:|---:|---:|
| mixed_coated_cube | CIE_sourcewise | 0.10565 | 0.32550 | 0.01006 | 0.000128339 | 4.63 |
| mixed_coated_cube | Dual_sourcewise | 0.10410 | 0.31209 | 0.01005 | 0.000128331 | 9.46 |
| mixed_coated_cube | CIE_global | 0.10502 | 0.32550 | 0.01005 | 0.000128339 | 9.18 |
| mixed_coated_cube | Aplus_global | 0.10537 | 0.32562 | 0.01005 | 0.000128338 | 17.45 |
| mixed_coated_cube | AP | 0.09252 | 0.09777 | 0.01005 | 0.000128297 | 67.99 |
| mixed_coated_cube | DR | 0.09456 | 0.11892 | 0.01004 | 0.000128299 | 60.75 |
| connected_inclusions | CIE_sourcewise | 0.13614 | 0.63772 | 0.01033 | 0.000133553 | 8.60 |
| connected_inclusions | Dual_sourcewise | 0.13234 | 0.62107 | 0.01033 | 0.000133649 | 8.27 |
| connected_inclusions | CIE_global | 0.13937 | 0.63772 | 0.01034 | 0.000133762 | 8.00 |
| connected_inclusions | Aplus_global | 0.13331 | 0.63510 | 0.01033 | 0.000133653 | 17.06 |
| connected_inclusions | AP | 0.08948 | 0.24311 | 0.01027 | 0.000132603 | 38.67 |
| connected_inclusions | DR | 0.08146 | 0.07712 | 0.01026 | 0.000132379 | 40.51 |
| mixed_coated_cube_crossgrid8 | CIE_sourcewise | 0.22367 | 0.30561 | 0.01031 | 0.000138152 | 8.37 |
| mixed_coated_cube_crossgrid8 | Dual_sourcewise | 0.22309 | 0.29119 | 0.01031 | 0.000138144 | 8.33 |
| mixed_coated_cube_crossgrid8 | CIE_global | 0.22350 | 0.30561 | 0.01031 | 0.000138155 | 6.57 |
| mixed_coated_cube_crossgrid8 | Aplus_global | 0.22352 | 0.30606 | 0.01031 | 0.000138158 | 10.72 |
| mixed_coated_cube_crossgrid8 | AP | 0.21923 | 0.14442 | 0.01031 | 0.000138208 | 37.70 |
| mixed_coated_cube_crossgrid8 | DR | 0.22462 | 0.27054 | 0.01031 | 0.000138152 | 40.02 |

All original120 common correctors hit their iteration cap. Before any GT score, root requested an identical additional max1000iterations/60CPU seconds physical continuation for all18 endpoints; CONTINUATION_CONTRACT.json records that common budget diagnostic. Original120 outputs and costs remain intact. Both18+18 endpoints were sealed before scoring.

| Case | Method | Longer-corrector NRMSE | Physical residual | Phi | Extra CPU s | Stop |
|---|---|---:|---:|---:|---:|---|
| mixed_coated_cube | CIE_sourcewise | 0.09764 | 0.01005 | 0.0001282785 | 3.92 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube | Dual_sourcewise | 0.09764 | 0.01005 | 0.0001282785 | 3.64 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube | CIE_global | 0.09765 | 0.01005 | 0.0001282785 | 3.63 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube | Aplus_global | 0.09765 | 0.01005 | 0.0001282785 | 3.72 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube | AP | 0.09763 | 0.01005 | 0.0001282785 | 3.83 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube | DR | 0.09766 | 0.01005 | 0.0001282785 | 3.25 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| connected_inclusions | CIE_sourcewise | 0.08239 | 0.01025 | 0.0001322307 | 6.83 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| connected_inclusions | Dual_sourcewise | 0.08237 | 0.01025 | 0.0001322307 | 7.25 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| connected_inclusions | CIE_global | 0.08241 | 0.01025 | 0.0001322307 | 6.88 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| connected_inclusions | Aplus_global | 0.08238 | 0.01025 | 0.0001322307 | 6.79 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| connected_inclusions | AP | 0.08231 | 0.01025 | 0.0001322307 | 5.79 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| connected_inclusions | DR | 0.08242 | 0.01025 | 0.0001322307 | 5.58 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube_crossgrid8 | CIE_sourcewise | 0.22557 | 0.01031 | 0.0001380929 | 4.34 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube_crossgrid8 | Dual_sourcewise | 0.22558 | 0.01031 | 0.0001380929 | 4.49 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube_crossgrid8 | CIE_global | 0.22556 | 0.01031 | 0.0001380929 | 4.37 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube_crossgrid8 | Aplus_global | 0.22556 | 0.01031 | 0.0001380929 | 4.49 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube_crossgrid8 | AP | 0.22554 | 0.01031 | 0.0001380929 | 4.56 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |
| mixed_coated_cube_crossgrid8 | DR | 0.22557 | 0.01031 | 0.0001380929 | 4.25 | CONVERGENCE: RELATIVE REDUCTION OF F <= FACTR*EPSMCH |

A relative-function optimizer stop is not a local-minimum certificate. The original120 stage and longer common continuation answer different budget questions. Neither extending the common corrector nor a lower objective establishes superior material recovery by itself.

CPU is measured work, not equal-compute evidence. All fixed-background Cholesky/setup/eigensystems, actions, projection multiplier evaluations, physical audits, common correctors and saved full audit-state archives are recorded separately in COST_LEDGER.json and receipts. The declared work differs substantially across method families; paid prefix histories permit physical-Phi comparisons but no claim of equal CPU superiority is made. Archive serialization is charged and can be substantial. Numerical recovery/residual and scientific mechanism support remain separate judgments.

Evidence: CONTRACT.json, PROJECTION_CONTRACT.json, GENERATION_RECEIPT.json, data_public/, sealed_truth/, RESULTS/{case}/{method}/prediction.npz and RECEIPT.json, projection/history.json and states.npz, PREDICTION_SEAL.json, SCORES.json, COST_LEDGER.json, MATERIAL_SLICES.png. No new physical solve occurs in this report compiler.

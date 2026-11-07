# Fixed absorptive DEV: factual result

All eight fixed runs finished with zero execution failures on one Mac BLAS1 lane. All eight public paid-incumbent predictions were sealed before this audit read scoring GT. The four objects were already exposed DEV, not independent confirmation. No phase/beta scan, extra solver, or GT candidate selection was performed.

|Exposed case|Green NRMSE / IoU|Green material gate|K0 NRMSE / IoU|K0 material gate|
|---|---:|---|---:|---|
|r10_twod_0001|0.436362 / 0.339869|FAIL|0.421783 / 0.350993|FAIL|
|r10_twod_0004|0.110632 / 0.748837|PASS|0.515137 / 0.086957|FAIL|
|r10_twod_0022|0.492486 / 0.256545|FAIL|0.097206 / 0.706349|PASS|
|r11_bothblocks16_0000|0.953935 / 0.051724|FAIL|0.928940 / 0.079038|FAIL|

Each arm passes material and CN99 on only one of four objects: Green on the shell, K0 on the bridge. The other six predictions fail CN99 as well as material recovery. Green corrector stops: ring PGTOL, other three max3000 iterations; K0 stops: shell/bridge PGTOL, ring/new0000 relative-F. None reaches the 480-CPU cap. These are complete trajectory interventions with identical frozen prior/initialization/physics settings, not same-state counterfactuals. Existing real-BothBlocks and original-method endpoints remain in `COMPARISON_WITH_EXISTING_REFERENCES.csv`; their costs are not charged again.

Complex-current adjoint wiring errors were 1.33e-15 (K19) and 1.66e-15 (K0). Real-material gradient finite-difference errors including the original prior were 1.78e-8 and 6.41e-9. Wiring used one public saved state and zero new FW. Production source matched the pre-run snapshot at final audit.

The imaginary shift itself is an established complex-CIE option; K0 uses the same eight-step material readout/real-CIE seed proposal accepted by the new objective and is not mature original CIE. No novelty, uniform recovery, or causal uniqueness claim follows from this DEV audit. The parent thread owns theoretical interpretation; no further parameter search is authorized here.

Exact recorded CPU: inversion children 759.281010s (Green 494.430131, K0 264.850879), controller 0.141694s, wiring 0.870579s, public seal+offline score 2.121386s; total 762.414669s. Counts: 17,522 LU/FW at576 order,17,487 adjoint calls;210,264 forward+209,844 adjoint=420,108 complex RHS. Small denominator LU/FFT counts are separate auxiliary costs, not physical LU. Mac shared-machine campaign wall is not a throughput guarantee or cross-host speed comparison.

The parent scientific review concludes that the bounded current transform changes the successful object from the ring to the shell rather than supplying a default replacement; K0 also gives no general improvement. No additional absorptive experiment is planned. The frozen48 and3D campaigns continue independently.

[Four-object material plate, shared0–20 PNG](../../../FIGURES/absorptive_dev/ABSORPTIVE_DEV_4x6.png) / [PDF](../../../FIGURES/absorptive_dev/ABSORPTIVE_DEV_4x6.pdf). All24 panels have saved materials; no reference was missing. Older references may come from different hosts, so CPU is not horizontally comparable with those columns. The eight new arms all ran on the same Mac and their recorded CPU is comparable within that fixed batch.

Raw material paths, availability, inherited metrics and exposure/CPU scope are recorded for every panel in the [source table](../../../FIGURES/absorptive_dev/MATERIAL_SOURCE_TABLE.csv). Production sources: [absorptive operator](../../../src/absorptive_denominator.py), [adapter](../../../src/absorptive_adapter.py), [isolated runner](../../../experiments/absorptive_dev_run.py), [plot script](../../../experiments/plot_absorptive_dev.py). The [frozen task manifest](TASKS.json), [solver config](SOLVER_CONFIG.json), [all8 pre-GT seal](ALL8_PREDICTIONS_SEALED_BEFORE_GT.json) and [source snapshot](source_snapshot/) retain executable provenance.

New saved-material rendering cost: external full-child CPU 1.699283s,0new solver/FW/adjoint,0fresh48GT reads. This is additional to the 762.414669s batch accounting above; total including this render is 764.113952s. Internal plot timer is nested, not added again.

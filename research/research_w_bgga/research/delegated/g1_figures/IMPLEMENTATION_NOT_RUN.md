# G1 factual postprocessor handoff

Status: **CODE_WRITTEN_NOT_RUN**. The assigned worker did not execute this script, physics, training, or a syntax/test runner. No source replay artifact, gate, main report, GT, historical array, or NPZ was modified or opened by this code task.

Parent entry point: `tools/g1_figures.py` from the campaign root, using the campaign's existing Python environment with NumPy and Matplotlib. The default run generates all tables and five PNG/PDF pairs. `--skip-figures` writes explicit NOT_RUN figure records. A later run replaces only the three named root result files and this tool's own generated files.

The root outputs are `results/G1_METHOD_TABLE.csv`, `results/G1_PARENT_TABLE.csv`, and `results/G1_FACTUAL_SUMMARY.json`. Detail CSVs, paired-parent bootstrap, figure manifest, and PNG/PDF files go in `research/delegated/g1_figures/`.

Per-case comparisons use GATES.G1 thresholds for joint TV, both coordinate CDF sup errors, both event probability errors, and important predeclared joint-cell errors. The actual important-cell reference-mass cutoff is applied to saved joint-cell scores. Regret uses the saved per-case limit after checking that its recorded formula agrees with GATES in the replay's saved dimensionless convention. Missing/undefined values remain explicit, and every method's pass fraction keeps the original 24-case denominator. Original exact-bound violations use strictly positive signed violations with no tolerance or clipping. Interval audit status remains separate from probability/risk accuracy.

Method mean/median statistics first average both noise repeats within each parent. A parent with either repeat missing is not silently reduced to one repeat. Complete-parent counts are exposed for every metric. The paired bootstrap requires all 12 parent pairs, uses the same 2000 resamples for every method/metric with seed 27102801, and compares each method minus M0. These are descriptive exposed-development intervals; scientific/publication approval and G2 are not assessed.

The five figure stems are:

- `G1_parent_joint_TV_and_regret`: all seven methods, parent means and full two-repeat ranges; exact zeros shown on symlog axes.
- `G1_M4_unknown_mass_queries`: all 24 unknown-mass upper traces at saved checkpoints, a log complement view through query 2949, and stop query counts. Saved U=1 entries are visibly marked without altering their CSV values.
- `G1_M5_cluster_count`: all 24 wall-cluster counts relative to the frozen 118 groups.
- `G1_cached_method_CPU`: cached-array method CPU excluding I/O, explicitly non-cold and not PDE-speed evidence.
- `G1_raw_bounds_and_normalization_scope`: the full positive original signed-component distribution, original and directed-audit case counts by checkpoint, and the conditional normalization scope. Positive likelihood-unit violations are counted separately from probability-unit histogram values.

Normalization repair is complete only when all 24 expected cases and seven checkpoints exist and the decision records completion. Original violations remain intact even when the directed audit observes none. The repair figure and summary explicitly exclude Maxwell/FEM forward error, float64 Gaussian-potential formation, and continuum PDE certification. Its recorded cost includes independent audit work and is not relabeled an isolated deployment certificate timing.

Known current factual inputs inspected during coding: the replay receipt records 168 case-method rows, 495600 online logical cached queries, 70800 offline-reference logical cached queries, actual new FE/RHS 0, and a raw interval violation. SCORE_TABLE has all seven methods completed on 24 exposed cases with 2950 logical cached queries each. The normalization decision records 24 completed cases and no corrected bound violations, conditional on the finite numerical log-potential model. These observations are inputs, not output-generation or scientific-acceptance claims.

The parent must run the script, inspect its input-problem list, and visually inspect the figures before calling their generation or review complete. This worker leaves final scientific interpretation to the parent.

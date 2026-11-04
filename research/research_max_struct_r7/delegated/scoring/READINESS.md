# Independent R7 scorer readiness

Status: code implemented, syntax parsed; primary GT scoring NOT_RUN. Awaiting parent message `ALL_PRIMARY_SEALED`.

No R7 sealed_truth file or private mapping was opened. Generator schema/source inspection exposed the deterministic family ordering to this scorer before sealing; this is disclosed rather than claiming scorer family blindness. The scorer has no inverse-provider or optimizer role. All truth arrays remain unopened until explicit authorization.

## Public source audit

- Public source count is variable: 48 receivers ×6 or12 sources. Scoring validates against the actual public Y shape, not the inherited fixed48×12 assertion.
- Final prediction is the paid original-Phi incumbent, with576 float64 materials and already-paid Y_prediction; scoring performs no physical solve.
- Budget materials are float32 incumbents. The four cost columns are CPU including process prefix, actual physical LU, actual complex RHS, original Phi. Last saved checkpoint <= shared budget is used without interpolation or GT iterate selection.
- `record` may be triggered on cached physical states as well as newly charged states; ten record events do not mean ten new LU. Actual recorded counters are retained and checked monotone. Equal LU/RHS excludes auxiliary FFT/readout work; actual CPU includes it.
- Direct physical max6000 iterations and joint corrector max3000 plus auxiliary stages are not equal computation. Endpoint results are therefore accompanied by matched CPU, LU and RHS material scores.
- A paid incumbent from an auxiliary stage may lack a saved gradient; projected-gradient availability is optional and must not be fabricated.
- Same-grid complex-Gaussian noise discrepancy and material success are separate fields. Six36-to24 mismatched cases retain raw noise-only discrepancy without interpreting sigma as a full error budget.
- The30 randomized object realizations, not180 runs or extra acquisition labels, define bootstrap units. Arms are paired on objects. Continuous complete-pair counts remain explicit if failures occur; failure denominators are retained for success rates.

## Commands after explicit parent authorization only

From repository root:

```sh
python research/research_max_struct_r7/experiments/score_primary.py seal --manifest research/research_max_struct_r7/configs/primary_final_manifest.json --output research/research_max_struct_r7/RESULTS/scoring/PRIMARY_PUBLIC_SEAL.json
python research/research_max_struct_r7/experiments/score_primary.py score --seal research/research_max_struct_r7/RESULTS/scoring/PRIMARY_PUBLIC_SEAL.json --output research/research_max_struct_r7/RESULTS/scoring/primary --authorization ALL_PRIMARY_SEALED
python research/research_max_struct_r7/experiments/analyse_primary.py --scores research/research_max_struct_r7/RESULTS/scoring/primary --output research/research_max_struct_r7/RESULTS/scoring/analysis
python research/research_max_struct_r7/experiments/score_rational.py --primary-seal research/research_max_struct_r7/RESULTS/scoring/PRIMARY_PUBLIC_SEAL.json --authorization ALL_PRIMARY_SEALED --output research/research_max_struct_r7/RESULTS/scoring/rational
```

Outputs are write-once. Seal checks use file size and modification time; no hashes. Scorer/analysis separate receipts report wall/CPU and zero new physical calls. All30×6 maps appear in six fixed groups of five; continuous material distributions, paired bootstrap95% intervals, condition/family summaries and matched-budget material scores are written.

Budget-matching caveat: auxiliary-stage physical validations are sparse; common CPU checkpoints can leave substantial unused budget. The table reports that gap. No interpolation or unrecorded intermediate auxiliary material is scored.

## Added post-seal diagnostics

- Primary sealing now includes every explicit per-task config, including the frozen Aplus config.
- Same-grid24 GT feasible-Phi witnesses reuse stored m/Y_clean/noise and the exact frozen smooth/TV prior; six cross-grid witnesses are NOT_AVAILABLE_DIFFERENT_OPERATOR. No forward is added.
- Pairwise A-vs-globalDual, Aplus-vs-A and secondary A-vs-sourcewiseDual cost frontiers use only the overlapping observed CPU/LU/RHS ranges, without stopped-method extrapolation. All-six common-budget results remain separately labeled.
- Extra scorer defaults to47 fixed endpoints: init18 repeated-object diagnostics, olddigit9 regression, and20 disclosed Aplus-development trajectories. It excludes primary180 task IDs, deduplicates shared manifests by output directory, never enlarges n30, and labels FFT-vs-old-dense comparisons as confounded.

After explicit parent authorization, primary scoring and analysis, run:

```sh
python research/research_max_struct_r7/experiments/score_extra.py --primary-seal research/research_max_struct_r7/RESULTS/scoring/PRIMARY_PUBLIC_SEAL.json --primary-scores research/research_max_struct_r7/RESULTS/scoring/primary --authorization ALL_PRIMARY_SEALED --output research/research_max_struct_r7/RESULTS/scoring/extra
```

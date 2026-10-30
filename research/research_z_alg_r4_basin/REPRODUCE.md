# Reproduction

Use the existing `migodam/b1-probabilistic-maxwell-imaging` repository. R4 imports the frozen R3-DIGITS physics and generic material regularizer directly. Python packages: NumPy, SciPy, Pillow, matplotlib, threadpoolctl; scikit-image is optional for the inherited SSIM metric. No GPU, NN weights, CUDA or additional service is needed.

```bash
git clone https://github.com/migodam/b1-probabilistic-maxwell-imaging.git
cd b1-probabilistic-maxwell-imaging/research/research_z_alg_r4_basin
python -m venv .venv
source .venv/bin/activate
pip install numpy scipy pillow matplotlib threadpoolctl scikit-image
```

## One fresh observation to one material image

The command receives an anonymous measurement and public operators. It never imports the GT generator or scorer.

```bash
python run_solver.py \
  --public public_data/fresh/fresh_37c60bb07853.npz \
  --operators ../research_z_alg_r3_digits/public_data/operators.npz \
  --method uniform_bridge --wall-seconds 60 --max-lu 2500 --max-rhs 60000 \
  --output /tmp/r4_uniform_bridge.json
```

Other final methods: `trf`, `data_lbfgs`, `globaldual`, `uniform_best_only`, `uniform_multistart`. All use the same final frozen defaults and physical model. Results include online top1, finite candidate archive, complete material trajectories, failed/truncated proposals and actual operation counts. A finite numerical result does not imply a data solution.

## Re-run the complete final confirmation

```bash
python run_experiments.py \
  --cases fresh_37c60bb07853 fresh_b840e9363b40 fresh_680c54b8f3ea \
  --methods trf data_lbfgs globaldual uniform_bridge uniform_best_only uniform_multistart \
  --wall-seconds 60 --max-lu 2500 --max-rhs 60000 --workers 2 \
  --run-name fresh_v3_reproduction
```

The original frozen predictions are never overwritten. Floating timing and resource-prefix endpoints can differ across hardware. This reproduction uses the released measurements, not a newly generated favorable target. New generator entropy is intentionally not published.

## Exact archived-score and figure replay

Download the existing-repository release [R4 full reproduction archive](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r4-basin-20261030), and extract it onto the repository root. It supplies all raw `results/` JSON and trajectory NPZ files omitted from the readable Git snapshot because of their size. Compact endpoint records are also readable in `published_results/`.

```bash
python score_fresh.py
python aggregate_costs.py
python prepare_public_artifacts.py
```

These commands read already-disclosed material arrays and saved predictions. They perform zero Maxwell solves. `score_fresh.py` regenerates the three-target/six-method figure and resource-prefix curves from the archived original predictions. Scoring is retrospective; it never changes online top1 or injects truth into an inverse run.

## Exposed-history mechanisms

```bash
python src/diagnostics.py
python src/step_qualification.py
python src/hessian_geometry.py
```

These are explicit offline development experiments, not blind confirmation. Their historical input predictions are in the earlier R3-DIGITS release. Geometry spectrum and certificate routines have separate costs. The Hessian script preserves its capped attempts; a dense-face followup is recorded in its JSON and source, not an invitation to increase precision or budget.

`prepare_cases.py fresh` creates a new random three-target study; it does not reproduce the original anonymous targets. Do not replace archived reported measurements with its output. The seedless truth arrays and public observations are sufficient to replay all reported metrics.

# Post-scoring mechanism analyzer: implemented, not executed

Dedicated new source: `experiments/analyse_warm_mechanisms.py`. Existing analyzers were not modified. Only public receipts and existing public diagnostic schemas were inspected while implementing this script. It has received a syntax-only compile check; it has not imported study scores, read test truth, or executed confirmation analysis.

After root explicitly declares readiness, run from the R8 worktree:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r8/experiments/analyse_warm_mechanisms.py --authorized-post-scoring
```

If root's all42 public projection diagnostics exist, append `--projection-diagnostics` followed by their exact `DIAGNOSTICS.json` path. The sibling receipt must declare `truth_access: NONE`. The script does not recompute these diagnostics and does not open validation-only diagnostic files by default.

The analyzer requires exact unchanged public file-stat seals for750 original and126 anchor tasks plus existing original750 scoring. Optional scored126 anchor results are included only with their completed scoring receipt. It opens no private truth path. Material metrics, including initial metrics, are inherited from existing scoring rows. Initial/pre-corrector Phi and selected stage come from public paid receipts/JSON. Saved projected-gradient algebra, when needed, uses only public prediction m/gradient and the inherited box; no LU or physics.

Outputs: all42 one-row-per-case wide table; original four-arm NN/uniform and selected anchor/direct/joint long table; separate24/12/6 cohort survival counts; all42 fixed factorial signed interactions; same-arm anchor and global nn_joint paired differences; direct extra-CPU summaries only within hardware; report; measured outer CPU/wall receipt with zero LU/RHS/new inverse solves.

Interpretation limits: stage−1 means the initial candidate entered correction, not an unchanged final image. The pre-Phi comparison explicitly distinguishes own initial Phi and NNdirect pre-corrector Phi. A lower pre-corrector Phi followed by worse final Phi/material is a trajectory observation, not proof of a basin/barrier/uniqueness claim. Zero is the sign of a recorded difference, not a tuned success threshold. No failed metric is replaced by zero. Global nn_joint is the strict truncated-physical-current control; Dual/sourcewise differences still combine current initialization and affine-space changes. Public projection covariates cannot establish causality from correlation and are not turned into an online selector.

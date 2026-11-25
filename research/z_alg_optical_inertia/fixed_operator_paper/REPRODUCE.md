# Reproduce the material and illumination tests

The complete [research archive](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/fixed-operator-material-sets-20261126/fixed_operator_reproduction.zip) contains the raw readings, material classes, source proposals, full-wave witnesses, one-sided calculations, failed attempts and the earlier operator files used by this study. Extract into a fresh directory, retaining the folder tree. The smaller IP_Anonymous_LaTeX.zip builds the article by itself. The online repository presents the article, code and compact results. Large replay arrays are in the research archive.

Calculations ran on an Apple M4 using Python 3.13.13. Numerical versions are in requirements.txt. All main stages ran on CPU. A and the moving-support comparison use single-threaded BLAS. Native source calculations record their thread counts in the result files. Wall-clock sums of concurrent jobs and CPU time have different meanings.

~~~sh
python3 -m venv .venv
. .venv/bin/activate
P=research/z_alg_optical_inertia/fixed_operator_paper
python -m pip install -r "$P/requirements.txt"
export OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1
python -B "$P/src/replay_package.py"
python -B "$P/src/p0.py"
python -B "$P/src/replay_native.py"
python -B "$P/src/release_p0.py"
python -B "$P/protocol/selfcheck.py"
python -B "$P/protocol/chart_check.py"
~~~

P0 includes all six supplied mechanism programs, full and reduced inertia, singular and null-space examples, coverage, retention of truth controls and native interval checks. The native disk had completed when the blind acquisition was sealed. The original matching record is results/p0_at_acquisition.json. The final P0 record adds the completed ring. Existing blind measurements can be replayed without knowing the hidden generator seed.

| Question | Computation entry and records |
|---|---|
| All 2688 data–class decisions | experiments/A/README.md. Use the saved acquisition. Run replay.py --begin 1 --end 24 and repeat with --audit-inertia. Score only after decisions are sealed. The first failed implementation's results are retained. |
| Cold-start convex optimization and positive source tails | experiments/A/psd_followup/README.md, summary.json, positive_tail_summary.json, positive_tail_upper_audit.json. Includes trial weights and one-sided proof frames. |
| Core location, single-cell peaks, 100/1000 moving boxes | experiments/B/README.md. run.py online reads public geometry and readings. replay.py checks all excluding map nodes and witnesses, and full-matrix reference subsamples. benchmark.py measures both implementations. |
| Prospective acquisition and paired sensitivity | experiments/C/README_ZH.md. Fresh replay starts with register.py and acquire.py old in a working copy. design.py writes proposal.json before acquire.py future. readout.py and menu_replay.py check the decisions. Repeat with register_diagnostic.py and the documented setting. |
| Receiver SOCP and indefinite forms | experiments/C/receive_socp/README_ZH.md and experiments/C/indefinite/README.md. Both are post-study checks of the old readings. |
| Source-side index, energy and coherent illumination | Run the native scripts below in order. source_index.json, energy.json, coherent_control.json and occupied_wide_control.json contain the signed spectral results, residuals, dimensions and costs. |
| Physical background compensation and a hypothetical beam | results/mechanisms/background_compensation_exact.json contains rational sign-changing brackets and a nonvanishing-denominator proof. experiments/virtual/verify.py and followup.py check an old-only intersection and physical negative controls. |
| Discretization, gain error and laboratory margins | experiments/D/error_audit.py, audit.json and ports.csv. These reproduce the discrete example and tabulate measurements from two facilities. Loaded-instrument error remains unknown. |
| Inputs, scope and independent review | SOURCE_AUDIT.txt, METHOD_LIMITS.md, protocol/API.md and reviews/RESPONSE_ZH.md. RUN_MANIFEST.json indexes evidence and costs. |

~~~sh
python -B "$P/experiments/native/source_index.py"
python -B "$P/experiments/native/energy.py" --precision 192
python -B "$P/experiments/native/compatible_receiver.py"
python -B "$P/experiments/native/coherent_control.py"
python -B "$P/experiments/native/plot.py"
python -B "$P/experiments/native/plot_illumination.py"
python -B "$P/src/plot_strengthened_phase.py"
python -B "$P/experiments/D/error_audit.py"
python -B "$P/experiments/virtual/verify.py"
python -B "$P/experiments/virtual/followup.py"
~~~

Source-index generation can resume from its existing prefix results. For a cold recomputation use a separate extracted tree and move source_index.json out of that tree's native experiment folder before starting. The 1153-state dense restriction uses exact exterior-vacuum equations on the 9216-cell stored operator. Older studies' files in the archive are read-only inputs.

## Article

Unzip IP_Anonymous_LaTeX.zip into an empty directory. Run pdflatex on manuscript_ip.tex until cross-references settle, or run tectonic manuscript_ip.tex. The IOP class and the figures actually cited are included. The distributed PDF was produced by pdfLaTeX. With both engines on the path, python src/package_manuscript.py repeats the independent builds from the manuscript archive. The result and page count are in results/typesetting.json.

Study-specific acquisition and proof interfaces are documented in protocol/API.md. Older source files keep their original paths. Each experiment records its finite search budgets, failed proposals and the information unavailable in older records.

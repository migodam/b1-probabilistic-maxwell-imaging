# Packaging implementation notice

Ready for parent review; packaging **NOT_RUN**. No archives created, arrays opened, physics, network, hashes, source deletion or active-job changes.

- `tools/package_release.py`: bounded metadata-only discovery (or supplied RELEASE_LAYOUT), retains all actual NPZ/NPY including late/failed/full-field arrays; PNGs included by layout/enumeration. Stage-relative ZIP members unzip directly into `research/research_w_slam/`. Groups E0E1, A0, controlled_E2E3, E4E5, inputs, auxiliary evidence. Streams bytes only during later execution. ZIP/member sizes and source paths recorded in `release_assets/ASSET_MANIFEST.json`; no scientific completeness claim.
- Root `.gitignore` excludes `/release_assets/`.
- `--mode final` requires terminal saved bank/cost/inference evidence, never scientific PASS. Failure-publication path: `--allow-not-run --status-manifest GATE_SUMMARY.json` accepts explicit matching terminal blocked/not-run/failed statuses or reasoned executed=false declarations. Missing arrays stay missing and are listed; no replacement solve or synthetic member.
- `REPRODUCE.md`: added exact cached-statistics commands, bank terminal-cost finalizer, paid independent bank validation, E2/E3/E4/E5 sequence, packaging and stage-root unzip commands. Removed extra approval phrasing; retained budget and exclusive-output-directory warnings, unmeasured historical cost limits and portability/licensing scope.

Source compiled in memory only, not executed. Suggested final wrapper: `python tools/run_job.py E0E1 release-package-001 python tools/package_release.py --mode final`. Review current terminal declarations first; partial mode remains labeled snapshot.

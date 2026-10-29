# R9 reproduction and portable evidence

This guide separates existing-evidence evaluation, archived prediction validation without truth, and separately charged scientific replication. Original frozen outputs and ledgers are never reset.

E2/E3 and weak-line REV1 were development on exposed R7 objects. E4 has limited conditional probability models. E5 is a separately approved32-parent possibility-volume campaign. Its before-truth seal and scoring receipt are authoritative; approval alone is not completion.

## Runtime and paths

Use the project virtual-environment executable directly. Do not resolve its symlink or replace it with the underlying base interpreter: that can lose the installed environment.

~~~sh
INV_ROOT="/Volumes/migodam's-external-brain/Research/Inv_SLAM"
PY="$INV_ROOT/Gaussian/.venv_nn/bin/python"
P="$INV_ROOT/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace"
~~~

Observed: Python3.13.13 on macOS arm64; NumPy2.5.3, SciPy1.18.1, Matplotlib3.11.1, threadpoolctl3.7.0. These are observed versions, not a new installation or a byte-identical cross-platform guarantee. Details: audit/reproducibility/DEPENDENCY_STATUS.json and PYTHON_REQUIREMENTS_OBSERVED.txt.

Numerical branches use one worker, BLAS1, float64/complex128. RHS256 and RHS576 remain different accounting units.

## Evidence entry points

| Branch | Authoritative output | Source/replay entry |
|---|---|---|
| E0 | audit/exploration/E0_SOURCE_AUDIT.json | exploration.py e0 |
| E1 | audit/exploration/E1_HISTORICAL_SUBSPACE_STABILITY.json | exploration.py e1 |
| Original E2/E3 | results/exploration/SUMMARY.json, E3_ANALYSIS.json, RESOURCE_LEDGER.jsonl | exploration.py; exploration_analysis.py |
| Separate REV1 | results/exploration_rev1/SUMMARY.json, REFERENCE_COMPARISON.json, RESOURCE_LEDGER.jsonl | exploration_rev1.py; exploration_rev1_analysis.py |
| Exact nonuniform ROI | results/stable_roi/RESULTS.json, SETUP.npz, COST_LEDGER.jsonl | stable_roi_operator.py |
| E4 | E4_PROBABILITY_REPORT_ZH.md; results/probability/STATUS.json, final_numerical_comparison_v3.json, original1d_results_v2.json, aggregate_costs.json | probability_transfer/integrate/envelopes/finalize/original_1d.py |
| E5 predictions | results/possibility_volume/PREDICTIONS_SEALED.json, parents, jobs, RESOURCE_LEDGER.jsonl | possibility_volume.py |
| E5 scores | results/possibility_volume/offline_statistics/SCORING_COMPLETE.json, SUMMARY.json | possibility_statistics.py |
| E5 scored plots | figures/e5, PLOT_PROVENANCE.json | plot_possibility_volume_scored.py |
| E5 galleries | figures/e5/gallery_*, GALLERY_POSTHOC_PROVENANCE.json | plot_possibility_gallery_offline.py |

Do not use the pre-E5 aggregate RESULTS_SUMMARY.json/RESOURCE_SUMMARY.json in place of E5's ledger and seal. Legacy assemble_evidence.py predates E5 and can emit older scope statements; do not rerun it over final combined reports without a reviewed reporting update.

## Dependencies and archive layout

audit/reproducibility/EXTERNAL_DEPENDENCY_INVENTORY.csv lists179 required/cited assets, original absolute paths, Gaussian-relative archive paths, usage and file metadata. The inventory did not open NPZ contents or private truth.

The existing raw shards are **repository-relative**, not Gaussian-relative. Extract the four E5 parent shards and one older-arrays shard over the same archived B1 repository checkout:

~~~text
<extracted repository root>/
  research/p3l_r9_stabilized_subspace/
    results/...
    audit/...
    private/possibility_volume/...
  src/...
~~~

Historical/source dependencies are a separate **Gaussian-relative** overlay:

~~~text
<dependency archive root>/
  Gaussian/
    B1/b1-p3l-r4-conditional-transport/...
    B1/b1-p3l-r6-task-optimal/...
    B1/b1-p3l-r7-material-fibers/...
    A30/a30-maxwell-material-manifold/...
~~~

The relocation helper combines these into a NEW runtime tree with Gaussian/B1/b1-p3l-r9-stabilized-subspace and its existing sibling names. It also accepts a previously assembled Gaussian-root archive. Do not extract repo-relative shards directly under Gaussian/ without the B1 repository prefix.

Only inventory-listed historical files/subtrees are required. These include R4 P1/core source, target/geometry and public weak Green/acquisition arrays; R7 source and exposed DEV anchors/paths; R6 controlled observations and paid quadrature imports; and A30 plus its A24/A20 and nested vendored A2 closure. Omitting vendor/a2/physics.py breaks the dynamically loaded historical solver.

The exact nonuniform ROI SETUP.npz and E4 scalar-family transfer_setup.npz have different scopes. They are not interchangeable. Formula-index source citations, including F/F3/A17X, are also listed.

The reading archive may include reports, atlas, figures, source, protocols, public score JSON and summary receipts. Full raw evidence also needs the listed dependencies, parent/trial arrays and ledgers. Proposed release paths are in audit/reproducibility/RELEASE_ARCHIVE_MAPPING.csv; they are plans, not assertions of publication.

Here **private/** means blinded synthetic experimental GT, not personal/private user data. After all32 prediction sealing and approved scoring, the full synthetic research archive includes E5 truth and generator states with role OFFLINE_ONLY. Their97-file metadata inventory is audit/reproducibility/OFFLINE_SYNTHETIC_EVIDENCE_INVENTORY.csv. Independent scoring replay must not require paid regeneration.

Online code isolation remains unchanged: the online engine never opens these files. Relocation includes them only with explicit --include-offline-synthetic-truth, after checking completed seal and scoring. Other unlisted blinded namespaces remain excluded. Virtual environments, Git internals and interpreter caches are excluded. Posthoc GT galleries are explanatory derived outputs with separate provenance.

## Transport-only relocation

Frozen protocols, seals and VISITED records contain original-host absolute paths. Layout preservation resolves relative imports, but textual paths still need relocation. probability_report.py also hardcodes the study path.

The helper reads completed archived evidence and writes a NEW directory. It never mutates archives, imports scientific runners or executes physics. Default copying skips blinded evidence. Explicit --include-offline-synthetic-truth permits byte-copying and metadata/digest recording for this declared E5 synthetic namespace only; it does not parse array values. This opt-in requires both the all32 prediction seal and completed scoring. Physics-enable files remain disabled.

~~~sh
REPLAY_PY="/path/to/a/compatible/venv/bin/python"
"$REPLAY_PY" "/path/to/extracted-repository/research/p3l_r9_stabilized_subspace/src/relocate_r9_replay.py" \
  --archive-root "/path/to/extracted-repository" \
  --dependency-root "/path/to/Gaussian-relative-dependencies" \
  --replay-root "/new/path/r9-replay" \
  --original-root "/Volumes/migodam's-external-brain/Research/Inv_SLAM"
~~~

The new tree records TRANSPORT_RELOCATION_MANIFEST.json and TRANSPORT_CHANGED_FILES.json, original/runtime content identities, prefix mapping and all changes. Runtime source identities and copied seal file paths/byte sizes are refreshed explicitly. Physics-enable files are disabled.

Scientific seeds, thresholds, model choices and original freeze times remain unchanged. The copied seal is a transport derivative of an existing before-truth seal, not a new scientific freeze. Original and runtime hashes intentionally differ where transport changes bytes. No original byte identity is falsely claimed.

A genuine completed-evidence replay has now passed in packaging/replay_audit/replay_attempt_1. The minimal repository-relative evidence and separate179-file dependency overlay were staged with APFS clones, relocated, publicly validated and rescored once. All seven canonical outputs and all four primary contrasts match exactly; new Maxwell RHS=0. See audit/reproducibility/RELOCATION_REPLAY_AUDIT.json/md. Originals remain read-only.

## Zero-physics evaluation in the NEW tree

~~~sh
RP="/new/path/r9-replay/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace"

"$REPLAY_PY" "$RP/src/possibility_statistics.py" validate_predictions \
  --output "$RP/results/possibility_volume" \
  --protocol "$RP/config/possibility_volume_protocol.json" \
  --statistics "$RP/config/possibility_statistics_frozen.json"

"$REPLAY_PY" "$RP/src/plot_possibility_volume_scored.py" \
  --statistics-dir "$RP/results/possibility_volume/offline_statistics" \
  --output "$RP/figures/e5_reproduced"
~~~

Validation reads public sealed predictions without GT. Plotting reads completed public score summaries only and makes no endpoint/model selection. Finite32-parent coverage, its uncertainty and approximate ensemble actions remain separate from posterior accuracy.

Historical zero-RHS analysis also belongs in a fresh derivative because it writes reports:

~~~sh
"$REPLAY_PY" "$RP/src/exploration.py" summarize
"$REPLAY_PY" "$RP/src/exploration_analysis.py"
"$REPLAY_PY" "$RP/src/exploration_rev1.py" summarize
"$REPLAY_PY" "$RP/src/exploration_rev1_analysis.py"
"$REPLAY_PY" "$RP/src/probability_qualification_v3.py"
"$REPLAY_PY" "$RP/src/probability_cost_from_existing.py"
"$REPLAY_PY" "$RP/src/probability_report.py"
~~~

REV1 analysis appends a reference section; run once in a fresh copy when comparing bytes. E0/e1 can be recomputed from bundled exposed arrays, but rewrite their audit products.

The separate gallery helper requires SCORING_COMPLETE and matching explicit offline approval. With the completed synthetic truth archive explicitly included, it can be reproduced without new physics or regeneration. Invoke it with --posthoc-approved only under the recorded posthoc authorization. First-parent/worst-failure choice is explanatory posthoc selection, never online selection.

## Independent score replay from archived synthetic GT (zero new physics)

Create a SECOND new derivative for a clean deterministic scoring replay. The original archive and original scoring stay intact:

~~~sh
"$REPLAY_PY" "/path/to/extracted-repository/research/p3l_r9_stabilized_subspace/src/relocate_r9_replay.py" \
  --archive-root "/path/to/extracted-repository" \
  --dependency-root "/path/to/Gaussian-relative-dependencies" \
  --replay-root "/new/path/r9-score-replay" \
  --original-root "/Volumes/migodam's-external-brain/Research/Inv_SLAM" \
  --include-offline-synthetic-truth \
  --prepare-score-replay

SP="/new/path/r9-score-replay/Gaussian/B1/b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace"
"$REPLAY_PY" "$SP/src/possibility_statistics.py" validate_predictions \
  --output "$SP/results/possibility_volume" \
  --protocol "$SP/config/possibility_volume_protocol.json" \
  --statistics "$SP/config/possibility_statistics_frozen.json"

"$REPLAY_PY" "$SP/src/possibility_statistics.py" score \
  --output "$SP/results/possibility_volume" \
  --protocol "$SP/config/possibility_volume_protocol.json" \
  --statistics "$SP/config/possibility_statistics_frozen.json" \
  --root-approval "$SP/config/possibility_volume_offline_score_REPLAY_approval.json"
~~~

The helper renames only the copied derivative's offline_statistics to offline_statistics_ARCHIVED_ORIGINAL_SCORING, preserving that copy too. It never deletes original scoring or original files. The ordinary score entry then creates a fresh offline_statistics directory. The new replay-approval record points to preserved original approval and records explicit transport/replay opt-in under the standing synthetic-evidence authorization; it is not a new model/protocol approval.

OFFLINE_ONLY truth files are byte-identical copies with archived/runtime identities. Online predictions, seeds, gates and model choices remain unchanged. Comparison with archived original scoring is numerical evaluation, not new physics, new calibration, policy tuning or new generalization data. Plot the fresh scores with the same scored-summary plotter. Physics remains disabled throughout.

## Separately charged scientific replication

Existing jobs and receipts cannot be automatically replayed. Root must approve a NEW output/private namespace, versioned replication protocol with runtime identities, and explicit enablement. Reusing the original seeds replicates the same32-parent experiment; it is not a new independent generalization sample.

~~~sh
"$REPLAY_PY" "$RP/src/possibility_volume_generate.py" --config "/path/to/ROOT_APPROVED_REPLICATION_PROTOCOL.json"
"$REPLAY_PY" "$RP/src/possibility_volume.py" all --config "/path/to/ROOT_APPROVED_REPLICATION_PROTOCOL.json"
"$REPLAY_PY" "$RP/src/possibility_volume.py" seal --config "/path/to/ROOT_APPROVED_REPLICATION_PROTOCOL.json"
~~~

Only after all32 seal does root issue a new offline-scoring approval:

~~~sh
"$REPLAY_PY" "$RP/src/possibility_statistics.py" score \
  --output "/new/approved/output/namespace" \
  --protocol "/path/to/ROOT_APPROVED_REPLICATION_PROTOCOL.json" \
  --statistics "$RP/config/possibility_statistics_frozen.json" \
  --root-approval "/path/to/NEW_AFTER_SEAL_SCORE_APPROVAL.json"
~~~

Statistical expected identities must match the declared contract; a changed scientific contract requires its own pre-data freeze. A path-only derivative keeps its relationship to the original explicit.

E5 worst bound:352512 RHS256,11296 LU256; caps400000/12000/7200 seconds. Cached setup264 RHS576 is historical and not recharged. Rebuilding it is new charged work.

Historical physical entry points, only for separately authorized replications:

~~~sh
"$REPLAY_PY" "$RP/src/exploration.py" costprobe
"$REPLAY_PY" "$RP/src/exploration.py" all
"$REPLAY_PY" "$RP/src/exploration_rev1.py" all
"$REPLAY_PY" "$RP/src/stable_roi_operator.py" run
"$REPLAY_PY" "$RP/src/probability_transfer.py" probe
"$REPLAY_PY" "$RP/src/probability_integrate.py"
"$REPLAY_PY" "$RP/src/probability_envelopes.py"
"$REPLAY_PY" "$RP/src/probability_finalize.py"
"$REPLAY_PY" "$RP/src/probability_original_1d.py"
~~~

These retain different branch-specific frozen budgets, ledgers and no-replay guards. Clearing a receipt is not replication. Preserve original failures, underflow/reporting repairs, E4 qualification revisions, and E5's pre-ledger duplicate-role failure plus narrow R2 fix.

No SHA check, extra micro test, posterior upgrade or continuous-fiber assertion is implied.

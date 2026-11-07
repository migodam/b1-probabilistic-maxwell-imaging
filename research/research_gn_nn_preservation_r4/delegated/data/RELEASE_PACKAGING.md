# Finite Release packaging entrypoint

`experiments/package_release.py` prepares assets for a GitHub Release on the **existing B1 repository**. It does not create a repository, upload, install anything, start a service, decode GT arrays, run physics, or score TEST. This delivery ran syntax compilation and metadata-only planning; no ZIP was created. Final reports/plots and actual publication remain the parent's work.

## Plan now, pack only after closure

From the R4 study directory, using the established Python:

```
python experiments/package_release.py plan
```

This writes `audit/RELEASE_PACKAGE_PLAN.json` using file names, sizes, exact nanosecond mtimes, public ID metadata and an AST-only Python import closure. It is a live inventory snapshot, not a claim that currently running work is frozen. Missing work is explicit NOT_RUN/missing metadata; empty/absent groups do not receive fabricated ZIPs. The plan records unrun sealed cases and current primary/traditional-control coverage. It does not read GT array values or evaluation scores.

Only the explicit `pack` action creates ZIPs. It requires `audit/METHOD_FREEZE.json` and a valid `evaluation/TEST_SEALED.json` with `ALL_PLANNED_PREDICTIONS_SEALED` and split TEST. It checks the seal's existing size/mtime stamps without hashes, then checks every inventoried source size/mtime before and after byte streaming. Changed files cause a preserved partial attempt, not a silent new scientific run. Call after relevant writers finish and after final outputs are ready.

The parent supplies actual verified repository/ref/old-Release values; no URL is invented:

```
python experiments/package_release.py pack \
  --output /absolute/path/to/fresh-release-directory \
  --repository-url "$B1_REPOSITORY_URL" \
  --publication-ref "$B1_RELEASE_REF" \
  --historical-e0-release-url "$HISTORICAL_E0_RELEASE_URL"
```

The output directory must be fresh/empty and should sit outside the scanned data/model/output subtrees. Existing assets and failed partial attempts are never overwritten or deleted. There is no automatic upload or automatic retry. A failed attempt leaves `.zip.partial` evidence and `PACKAGING_RECEIPT.json`; retry into a new directory after inspecting the cause.

## Functional assets and limits

All archive paths preserve the existing repository-relative layout.

| ZIP name prefix | Actual included payload |
|---|---|
| `gn-nn-r4-data-public` | R4 operators, reference Jacobian/receipt, all public observations, online manifest and generation receipt |
| `gn-nn-r4-data-evaluator` | Actual R4 truth NPZs and evaluator manifest, copied as opaque bytes only after seal gates |
| `gn-nn-r4-anchors` | Primary anchors, old/new traditional controls, misconfigured-attempt evidence and existing batch receipts |
| `gn-nn-r4-models` | All actual model-directory weights, configs, histories and summaries, including git-ignored `best.pt` |
| `gn-nn-r4-predictions` | Existing prediction NPZ/JSON files, including supplementary/failed work |
| `gn-nn-r4-evaluation` | Existing canonical evaluation NPZ/JSON, seals/manifests and batch receipts |
| `gn-nn-r4-fixed-projectors` | Study projector cache/receipts plus the explicitly referenced old `fixed_spectra.npz` |
| `gn-nn-r4-analysis-tables` | Actual E0/test_energy and other saved result CSV/TSV/JSON/NPZ, including NOT_RUN tables and energy receipts |
| `gn-nn-r4-canonical-receipts` | Actual audit/delegated JSON/JSONL/CSV/log/error/text receipts; transport ZIPs and temporary files are excluded |
| `gn-nn-r4-frozen-python-dependencies` | Minimal AST-referenced old R3/R2/R4 Python closure, R2's explicit importlib physics file and historical E0 seal/plan metadata |
| `gn-nn-r4-atlas` | Optional existing atlas directory, only with `--include-atlas --atlas-dir figures/atlas` |

Each created asset is named `<prefix>.part001.zip`, etc. Default maximum is **1,500,000,000 bytes**, a decimal1.5GB. ZIP_STORED, fixed1980 ZIP metadata, sorted paths and exact header/payload accounting make partitioning deterministic. Members are streamed in8MiB blocks. If one file itself exceeds an asset, deterministic opaque chunks are indexed with original path, offset and length. No fake empty archive is created. A completed asset's actual size must equal its plan and stay within the requested cap.

`gn-nn-r4-release-index.json` lists actual created asset sizes, member/logical-file counts, every source size/mtime, chunk reconstruction map, provenance, missing cases and recorded publication ref. `PACKAGING_RECEIPT.json` records elapsed packaging work and preserved failures. Neither claims that an asset was uploaded or publicly accessible. Publish only CREATED ZIPs and the final index through the existing repository's Release; the parent performs access verification independently.

R4 source, reports, Markdown instructions, vendor source/provenance and ordinary figures remain normal git content. The optional large atlas can be an asset. No virtual environment, site-packages directory, credentials, unrelated private tree, or old-study directory-wide sweep is included. Old Python dependencies are resolved from source AST without importing or running them.

Historical E0 CSVs, receipts and frozen basis are included. Rerunning the original historical E0 additionally needs the old R3 dataset/evaluation arrays from its **existing public study**. Those old raw arrays are deliberately not silently swept into this package. Supply the verified direct existing Release page/asset URL using `--historical-e0-release-url`; an absent link stays `NOT_PROVIDED_PARENT_MUST_SUPPLY_VERIFIED_DIRECT_RELEASE_LINK`. This script makes no network request to attest access.

## Restoring bytes and preserving provenance

Place the downloaded index and every indexed ZIP together. The source entrypoint is retained in normal git. Restore to an empty checkout-layout staging directory:

```
python experiments/package_release.py restore \
  --index /downloaded/assets/gn-nn-r4-release-index.json \
  --target /absolute/path/to/restored-repository-layout
```

It reconstructs original files/chunks without decoding arrays. It preserves exact recorded `mtime_ns` with `os.utime`; ZIP timestamps alone have insufficient precision. Filesystem timestamp rounding/failures appear in a separate restore receipt. It never rewrites any seal. The source archive's original seals retain their original absolute machine paths as historical provenance.

By default, existing destination files are preserved and cause a hold. For a clone containing identical source/receipt files, explicit `--reuse-identical` allows a streaming exact-byte comparison and reuse without digests; differing existing bytes still hold and are never overwritten. Empty staging is the simplest way to preserve both versions when git and asset metadata differ.

**An archived original seal is not a new-machine execution seal.** Data/evaluation manifest public/truth paths are relative, but the original sealed provenance paths may be absolute. Do not edit old seal paths or timestamps to pretend continuity. Either perform a separately documented fresh public inference/correction replay, or create a separate new local seal from the restored canonical pairs before offline scoring. For the latter, from the restored R4 study directory, use the exact frozen run names from FINAL_RUNS:

```bash
python - <<'PY'
import json, subprocess, sys
runs = json.load(open('audit/FINAL_RUNS.json'))['runs']
seal = 'evaluation/TEST_RESTORED_LOCAL_SEALED.json'
subprocess.run([sys.executable, 'experiments/score_predictions.py', 'seal',
                '--split', 'test', '--runs', *runs, '--seal', seal], check=True)
subprocess.run([sys.executable, 'experiments/score_predictions.py', 'score',
                '--split', 'test', '--seal', seal,
                '--output', 'results/RESTORED_TEST_SCORES.json'], check=True)
PY
```

These are future reproduction instructions, not work performed by the packager. The original TEST seal and scientific outcomes remain unchanged. Missing planned cases must retain NOT_RUN in the new local seal as well.

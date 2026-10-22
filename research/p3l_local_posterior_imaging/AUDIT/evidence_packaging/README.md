# P3-L existing-B1 Release evidence packaging

This tool only packages existing evidence for the existing B1 repository. It
does not upload or publish a Release, create a repository, calculate hashes,
interpret arrays, alter source artifacts, or make scientific judgments.

From the P3-L worktree, metadata-only size planning is:

```
../../.venv_nn/bin/python tools/package_p3l_evidence.py
```

After all paid jobs close, explicitly create or resume the package:

```
../../.venv_nn/bin/python tools/package_p3l_evidence.py --finalize
```

Finalize requires the parent's final `CAMPAIGN_CLOSURE.json` before freezing
the source inventory, and refuses active resource reservations or held worker/job locks. It holds
the worker lock while streaming. It checks source path, size and modification
time before/after reads and on resume. Existing completed parts are reused only
when their registered path and compressed byte size agree. This is a resumable
path/size inventory, not a cryptographic content-integrity certificate.

Independent tar.gz parts are limited to 700 MiB each, with a conservative
680 MiB uncompressed tar budget. Files exceeding that payload budget are split
into 64 MiB fragments inside archives. Every fragment maps to its original
root-relative path and byte offset through INVENTORY.jsonl and the fixed
partition rule in the included source script. The complete append-only
RESOURCE_LEDGER.jsonl is included verbatim, without parsing it into memory.

Keep MANIFEST.json, INVENTORY.jsonl, this README, and every part together under
`research/p3l_local_posterior_imaging/AUDIT/evidence_packaging/` in an evidence
bundle that includes the source script. Reconstruct all original logical files
into a new directory with:

```
../../.venv_nn/bin/python tools/package_p3l_evidence.py --restore /absolute/new/restored-evidence
```

For a standalone downloaded Release, place this source script beside the
downloaded manifest/inventory/parts (or fetch it from the existing repository),
then use ordinary Python without an installed project:

```
python3 package_p3l_evidence.py --bundle-dir /absolute/downloaded-parts --restore /absolute/new/restored-evidence
```

Restoration refuses an existing destination and checks archive member paths,
sizes, regular-file types and fragment offsets. It never extracts symlinks.
Gzip/tar decoding supplies ordinary format checks; no SHA or other hashes are
computed. Paths are relative to the P3-L worktree, including the campaign prefix
and project-owned src/tools/configs snapshots. External dependency trees are
excluded, and all included symlinks are refused.

Included: campaign raw NPZs, weights/checkpoints, every configuration and saved
summary, source snapshots, append-only journals and diagnostic logs. Excluded:
transient `.lock`, `.tmp`, `.pyc`, `.DS_Store`, `__pycache__`, `.git`, `.venv`,
`node_modules`, and generated packaging outputs other than this README. Logs
are retained because they may contain failure evidence. No parent reports,
frozen estimator code or source arrays are modified by packaging.

The metadata-only plan reports logical bytes and a conservative additional-disk
projection. Existing compressed NPZs/weights may compress little; reserve roughly
the full inventoried logical size again, plus archive overhead. Final measured
compressed bytes appear in MANIFEST.json. Memory use is bounded to one archive
partition's metadata and small streaming buffers, not the artifact byte volume.

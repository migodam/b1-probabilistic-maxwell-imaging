# Full raw evidence restore

These manifests describe ORIGINAL full-field files, all candidates/checkpoints and queue receipts, not export_scoring material-only derivatives. Plans contain no release URLs until actual upload/readback.

Text JSON/CSV/MD/source/log uses lossless ZIP_DEFLATED level1; NPZ/NPY/PNG/PDF and other binary assets use ZIP_STORED. All original file bytes restore without numerical transformation or lossy reduction. Manifest records each file compression and restore checks it. Volumes target600MB by raw input size; a single larger text file occupies its own compressed volume, whose actual ZIP must still stay≤650MB or execution fails. No SHA operations.

After publication download every listed part plus PACKAGE_PLAN.json to one external-disk folder. Run:

```text
python package_evidence.py restore --plan /external/PACKAGE_PLAN.json --archive-dir /external/downloads --destination /external/restored_r14
```

Restore refuses an existing destination, checks exact member path sets and byte sizes, and opens each ZIP member as a stream. Python zipfile performs ordinary ZIP CRC checks; no scientific rerun, SHA, optimizer, scorer or pickle is used. Repeated config files across source/raw packages must compare byte-identical without hashes. Paths reconstruct RESULTS/<campaign>/..., queues/<campaign>/..., configs/..., data_public/... and original source layout. Missing parts are fatal; originals are never deleted.

To produce actual archives later:

```text
python package_evidence.py execute --root /external/r14 --plan /external/PACKAGE_PLAN.json --output /external/release_shards
```

Before execution terminal campaign status and source sizes/mtime must still match the plan; stale plans must be regenerated. Successful creation reports archive actual sizes and per-file membership. Upload, anonymous access verification, data closure and scientific acceptance are parent tasks. Never substitute *_scoring_only.zip for these full raw archives.

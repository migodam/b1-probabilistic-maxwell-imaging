# Post-seal stationarity diagnostic: source prepared, not executed

New source: `src/analyze_stationarity.py`. This worker performed syntax compilation only. No public numeric arrays, primary materials, gradients, fresh quality or GT were read/evaluated, and no diagnostic physics was executed.

```text
python src/analyze_stationarity.py \
  --root /ABSOLUTE/PATH/TO/research_max_struct_r13 \
  --authorization /ABSOLUTE/PATH/TO/audit/JOINT_OFFLINE_AUTHORIZATION.json \
  --legacy-repo /ABSOLUTE/PATH/TO/FROZEN_LEGACY_REPOSITORY \
  --out /ABSOLUTE/PATH/TO/NEW_ANALYSIS_DIRECTORY
```

## Exact joint authorization guard

The script reuses the parent's existing fields: `authorized_joint_offline_scoring=true`, `required_campaigns`, `score_campaigns`, `public_selection_files`, and optional per-campaign `seals` mapping. It validates every required authorized seal, campaign/root correspondence and each frozen public-selection file's bytes/mtime before importing NumPy or the Maxwell model or reading any material/operator array. The required joint set includes core48, mechanism16, physics_v1, physics_recovery_v1, online controls16, rescue48, local reference16 and local TRF16. It does not authorize anything or silently accept a subset.

Only after that guard does it import the requested `--legacy-repo` MaxwellModel/Regularizer, checking the imported module's real path. Output must be a new analysis directory outside all frozen prediction/queue and private truth directories. Backslash relative paths are normalized at stored-path consumers. No SHA check is introduced.

## Fixed denominator and measurements

The requested scope is fixed before quality access:

- core48 × `{cie,dual_cie,g19,dr,event_integer,direct}`:288 rows;
- the same predeclared every-third-ID16 subset in local reference × `{cie,dual_cie,g19,dr,direct,physical_globaldual}`:96 rows;
- the local TRF16 fresh subset only:16 rows; exposed B is excluded.

All400 requested rows remain, including missing, ambiguous, invalid or failed points. NN variants are outside this declared scope. Case/material selection is never based on GT, residual quality or an observed gradient.

At each frozen primary material m, using its original public G/E/S/Y, box and exact scalar Regularizer parameters, the script computes the **original unscaled Phi** and gradient, a separate original prior gradient gR, data gradient gD=g−gR, original data residual, g/gD/gR infinity/L2 norms, and the gD/gR cosine. A zero norm makes the cosine null, not zero. It records the exact unit-step box gradient mapping

`m - clip(m - grad Phi(m), lower, upper)`

and its infinity/L2 norms. It does not use the optimizer's1e6 scaling, optimize a new point, compute a Hessian, clip/replace the frozen material, or create a candidate. Exact bound activity and actual box violations are reported without a near-bound numerical pass tolerance. `computed Phi - saved primary Phi` is a diagnostic difference only, without an error threshold.

These measurements are **not a stationarity pass gate, local-minimum certificate, recovery certificate or independent physical validation**. They address the distinction between budget/work truncation and the local behavior of the original finite objective at already selected points. Their scientific interpretation remains with the parent.

## Deduplication and separate analysis cost

Within each case, exact identical material arrays under exact identical public G/E/S/Y/grid, box and prior share one measured evaluation, including across campaigns. Each alias retains its own run/candidate/saved-Phi reference. No approximate deduplication or material-quality sorting is used. Public arrays/models/material caches are discarded between cases.

The existing Maxwell forward and analytic adjoint are charged once per unique point. `UNIQUE_FORWARD_COSTS.json/csv` records actual model counter deltas, LU, source/adjoint/aggregate complex RHS, CPU/wall and every referencing run. Failed evaluations retain their spent counters/errors and are not automatically retried. Setup/imports/reads are covered by the total receipt; model setup also has its own table and is not added twice to per-evaluation deltas. Missing counters remain null, with known sums labeled as such.

This uses the same finite public physical discretization on the parent-declared Apple M4 host; platform metadata and the real imported model path are recorded. It is additional post-seal analysis work, not free solver work or independent discretization validation.

Outputs are compact `STATIONARITY.json/csv`, `UNIQUE_FORWARD_COSTS.json/csv`, `MODEL_SETUP_COSTS.json`, and `RECEIPT.json`, all outside RESULTS. Execution is deferred until the parent writes joint authorization and chooses to freeze/run this source.

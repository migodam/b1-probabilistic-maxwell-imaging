# Zero-new-physics cache audit for finite feedback

2026-10-16. Read-only audit of historical source and NPZ arrays. No solver/backend imported or run; no forward states, LU, adjoints, JVPs, training, or new noise. Only /tmp/fff_cache_audit written. All development cases discussed are historically exposed; GT use is offline diagnostic only.

## Frozen first case and usable comparison

First case is **p250065_weak_a1_snr20**, fixed by A17/N3_APPLICATION_R1/configs/RESOLVED_CONFIG.json: first parent 250065, first background weak, amplitude 1.2, noise20. This is family shell, not a quality-based choice. It is row **7** in both A26 PREDICTIONS.npz and A25 validation.npz (zero-based). Historical noise realization lies outside tau=23.785449 even at truth (24.090925): retain original sigma/tau and report this fact, not model mismatch.

Saved physics-only images on this exact case exist in:
`Gaussian/A26/a26-om-feedback-prior-completion/results/a26/evaluation/261009/PREDICTIONS.npz`

- COMMON_BP_INITIAL: known-background adjoint initialization, projected box; no GT input.
- BASE_N1 and BASE_N3: common BP initialization followed by one L material-gradient update, projected box; no network/GT input.
- BASE_N2: same initialization followed by one R closure-gradient update, projected box; no network/GT input. This is a physics-only finite ROM-based **one-step base**, not a converged reconstruction.
- A25_CERTIFIED_TV_INITIAL_OFFLINE_CACHE: previously computed background-started TV initializer, not GT initialized. Provenance: A25 src/a25/dataset.py 413–424 passes only cache, observation, sigma and case ID to initialize; src/a25/campaign.py 38–49 starts xb=mass_pack(chi_b).
- N3: frozen learned output on same observation; trained on TRAIN and checkpoint-selected on exposed validation.

A26 source encoder.py 425–436 establishes BP and BASE arms. evaluation.py 154–160 exports all above images. There is no cached current/field/full forward prediction for these physics-only base images in PREDICTIONS; evaluation metrics exist but do not replace arrays. An independent background-started fullGN or iterative TSOM solution on this A26 case was **not located** in requested scope. Older A25 B0 image files have parents 2401–2404, hence are not a same-case N3 comparison. A24/A26 anatomy TSOM routines explicitly say “not a TSOM reproduction” and “no reliable iterative SOM/FFT-TSOM baseline claimed.”

## N3-refined fullGN/ROM are not independent physics-only

`Gaussian/A17/N3_APPLICATION_R1/src/pilot.py` 213–235 loads existing N3 exact state, builds J and gradient at N3, and applies full_GN_anchored and current-ROM steps. Thus these are **PHYSICS_REFINED_N3**, no GT initialization, but depend on learned N3. Never relabel them independent physics-only reconstructions. The method outputs are suitable same-case local diagnostics.

First-case persisted files under `Gaussian/A17/N3_APPLICATION_R1/raw/`:

- `<case>_shared.npz`: chi_N3[576], observation[32,8], sigma, r[512], J_full[512,1152], full_step[1152], receiver[576,32], candidates[576,12], anchor[576,64], anchor_J[512,1152], gradient[1152]. J_full is at N3, not at fullGN or ROM endpoint, and was historically paid for.
- `<case>_online_full_GN_anchored_k576.npz`: chi, prediction, proposed_step. Same for online_receiver_ROM_k8 / k16.
- Non-online full_GN_anchored and receiver_ROM files additionally save truth/background and proposed/executed steps; keep these offline.

## Exact N3 state, GT and operators

Portable root:
`Gaussian/A30/a30-maxwell-material-manifold/inputs/a29_release/`

- `inputs/a28/results/a28/exact_states/<case>.npz`: chi[1,576], current/field/w[1,576,8], data[1,32,8], LU[1,576,576], pivots[1,576], residual, background. This is N3 state, not GT state.
- `inputs/a28/results/a28/anatomy/<case>.npz`: separately labelled TRUTH_CENTERED_ORACLE material/current/field/F/L arrays (5 endpoints), CNN_CENTERED_DATA_ONLY arrays (4 endpoints), names, observation, truth/background/sigma. Exposed offline benchmark only. Do not initialize from the TRUTH_CENTERED_ORACLE entries or imply they are recovered solutions.
- `inputs/a25/validation.npz`: row7 chi GT[576], clean[32,8], current[576,8], noisy data, sigma, noise seed, initial_mass, split/family/parent/case metadata.
- `inputs/a25/background/weak.npz`: G[576,576], S[32,576], incident/eb[576,8], yb[32,8], chi_b[576], Db, L_b, background real-mass J[512,1152], volume/grid/points/packing provenance.

Full inventory, absolute paths, array shapes and dtypes: ARRAY_INVENTORY.json. This inspection loaded existing arrays only; it did not recompute responses.

## Missing assets and boundaries

1. N3-refined fullGN/ROM endpoint currents, total fields, exact endpoint Jacobians and LU are not in their saved endpoint files. Do not infer availability from saved prediction.
2. Independent physics-only BP/L/R/TV image arrays exist but their exact endpoint state/prediction/J caches were not found in this scoped audit. Zero-new-physics work must stop at saved-image or background-operator diagnostics for these arms.
3. No original-P0 phase/range certificate follows from these images: LXZ E0 established 0/32 N3 original-phase passes; N3 box validity is weaker than original P0 membership.
4. LXZ's “only background J available” claim is specific to its qualified original-P0 anchor audit. A17 has a historical exact J at N3, but N3 remains non-P0. Do not turn either statement into a contradiction or a P0 certificate.
5. `A30/src/a30/physics.py:cached_anchor` must not run under zero-new-physics: it performs 8 charged audit RHS even with saved LU. `receiver_transfer` costs 32 adjoint RHS. Existing LU does not authorize free new solves.
6. T research's MISSING_CACHE_MANIFEST.json likewise records missing Research-A endpoint current/field and rigorous output/stability bounds; sparse64 capacity has no saved full-forward endpoint predictions.

## Source records read

- B1/b1-research-lxz/research/research_lxz_20261016/E0_REPORT_ZH.md and delegated/e0_2d/assets.json, N3_ONLINE_CACHE_AUDIT.csv.
- B1/b1-research-lxz/research/research_t_20261015/E0_AUDIT.md, RESEARCH_REPORT_ZH.md, METHOD_COMPARISON.md, e0/MISSING_CACHE_MANIFEST.json.
- A17/N3_APPLICATION_R1/research/delegated/e0_n3/E0_N3_SOURCE_MAP.md, configs/RESOLVED_CONFIG.json, src/pilot.py, matched comparison and adaptation reports.
- A26 encoder.py/evaluation.py and A25 dataset.py/campaign.py as cited above.

Scientific status: cache availability established; no new physical acceptance, novelty, robustness, improvement, or finite-feedback success established. Independent physics-only fullGN/TSOM matched reconstruction remains missing.

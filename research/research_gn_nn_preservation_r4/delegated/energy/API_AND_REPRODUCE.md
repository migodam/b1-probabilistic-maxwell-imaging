# Fixed-energy API and E0 reproduction

```python
from research.research_gn_nn_preservation_r4.src.energy import FixedEnergyGeometry
geometry = FixedEnergyGeometry.from_cache(frozen_spectra_path, G=G, S=S)
energy = geometry.energies(error)       # (576,) or (576, source_count)
projected = geometry.project(error, 'D509')
change = geometry.transition(error_before, error_after)
# change['P_in_strong19']: before, after, delta, update, signed_cross
```

The fixed spectra file is `research/research_gn_green_error_20261106/results/fixed_spectra.npz`. Its right-domain complex bases define all Green projections. It contains the complete receiver basis, including structural nulls. The class checks the frozen 19/557, 9/39/528 and 67/509 rank contract. `energies(..., twofold=False)` and `transition(..., twofold=False)` omit sequential and intersection columns if the consumer only needs single projections.

G/S field energies are included when the corresponding matrices are passed. Current errors require actual nonlinear Maxwell currents computed independently at each material. `j_aux` is an auxiliary-head discrepancy and cannot stand in for that current. For a real material vector, complex Green-domain embedding is only a geometric diagnostic, distinct from the real physical material Jacobian.

Sequential `I_to_O__...` means P_O(P_I x), and `O_to_I__...` means P_I(P_O x). These families are kept separate from `intersection__...`. The numerical intersections have weak/weak rank 20 and weak/null rank 509; the remaining intersection cells have rank zero. `intersection_remainder` retains the orthogonal remainder. No projected-span/TSOM energy is substituted for an intersection.

The signed cross field is 2 Re <T error_before, T update> for each fixed linear projection or field map T, computed by polarization. Negative delta is an improvement in absolute squared error. No normalization by total error is performed.

E0 is an exposed historical R3 diagnostic. Its runner never reads R4 test truth. From the repository root, use the existing NN environment, one CPU process, BLAS=1:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_gn_nn_preservation_r4/experiments/e0_energy.py --max-minutes 30
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 "/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_gn_nn_preservation_r4/delegated/energy/finalize.py
```

The second command only reads saved energy tables and frozen physical scalar receipts. It adds the decisive true-current table, residual/noise relation, real-J summary and the true-current strong-error figure; it performs zero LU and no new optimization. Raw means the historical unclipped NN material. All 1278 differences from the verified clipped raw_1 are lower-bound float32 roundoff only (maximum Linf4.47035e-10, RMS3.56345e-10); there is no material clipping confound. Exact unclipped-raw observed-misfit relations remain formally NOT_RUN for these cases, with the actual verified raw_1 residual preserved separately. Absolute true-current energies are fully evaluated.

## New TEST runner (prepared, NOT_RUN)

`experiments/test_energy.py --seal PREDICTIONS_SEALED.json --runs-json FROZEN_RUNS.json` defaults to `results/test_energy`. The explicit run JSON is a list of names or `{plans:[{run,arm,seed,passes},...]}`. Metadata dictionaries are preferred so seeds/cascades remain explicit.

Before the first truth read it verifies the new R4 study identity, TEST-only full observation scope, exact explicit run list, complete finite closure inventory, and every sealed JSON/NPZ size/mtime stamp. It compares the sealed test identities with the entire public online TEST manifest. It performs no hash checks, no new fitting and no method selection. Sealed NOT_RUN predictions stay NOT_RUN. Only after that boundary may it read GT and recompute actual nonlinear Maxwell currents.

Unlike the historical E0 formal raw definition, this runner uses the **actual physically evaluated raw_1 candidate material**, including its original box clipping. It records actual-current and distinct material-coordinate errors, auxiliary-head discrepancies separately, absolute before/after/delta/update/signed cross, all fixed single/twofold/intersection projections, grouped noise/amplitude/arm/seed/cohort means and mother-paired bootstrap intervals. Exact observed-misfit/noise/true-S-field relations are calculated from the same real physical states. No new TEST truth has been read to prepare or syntax-check it.

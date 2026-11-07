# E0 mechanism figures

**Scope:** read-only plotting of completed historical R3 E0 CSVs. No new TEST truth, raw GT file, LU, Jacobian, SVD, optimizer or training call is made. The plot script imports no project physics or network module.

Run `python research/research_gn_nn_preservation_r4/experiments/plot_energy_mechanism.py` from the repository. It writes only the two requested PNGs and this note. Inputs are `results/e0/grouped_absolute_energy.csv` and `results/e0/per_case_energy.csv`. The former supplies saved energy sums; the latter is read only to confirm matching observation/mother identities before connecting stages.

Both figures use seed7, one pass, low amplitude, clean/noisy separately, with all 64 mothers per noise. The selected arms are J-full, J-null and voxel physics (`v_physics`), the voxel arm with physical training loss. `v_image`, other seeds, cascades and A2.5 are outside these figures, not removed from E0. Amplitude groups are pooled with their saved case counts; means are not unweighted means of amplitude means.

The unit is the discrete squared current norm Σ over 576 cells and 12 sources of |error|², averaged across the 64 observations. It is not physical power, flux, a continuum integral or an error-energy share. All signs and numeric magnitudes come from saved CSV fields. No new confidence intervals or significance claim are produced.

## Twofold absolute delta

`figures/e0/twofold_absolute_delta.png` has six rows (noise × arm) and three strictly separate definition columns. Each small matrix is 2×3: internal strong19/weak557 by receiver strong9/weak39/null528. Cell text gives original GN energy, Δraw versus GN and Δfinal versus GN. Cell color shows signed Δfinal on a shared symmetric-log color scale; negative decreases remain negative and blue. Color scaling changes display only, not the saved values.

- I→O means ‖P_O P_I e‖². O→I means ‖P_I P_O e‖². They are sequential fixed-operator energies. The resulting joint cells are not labelled mutually orthogonal subspaces or a true intersection partition; the two orders must not be interchanged.
- The numerical-intersection column is a separate definition. Nonzero intersection ranks are weak/weak20 and weak/null509 (D509); the other four cells have rank0. The rank47 orthogonal remainder is printed explicitly below every matrix. Its GN value and both signed deltas are retained. This column does not substitute TSOM projected spans or sequential products for intersections.
- Errors are the actual nonlinear Maxwell-current errors already computed in E0. Raw is the frozen unclipped proposal convention used by that historical audit. All differences between that proposal and the verified raw_1 material were previously documented as float32 lower-box roundoff only (max Linf4.47035e-10), not substantial clipping. No new solve is made here.

## Auxiliary-to-physical chain

`figures/e0/auxiliary_to_physical_chain.png` shows total, internal strong19, receiver strong9 and D509 energies, separately for each noise level. The four positions are GN physical current, auxiliary head, actual Maxwell current at the raw material, and actual Maxwell current at the selected final top1.

The scalar energies are comparable as discrepancies to the same stored GT current under the same fixed projections, same discrete norm and exactly matched 64 observation identities. J-full/J-null have all four saved positions. Voxel physics has no auxiliary-head rows: its missing point is N/A and is neither zero-filled nor joined by an invented intermediate value. The auxiliary point is always labelled a head discrepancy; only the GN/raw/final points are true nonlinear Maxwell-solution current errors.

**Limits/gaps:** the saved E0 CSVs do not contain the full head→raw-physical vector displacement or its signed projected cross term, a constitutive-defect table, or an isolated readout-only intervention. Final top1 also includes finite-candidate physical correction and selection. Consequently these connecting lines do not establish a causal readout contraction, explain the source of an increase, or isolate a TRF effect. Only the available absolute stage energies are plotted; no missing vector or defect value is reconstructed.

## Saved scalar anchors for interpreting the chain

|noise|arm|band|GN physical|auxiliary head|raw physical|final physical|
|---:|---|---|---:|---:|---:|---:|
|0.0|J-full|S_strong9|2.083934e-13|0.03268086|7.495889e-05|3.685732e-08|
|0.0|J-full|D509|0.01780531|0.04006358|0.01095182|0.0131727|
|0.0|J-null|S_strong9|2.083934e-13|2.08675e-13|3.472778e-05|3.631096e-08|
|0.0|J-null|D509|0.01780531|0.04041752|0.01018076|0.01134168|
|0.0|voxel physics|S_strong9|2.083934e-13|N/A|3.213617e-05|7.533296e-08|
|0.0|voxel physics|D509|0.01780531|N/A|0.008102429|0.009004057|
|0.001|J-full|S_strong9|1.258539e-07|40.82543|0.001368758|2.49042e-06|
|0.001|J-full|D509|0.8732758|1.894103|0.02274629|0.03963404|
|0.001|J-null|S_strong9|1.258539e-07|1.258538e-07|0.00196115|6.477354e-06|
|0.001|J-null|D509|0.8732758|2.7338|0.03229202|0.08855825|
|0.001|voxel physics|S_strong9|1.258539e-07|N/A|0.000144869|1.523893e-06|
|0.001|voxel physics|D509|0.8732758|N/A|0.02303019|0.02864296|

For example, clean J-null receiver-strong9 head discrepancy stays near the GN baseline, while the corresponding raw-material real Maxwell-current error is much larger; physical correction lowers that raw error but leaves it above the extremely small GN baseline. This is a descriptive historical stage association, not a new TEST result or a theorem about arbitrary readouts. D509 can improve at the physical endpoints while these strong-direction comparisons differ. The parent retains mechanism interpretation and final scientific judgment.

## Coverage and writes

Every plotted physical raw/final group and available auxiliary group has 64 identical observation/mother pairs. Voxel auxiliary data are intentionally absent. No E0 source table or pre-existing source file is edited. The only new outputs are this note, `twofold_absolute_delta.png`, and `auxiliary_to_physical_chain.png`; the new plotting source is `experiments/plot_energy_mechanism.py`.

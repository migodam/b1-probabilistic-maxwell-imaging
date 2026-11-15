# Fixed-state explanatory diagnostic implementation (not executed)

Only `src/fixed_state_analysis.py` is added. No frozen main factory/wrapper/config/runtime/theory source is changed for this diagnostic. It has been syntax-compiled and its `--help` inspected; no diagnostic physics, optimization or main data analysis has run.

The CLI requires `--root`, `--legacy-repo`, `--seal`, `--out`. The seal is the parent-approved terminal `observable19_prediction_seal_v1` from the existing `analyze_results.py`; its exact `authorized=true`, frozen144-task configuration, complete terminal queues and result inventory are validated before even creating output. Output must be a new supplement audit/ or analysis/ directory. There is no GT reader or GT helper import.

The fixed21 source slots are cases0000/0002/0009/0012/0015/0020/0021 and cumulative readout1/128/512 of internal_galerkin. Actual existing RETURNED inputs are used; absent/INTERRUPTED/incompatible inputs generate six missing/null rows without replacement runs. Each arm's U/B is read from its own sealed `DENOMINATOR_BASIS.npz`. A thin subclass only sets U/B/n/rank/projection/cost and inherits all LowRankDenominator methods; it does not call its G@U constructor or compute a basis SVD.

Per valid common source, one original `MaxwellModel.data_value_gradient` pays the common original state/LU plus ns source-adjoint RHS. The exact retained LU then pays one nr adjoint-RHS solve for `O_m=S(I-D_mG)^-1`. The original state, LU factors/pivots, physical gradient and O_m are retained in the common NPZ before six proposals change the model cache. Repeated identical proposal material can reuse an already-paid state; it never forces another common-m forward. The source DPCS rank/Vd are loaded from baseline DPCS_RULE/RECEIVER_DPCS; missing basis overlap remains null, no new DPCS SVD.

Recorded current-space definitions are explicit: `Pe=U(U*e)`, raw energy ratio `||Pe||²/||e||²`, receiver transport omission `||O_m(I-P)e||/||O_m e||`, Green input omission `G(I-P)e` relative to `Ge`, and PGP extra output trim `(I-P)GPe` relative to `GPe`. The last quantity is still recorded counterfactually for GP with `PGP_output_trim_applied_to_this_arm=false`. The transported ratio may exceed1 and is never clipped or called a material-information percentage. DPCS subspace trace overlap and projections of e/Pe/(I-P)e are saved, with zero-denominator reasons.

Full/constitutive/metric/prior gradients use original R11 `gradient_parts`. The original R10 `matrix_material_readout` receives the same J/Q/m/beta/rho/Os/Ls for all arms, original real box/prior and reconstructed R8 arm_config('a_plus') material flags, matrix_inner_maxiter8, original seed/old comparison. One call creates each OFFLINE_DIAGNOSTIC proposal, then original full Maxwell Phi is paid and recorded. No proposal enters any main pool. The19×19 Woodbury T spectrum is a separately paid `svdvals` diagnostic only; failure/near-singularity does not gate readout or alter rank.

Outputs include common states/retained LUs, per-arm projection and gradient/proposal arrays, all warnings/failures/null reasons, common and arm tables, event-level elapsedCPU plus original physical counters, separate lowrank/dense projection costs, and a no-GT diagnostic output inventory seal. Actual aggregate physical RHS takes precedence over role counters; small19-LU/RHS are not added to576-dimensional LU/RHS. Inherited logical prior/helper counters are not summed again as disjoint evaluations. No scientific result is claimed because this implementation has not run.

After the parent has reviewed144 terminal outputs and written an authorized seal, a real Windows command in the existing numerical Python environment is:

```powershell
$r14Root = 'D:/AI/MAX_STRUCT_R14_DUALMETRIC_20261115'
$supplementRoot = "$r14Root/supplement_observable19"
$legacyRoot = 'D:/AI/MAX_STRUCT_R11_DEV_20261107'
python "$supplementRoot/src/fixed_state_analysis.py" --root $supplementRoot --legacy-repo $legacyRoot --seal "$supplementRoot/audit/PARENT_AUTHORIZED_TERMINAL_SEAL.json" --out "$supplementRoot/audit/fixed_state_analysis_v1"
```

`PARENT_AUTHORIZED_TERMINAL_SEAL.json` is a placeholder filename for the actual parent-approved observable19 seal, not a file claimed to exist now. The script does not create or approve that authorization. No command above has been executed.

# R15 denominator wiring receipt

Implementation: `research/research_max_struct_r15/src/twofold_denominator.py`.

`TwofoldDenominator(G,U,shape,one_hop=False,cost=None)` takes the externally supplied fixed orthonormal basis; no rank selection/SVD is performed. The parent-supplied `U11` was loaded from the existing current-complex characterization for public r14_new_0000. G19/S19 controls use the already-saved right frames. Rank0 has a diagonal factor/solve with no small LU and no padding.

The original R10 low-rank constructor/factor/solve/action and `matrix_material_readout` are reused. `twofold_material_readout` is an alias to that original function, so the LBFGS8, scalar seed, old-point comparison, bounds and prior are not reimplemented. Supply the operator as `config['matrix_denominator']`; the inherited helper's reassignment of the shared counter also updates the original OneHop block counter. No optimizer was run in this check.

Optional hop uses the original fixed disjoint4cell off-diagonal block N and exactly R=N−PNP. Setup retains U*NU once, so each R or R* action uses one original block action and the low-rank retained-neighbor subtraction; the retained PNP component is not added twice. `transform(m,beta,rhs)` and `adjoint(m,beta,rhs)` (also `transform_adjoint`) provide the exact adjoint pair of the truncated map. `transform(...,adjoint=True)` is equivalent.

No-hop `value_gradient` delegates directly to the original helper. The hop branch implements the parent-supplied v,w,p,q,H,T equations and uses the original prior callback once per evaluation. This is a finite truncated residual correction: **no all-box contraction guarantee is claimed**.

Required local check: public r14_new_0000 G/E0/grid_shape only, fixed auxiliary J=.1 E0 and Q=E0+GJ (not a physical PDE solve), fixed source-normalized diagonal Os. Two legal materials: uniform .1 and deterministic .1+.035 sin(.17 index)+.02 cos(.31 index). Two fixed unit directions, central difference epsilon1e-4, beta=rho=1. Original smooth/TV regularizer defaults were retained. The same public diagnostic point was used for U11 and the legitimate rank0 branch; existing G19/S19 frames are only no-hop equality controls. These are wiring inputs, not selected reconstruction candidates or new parameter recommendations.

Results:

- Eight central directional comparisons: max relative error 2.514068112894441e-08.
- Six no-hop comparisons (U11/G19/S19 × two materials): value exactly equal, gradient arrays bitwise equal to original R10 helper.
- Four transform/adjoint pairings (U11/rank0 × two materials): max relative difference 2.559663404696728e-16.
- Shared-counter reassignment and inherited readout function identity passed.
- No truth access, physical LU/RHS, new SVD, online inversion or optimizer execution.

`DERIVATIVE_VALUES.csv`, `VALUES.json` and `FIXED_WIRING_INPUTS.npz` retain actual inputs/results. `RECEIPT.json` records setup/action/solve counts, CPU and wall. Original small-LU/block counters and newly added request/action labels are overlapping descriptions of the same work and must not be summed as additional operations. Reference setup and the separate counter-reassignment construction are explicitly counted.

The first required check passed. A second identical check was run only after correcting diagnostic cost reporting for the original reference constructor and the separate counter-reassignment setup. First receipt is preserved as ATTEMPT01_RECEIPT.json; there were no numerical failures. Combined check costs are recorded in TOTAL_COST.json. Source/derivative checks do not provide physical qualification or final proof acceptance.

Actual command, from r14_worktree:

```sh
"/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python" research/research_max_struct_r15/audit/denominator_wiring/check_wiring.py
```

The check sets all supported BLAS thread environment variables to1 before NumPy imports. No installation or remote action was performed.

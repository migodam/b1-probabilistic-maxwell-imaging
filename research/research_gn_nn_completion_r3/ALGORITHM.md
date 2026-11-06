# Algorithm and attribution boundaries

All methods use the frozen nonlinear scalar-TM VIE, 576 independent real materials, 12 source columns and 48 receivers. Public inputs are observations, operators, noise level and the physical GN anchor. Truth is read only by generation, TRAIN/VAL supervision, and sealed offline scoring.

## Representation arms

The common residual U-Net has channels 16/32/64 and a 24→12→6 spatial hierarchy, with a zero-initialized output head. AdamW uses learning rate 1e-3, batch 64, at most 100 epochs and validation patience 15.

- V-image: material image and known noise; physics input channels are zero.
- V-physics: full material residual with background-Jacobian residual backprojection and complex total-field inputs.
- J-full: full complex auxiliary-current residual, followed by shared-real material readout.
- J-null: the same head, projected onto the entire receiver nullspace before readout.
- J-null-current: the nullspace head with projected exact-current supervision instead of material supervision.

The first four arms share material L1 plus 0.1 edge L1. Current-space and material-space DOF are not equal: 576 material outputs versus 13,824 real current outputs, or 12,672 receiver-null current coordinates. The two output heads have 121,137 and 121,528 network parameters respectively. Parameter counts do not establish equal admissible physical families.

## Online inference and correction

1. Generate a public-only proposal. Current proposals recompute E+G j and use the shared-real rowwise least-squares material readout.
2. Form material proposals at scales 0.5 and 1, under the common material box.
3. Evaluate the full nonlinear Maxwell model; correct each proposal with bounded analytic-J TRF for at most 3 seconds or 8 function evaluations, stopping at the declared data gate.
4. Keep at most five slots: anchor, two proposals, two corrected proposals.
5. Among physically qualified candidates choose the material closest to the full learned proposal. Otherwise retain a qualified anchor; if none is qualified report FAIL. Also record the old minimum data-plus-TV selection.

Each arm has a 10-second / 100-LU / 6,000-complex-RHS incremental budget, including failed attempts and preparation. Anchor construction and shared reference setup are charged separately. Continued and weakly regularized ordinary TRF receive the same incremental budget.

No primary arm truncates a weak spectrum. Full-operator soft filtering is a secondary ablation. Shared-weight two-stage cascades require validation evidence before activation. Availability in source code is not evidence that an optional arm ran; consult actual receipts and the final summary.

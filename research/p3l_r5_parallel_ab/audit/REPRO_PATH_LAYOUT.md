# Reproduction path layout

Keep the existing sibling layout under Gaussian/: B1/b1-p3l-r5-parallel-ab (newR5), B1/b1-p3l-r4-conditional-transport (historical R4), and A30/a30-maxwell-material-manifold (frozen paid solver). The R5 shared adapter resolves these sibling locations directly. Exact files, byte sizes and source Git identities are in EVIDENCE_INPUT_REFERENCES.json. Historical R4 data are referenced, not copied wholesale.

The original P1 decoder/config/points remain inside the R5 repository at src/p1_material/ and research/p1_material_target_20261018/. A30 requires its src/ and vendor/frozen_a29/src/ code trees plus the named paid weak.npz background. Importing background arrays is not a new solve. Any fresh forward/inverse action must use the shared cumulative R5 ledger; published reproduction is not authorization to reset or silently replay paid physics.

Use the existing Gaussian/.venv_nn/bin/python with NumPy/SciPy/Torch and VECLIB_MAXIMUM_THREADS=1, OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1. B's final observation-only entry is stage_b_material_risk/inference.py:Action. y is finite real source-major512; online keys are observation,residual,q0,chi0,J from the unchanged historical96+64 strong anchor. Original K32 additionally requires parent IDs for its frozen sample-seed law. Student inference needs only its R5 PROVIDER.pt/HEAD.pt plus historical feature-source code and online anchor input.

The executable replay script uses retained VAL data and therefore tests serialization/device replay, not fresh independent generalization. Preserve INFERENCE_REPLAY.json's declared tolerance and actual status. It may contain small CPU/MPS float32 drift while material risk agrees.

A5 labels are directional: global decoder constants are along one fixed unit strong W segment at a fixed weak center. The legacy RESULTS key L_on_rho_ball does not establish a full2D/66-dimensional ball. Norm values are Float64 evaluations of exact-arithmetic inequalities, without directed-rounding enclosure. Parent captions must state this scope.

# R4 representation/training implementation

Status: six representations and three loss configurations implemented; initial five seed7 training waits for the data worker's strict full TRAIN/VAL release. No historical weights loaded. Frozen R3 files untouched.

- `src/networks.py`: shared-weight two-pass U-Net widths16/32/64,27input channels; current heads all24outputs/121528parameters; voxel1output/121137parameters.
- Fixed projections: `j_full=I`, `j_null=I-V48 V48*`, `g_weak=I-V19 V19*`, `g_dark=I-WW*`, `soft_dark=I-(1-alpha)WW*`. G-right top19; S-right all48; W numerical unionrank67. Complex128 CPU SVD cached under `cache/projectors.npz`; complex64/pairfloat32 application. No Weak64 or receiver truncation.
- Both stages update material. Current stages add projected increments to cumulative current. gb/Q/noise and all normalization scales remain at the public anchor. Only material channel0 changes on pass2.
- Final material recovered by shared-real row least squares and box[.02,20]. Material loss is normalizedL1+.1edge. `const` adds .01*||j-Dm(E+Gj)||²/||j0||²; `internal` adds .01*||G(j-j0)||²/||Q0||². Anchor denominators fixed. No current truth used.
- `train.train_model(...,cascade_passes=2,alpha=.3,loss_name='material'|'const'|'internal',operators=shared,feature_stats=stats,cohort=released)` fits fromscratch. AdamWlr.001,batch64,max100epochs,patience15. VAL materialL1+.1edge selects all configurations consistently; regularizer excluded from selection. Seeds7/17/27 supported; initial authorized runs seed7 only.
- `strict_training_release(study)` is the central readiness predicate. `released_training_cohort(...)` performs count/split/provenance check once per multiarm launch and passes one cohort to each arm. Requires1008TRAIN/156VALobs from504/78parents and strict full24 anchors; withheld4/7 excluded fromTRAIN.
- `predict.predict_anchors(checkpoint,anchor_paths,operators_path,out,device,operators=shared)` opens only explicit public anchors and checkpoint; no manifest or truths. Old-compatible checkpointformat`gn_nn_r3_v1` plus mandatory`r4_format=gn_nn_r4_preservation_v1` refuses historicalR3weights. Outputs `m_proposal`,currentarms`j_aux,q_aux,delta_current`, plus `stage_materials`,currentarms`stage_currents,stage_fields`, and`preparation_cost_json`.
- Voxel stage materials are saved; it has no learned auxiliary current and does not fabricate stage currents.
- `experiments/train_models.py --data-ready --device mps --arms v_physics j_full j_null g_weak g_dark --seeds7 --predict-val` writes `models/{arm}_seed7_cascade2/best.pt` and `predictions/{run}/{obs_id}.npz` immediately after each trained arm. Shared public operators and TRAINstats uploaded/fitted once.
- `experiments/predict_models.py` takes explicitJSONlistofpublicanchorpaths and checkpoints, sharing one operator upload.
- Finite waiting launcher `delegated/networks/run_initial_five.py` polls readiness every30s and stops on protocoldeadline or explicit anchorfailure metadata. Logs `initial5_seed7.log`; processreceipt`audit/training_initial5_process.json`.

One meaningful MPS synthetic wiring check passed for allsixarms: zero-head firststep gives positivehead/zerostem gradients; secondstep gives finitepositivehead/stem gradients. Bothstageoutputs and const/internal terms finite. Diagnosticweightsdiscarded; no dataset/GT read. Receipt `audit/training_wiring.json`.

## Optional full-VIE stage (parent gate required)

Implemented but not executed: `train.FullVIEBridge(operators_path)(m_B576)->Y_B48x12x2`. CPUcomplex128 SciPy fullLU retained for analytic adjoint; each forwardcase1LU+12sourceRHS, backwardcase12adjointRHS. VJP useslambda=L^-*S*gradY andRe sumconj(lambda)Q without anextra2. `bridge.receipt()` records eachcall/cost/transfer/setup; `full_maxwell_hinge(m,publicYpair,bridge,tolerance)` requires explicit public discrepancy threshold. Base training never invokes this backend. Parent controls5%VALgate, exact128TRAINparentset,batch8,10epoch matched purematerial/hingestage. No finite-difference physics or all576 checks run.

Softalpha training, extra seeds, E2 ablations and optional full-VIE fine-tuning remain WAIT_PARENT. Posthoc JfullVAL projection may selectalpha; it does not replace fromscratch SoftDark training.

# R5 Stage I：只读状态库存

模型配置：Sol medium。仅读取公开输入、预测、receipt/history 和源码；物理求解 0、测试 0、GT 读取 NONE。原结果未修改。数组仅作键、shape、dtype 库存，没有求解或数值验证。

工作树：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree`

## 结论

- 旧开发 b20d 的 global/onehop/phase_pi/self 各有三个完整阶段端点 m/J/Q/A，均 mf3/6/10、500/500/500 步；PR 标量历史存在，但 PR 向量状态和逐步材料轨迹未保存。
- 三个确认对象的 global/onehop 只有四行 candidate_m、1503 行标量历史、corrector 标量历史和最终精确物理 m/j/Q。confirmation_solver 没有输出 aux.stage_states，无法从阶段材料重建当时的辅助 J。
- 三个确认对象没有 phase_pi/self 历史结果；不能把旧开发控制移植为已完成确认。
- 最终 prediction 的小写 j 是材料对应精确物理解，stage 的大写 J 是优化中的辅助电流。二者不可互换。

## 可重建 J 与状态语义

S=Us diag(sig) Vh; r=count(sig>0.001*sig[0]); Vd=Vh[:r].conj().T; JD=Vd@((Us[:,:r].conj().T@Y)/sig[:r,None]); construct U with original fourier_apcs(shape,mf,Vd); J=JD+U@A; Q=E0+G@J. This requires SVD/basis arithmetic, not a Maxwell material solve; not executed in this inventory.

fourier_apcs 的 U 使用投影后 Fourier 原系数坐标：F 的角频率枚举是每轴 `[0,…,mf−1] ∪ [size−mf,…,size−1]`，`indexing="ij"`，再用原 SVD 截断舍入级 APCS 零动作，不能随意正交化后复用 A。阶段 A shape 为 (36,12)、(144,12)、(400,12)。本库存不运行这些重建。

阶段 NPZ 的 m 是当前更新后的材料，J/Q 是当次 current 更新的辅助状态；m 更新后不会重新求精确物理解。阶段 metadata（β、mf、rank、stage stop）在 AUXILIARY 的 stage_end，不在 NPZ 内。

## PR 与最小重放边界

两个主源码在每阶段初始化 `old_grad=old_direction=None`；阶段间按频率零填充旧 Fourier A。PR+ 每步用当前梯度和上一梯度/方向，history 仅保留 PR_plus、line_alpha、direction_restarted。已有 iteration_callback 是 post-readout m/J/A；wrapper 保留上一次 callback m（首步用 host）即可获得 pre-readout m。Q=E0+G@J 只在预定 checkpoint 计算，避免每步额外 G 动作。stage_callback 已直接提供 m/J/Q/A。以上无需修改数值源码；精确内部 PR 连续重启仍缺 old_grad/old_direction。

- **old_stage_diagnostics**：No replay needed for m/J/Q at saved stage ends; use metadata from matching AUXILIARY history stage_end.
- **old_boundary_restart**：Next stage can start from saved m and A using original frequency zero-padding; PR is reset at stage boundaries. Exact numerical reproduction remains environment-dependent.
- **old_interior_trajectory**：Replay only auxiliary stages to requested point; no corrector needed. Save pre-current m/A plus PR old_grad/old_direction if exact interior restart is required.
- **confirmation_stage_states**：Replay auxiliary only from public host m0 and original frozen config because only stage m survived. Stop at needed stage; do not rerun corrector to obtain auxiliary state.
- **callback_gap**：Existing iteration_callback exposes post-readout m/J/A; stage_callback exposes m/J/Q/A at stage ends. A wrapper can retain previous callback m (initial host before first call) for m_before and reconstruct Q=E0+G@J at predetermined checkpoints only. No solver edit needed for these diagnostics. old_grad/old_direction remain unavailable for exact interior PR restart.
- **absent_controls**：phase_pi/self on confirmation cases have no saved run; cannot be called replay. A new control execution requires explicit separate scope and cost.

## 逐例逐臂与真实路径

### obs_b20d5a55677d4eca / global

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r2/results/blind/obs_b20d5a55677d4eca/global/stage_0.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [36, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r2/results/blind/obs_b20d5a55677d4eca/global/stage_1.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [144, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r2/results/blind/obs_b20d5a55677d4eca/global/stage_2.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [400, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/global_common_corrector/SELECTED_BEFORE_CORRECTOR.npz`
  - `m` [576] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/global_common_corrector/prediction.npz`
  - `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r2/results/blind/obs_b20d5a55677d4eca/global/AUXILIARY.json`
  - history 1503 行 / 1500 更新行；stage_end：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 1.2875619158148766}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 2.812621207907796}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 5.181019040755928}]
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/global_common_corrector/RECEIPT.json`
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/global_common_corrector/CORRECTOR.json`
- 缺失：per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_b20d5a55677d4eca / onehop

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/ALL_PHYSICAL_CANDIDATES.npz`
  - `m` [4, 576] float64; `original_Phi` [4] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/SELECTED_BEFORE_CORRECTOR.npz`
  - `m` [576] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/prediction.npz`
  - `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/stage_0.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [36, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/stage_1.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [144, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/stage_2.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [400, 12] complex128
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/AUXILIARY.json`
  - history 1503 行 / 1500 更新行；stage_end：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 14.884705666918308}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 23.56946129212156}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 31.14269237499684}]
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/RECEIPT.json`
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_absorption/CORRECTOR.json`
- 缺失：per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_b20d5a55677d4eca / phase_pi

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/ALL_PHYSICAL_CANDIDATES.npz`
  - `m` [4, 576] float64; `original_Phi` [4] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/SELECTED_BEFORE_CORRECTOR.npz`
  - `m` [576] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/prediction.npz`
  - `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/stage_0.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [36, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/stage_1.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [144, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/stage_2.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [400, 12] complex128
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/AUXILIARY.json`
  - history 1503 行 / 1500 更新行；stage_end：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 18.76197945792228}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 37.548972375225276}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 56.02026104182005}]
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/RECEIPT.json`
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/onehop_phase_pi/CORRECTOR.json`
- 缺失：per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_b20d5a55677d4eca / self

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/ALL_PHYSICAL_CANDIDATES.npz`
  - `m` [4, 576] float64; `original_Phi` [4] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/SELECTED_BEFORE_CORRECTOR.npz`
  - `m` [576] float64
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/prediction.npz`
  - `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/stage_0.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [36, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/stage_1.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [144, 12] complex128
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/stage_2.npz`
  - `m` [576] float64; `J` [576, 12] complex128; `Q` [576, 12] complex128; `coefficients` [400, 12] complex128
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/AUXILIARY.json`
  - history 1503 行 / 1500 更新行；stage_end：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 11.506688167341053}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 19.597589417360723}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 26.0613762922585}]
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/RECEIPT.json`
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/self_absorption/CORRECTOR.json`
- 缺失：per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_42cb233e4f224bb3 / global

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/predictions/obs_42cb233e4f224bb3/global.npz`
  - `candidate_m` [4, 576] float64; `candidate_original_Phi` [4] float64; `selected_before_corrector_m` [576] float64; `corrector_history` [1019, 2] float64; `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
  - 标量列 `history_PR_plus`, `history_apcs_rank`, `history_beta_current`, `history_beta_material`, `history_coefficient_relative_change_squared`, `history_current_solver`, `history_current_step_objective_change`, `history_data_after_current`, `history_data_before_current`, `history_data_term`, `history_direction_restarted`, `history_event`, `history_fourier_columns`, `history_global_iteration`, `history_iteration`, `history_line_alpha`, `history_material_inner_message`, `history_material_inner_nit`, `history_material_inner_success`, `history_mf`, `history_readout_material_change`, `history_receiver_rank`, `history_regularizer`, `history_sourcewise_cg`, `history_stage`, `history_state_after_current`, `history_state_before_current`, `history_state_term_current`, `history_state_term_material`, `history_status`, `history_wall`；shape [1503]。完整 dtype 见 JSON。
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/results/confirmation/obs_42cb233e4f224bb3__global.json`
  - result.stage_history：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 5.21257669996703}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 15.26454989996273}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 34.2877196999616}]
- 缺失：stage auxiliary J/Q/A, per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_42cb233e4f224bb3 / onehop

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/predictions/obs_42cb233e4f224bb3/onehop.npz`
  - `candidate_m` [4, 576] float64; `candidate_original_Phi` [4] float64; `selected_before_corrector_m` [576] float64; `corrector_history` [1021, 2] float64; `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
  - 标量列 `history_PR_plus`, `history_apcs_rank`, `history_beta_current`, `history_beta_material`, `history_coefficient_relative_change_squared`, `history_current_solver`, `history_current_step_objective_change`, `history_data_after_current`, `history_data_before_current`, `history_data_term`, `history_direction_restarted`, `history_event`, `history_fourier_columns`, `history_global_iteration`, `history_iteration`, `history_line_alpha`, `history_material_fallback`, `history_material_final_block_loss`, `history_material_inner_message`, `history_material_inner_nfev`, `history_material_inner_nit`, `history_material_inner_success`, `history_material_inner_total_evaluations`, `history_material_old_block_loss`, `history_material_scalar_seed_block_loss`, `history_mf`, `history_readout_material_change`, `history_receiver_rank`, `history_regularizer`, `history_sourcewise_cg`, `history_stage`, `history_state_after_current`, `history_state_before_current`, `history_state_term_current`, `history_state_term_material`, `history_status`, `history_wall`；shape [1503]。完整 dtype 见 JSON。
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/results/confirmation/obs_42cb233e4f224bb3__onehop.json`
  - result.stage_history：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 53.10481630003778}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 110.71455740003148}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 167.54560990002938}]
- 缺失：stage auxiliary J/Q/A, per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_42cb233e4f224bb3 / phase_pi

NOT_RUN / 无匹配已保存控制。

### obs_42cb233e4f224bb3 / self

NOT_RUN / 无匹配已保存控制。

### obs_63abafa750fc483c / global

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/predictions/obs_63abafa750fc483c/global.npz`
  - `candidate_m` [4, 576] float64; `candidate_original_Phi` [4] float64; `selected_before_corrector_m` [576] float64; `corrector_history` [1029, 2] float64; `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
  - 标量列 `history_PR_plus`, `history_apcs_rank`, `history_beta_current`, `history_beta_material`, `history_coefficient_relative_change_squared`, `history_current_solver`, `history_current_step_objective_change`, `history_data_after_current`, `history_data_before_current`, `history_data_term`, `history_direction_restarted`, `history_event`, `history_fourier_columns`, `history_global_iteration`, `history_iteration`, `history_line_alpha`, `history_material_inner_message`, `history_material_inner_nit`, `history_material_inner_success`, `history_mf`, `history_readout_material_change`, `history_receiver_rank`, `history_regularizer`, `history_sourcewise_cg`, `history_stage`, `history_state_after_current`, `history_state_before_current`, `history_state_term_current`, `history_state_term_material`, `history_status`, `history_wall`；shape [1503]。完整 dtype 见 JSON。
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/results/confirmation/obs_63abafa750fc483c__global.json`
  - result.stage_history：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 5.472460399963893}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 15.814570599992294}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 35.67278159997659}]
- 缺失：stage auxiliary J/Q/A, per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_63abafa750fc483c / onehop

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/predictions/obs_63abafa750fc483c/onehop.npz`
  - `candidate_m` [4, 576] float64; `candidate_original_Phi` [4] float64; `selected_before_corrector_m` [576] float64; `corrector_history` [1008, 2] float64; `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
  - 标量列 `history_PR_plus`, `history_apcs_rank`, `history_beta_current`, `history_beta_material`, `history_coefficient_relative_change_squared`, `history_current_solver`, `history_current_step_objective_change`, `history_data_after_current`, `history_data_before_current`, `history_data_term`, `history_direction_restarted`, `history_event`, `history_fourier_columns`, `history_global_iteration`, `history_iteration`, `history_line_alpha`, `history_material_fallback`, `history_material_final_block_loss`, `history_material_inner_message`, `history_material_inner_nfev`, `history_material_inner_nit`, `history_material_inner_success`, `history_material_inner_total_evaluations`, `history_material_old_block_loss`, `history_material_scalar_seed_block_loss`, `history_mf`, `history_readout_material_change`, `history_receiver_rank`, `history_regularizer`, `history_sourcewise_cg`, `history_stage`, `history_state_after_current`, `history_state_before_current`, `history_state_term_current`, `history_state_term_material`, `history_status`, `history_wall`；shape [1503]。完整 dtype 见 JSON。
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/results/confirmation/obs_63abafa750fc483c__onehop.json`
  - result.stage_history：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 53.15994740003953}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 97.91901800001506}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 155.0194513000315}]
- 缺失：stage auxiliary J/Q/A, per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_63abafa750fc483c / phase_pi

NOT_RUN / 无匹配已保存控制。

### obs_63abafa750fc483c / self

NOT_RUN / 无匹配已保存控制。

### obs_100bd65498644322 / global

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/predictions/obs_100bd65498644322/global.npz`
  - `candidate_m` [4, 576] float64; `candidate_original_Phi` [4] float64; `selected_before_corrector_m` [576] float64; `corrector_history` [1046, 2] float64; `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
  - 标量列 `history_PR_plus`, `history_apcs_rank`, `history_beta_current`, `history_beta_material`, `history_coefficient_relative_change_squared`, `history_current_solver`, `history_current_step_objective_change`, `history_data_after_current`, `history_data_before_current`, `history_data_term`, `history_direction_restarted`, `history_event`, `history_fourier_columns`, `history_global_iteration`, `history_iteration`, `history_line_alpha`, `history_material_inner_message`, `history_material_inner_nit`, `history_material_inner_success`, `history_mf`, `history_readout_material_change`, `history_receiver_rank`, `history_regularizer`, `history_sourcewise_cg`, `history_stage`, `history_state_after_current`, `history_state_before_current`, `history_state_term_current`, `history_state_term_material`, `history_status`, `history_wall`；shape [1503]。完整 dtype 见 JSON。
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/results/confirmation/obs_100bd65498644322__global.json`
  - result.stage_history：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 5.321015399997123}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 15.17844490002608}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 33.732301100040786}]
- 缺失：stage auxiliary J/Q/A, per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_100bd65498644322 / onehop

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/predictions/obs_100bd65498644322/onehop.npz`
  - `candidate_m` [4, 576] float64; `candidate_original_Phi` [4] float64; `selected_before_corrector_m` [576] float64; `corrector_history` [1052, 2] float64; `m` [576] float64; `Y_prediction` [48, 12] complex128; `j` [576, 12] complex128; `Q` [576, 12] complex128
  - 标量列 `history_PR_plus`, `history_apcs_rank`, `history_beta_current`, `history_beta_material`, `history_coefficient_relative_change_squared`, `history_current_solver`, `history_current_step_objective_change`, `history_data_after_current`, `history_data_before_current`, `history_data_term`, `history_direction_restarted`, `history_event`, `history_fourier_columns`, `history_global_iteration`, `history_iteration`, `history_line_alpha`, `history_material_fallback`, `history_material_final_block_loss`, `history_material_inner_message`, `history_material_inner_nfev`, `history_material_inner_nit`, `history_material_inner_success`, `history_material_inner_total_evaluations`, `history_material_old_block_loss`, `history_material_scalar_seed_block_loss`, `history_mf`, `history_readout_material_change`, `history_receiver_rank`, `history_regularizer`, `history_sourcewise_cg`, `history_stage`, `history_state_after_current`, `history_state_before_current`, `history_state_term_current`, `history_state_term_material`, `history_status`, `history_wall`；shape [1503]。完整 dtype 见 JSON。
- JSON：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/results/confirmation_remote/results/confirmation/obs_100bd65498644322__onehop.json`
  - result.stage_history：[{"stage": 0, "event": "stage_end", "status": "MAXITER", "mf": 3, "fourier_columns": 36, "apcs_rank": 36, "beta_current": 1.0, "beta_material": 1.0, "wall": 37.84470750001492}, {"stage": 1, "event": "stage_end", "status": "MAXITER", "mf": 6, "fourier_columns": 144, "apcs_rank": 144, "beta_current": 1.0, "beta_material": 0.2, "wall": 63.009114200016484}, {"stage": 2, "event": "stage_end", "status": "MAXITER", "mf": 10, "fourier_columns": 400, "apcs_rank": 400, "beta_current": 1.0, "beta_material": 0.2, "wall": 101.1219918000279}]
- 缺失：stage auxiliary J/Q/A, per-iteration m/J/Q/A, old_grad, old_direction, material-inner L-BFGS memory, corrector material iterates。

### obs_100bd65498644322 / phase_pi

NOT_RUN / 无匹配已保存控制。

### obs_100bd65498644322 / self

NOT_RUN / 无匹配已保存控制。

## JSON 结构、公开载荷和来源

`STATE_INVENTORY.json` 保存每个已检查 NPZ 的所有键/shape/dtype、相关 JSON 完整递归结构、阶段 metadata、配置与绝对路径。确认 JSON 外层是 event/slot/pid/task_id/status/result/task；科学 receipt 在 result，配置在 result.configuration，阶段事件在 result.stage_history，candidate_scores 与 candidate_m 顺序一致。corrector_history 是每次目标梯度评估的 `[原 Phi, ||grad Phi||∞]` 对，不是物理残差或材料向量；行数可大于 optimizer iteration 数。

公开输入：

- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r2/data_public/obs_b20d5a55677d4eca.npz`：G[576, 576], E0[576, 12], S[48, 576], shape[2], grid_shape[2], coords[576, 2], sources[12, 2], receivers[48, 2], Y[48, 12], k[], radius[], sigma[], noise_radius[], material_bounds[2], amplitude[], host[], smooth_weight[], tv_weight[], tv_epsilon[], prior_width[]
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/data_public/confirmation/operators.npz`：G[576, 576], E0[576, 12], S[48, 576], sources[12, 2], receivers[48, 2], k[], radius[], coords[576, 2], shape[2], grid_shape[2]
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/data_public/confirmation/obs_42cb233e4f224bb3.npz`：Y[48, 12], sigma[], noise_radius[], host[], amplitude[], material_bounds[2], smooth_weight[], tv_weight[], tv_epsilon[], prior_width[]
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/data_public/confirmation/obs_63abafa750fc483c.npz`：Y[48, 12], sigma[], noise_radius[], host[], amplitude[], material_bounds[2], smooth_weight[], tv_weight[], tv_epsilon[], prior_width[]
- `/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r5_worktree/research/research_max_struct_r4/data_public/confirmation/obs_100bd65498644322.npz`：Y[48, 12], sigma[], noise_radius[], host[], amplitude[], material_bounds[2], smooth_weight[], tv_weight[], tv_epsilon[], prior_width[]

来源代码：R4 src/cie_onehop.py、src/confirmation_solver.py；R2 src/cie_global.py、experiments/run_blind.py；成熟基线 src/cie_tsom_baseline.py 的 fourier_apcs；R4 onehop/phase/self 的 executed_source_snapshot 与 FROZEN_CONFIG。global_common_corrector 复用 R2 阶段候选并保留历史 incumbent，故其最终选择还需查 RECEIPT.selected_candidate。

只读脚本仅依赖 pathlib/json/numpy，不导入项目求解模块；只枚举明确 allowlist 的结果和公开载荷，不递归扫描 sealed_truth/private。执行环境按任务指定为 Gaussian/.venv_nn/bin/python。该库存不构成数值重放、科学确认或新结果发布。

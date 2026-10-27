# 共享物理入口、在线隔离与合同静态审查

只读审查 `shared/physics.py`、`configs/P0_FREEZE.json`、`src/generate.py`、`src/run_c.py`、`track_c/anchors.py`及导入的原R4核心/优化源码。未导入Target、未求解、未打开GT数组；本审计没有改共享源码。

- 现行生成器单独写 `raw/offline_truth/<parent>.npz` 与仅观测的 `raw/online/<parent>.npz`。其先抽Gaussianized X66再Phi变换后两维，与Gaussian64×Uniform2的分布合同一致；不要求与旧MaterialPrior.draw在同seed下逐样本相同。
- `run_c`与`anchors`无offline_truth读取或generator导入。callback只读取已缓存Point、当前y、已存J，记录`GT_used=false/future_endpoint_used=false`；当前未传decision_predictor。最终context所需J可能追加收费，已明确标final_context_J_required。
- `Target`导入原prior/config及96/64优化函数；native默认polish=False。forward在backend full_state内部先`reserve_state`付8RHS，随后consume credit；receiver/source adjoint经reserve_rhs另付。失败无refund，状态、LU按attempt计。
- `Book`独占physics_worker锁，所有ledger写另持ledger锁。job/action开始结束显式登记；未协调中断阻止再开job。大账本的只读snapshot prefix复核结果另见ONLINE_LEDGER_AUDIT.json；运行中的active job不能当失败或封账。
- `CONTEXT.npz`路径：`track_c/contexts/<group>/<parent>/<config>/CONTEXT.npz`。键为observation/prediction/residual (32,8 complex)、q0/native_q (66)、chi0/native_chi (576 complex)、J (512,66，完整采集)、x (66)、H (66,66)、phi scalar、parent scalar。上游head可直接用q0/chi0/J/x/H；report的output_context指向同文件。原optimizer另写best_cache与RAW_MODE_CACHE，保存current、field、receiver_T等。

## 必须保留的限定与修复建议

1. **source-subset分支尚不可借用完整512资格。** Target支持mask并正确限制white/phi；但原reconstruct的noise-ball自由度仍固定512，若启动subset，须调整或标不适用该报告项。当前CDEV完整32×8无此问题。
2. **冻结存在性不等于冻结内容有效性。** Book只检查P0_FREEZE存在，不检查其状态、源文件身份或研究门。必须由顶层执行前确认合同/代码身份；不要把可创建Book当P0全PASS。
3. **19,000 A states与4个live LU限制不由Book自动硬限。** 当前Book硬限stage RHS/总RHS和job wall；A首轮states需caller另核对。原native缓存4、probit当前/最佳/native及backend transfer缓存2分别有界，但当前代码未给“union最多4个live LU”的证明。建议将max_live_LU_points写成实测/已证上限前保持未验证，或另审对象联合所有权，勿以单个缓存长度代替全局峰值。
4. **历史仍CLOSED。** 任何新manifest非重叠检查只支持其已扫描metadata范围；不等于历史真值重新盲化。最终评估要先seal预测，再独立offline评分；当前run_c自身不执行封存或GT评估，不能由GT_read=false元数据独自证明未来整个流程隔离。

## 实际只读核验结果

- 五组新manifest合计736 unique parents、1472 unique row prior/noise seeds；内部重复均0。与扫描到的历史R5B/R4/material-P2/joint-P2 config/manifest中的parent标识及seed值交集均为空。R4有全局generator seed而非逐parent seed的配置；本结论限于实际已声明的随机流身份，不声称检查了未知未落盘seed或逐truth内容。
- 账本snapshot序号31,200对应一致前缀：195,736 RHS=45,112 forward+150,624 derivative；5639 state/LU reservations；顺序、行累计、snapshot累计和active集合无差异。该时刻active为空；不是整个R6封账证明。
- 80个已完成CDEV CONTEXT只读取ZIP/NPY头，shape/dtype和报告路径已保存；没有读取truth/评分数组。
- `generate`与`run_c`目前stage路由均为`C if group=='CDEV' else D`。若后续运行C_TRAIN/ACQ_DEV，须确认这就是冻结费用归属；否则它们会消耗D cap，而不会自动归到C/B。当前CDEV路由正确。

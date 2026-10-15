# Research T → G 物理交接接口

合同版本`research-t-g/1.0`，严格保留原始`PRO_G_INTERFACE.md`的candidate-only语义。
此目录中的接口代码与纯数组测试不是已发现候选；正式候选只能由父实验运行器将真实端点写入候选子目录。
读取或校验交接包不会启动任何Maxwell求解。G目录未被访问或修改。

固定接口预期为24×24二维TM/VIE、所有8源共享单个复材料场、已知共享实增益1。
保存的原始complex channels为`(receivers,8)`。每源先receiver Re再receiver Im，再拼8源；
`sigma`是complex RMS，故`W=sqrt(2)/sigma`。具体solver、source、mesh、自项、时间谐波、
source revision、noise realization和tau由父实验的冻结target原样写入manifest，不能由本接口猜测。

实现：

- `src/research_t/structure.py`：绝对空气/材料相、二维Betti、阈值敏感性、固定材料距离、protected cavity edit和core/guard约束。
- `src/research_t/handoff.py`：真实保存端点的`write_bundle`与零求解`validate_bundle`。
- `new_measurement(target, chi, budget=..., acquisition=...)`、`correct_candidate`、`validate_candidate`只转交显式budget到父实验的计费Target。
- `FrozenPhysicsInterface`将实际A30 model、唯一Book、原观测、sigma、tau和固定CorrectorConfig封装。
  `predict`只调用`model.full_state`；`correct/compensation_score`复用Research T common corrector；
  `validate`对真实owned endpoint执行代数合法性、constitutive和残差检查，材料数组输入则先支付一次full-state。
  外来state、其他budget或改变冻结acquisition均被拒绝。`last_cost_receipt`保留原始RHS role和ACTION counts。

正式bundle包含`manifest.json/material.npz/observation.npz/prediction.npz/structure.npz/lineage.json/cost.json/legal_domain.json/distance_definition.json`。
`optional_state.npz`仅在后端实际暴露内部场/电流时保存；缺失时manifest为null并注明原因，不创建零数组。
预测缺失时只能NOT_VALIDATED或REJECTED；不能把低残差或box feasibility自动升级为full-validation。
`FULL_VALIDATED_WITH_BOUND`必须额外具有被论证的output error bound，backward residual不替代此界。

本地验证（不接触B/G或新增物理）：

```sh
PYTHONPATH=src python research/research_t_20261015/T_PHYSICS_HANDOFF/validate_handoff.py PATH_TO_REAL_BUNDLE
PYTHONPATH=src python research/research_t_20261015/T_PHYSICS_HANDOFF/test_semantics.py
```

校验器会复算已保存数组的白化残差，检查单一材料、源数、绝对空气mask、protected guard、概率null和路径存在。
它证明文件与语义一致，不是第二个独立Maxwell求解器。所有新物理必须经唯一Target的状态/RHS/wall三重预算执行。
新measurement action若不在原冻结target支持范围内，不能靠修改acquisition字典静默扩展物理，须由Target明确拒绝。

必要的7项语义检查已运行通过，新增full-state和full-solver RHS均为0，详见`INTERFACE_TEST_RESULTS.json`。
这些检查不构成真实基础解qualification或compensated nucleation证据。

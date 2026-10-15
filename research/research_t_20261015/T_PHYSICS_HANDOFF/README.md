# Research T → G 物理交接接口

合同版本`research-t-g/1.0`，严格保留原始`PRO_G_INTERFACE.md`的candidate-only语义。
本目录已交付真实保存端点，入口为`CANDIDATE_INDEX.json`。
8个最低残差历史base、8个minimum-change默认base及三轮32个合格补偿端点，共48个真实bundle。
它们是8份观测下的材料端点，不是48个已证明独立的世界或后验模式。
`MATERIAL_ARCHIVE.json`和`PAIRWISE_DISTANCES.csv`在同一观测内按冻结材料与相位门槛去重；原始端点全部保留。
接口代码的合成测试仍不是候选，测试fixture只存在于临时目录。
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

必要的8项语义检查已运行通过，新增full-state和full-solver RHS均为0，详见`INTERFACE_TEST_RESULTS.json`。
这些检查不构成真实基础解qualification或compensated nucleation证据。

G消费端必须新建其真实授权的Book/model，构造时传`consumer_role="G"`、
`authorization_identity=其授权记录`与`producer_ledger_root=T历史账本根目录`。
接口拒绝把T的Book作为G预算；关闭后的T run不能重新支付物理。读取旧bundle和保存数组重放仍为零新增物理。

## G应使用哪个base

使用`CANDIDATE_INDEX.json["default_bases_for_G"][observation_id]`，对应`__minchange_base`。
选择规则为：在已支付、合法且观测相容的原N3与两qualification端点中，选相对原N3材料体积L2变化最小者。
该规则不读取真值。它在本轮失败诊断后提出，属于exposed development修订，不是独立盲验证。
旧`__base`仍完整保留，代表lowest-observed-residual选择；它的数据拟合更好，
但离线材料MAE更差。材料MAE只用于报告，不作为在线排序输入。

## 零求解加载

```python
import json
from pathlib import Path
import numpy as np
from research_t.handoff import validate_bundle

handoff = Path("research/research_t_20261015/T_PHYSICS_HANDOFF")
catalog = json.loads((handoff / "CANDIDATE_INDEX.json").read_text())
observation_id = "p250065_ring_a1_snr20"
base_id = catalog["default_bases_for_G"][observation_id]
entry = next(r for r in catalog["candidate_inventory"] if r["candidate_id"] == base_id)
bundle = handoff / entry["bundle_ref"]
assert validate_bundle(bundle)["valid_bundle_semantics"]
manifest = json.loads((bundle / "manifest.json").read_text())
with np.load(bundle / "material.npz", allow_pickle=False) as material:
    chi = material["chi"].reshape(-1).copy()
with np.load(bundle / "observation.npz", allow_pickle=False) as observation:
    y = observation["y"].copy()
    sigma = float(observation["sigma_complex_rms"])
with np.load(bundle / "prediction.npz", allow_pickle=False) as prediction:
    saved_channels = prediction["channels"].copy()
```

这一步不实例化solver、不支付新状态或RHS，也不读取任何truth。
`optional_state.npz`包含实际保存的mass-normalized current/total field，不包含LU因子。
以后需要新修正时，不能免费重造一个LU；必须在G自己的新授权下支付matching full endpoint/factorization。

## 新物理与common corrector

`backend_assets/background/{weak,ring}.npz`只含已支付、已知背景的物理数组和几何；
不含对象真值、teacher J或oracle step。`BACKEND_ASSETS_PROVENANCE.json`明确历史setup成本未量化时为null。
复用原A30 revision `717d3d26a290248bd2b9aeb55f3d03694b407b84` 的
`a30.physics.FrozenMaxwellAdapter.from_background_npz`，传入G独立授权的budget book。
该Book须实现`counts/span/reserve_state/reserve_rhs/rows/check`，并在solve前同时执行其状态、RHS及wall上限。

```python
from a30.physics import FrozenMaxwellAdapter  # 原A30/src在环境路径中
from research_t.corrector import CorrectorConfig
from research_t.handoff import FrozenPhysicsInterface

# g_book由G在其真实独立授权下创建；此处不创建、不复制T的Book。
model = FrozenMaxwellAdapter.from_background_npz(
    handoff / "backend_assets/background/ring.npz", g_book,
    allowed_cache_paths=(), purpose_prefix="G_AUTHORIZED_PHYSICS")
config = CorrectorConfig(**json.loads((handoff / "DEFAULT_CORRECTOR_CONFIG.json").read_text()))
physics = FrozenPhysicsInterface(
    model, g_book, y, sigma,
    target_identity="A30-original-24x24-scalarTM-v1",
    observation_id=observation_id, tau=manifest["data"]["tau"],
    fixed_corrector_config=config,
    authorization_identity=g_authorization_identity, consumer_role="G",
    producer_ledger_root=handoff.parent)
state = physics.predict(chi, budget=g_book)  # 明确新增一个full-state与全部8源RHS
receipt = physics.validate(state, budget=g_book)  # 同一个owned endpoint不重复solve
```

protected edit的`correction_constraints.npz`保存实际`lower/upper/core/guard`；V2还保留其额外phase-lock。
将`{"lower": lower, "upper": upper}`传给`physics.correct(..., protected_geometry=..., budget=g_book)`，
禁止悄悄去掉guard或只对某个arm加入更强full-space fallback。
`physics.compensation_score`使用完全相同物理Jacobian和约束；无uniform nonlinear error bound时只给局部线性诊断。

这里的acquisition只支持原固定1频率、8源、32receiver配置。
`new_measurement`函数仅计算该冻结配置下的新材料prediction，不会产生真实新观测或新噪声样本。
新频率、新源/receiver或新传感数据需要另一个明确授权且冻结的物理合同；当前facade会拒绝静默更改。

## 成本与隔离

每个bundle的`cost.json`列dedicated trial费用、共享setup/parent replay、其它base-selection arm的开销及生产ledger行范围。
`setup_seconds`未单独隔离时为null，不能理解为0；所有失败与未导出的失败尝试仍在原始`COST_LEDGER.jsonl`及各轮SUMMARY中。
同一shared receipt可被多个候选引用，不可简单逐候选相加后重复收费；也不能只算成功端点而忽略失败trial。
`EXPORT_RECEIPTS.jsonl`保存每次零求解导出本身的wall和缓存矩阵应用，既不重置也不借用历史预算。

`e0/OFFLINE_AUDIT_MATERIAL_ARRAYS.npz`、离线评分脚本和truth图属于隔离的审计材料，
不是G在线接口输入。G在线加载应仅使用本目录candidate、固定known-background assets与原观测。
所有权重为`candidate_only/null`，没有posterior mass或complete coverage声称。

## 当前接口代码版本

G应使用仓库当前src/research_t。主证据zip保留旧执行快照，若直接提取全部源文件，请最后应用同Release的interface-overlay.zip；或用unzip -n不覆盖已有checkout文件。详见[交付说明](../DELIVERY_NOTES.md)。NumPy材料数组的.data buffer不被当作full endpoint；这一纯类型修复未重跑任何物理。

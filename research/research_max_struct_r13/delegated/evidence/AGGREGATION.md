# R13 封存后研究聚合：接口与边界

入口为 `src/aggregate_study.py`。它只消费已有离线评分JSON/CSV，不打开材料NPZ、预测NPZ或真值NPZ，不调用求解器，也不修改运行代码/配置/队列。当前仅执行已暴露的 `development_v1`、`development_v2`、`representation_v1`；产物在 `analysis/delegated/development_aggregation_v2/`。

## 已实现内容

- 合并所有估计器行，保留raw、corrected、primary、secondary、feasible_raw_retain及原评分表的initial/secondary_only。缺失结果不从预期分母删除。
- 保留每物理域、host、campaign独立统计；核心48可按matched/cross-grid、照明/孔径/噪声合同分层。不得把900运行之类数量写成900个独立材料。
- 配对表逐对象保留双方实际值与差值。bootstrap单位为材料对象，重复条件在对象内先平均；固定种子131313、4000重复，区间与Wilcoxon为未校正描述性统计，不新增成功门。
- NN固定16对子：`nn_dr↔dr`、`nn_g19↔g19`、`nn_direct↔direct`，对应同一个原观测和预先固定ID。Apple上的NN前处理费用单列，历史训练费用未知；不把Apple与i9秒数合成同host速度比。
- 在线控制可对实际存在的online_event/random/double/noevent/reset方法与相同16对象DR配对。若某控制没有实际评分记录，会明确列入 `ONLINE_CONTROL_COVERAGE`，不虚构执行结果或扩大样本。
- 条件策略按**全48**原G19对象构造：未触发沿用原G19 primary；触发使用该rescue方法的实际估计器结果。缺失触发结果保持不可用，不能悄悄改回parent掩盖缺失。每策略费用为全48parent+实际触发rescue；竞争策略总费用不能相加冒充研究总耗用。
- Apple `confirmation_local_reference16`是独立host/工作完成度诊断，始终与i9分开；相同16对象不能再计作新的16个独立材料。
- 图表只用真实评分表和费用记录，提供配对误差差值/区间、正负变化对象数及同host费用—材料散点；不制造材料图片、插值路径或GT选根。

## 独立数据授权

任何非development_/representation_前缀的campaign，必须先提供父线程授权JSON。授权检查在读取任何独立评分表之前发生；每个source还必须有 `SCORING_RECEIPT.json`，且其指向的预测seal仍通过授权、终态队列和mtime/字节检查。

授权文件示例只说明schema，不是已经授权：

```json
{
  "authorized_aggregation": true,
  "campaigns": [
    "confirmation_core48",
    "confirmation_mechanism16",
    "confirmation_rescue48",
    "physics_v1",
    "confirmation_local_reference16"
  ],
  "family_labels_authorized": false
}
```

父线程在所有必需预测封存、评分完成后才创建真实授权文件。需要家族分层时，还必须显式授权 `family_labels_authorized=true`、包含已评分core48，并通过 `--family-metadata`传入评分后允许读取的家族表/生成manifest；默认不读家族标签。生成manifest中的其他私有字段不会导入统计表。

后续调用形式：

```text
python src/aggregate_study.py --root R13_ROOT \
  --campaign confirmation_core48 \
  --campaign confirmation_mechanism16 \
  --campaign confirmation_rescue48 \
  --campaign physics_v1 \
  --campaign confirmation_local_reference16 \
  --authorization PARENT_AUTHORIZATION.json \
  --out R13_ROOT/analysis/delegated/study_aggregation_v1
```

输出目录必须为新的目录。缺失source评分或失效seal会停止聚合并保留原资料，不自动开真值、补跑或改算法。

## 主要输出

`ALL_ESTIMATOR_ROWS`、`STRATIFIED_SUMMARY`、`PAIRED_OBJECT_CHANGES`、`PAIRED_OBJECT_BOOTSTRAP`、`CONDITIONAL_DEPLOYMENT_ALL48`、`CONDITIONAL_DEPLOYMENT_COSTS`、`SOURCE_RUN_COSTS`、`SOURCE_COST_CONTEXT`、`SOURCE_QUEUE_ATTEMPT_HISTORY`、`SOURCE_ANALYSES`均提供JSON/CSV。另有生成范围回执与机械说明报告；独立NN条件可追加 `NN_PREPROCESS_COST_CONTEXT.json`。图片为实际表中数据的PNG/PDF。

费用来源保留原评分COST_LEDGER；不把内部CPU、外部CPU、累计history或辅助delta重新相加。中断/失败尝试的完整账仍在source评分产物；聚合不把未知项补零。条件部署的parent成本是策略成本重用，不能误称每条rescue都重新执行parent。

## 术语修正

R8 `global_dual`是global-PR的Dual-CIE，不是历史物理total-field-Gram GlobalDual。父线程已核对原 `full_maxwell_corrector.py` 的物理GlobalDual使用 **beta=0**；beta=.2属于feedback形式。此次local96 reference采用父线程冻结的历史物理GlobalDual补丁（来源 `a1317c738`，ready配置 `749703175`）。本聚合脚本没有修改该活动代码/配置，不引入新参数搜索。

科学结论、负例解释、图像评价及发表由父线程复核；聚合成功只说明数据整理与配对口径已经落实。


## 父线程验收修正：跨campaign共同16对象

源码现已将 `confirmation_core48` + `confirmation_mechanism16` + `confirmation_online_controls16` 合并为在线控制配对集合。C3配置的真实方法名为 `online_double`、`online_noevent`、`online_reset`；它们连同mechanism中的 `online_event`、`online_random`分别对普通DR配对，四个非event方法另对online_event配对。每一比较固定原16个ID，保留全部5个stage/deployment估计器视图；缺失保持不可用，不因campaign不同被丢掉。

源码另将 `confirmation_local_reference16` + `confirmation_local_trf16`组成 `LOCAL_REFERENCE_COMMON16_ROWS/SUMMARY`：cie、dual_cie、g19、dr、direct、physical_globaldual、physical_trf共7方法，全21配对；同为Apple、360CPU上限、同原始16对象。TRF额外的已暴露B对象明确排除出共同16分母，仍保留在自己的source campaign。派生比较集合只复用来源行，不新增研究运行或费用。

此次仅修正源码与说明，核对的是公开冻结config的方法名/ID/预算；没有读取新评分或执行独立聚合。后续父线程调用应同时列入C3与local_trf16两个source campaign，并在授权文件中覆盖它们。

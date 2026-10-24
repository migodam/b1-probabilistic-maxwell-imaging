# 冻结概率模型的新联合样本检验

本表来自已关闭的独立评分进程：L1为512个parent，L2为256个parent。每个family先在同一parent内平均训练seed的损失，再以parent为统计单位；不是把三个seed当作三倍独立对象，也不是先做图像ensemble。模型顺序、K128、CFM NFE32及展示的前两个L1对象/样本索引0..3均保持冻结规则。

## 完整66维密度的相对对数分数

26个模型中，14个显式密度有proper NLL，12个CFM的NLL为NOT_RUN，而非零。下列差值是冻结比较中的平均NLL减少量，单位nat：

| 冻结比较 | L1 | L2 |
|---|---:|---:|
| Raw spline相对容量匹配diagonal Gaussian |5.568208 |6.171013 |
| TRAIN Fisher frame triangular相对固定random frame triangular |1.357982 |1.601886 |

依主线程的条件proper-score解释，在相应可积条件下，这些绝对NLL差估计了平均conditional forward-KL的相对减少。真正条件熵在比较中相消；绝对KL仍未知，NLL减少百分比也不是未知绝对KL减少百分比。完整模型、seed及原有bootstrap区间见`MODEL26_BY_COHORT.csv`和`FROZEN_PAIRED_COMPARISONS.csv`。

Frame/data信息单列于`FRAME_AND_DATA_INFO.csv`。Raw family使用TRAIN32768/VAL2048；anchor CFM使用TRAIN1024/VAL256。TRAIN Fisher frame取自256个TRAIN anchor的全66维likelihood Gram，不是默认坐标前16维，也不是逐测试对象的新frame。Triangular仍有全部66维随机innovation；不能把strong block解释为独立先验填充。CFM bounded metric的random alignment是另一种control，不是该global random frame。

## 材料生成law的proper sample score

| CFM family | L1 material energy | L2 material energy |
|---|---:|---:|
| Anchor baseline |0.15893209 |0.16803881 |
| Bounded geometry |0.13066100 |0.13882215 |
| Random alignment |0.16530965 |0.17458628 |

Geometry相对baseline的绝对ES减少为0.02827109/0.02921667。在主线程列明的energy-distance恒等式下，对应的平均conditional material energy-distance减少量估计为其两倍，即0.05654218/0.05843333。这里不把score降低百分比解释为未知posterior-distance降低百分比。材料plateau原子不妨碍该proper-score比较；它也不确定latent fiber内的law是否正确。

`FAMILY_SEEDLOSS_BY_COHORT.csv`保留所有family的材料ES、固定K样本均值风险、相事件及覆盖率。Raw、anchor及single-seed refinement的数据/信息预算分开列示；没有根据GT改选模型或样本。

## 必要的Maxwell校准检查仍严重不相容

下表仅为预先指定四个S1生成器的fixed sample0 Full Phi均值；Gamma256参考均值为256：

| S1 q-model | L1 | L2 |
|---|---:|---:|
| Raw spline |8371.95 |8662.42 |
| Anchor CFM baseline |36282.18 |40924.89 |
| Geometry CFM |22139.70 |23557.96 |
| Random CFM |47633.59 |49493.64 |

PIT与固定前16个L1对象、前64 draws的数据依赖likelihood-rank图保留在`figures/final_probability/`。未做该扩展检查的L2或其他模型标为NOT_RUN；不补零、不删失败、无残差筛选。相对NLL或ES收益与这些绝对校准不相容同时报告，不能升级为FULL posterior合格。图中样本明确标注q-model/UNQUALIFIED。

闭合的三组paired VI验证曲线显示Phi减少同时proper NLL及熵发生明显tradeoff，三组预注册replication gate均失败。这里只读取既存历史并画图，没有重训或改变选择准则。自适应CFM诊断的原始partial结果及报告以`ADAPTIVE_CFM_*`附录保留：895s预算内仅完成3个solve/cycle；已完成的两容差NLL差0.13599nat，未达到0.05阈值；未完成的anchor模型不冒充已测试。

## 来源顺序与本次成本

`CHRONOLOGY.json`记录两cohort的density freeze均早于首份新observation；两cohort全部online输出均早于第一次GT评分gate；inference/score PID不同；所需三个density/Full评分Book均CLOSED COMPLETE。本次只在全局隔离评分COMPLETE后读取结果和固定展示GT。最后一个representative离线schema修复由主线程另行记录于`AUDIT/FRESH_REPRESENTATIVE_SCORE_SCHEMA_RECOVERY.json`，原失败及预测保存；没有影响这些已关闭概率分数。

本次表格/作图Book IV：4.553736s controlled wall，3.462356s CPU，peakRSS588595200B；新增Full RHS0、训练更新0、重复评分0。26个checkpoint的byte/mtime身份保持不变。主线程proper-score说明原文保留于`PARENT_PROPER_SCORE_INTERPRETATION.md`，最终科学判断归主线程。

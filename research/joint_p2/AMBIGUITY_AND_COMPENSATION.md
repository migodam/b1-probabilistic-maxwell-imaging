# Maxwell 材料歧义与可允许补偿

实验阶段20261013；公开交付2026-10-14。Fresh主概率参考与范围以 [ACCEPTED_REFERENCE_AUDIT](results/ACCEPTED_REFERENCE_AUDIT.json) 和 [全部fresh细化记录](results/ALL_FRESH_REFERENCE_REFINEMENT.json) 为准；旧13阶与中间结果保留为历史近似。

结构不同不必然意味着散射可区分。P2固定参考实心球 ε=2.5+.1i、gain=1，用声明范围内的另一个几何世界进行ε和gain补偿。最终distance都由full vector-Mie endpoint算得；proxy仅指导小维度优化。噪声由固定实心响应RMS的.12倍冻结，没有随候选变动。

| 空腔 | 频率kR | raw noise-distance | 仅ε补偿 | ε+gain补偿 |
|---|---|---:|---:|---:|
| 中心r=.08R | 1.6 | .056838 | .049314 | .035593 |
| 中心r=.08R | 1.6+2.4 | .125981 | .108569 | .107504 |
| 中心r=.25R | 1.6 | 1.529652 | 1.295221 | .900170 |
| 中心r=.25R | 1.6+2.4 | 3.002910 | 2.548344 | 2.467944 |
| 环形halfwidth=.08R | 1.6 | 13.600991 | 7.285575 | 4.450536 |
| 环形halfwidth=.08R | 1.6+2.4 | 19.924801 | 10.562367 | 8.022080 |
| 中心r=.20R+环形h=.06R | 1.6 | 11.076859 | 5.107862 | 3.725769 |
| 同上 | 1.6+2.4 | 16.414256 | 7.568978 | 6.570264 |

最清楚的可行witness是中心r=.25R：单频下补偿ε≈2.466617+.130009i、gain≈1.053221，gain prior z≈.355，处于prior内部。原whitened差1.52965下降为.90017，吸收约41.15%；归一化体积材料距离仍为.12846。候选具有实质材料差别，但在冻结联合实通道噪声尺度下response接近。它不是严格相等，也不是全局最小值。

增加预先固定的2.4频率，且坚持**同一个gain跨频共享**后，最佳已找到补偿距离升为2.46794。模型因此给active sensing方向一个实证起点。小空腔在两频依然很难区分，说明新增频率并非普遍解决全部歧义。

环形与prior-boundary stress cases是负结果。很多优化endpoint碰到ε或gain边界，即使允许补偿仍留下多个噪声尺度差异。不能因此称所有不同结构都会形成两个高概率世界；不能把达到边界的最佳residual解释成高evidence，更不能通过放大noise或放宽prior把它改成成功。所有参数位置、边界距离、材料距离和proxy/full差保留于 [COMPENSATION.json](results/COMPENSATION.json)。

`gain_negative_control`的truthgain=1.2，而A/B仍假设校准gain=1，因此它们在该case是**故意错误的known-calibration模型**，不是知道truthgain的公平informationablations。C/D的properboundedprior包含1.2，用于明确nuisance误设的影响；negativecontrol不混入普通case的known-versus-unknown优越性比较。

## 从witness到材料事件posterior

独立reference处理的是完整 prior-weighted likelihood，而不是上述最佳endpoint。预注册large-cavity observation，在未知ε与bounded gain的C条件下，full-Mie residual reference保留：

| 材料事件 | prior mass | posterior mass | residualMC SE |
|---|---:|---:|---:|
| 实心，event0 | .250000 | .162540 | .00002626 |
| 中心半径(.10,.20]且无环形空腔，event6 | .055788 | .071828 | .00001224 |
| 中心半径(.20,rmax]且无环形空腔，event9 | .180888 | .677226 | .00004582 |
| 无中心空腔、环形halfwidth(0,.06]，event1 | .033203 | .041259 | .00003132 |

前三个为预注册.05重要质量阈值以上的材料区域。它们不是三张精确图像的点概率，也未被认定为三个分离modes；其中可能有连续ridge和相邻空腔族。Top4累计.952852，剩余已建模质量.047148，未重新归一化。

在真实材料为solid的另一预注册observation中，event0质量.407667，event9质量.410889，event6质量.116300。单MAP无论选哪一个，都无法传递另一显著材料解释。truth仅作模拟与评价，不进入posterior feature。

## 信息消融不是算法排名

下表取同一large-cavity observation的连续proxy积分；不同列有不同target identity，不能互作方法优越性排名。

| 条件 | 实心event0 | 中等中心空腔event6 | 大中心空腔event9 |
|---|---:|---:|---:|
| A：ε/gain已知 | .257488 | .115363 | .611768 |
| B：ε未知、gain已知 | .048169 | .034325 | .775068 |
| C：ε+bounded gain未知 | .162591 | .071867 | .677377 |
| D：C加第二频率 | .000431 | .003754 | .966140 |

ε自由后实心质量下降、gain自由后又上升，显示信息条件变化并非“参数越多每个事件质量越大”。posterior由观测、prior体积与proper nuisance共同决定。D中重要质量更集中，大空腔仍被支持，但D此表是proxy numerical target；不是对所有条件完成full-Mie probability calibration。raw secant和full corrected C mass分别给物理与概率证据，二者各有明确范围。

局部projected Jacobian/Fisher诊断保留五个半径的值，用来判断一阶response被ε/gain方向吸收的程度；没有仅凭局部投影宣称大几何非线性等价，也没有在本实验建立普遍r⁵或其他幂次定律。

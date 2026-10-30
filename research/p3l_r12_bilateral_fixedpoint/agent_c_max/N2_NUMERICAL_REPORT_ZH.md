# N2 数字离线数值与图像交付

N2 的 9 个设计数字、6 条方法共 54 组封存预测已经离线评分。**主材料覆盖全部为 0/9；数字 6 在每条方法中都是空合格云，仍计入覆盖分母。** 本工作没有新增 Maxwell 求解，也没有改变在线选点、模型、方向或停止规则。最终视觉解读和科学结论由主线程负责。

## 阅读与评分顺序

1. 原有 `agent_c/offline_score.py blind` 使用封存预测生成 12 张匿名材料/occupancy PNG；没有接收 GT。
2. 主线程实际打开全部 12 张图，并写下 [`N2_MAIN_BLIND_READ.json`](../audit/N2_MAIN_BLIND_READ.json)。
3. 主线程显式授权后，原评分器执行一次 `score`，输出至 [`n2_score`](n2_score/)。其 [`TRUTH_ACCESS_BEGIN.json`](n2_score/TRUTH_ACCESS_BEGIN.json) 与 [`SCORING_COMPLETE.json`](n2_score/SCORING_COMPLETE.json) 保存访问及完成回执。
4. 最终图从同一封存数组展示 GT、共同 anchor 和六种方法的在线输出；GT 最近成员只用于离线云容量统计，从未替代部署输出。

第一次用系统 Python 启动时缺少 numpy，运行在读取预测前退出；随后使用既有 `Gaussian/.venv_nn/bin/python` 成功执行。原失败保存在 [`BLIND_STARTUP_FAILURE.txt`](BLIND_STARTUP_FAILURE.txt)，没有安装环境或更改数值程序。

## 主材料结果与在线失败

|方法|ε=.075 覆盖|合格云可用|选点真实 noise 通过|选点引擎集合资格 Q|平均选点 MSE，全部 9 个对象|材料拓扑 exact|occupancy 拓扑 exact|
|---|---:|---:|---:|---:|---:|---:|---:|
|H0_GN|0/9|8/9|8/9|8/9|0.02880951|5/9|3/9|
|H0_DATA_RML|0/9|8/9|8/9|6/9|0.03036269|5/9|4/9|
|R9_STYLE_FIXED8|0/9|8/9|8/9|8/9|0.03171430|5/9|4/9|
|H1_FULL16|0/9|8/9|8/9|8/9|0.03171430|5/9|4/9|
|H2_FULL16|0/9|8/9|8/9|8/9|0.03171430|5/9|4/9|
|HA_FULL16|0/9|8/9|8/9|8/9|0.03171430|5/9|4/9|

MSE 使用原始 576 像素的复材料均方误差。N_005 对应数字 6；该对象所有合格集合为空，真实 noise 门也失败。原评分摘要把包含此空集的整体最近成员 MSE 保持为 null，没有用可用的 8 个对象平均值替换。所有在线输出都保留在图和选点误差中。

`Q` 是引擎的 noise、合法性和贪心材料 RMS .002 去重的合取；`N` 才是实际 noise ball。`Q0` 可以只是 `GREEDY_DUPLICATE`，不能直接叫作 noise 失败。H0_DATA_RML 的选点 Q 通过 6/9，但真实 noise 门通过 8/9；其他方法真实 noise 门同样为 8/9。该区分保存在 [`N2_METHOD_SUMMARY.csv`](N2_METHOD_SUMMARY.csv) 及选点表。表内拓扑是全部 9 个在线输出的描述，包含真实 noise 失败的数字 6。

材料拓扑阈值为 Reχ≥.75；occupancy 阈值为 .5。前景 4 邻接、背景 8 邻接，计入所有孔洞，不平滑、不膨胀、不新增小孔过滤规则。拓扑 exact 只判断连通分量数与孔洞数同时相等，不等于像素级或材料恢复。

逐对象数值在 [`GROUP_SCORES.csv`](n2_score/GROUP_SCORES.csv)、[`TOPOLOGY_SCORES.csv`](n2_score/TOPOLOGY_SCORES.csv)；补充的结构与原失败旗表在 [`N2_SELECTED_STRUCTURE.csv`](../figures/n2_online_and_clouds/N2_SELECTED_STRUCTURE.csv)、[`N2_CLOUD_STRUCTURE.csv`](../figures/n2_online_and_clouds/N2_CLOUD_STRUCTURE.csv) 和 [`N2_CANDIDATE_STRUCTURE.csv`](../figures/n2_online_and_clouds/N2_CANDIDATE_STRUCTURE.csv)。

## 数字 8：已评估有限候选没有恢复双孔

GT 在两个冻结阈值下均为 1 个分量、2 个孔洞。每条方法保存 7 个 raw 候选，共 42 个；**这 42 个候选在材料 .75 或 occupancy .5 下都没有出现 2 孔。** 此结论只描述这些已评估点与阈值，不证明 decoder 或整个可行集合无法表示数字 8。

|方法|raw 候选|合格 raw|noise 失败 raw|合格材料孔洞数|合格 occupancy 孔洞数|离线最近材料 MSE|
|---|---:|---:|---:|---|---|---:|
|H0_GN|7|4|0|0 或 1|0 或 1|0.03193057|
|H0_DATA_RML|7|5|0|0|0|0.03521366|
|R9_STYLE_FIXED8|7|3|4|0|0|0.03726738|
|H1_FULL16|7|3|4|0|0|0.03632093|
|H2_FULL16|7|3|4|0|0|0.03632093|
|HA_FULL16|7|3|4|0|0|0.03632093|

raw 合格计数保留引擎资格旗，不能与有限集合去重计数混称。raw 序号、Q/noise/合法旗均保存；Q0 图框标红，未删除失败点。四条 weak 方法的数字 8 在线输出仍为原 anchor；H1/H2/HA 的较近离线成员没有被偷换成输出。

## 原盲图材料视图的显示限制

主线程原盲读记录中 A05/A06 的 real-material 图近乎全紫，而 occupancy 可辨。解盲后 A05 为 R9_STYLE_FIXED8、A06 为 H2_FULL16。为解释该现象，仅检查了现存选点数组和 decoder 公式，没有执行新物理。

decoder 为 χ=host+δa。四条 weak 方法所选材料的最大实部按对象落在 1.144–1.464，最大 `|χ−host−δa|` 为 1.12e−16；N_007 的最大实部为 1.196，occupancy 最大 .656，δ≈1.671+.139i。因此原图近空白不能解释为存储材料幅度塌缩。本子线程打开原 A05 PNG 时能看到数字，与主线程原观察不一致，归为原盲读显示限制，未改写主线程盲读结论或补记视觉 PASS。原始盲图仍保留。

最终比较图直接从封存的 `materials` 和 `occupancy` 数组作原 24×24 显示，统一材料范围 .1–1.9、occupancy 范围 0–1，没有按图自动拉伸或插值。数字轮廓可辨不能覆盖主材料门失败及孔洞错误。

诊断数值保存于 [`SELECTED_MATERIAL_AMPLITUDE_DIAGNOSIS.json`](SELECTED_MATERIAL_AMPLITUDE_DIAGNOSIS.json)；decoder 来源为 `shared/physics.py` 的 `decode`，本工作只读取其源码。

## 原始主线程盲读与解盲映射

[`N2_BLIND_READ_MAPPING.csv`](N2_BLIND_READ_MAPPING.csv) 将主线程回执中的 54 条原读数、0–5 分、原备注逐条复制，再附加解盲后的方法和真实数字。原始分数不改写；A05/A06 材料显示限制逐条保留。六条臂的同一 panel 具有相同原读数，不能补记为新的独立读者评测。

|真实数字|匿名 panel 后缀|原始读数|原始分数|
|---:|---|---|---:|
|1|01|1|4|
|2|02|2|3|
|3|04|5or3|2|
|4|05|4|3|
|5|08|3|3|
|6|06|6|3|
|7|09|7|4|
|8|07|0or8|1|
|9|03|9or3|2|

因此真实 3 仍是 5/3 含混，真实 5 被读为 3；不得在解盲后将这两项改成识别正确。原回执还明确说明：这是此前已知的 R11 benchmark，本轮匿名输出阅读，不是新的人类留出分类实验。

## 交付图

- [`DIGITS_123456789_RECONSTRUCTION.png`](../figures/n2_online_and_clouds/DIGITS_123456789_RECONSTRUCTION.png)：9 行数字 × GT、共同 anchor、六种封存在线输出。
- [`DIGITS_123456789_OCCUPANCY.png`](../figures/n2_online_and_clouds/DIGITS_123456789_OCCUPANCY.png)：相同布局的 occupancy。
- [`DIGITS_123456789_MATERIAL_MASK.png`](../figures/n2_online_and_clouds/DIGITS_123456789_MATERIAL_MASK.png)、[`DIGITS_123456789_OCCUPANCY_MASK.png`](../figures/n2_online_and_clouds/DIGITS_123456789_OCCUPANCY_MASK.png)：两个冻结阈值的原像素结构。
- 同目录 `DIGIT8_CLOUD_<METHOD>_<VIEW>.png`：每条方法的 GT、anchor、全部 7 个实际候选，四个视图分别为 real material、occupancy、材料 mask、occupancy mask。

全部 28 张图路径见 [`FIGURE_CONTENTS.json`](../figures/n2_online_and_clouds/FIGURE_CONTENTS.json)。这些图和数值是有限已评估集合的描述，不建立后验、总体高概率保证或独立确认结果；N1 原失败与 N2 算法预声明边界由主线程原样保存。

# R14 数据、物理和评价合同

## 分母与数据隔离

主实验为既有已曝光48材料×四metric=192运行。它是后续机制研究，不是48个新盲测。另有曝光16对象的self/global-Dual48运行、stage次序32运行、曝光四对象共享material solver16运行。次数不作为独立材料数。

新独立组为24个新seed材料、11方法、264运行；在算法/参数/次序规则/预算与公开选择器冻结后生成并执行。8种既有几何家族×low/high/mixed三个contrast regime；它检验新对象和物理条件，不宣称未见家族分布外。几何家族与contrast标签在预测闭合前保存在独立生成/评分区。4个cross-grid条件保留，不因模型误差而删掉。

公开输入只含case_id、原离散模型/几何、合法源接收器、Y、噪声与CN99半径、box/prior和初始材料。预测包不含GT材料、support、shape标签或生成私有assignment。在线runner不接受GT参数。先建立公开saved-view/prefix/portfolio选择，再读评分真值。

## 二维物理

exp(−iωt)，固定μ，实正介质χ=εr−1；576独立χ，box[.02,20]。生成材料的ROI host为.1，公开初值也为.1，但整个ROI仍然未知，没有按GT固定host区域或support。外部物理背景为χ=0。二维scalar-TM使用同原波数和有限disk-pulse Hankel–Green核；本文J=χQ是归一化对比极化量，物理对比电流另含−iωεb。源/接收器在未知ROI外；完整观测保持原near-field定义，未改称far-field证书。μ未知、一般tensor、连续PDE等不在本轮二维结论内。

新数据复用原算子与生成器：内部高分辨率几何216²做cell平均，inverse24²；四cross-grid由36²正演产生，inverse仍24²。没有用GT的support降低576自由度。孔径/6或12源/噪声条件按预先固定协议分配。实际逐对象网格、源、频率和noise以公开NPZ/manifest为准。

原物理目标为Φ=D+R，D=.5||F(m)−Y||F²/||Y||F²；R为原edge smooth/TV，参数不改。辅助source normalization与D的归一化分开。复噪声σ表示complex RMS，实/虚各σ²/2；CN99是原有限维复噪声球。cross-grid/实测误差不包含在这个噪声球中，不能将球内标签升级为模型误差证书。

## 三维物理

真实CM+RR TensorDDA/vector VIE，6³=216实isotropic χ、3N=648三分量current，自项s=−1/3+i k³V/(6π)，dyadic offdiag保留跨分量耦合。60个外部合法dipole源、180接收分量；box[0,8]、原3D prior/source metric/rho。非局部块沿每维两cell，即2×2×2体素、24分量块，不能套成2D四scalar行。

四公开旧配置是bracket、frame、steps三个独立几何，加同frame的8³→6³ cross-grid条件。6方法×4配置=24运行，不是24独立三维对象。三维CIE/Dual参考与四pair同host、同120总CPU、严格最多48corrector；和二维总120余额重分配规则分别报告。

继承物理核及既有Mie误差/网格限制，不将matched低残差称连续Maxwell成功。frame cross-grid clean模型差约8.0%，超过1%噪声，噪声球失败需按模型失配解释。central三切片用于展示，不是完整体积验证。

## 评分与选择

Primary为合法、已完成且已付费候选中的最小原Φ；secondary为原CN99球内最小R，无候选时保留primary回退并标注。raw保留原proposal，同时在解释层按公开Φ与initial比较；原数组不覆写。Corrected与primary分别记录，corrector的Φ下降不被标成材料改善。

离线材料RMSE使用全部未知量；amplitude NRMSE除以max(coarse truth)−host，另给relative L2与boxrange归一化。固定χ>1的IoU仅作该合成任务的支持描述；连续/多相材料仍看连续误差、幅度与整体结构。沿用NRMSE≤.15且IoU≥.70作后验联合描述，不调旧标准、不作为在线算法gate。IoU .5或NRMSE .156等边界例可保留主要结构，不能仅凭这个二元表宣称全结构错误。

统计单位是独立材料。新24的按对象配对差与bootstrap描述性区间不把迭代、source或多个算法当独立样本；同family/同三维几何重复条件明确标记。曝光48只作机制/开发证据，不能与新24混合后宣布确认显著。

## 来源

- 旧48与强基线：R13及`supplement_aplus48/`实际数组、封存和stack_comparison。
- 新24：本轮`generate_independent24.py`、公开manifest、独立封存与评分回执。
- 三维：R1/R7/R8真实TensorDDA和R13 vector adapter；本轮`audit/CROSS_PHYSICS_REUSE.md`逐API审计。
- 原ZIP全部文本读取范围与有限Chat读取：`sources/BUNDLE_*`；未取得的完整Pro对话不冒称已读。
- 新Fresnel/Austria/NN接口实验：NOT_RUN。本轮没有伪装外部实测或原论文作者复现。

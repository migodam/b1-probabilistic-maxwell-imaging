# Research W-GQ-R3-JI：未知墙–未知目标联合全波成像

[中文网页与真实图像](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r3_joint/) · [41组推导公式](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_gq_r3_joint/formulas.html)

已实现受限二维 scalar-TM 联合成像：六个墙宏格的复材料未知，目标为七个连续参数的椭圆／圆，不读取真值或旧标签库。三个独立 parent 首噪声 Full-GN 位置 RMSE 0.03234λ₀；一个 parent 同时通过位置、形状、实介电率门。目标损耗和其他形状／材料失败保留。Adaptive chart 没有超过普通 Full-GN 的同质量费用优势。

新增数字8为110变量真实材料反演，不使用数字模板；两孔未恢复，两起点未收敛。1–9仅8实际反演。不能将较小数据残差写成成像成功或物理不可辨识证明。

## 下载内容

- **RESEARCH_W_GQ_R3_OFFLINE_AUDIT.zip**：193,037,632 bytes，7552项。完整源码、必要历史后端、原始复观测、材料数组、全部分支、图件、失败和回执。解压目录为 `research_w_gq_r3_joint/`。
- 外部最终费用、打包和发布记录：保留科学封账与归档时间差，不回写原始历史状态。
- 主 `delivery/portable_runtime/tools/deploy_r3.py` 已实际执行测量到图像的独立运行；全新依赖安装未运行。数字8新预算独立部署封装未完成，原实验入口不能解除旧时间门后覆盖账本。

12,440 FE / 187,732 complex RHS；两个早期CPU/RSS记录缺失，未补造。已知process CPU 6468.60秒，打包另列。发布新增物理调用为0。

## 科学边界

有限开发样本、受限未知墙与椭圆目标；不是任意拓扑／一般未知墙成功，不是连续PDE认证。A17X继承和标准SOM/TSOM、Schur、resolvent归属明确。论文级新算法优势未成立。

原报告和机器裁决中的“未发布/未授权”字段描述科研封账时的状态。2026-10-30用户明确请求发布，本Release是后续公开记录，不修改原科学裁决。


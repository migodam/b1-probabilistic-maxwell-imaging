# Research T 的材料、结构与距离语义

本接口只讨论24×24二维栅格。`chi=epsilon_r-1`；严格真空编辑的core为`chi=0`。
绝对空气mask为`abs(chi)<=0.10`，并在`.075/.10/.125`三个固定阈值下报告拓扑敏感性。
`abs(chi-background)`只能描述相对背景的增量，不能作为空气或孔洞mask。

材料前景采用4连通，空气/背景采用8连通。二维`beta_0`为材料连通分量数，
`beta_1`为空气中不接触栅格边缘的连通分量数，即栅格上的封闭空腔。
这不是连续拓扑或网格无关证书；二维孔洞不能被称为三维隧道/空腔。

材料距离是`d=sqrt(sum(mass*abs(chi_a-chi_b)^2)/sum(mass))/1.2`。
归一化1.2固定，不依赖真值或候选。另列复材料L1、Re/Im MAE、绝对相位mask对称差、
封闭空腔mask对称差、Betti差、材料边界Hausdorff距离（cells）和空腔质心匹配距离。
相位差至少6个cell才满足显著栅格编辑的必要条件；调用方还必须执行其冻结的材料距离阈值。
单个像素或纯latent变化不能成为新的结构族。

`cavity_birth`的最小core半径由调用方显式冻结，函数不推断物理尺寸。
core固定为零；周围至少一圈完整guard要求`Re chi>=.30, Im chi>=0`，
core和guard均不得触及栅格边界。guard修复属于提出的材料编辑本身，须计入端点物理评估。
corrector必须使用`protected_bounds`，并在最终端点重新执行`check_protected`；
只固定core而不固定guard，无法排除把洞打通到外部的伪成功。

输入只有在线候选材料、网格和冻结阈值。任何truth-relative阈值（例如`.25*abs(truth increment)`）
均不进入本模块。真值只能由隔离的离线评分过程读取。

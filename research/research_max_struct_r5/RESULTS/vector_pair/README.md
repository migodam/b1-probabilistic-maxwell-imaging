# 固定两体素 vector Green 残差谱例

使用R1原TensorDDA：lambda=1、h=.04、epsilon=3、mu=1；两相邻立方体素中心沿x。电场自由度6，已知各向同性材料标量2；没有反演。两个+y传播平面波分别x/z极化，一次6×6 LU共享2 RHS。外部电场仅取沿三轴±.5的六个点。

读取的self=-0.333333333333+0.000842206242226i。密度Green offdiag纵向x=0.164102393376+0.000836898394107i；横向y/z=-0.0771825666985+0.000831602519009i。实际retarded Green的纵/横实部分别正/负；静电diag(2,-1,-1)仅作符号解释，未用于场解或权重。

m=epsilon−1=2；a=m/[1+(beta−gself)m]；T=I+a Goff，其中Goff是乘过体积的密度算子。以下权重均为实际复Green的 |1+a lambda|²。

| beta | 极化 | coherent | antisymmetric |
|---:|---|---:|---:|
| 0.2 | x_longitudinal | 1.34283686829 | 0.707604859536 |
| 0.2 | y_transverse | 0.856192833723 | 1.15496648488 |
| 0.2 | z_transverse | 0.856192833723 | 1.15496648488 |
| 1.0 | x_longitudinal | 1.18703265334 | 0.828991982169 |
| 1.0 | y_transverse | 0.917573017513 | 1.08597214283 |
| 1.0 | z_transverse | 0.917573017513 | 1.08597214283 |

两固定平面波只激励coherent x与coherent z；其余mode是完整残差空间的代数方向。纵向coherent与横向coherent的权重方向不同，不能据2D scalar的同相解释直接统一描述该vector算子。

kh=0.25132741，材料内部kh=0.43531185。这是粗离散两cell模型，不是解析球或连续PDE。继承R1 bounded-DDA适用界；历史Mie球结果不认证这两个体素。没有新增恒等式/导数测试，没有GT评分，不宣布3D imaging PASS。科学解释与新颖性由父线程审查。

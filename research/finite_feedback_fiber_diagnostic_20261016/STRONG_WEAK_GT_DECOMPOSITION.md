# GT差量的强/弱材料分解

仅ORACLE_GEOMETRY_DIAGNOSTIC。GT不进入任何在线方向、corrector或候选选择。

| Anchor | GT材料RMS | 白化一阶响应 | 实际full差异 | 弱能量@.001/.01/.1 | 原观测残差 |
|---|---:|---:|---:|---|---:|
| physics_R_base | 0.204978 | 115.283215 | 117.714236 | 71.91%/73.31%/76.59% | 120.432293 |
| N3 | 0.143047 | 16.401827 | 15.552557 | 98.09%/98.75%/99.06% | 28.034654 |


## 测度与阈值

材料坐标m=sqrt(v)[Re h;Im h]，v=1/576，因此范数为复材料RMS。J=W DF/sqrt(v)，W使用sqrt(2)/原复sigma。512行、1152列；数据packing每source的32Re/32Im。对J=U diag(s) Vᵀ，强投影V_s V_sᵀm，强能量为其平方范数/||m||²，弱能量为补集，包含矩形nullspace。相对cutoff .001/.01/.1预先冻结。不是参考后验measure，也不是“可自由生成”的方向数。

N3和physics基线的原P0 chart缺失且phase必要条件失败，故真正prior tangent占比NOT_ESTABLISHED。固定增量相位子空间只是必要条件：N3差量99.8896%、physics差量99.2166%的能量在该线性空间内，并不解决有限先验支持/峰值/rasterizer约束。TRAIN chart存在、rank7，但其hGT=0，能量比例null。

各切点的谱隙和可观测condition在METRICS.json逐项保存；不存在明显唯一cutoff。完整谱、V和GT系数在results/*_DECOMPOSITION.npz。Full J与旧N3独立缓存相对差6.75e-16。不能因全空间弱比例高就推定原P0里的真实等数据流形。

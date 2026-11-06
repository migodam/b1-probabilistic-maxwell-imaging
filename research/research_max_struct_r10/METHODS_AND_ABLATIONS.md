# R10 方法与对照

本轮保留充分运行的成熟CIE/FFT-TSOM，实现来源是R8/R9适配代码，而非作者原始代码。NIE未新增。下表列真实算法差别；Global-Dual只是内部机制控制。

|名称|辅助current|材料读出|角色|
|---|---|---|---|
|CIE-single|sourcewise PR，current β随材料|普通CIE|主要文献方法适配|
|Dual-CIE|sourcewise PR，current β=1|普通CIE|主要成熟强对照|
|Global-Dual|全源global PR，current β=1|普通CIE|同current消融|
|A+sourcewise|与Dual-CIE相同|局部one-hop Möbius/MM|隔离材料规则|
|A+global|与Global-Dual相同|局部one-hop Möbius/MM|R9继承方法|
|Green19矩阵分母|与Global-Dual相同|R_K，K=P_G G P_G，rank19|新增简单候选|
|Two-fold28矩阵分母|与Global-Dual相同|R_K，K=P_T G P_T，rank28|新增接收子空间候选|

两个P之间保留完整G；投影不修改原Maxwell forward。所有新材料仍576维，19/28只是分母低秩小系统的维数，不是材料参数化维数。Two-fold是公开G/S强右子空间构造，不能把其秩称为实际DPCS秩。

## 共同预算与选择

二维mf3/6/10各最多500步，mf12最多1500步；材料β=1/.2/.2/.2。现有系数早停设为0，至少50步，不设置45秒等历史wall截断。所有方法共享完整原物理L-BFGS-B corrector，最多3000接受迭代/12000函数评估，原停止准则保持。所有576个材料自由度、box、smooth/TV权重、观测与初始化相同。

在线选点只看已付费原Phi：初始点、阶段端点及corrector保留的Phi incumbent。GT只在封存后评分，不能为了表格择优改用材料最优迭代。辅助objective下降不等于Phi下降；sourcewise/global、dualβ、material metric通常没有统一全程下降目标。

固定迭代上限是共同允许量，不是等实际算量。报告完整child CPU、wall、576阶物理LU/正伴随RHS、FFT/current/材料内循环和19/28阶小LU。实际已保存预算点用当时已付费Phi incumbent，不能插值伪造算法达到的质量；未覆盖预算记缺失。外部完整进程CPU与内部prefix时间分别标记。

## 新分母与消融的边界

新增R_K=β[I+βD_m−D_mK]⁻¹[J−D_m(Einc+GJ)]，包括分母随m变化的精确导数。固定current时最多8次bounded L-BFGS-B材料迭代，以同一新辅助objective比较旧点、普通读出种子与优化点；失败候选保留并收费。原A+四颜色条件独立性不适用于非局部K，因此未复用该假设。

开发控制包括等秩Green28、K=GP固定状态控制、保留全部原Gii的Green19自项控制。它们只回答机制问题，没有根据fresh评分追加进确认集。矩阵分母不保证全局凸化，不保证整个材料box可逆；未证明部分不称证书。具体推导、反例及误差式见[理论](THEORY_DERIVATIONS.md)。

## 因果比较的层级

1. 相同(J,m,Y)下切换材料metric：隔离一步方向，不能预言全程basin。
2. 同current更新规则、不同材料读出：隔离算法规则，但后续J轨迹自然不同。
3. 同FFT/APCS空间和容量：排除把多加Fourier自由度当成材料规则收益；不能排除实际求解精度差异。
4. 同物理corrector：它的贡献必须单独报告，不能将最终全部改善归给辅助阶段。
5. 同新对象与真实费用：检验哪类材料值得采用，不能要求所有对象稳定获胜，也不能只展示获胜对象。

三维三臂分别sourceCIE、GlobalDual和AplusGlobal。GlobalDual/AplusGlobal为材料读出匹配对照；sourceCIE比较同时包含sourcewise与β差异。新Green19/Two-fold矩阵分母尚未移植三维，不能把二维新候选的结论转移给三维A+。

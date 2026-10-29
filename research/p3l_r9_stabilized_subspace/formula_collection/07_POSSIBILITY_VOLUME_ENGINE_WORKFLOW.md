# 07 可能材料云：新目标、算法与证据流

这是新 possibility-volume namespace 的源码流程说明，写作时没有读取新 truth 或运行新物理。它不改变03章保留的原 E2/E3与REV1。执行批准、最终协议冻结、生成、在线运行与封存由 root负责；统计配置独立冻结。

## 目标变化必须明示

原探索主要使用相对锚点 DeltaPhi≤2与DeltaPrior≤10评估局部候选。新任务的主目标是实际材料可能性覆盖：候选满足原绝对95%noise ball、原decoder合法性、指定有限材料半径和signed progress。旧相对门槛保留为secondary fields，不参与新primary reject。

\[
\Phi_y(x)=\|F(\chi(x))-y\|_F^2/\sigma^2
\le\tfrac12\chi^2_{512,.95},
\quad |\|\mathcal R_m(\chi(x)-\chi(x_0))\|-\rho|\le.02\rho,
\]
\[
u_0^T\mathcal R_m(\chi(x)-\chi(x_0))\ge.5\rho.
\]

校正器的建议步仍使用原MAP二次模型
\[
m_k(d)=\tfrac12\|r_k+J_kd\|^2+\tfrac12\|x_k+d\|^2.
\]
但每个proposal的最终选择是已访问、合法且有限半径/progress合格状态中的最小真实dataPhi。模型目标和选择目标不同，必须同时记录。原raw proposal始终保留；若没有更小dataPhi的合格访问点，不强制用修正点替代raw。最小dataPhi仍不在noise ball时，点保存但primary_qualified=false。

## 固定8维弱空间、59维允许family和更强随机对照

用03章相同的material-effective阈值、prior-aware广义谱，取最弱8维，构造Gaussian正交frame \(Q_w\in\mathbb R^{66\times8}\)，其补空间 \(Q_s\in\mathbb R^{66\times58}\)。初始方向由8个标准Gaussian系数作用于材料单位weak lift \(TW_8\)，B/C完全共用。

固定方法允许
\[
x-x_0=\alpha v_0+Q_s s,\qquad
(\alpha,s)\in\mathbb R^{1+58}.
\]
其GN允许空间由 \([Q_s,v_0/\|v_0\|]\) 张成，再限制到材料有效空间并施加原signed-progress tangent row；通常维数至多58。随后使用原R7 radial decoder root匹配有限材料半径。

这与REV1“只在Qs中做GN，再由radial root调整幅度”不同：新runner的GN本身包含所选weak line的幅度方向。也不能称为完整66维GN。Gaussian Qs与材料方向并非材料度量正交，所以signed material progress不自动保证native weak coefficient为正；冻结合同允许signed amplitude，历史名称中的RAY应理解为weak line，实际系数/符号另存。

动态方法每次有已付费当前J时刷新rank8 weak frame，取当前位移在该frame中的weak projection作为本步weak line；不声称跨frame平行运输。随机方法用独立nativeGaussian66 antithetic方向与原ordinary full-effective66 signed-row校正器，是更强几何class，不能人为削弱成59维后再声称公平。

## K16不是每半径16

每个arm/parent总共16个非anchor proposal：两个antithetic pairs乘四个半径 \(.02,.05,.10,.15\)，每个pair有+/-两个方向。顺序为pair0的全部半径+/-，再pair1的全部半径+/-。因此prefix K=1/4/8/16是确定且公开的。

每个半径独立从同一anchor出发，没有从较小半径的已选点warm-start。latent距离上限固定4，不自动放宽；返回RADIUS_BRACKET_NOT_FOUND_WITHIN_LATENT_CAP只表示给定射线搜索和cap未找到bracket，不证明非单调decoder距离下的全局不可行性。

共享anchor加16个slot至多17个点，可能重复、可能不合格、可能slot unavailable。PRIMARY_CLOUD.npz存放所有selected states及资格flags，并非“仅合格点文件”，更不是posterior samples。anchor也只有原decoder合法且在绝对noise ball内才进入qualified cloud。

## 真物理与成本分开写

固定host exact ROI elimination对所有原P1非均匀ROI材料有效：
\[
L_R=I_{256}-\operatorname{diag}(\delta)C,\quad
j_R=L_R^{-1}\operatorname{diag}(\delta)E_R^h,\quad
F=F_h+Oj_R.
\]
\[
E_R=E_R^h+Cj_R,\qquad
J_x=\mathcal W\{O L_R^{-1}\operatorname{diag}(E_R)D_x\chi_R\}.
\]

一次forward实际支付8个RHS256和一次256LU；full66 J用同一LU的32个reduced adjoint RHS，不是528个forward directional solves。只保留一个live LU。这个精确块消元是共同的普通physics backend，不是分配给fixed Q的独占科学优势。

共享SETUP缓存已经支付264个full576 RHS；此前完整stable ROI timing study共2184个full576 RHS。两者分别注明历史身份，不在新campaign再次收费。当前协议没有新full576审计或fallback。

上界是1536个proposal，各最多6次forward和4个J-bearing states：
\[
1536(6\cdot8+4\cdot32)=270336\ {\rm RHS256}.
\]
加32个cap64 anchor的81920 RHS256及生成256 RHS256，总352512 RHS256；最坏11296次LU256。硬cap为400000 RHS256、12000 LU256、7200 controlled seconds、8GB。不同维度RHS不相互冒充。拒绝、失败、cache恢复、源几何/SVD、decoder radius search与I/O都保留成本；sourcegeometry/SVD通常不增加RHS但增加CPU/wall。

## 数据流和seal边界

1. 独立generator使用20261030_R9_POSSIBILITY_VOLUME_V1 SeedSequence命名空间，prior与noise分角色。所有truth、generator预测及generator求解attempt arrays写private output；online文件只有observation/sigma/parent。
2. online anchor从zeros66开始，原R7 direct-probit TRF设置不变；online代码不导入generator，不打开truth。
3. 每个parent记录48行（3arm×16slot），包括不可用与失败slot；保存raw、每个实际visited trial、接受标志、最终selected数组、J是否真正计算、timing window/arm order和所有计费。
4. 每个parent的seal同时包含parents/<id>和jobs/<id>，因为VISITED.x_path指向后者；绝不加入private generator目录。
5. 只有32个parent全部完成并产生PREDICTIONS_SEALED.json后，下游统计/离线truth评分才能读相应数组。运行器不自动把候选云变成概率参考或posterior结论。

源冻结记录BLAKE2b与文件stat，包含实际engine/generator、原P1源码/target、共享SETUP缓存及root批准的statistics frozen config。没有SHA要求，也没有额外微型数值测试来替代最终32-parent语义验证。

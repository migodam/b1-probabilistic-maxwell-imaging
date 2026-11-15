# 下一轮Pro：优先读这里

## 本轮已经回答

用户提出按内部G与接收S的四种谱强弱重写G19，并周期融合DR。实际完成264主反演、56完整state干预、16材料读出反事实、44终点梯度验证及旧48端点分解。全部代码、字段、成本和反例留存。

**分母对齐观测有部分实质收益，但本轮DR融合没有净增加恢复。** primary CIE/Dual/G19/S19/U11/hop/routedDR/AP/full-cyclic/reset/standaloneDR分别9/9/11/12/12/12/12/12/8/12/13，均/24。原secondary球内最小R加回退分别10/9/10/14/13/12/12/11/7/13/14。不得按GT逐对象混选primary/secondary。

只涉及已曝光r14_new_0000–0023；不要与r13_fresh_XXXX或R11/R12同尾号混用。没有新的盲材料、矢量3D或实测证据。

## 最有用的公式

1. 原material Jacobian是 h -> [WS(I−DmG)^−1 D_Qs h]_s。背景Sh的核不是这个算子的核。旧G19误差在S-null的pooled比例55.45%，但经真实场/传播有显著可见作用；交叉项强烈抵消，不能把孤立端口能量强行归一成100%。
2. 新四空间按S条件下G*G分解；不是共同SVD。complex ranks全孔径11/0/9/556、半孔径6/0/13/557，U10为空，所以U11=P_Sstrong。此批不能证明G进一步筛选S强方向的价值。
3. 新变换T=β[A^-1+A^-1Dm(N−PNP)A^-1]。完整材料导数已推导/接线，只有一次余项插入；T*T未自动取得全box正定资格。物理G从未被它替换。
4. Dm与非坐标谱投影不交换，因此把S弱G强叫作另一个独立Maxwell材料集合不成立。路由的是复current/field增量，不是材料坐标。
5. 固定可行anchor及正交切投影假设下，路由DR切映射非扩张是标准后果。实际移动G19、非根曲率、PR和注入均不受此定理覆盖。传播G的谱上界不约束它们的闭环谱。

## 三个最关键的反例

- **0020对0021：** 路由DR相对hop救回0020（NRMSE .4565→.0667），却丢0021（.0489→.4723）。整体12不变是真实交换，不是“都一样”。全空间周期DR更差；独立DR仍13。
- **即时下降不预测完整恢复：** 同0021 checkpoint512，G19单次原Phi下降.01811，S19上升.003416，hop上升.000616；完整程序却是后两者恢复。0020的G19有限步线性预测+.03609，实际−.05852。不要再用局部梯度的一个符号给未来搜索设接受门。
- **corrector与选择：** 0013原G19 raw NRMSE .0915/IoU .962，primary .1856/.832；secondary .0912/.965。D下降而R上升，不能都怪过强平滑。0019许多方法已恢复两岛主体，IoU .67与.72不代表截然不同的材料结局。

## 真实动力学材料

早期4cycle：路由DR材料变化在real U01之外的中位能量81.985%；current注入相对L2误差中位23.896%，AP17.346%、full32.796%。不是信息损失率，也不是证明每个后期周期都如此。

四对象共享checkpoint512，五臂64step；0021/23三臂三band±扰动16step，共56。全部stage/beta不混淆。16步材料增益3.7–50.5且ordinary也出现；不是DR特有，不是局部谱/chaos。PR-only在第16步事件后才影响第17步材料，解释其16步与hop相同。full初末PR状态保留，窄STATE016不能冒称包含完整PR历史。

0021的G19坏primary lower/upper为515/31，route518/29，S19正确为96/0；弱投影梯度不证明全Hessian局部极小，却提示活动面的可访问性比原G谱更接近实际瓶颈。

## 文件入口

- RESEARCH_REPORT_ZH、THEORY_DERIVATIONS、TRAJECTORY_ANALYSIS、CONFIG_PRECEDENCE。
- analysis/twofold_exposed24_v1：ALL_VIEWS、RECOVERY_GAIN_LOSS、STRATIFIED_SUMMARY、公共预算前缀及bootstrap；全部选择先于GT。
- analysis/fullstate_v1：480cycle pairs、36finite responses、same-state比较、真实counter。
- audit/readout_counterfactual_v1：4×4同状态材料步、grad方向、完整forward变化。
- audit/paid_endpoint_diagnostic_v1：4×11 primary的physical gradient/active box/full fields。
- analysis/intervention_scoring_v1：4对象中间state和材料步的封存后GT诊断，不是新完整反演。
- audit/old_endpoint_geometry_v1：48旧G19/DR材料/current误差、chord/local动作与有符号交叉能量。
- analysis/selection_geometry_v1：24粗GT物理参考、全部已选视图与球内候选对。
- release full raw：全部264/56 current/field、PR、候选和checkpoint。KEY_FULLSTATES选集只提供快捷入口，不冒充全量。

原生DR更早进入corrector，中位首完成审计43.25CPU秒，stateful约73秒。最高恢复数属于整个付费流程，不宜只归功于反射。共同LU/RHS池中G19/standaloneDR为12，S/TF/hop/route为11；等CPU与等物理预算不能混称。

## 下一轮值得推导，而不是现在已证明

1. 用真实材料transfer/Jacobian而不是裸S/G定义可观测商，研究source场、resolvent与材料normal cone怎样耦合。
2. 给出能解释好/坏活动面的材料有限动作判别。不能读取GT，也不能把即时Phi下降当充分条件；当前反例已排除这些简单故事。
3. 研究表示注入能否保护material/normal memory，区分joint P_C接口与路由增量；本轮还没有H=0的严格完整24对象控制，不要把所有作用都归给所选H。
4. 噪声球内保留多个材料与搜索能力分开；原Phi低于GT参考的错误材料不是靠更精确优化同Phi自动修好。

本轮没有新算法循环等待执行。任何后续主张须面对0020/0021反向结果、0013/0019选择反例，以及独立DR的效率。SOM/twofold/DR本身有先例；新颖性仍需实质理论或新的独立物理恢复证据支撑。

# Full-J 跨会话源码审查：基线充分性重新打开

最终状态：`TARGETED_TRF_DIAGNOSTICS_COMPLETE / MATURE_BASELINE_DOMINANCE_NOT_ESTABLISHED`。第1–5节为质疑时审查，第6节为科学冻结后的更新。本次只读源码、公开配置及已冻结结果事件；不运行新物理，不访问private GT、种子或未揭晓材料。审查者为Sol/high，没有Astra或子代理。原主与补充预测保持原样。

**当前判断：不能用“本研究对比度更高”替代对普通求解器实现和参数的核查。** 本轮Full-J的解析导数、共同数据和参数空间有审查依据，但自编投影LM的全局化、cap和停机不等于成熟有界最小二乘；因此“Gram在充分强普通基线上有独立增量”的充分性仍待审查。数字9在旧冻结方法间的观察成立，但不能升级为最优普通方法不可达。

## 1. 真实文件与容易混淆的命名

本研究：`src/solver.py::Run.gn`（122–189行）、`src/physics.py::Model.jacobian`、`configs/main_methods.json`，根为本文件上两级。

Route B实际R11根：

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-p3l-r11-stable-digits/research/p3l_r11_stable_structure_digits/`

读取 `shared/invert_digits.py`、`shared/core.py`、`INVERSE_CONFIG.json`、`inputs/PRIOR_CONFIG.json`、`data/GENERATION_CONFIG.json`及公开 `METHOD_COMPARISON.csv`。**R11实际逆解是自编投影GN，不是SciPy `least_squares`或TRF。** `invert_digits.py` 导入的是 `scipy.linalg.solve`，用于其阻尼法方程。不能把跨会话“Route B”标签自动等同于成熟TRF。若另有一个明确使用 `least_squares` 的Route B版本，必须另给源码路径核对，不能以当前R11替代。

Route A当前FWI根：

`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/Research_W/publish-w-p0/research/research_w_gsda_r4/`

读取 `src/gsda/fwi.py`、`src/gsda/fem_model.py`、`configs/main.json`、`configs/base.json`、`METHOD_FREEZE.json`。实际是有尺度变换的有界L-BFGS-B、两公开起点、三频延续；不是Full-J/TRF。它的当前有限结果不能未经审阅就引用为全数字已成功的反例。

## 2. 优化器的逐项比较

|项目|本轮R3-DIGITS|Route B R11真实实现|Route A R4 FWI真实实现|
|---|---|---|---|
|优化变量|576独立实材料；8²/12²/24²连续释放|Native N为256 occupancy＋2全局材料属性=258；P为64latent＋2属性=66|49target contrast＋49target loss＋6wall epsilon＋6wall loss=110|
|主优化算法|自编投影GN/LM，GN方向后逐分量clip到box|自编投影GN，方向norm受trust cap，occupancy clip；两属性经Gaussian CDF进入区间|SciPy L-BFGS-B，q=z/scale后box；解析伴随gradient，无显式Full-J|
|方向缩放|Hessian对角damping，floor为maxdiag×1e-8或1e-12；无独立x_scale|Hessian对角floor1，初始damping .01；trust=3，最大6|wall/target epsilon、loss分别用显式scale；target contrast2、loss.1等|
|全局化|8次阻尼试探，投影后只要原目标下降即接受；无alpha回溯或反射/active-set子问题|同一方向alpha=[1,.5,.25,.125]；方向norm受trust限制；拒绝后damping×10、trust×.5并继续外迭代|L-BFGS-B库line-search maxls8；有界准Newton；maxfun/maxiter按阶段计数|
|拒绝/停机|8次未接受即`GN_NO_ACCEPTED_STEP`且break；过小step继续增damping可能令后续更小；无gradient备用步/TRFfallback|拒绝累计patience；相对目标小或拒绝累计4次才停；cap为80Jac/240forward/120秒|库停止条件ftol1e-9、gtol1e-4；calls[18,25,110]、iterations[12,18,90]，每起点三阶段；保留未收敛状态|
|初值/延续|共同常数1；所有臂共享native8/12各320iter后fine分流；原fine cap1000|occupancy .2；属性Gaussianized0，即delta1.3+.105i；无本轮多尺度或双起点|两个公开uniform初值；频率1→2→3；main spatial_shapes=null，即7×7从始至终|
|数据接受/选择|共同全观测相对半平方误差＋平滑TV/二次项|共同噪声白化半平方误差＋occupancy TV/L1/ridge和属性prior|共同full3freq白化误差＋TV/L2/wall Gaussian/spatial prior；低频阶段按nf/nactive缩放|

这些差异不意味着R11一定优于本轮，也不能从一条line-search技巧直接断言本轮失败根因。它们足以说明“都有解析Full-J”不等于“优化器相当”。

## 3. 本轮停机的实际证据，不把潜在缺陷冒充已定位原因

直接读取原27份冻结事件：native的1、2、3、5、6、7、8、9均完成fine1000迭代到cap，只有4在投影梯度门停止；这些原native没有`GN_NO_ACCEPTED_STEP`。reciprocal为6例到1000cap、1/9到wall、4投影梯度停止。Gram原native-polish有3和8到1000cap、4梯度停、其余wall。

所以“8次阻尼失败导致原native六个图像失败”**没有被这些原运行事件支持**。更具体的已知不足是：多数原方法在远离准确材料时到达1000cap，且最短45.95秒就提前结束90秒额度；后加budget-fill在剩余时间内继续，但没有替换投影GN的全局化策略。普通data-adjoint补充另有六例fine GN没有接受步；这才是8次有限试探/无备用步可能直接影响的运行集合，仍需更强求解器对照来辨别局部最小与实现停滞。

逐坐标clip无约束GN方向不等于解box约束GN子问题；它不继承无约束方向的下降保证。有限阻尼次数后停机，不构成满足KKT或普通最小二乘已充分收敛的证据。本轮虽记录投影梯度门，`COMPLETED_NUMERICAL`可表示到迭代cap，不能读作收敛。阶段热启动从已评分最优合法候选出发是合理保护，但不补足全局化。

## 4. 物理、自由度与先验必须量化，但不是免责理由

|合同|本轮|R11 Native N|Route A R4|
|---|---|---|---|
|物理/观测|24² scalar-TM disk-pulse VIE，单频k=2pi，12源×48接收=1152实观测|同24²类TM网格，ROI精确宿主消元，8源×32接收=512实观测|有限元散射/穿墙模型，4源×16接收×3频=.85/1/1.15，共384实观测|
|未知材料支持|576全网格均未知，box[.5,20]|ROI256 occupancy未知，320外部单元固定已知宿主；另2个全局增量属性|49target单元复材料＋6wall单元复材料=110实参数；已知目标/wall参数分区|
|生成材料公开范围|背景.8–1.5、笔画8–18，实数，无loss|HOST=.1+.04i；公开固定生成增量1.3+.105i；满occupancy foreground1.4+.145i|target contrast box[0,8]，loss[0,.25]；wall epsilon[1,6]、loss[0,.25]；实际每例GT不在本静态审查访问|
|材料族结构|各cell实值独立＋generic光滑先验，无已知外部固定区|各cell occupancy[0,1]、所有cell共享一个复增量；属性区间real[.8,1.8]/imag[.03,.18]|相同cell的real/loss独立，但使用明确ROI/wall分区和参数尺度/先验|

R11已知外部320/576=55.6%，Native未知258与本轮576相差2.23倍；本轮笔画real8–18约为R11满occupancy real1.4的5.71–12.86倍。以笔画–背景差6.5–17.2对R11增量模长约1.304计，振幅比约4.98–13.19。它们是材料值/生成合同的比较，**不能替代实际Born非线性、resolvent条件或优化难度的测量**，也不能证明增加loss或减少对比度必然足以解释成败。

R11 Native TV_Huber_lambda=.4、epsilon=.05、ridge=.02、L1=.05，作用于occupancy的边差之和，并有2属性Gaussian prior；本轮lambda=1e-4作用于材料差/19.5的边均值，无L1与固定外区；Route A target isotropicTV=.2、scaledL2=.01且wall Gaussian/spatial prior。数据白化也不同，不能用原始lambda比值声称“正则强4000倍”或直接归一化成绩。R11 full-grid误差和SSIM还包含已知宿主外区，须同时看ROI误差；本轮全图均未知，原始RMSE/IoU不能当同难度排行榜。

R11公开表记录Native Full256平均IoU .6382、full-grid complex MSE .031714，8/9通过其noise-ball，但该研究自身`EXPERIMENT_RESULTS`称资源protocol FAIL，不能选择性引用成功而抹去资源资格。无论这些差异有多大，本轮仍有义务在自身相同Y、box、prior与时间下检查成熟有界基线。

## 5. 当前应采取的科学状态

1. 原冻结原始结果全部保留，数字9“旧Gram方法优于旧native/reciprocal/ordinary组合”的观察不删除；升格为“充分强普通基线后仍成立”的结论改为 `UNDER_REVIEW`。
2. 父线程已授权两公开失败状态上的独立SciPy TRF诊断，源码/配置后GT隔离冻结。它可以检验该基线是否被改进，但属于定向事后诊断，不是补回原预注册盲成功率，也不能将两例结果泛化为九例。
3. 若TRF显著改善普通解，本轮方法比较须相应降级，不能把收益继续全归Gram；若TRF两例未改善，也不能据此宣称普通求解器已经最优。需要实际停机、预算、原目标与材料指标，不能只看库success或高对比度口头解释。
4. 本次审查新增物理RHS=0、LU=0；对RouteA/B只读源码/公开配置和既有汇总，没有碰private GT/seed。最终核心方向由父线程裁决。


## 6. 科学冻结后的TRF、幅度控制与最终边界

截至 `SCIENCE_FREEZE.json`（20:41:43 UTC），不再允许物理运行。本审查只读最终收据：两例重走全路径的SciPy TRF残差0.8901475和1.0582021，没有恢复高对比8/9形。对9形另以原native/Gram完全相同的公开fine初值做TRF，新fine62.86499秒＋原coarse参照27.36016秒＝90.22515秒，残差0.3909565；自编补齐native为0.4217865，Gram为0.000228314。新TRF确有适度拟合改善，仍失败；不能据两定向事后case宣布成熟普通基线最佳或Gram普遍支配。coarse参照不是新计算，成本不可重复累加。

坏端点单方向Jv相对误差1.3570e-9（3LU/84RHS），三TRF端点独立总场验证再3LU/36显式RHS并一致，分别见 `FULLJ_DIAGNOSTIC_ENDPOINT_DIRECTION.json`、`FULLJ_DIAGNOSTIC_INDEPENDENT_VALIDATION.json`。这些检查反对该范围内明显导数/物理不一致的解释，但不证明所有参数/全域实现或优化充分性。

独立生成的同8/9形态低幅度控制 `m=1+1.3alpha` 保持算子、box、初值、正则和native算法；45秒cap实际早停。8/9的残差2.12502e-4/1.88303e-4，RMSE0.231967/0.231460，支持同一实现对这些低幅度图像可工作。fine576各1J、0接受更新，粗阶段已拟合，不能虚称fine充分迭代。它改变物理材料/Y，是主GT后灵敏度诊断，不是原高对比失败被修复，也不等于复现RouteB全部ROI/loss/prior。其新support阈值1.65不得混入主阈值5成绩。

结论收口：用户对基线质量的质疑合理，已补源码对照、成熟TRF及匹配fine初值诊断；没有发现足以解释全部失败的简单Jacobian bug，也没有得到高对比普通法或全九图像恢复。原数字9方法间观察保留，成熟基线支配性和一般算法归因不成立，全部强主张仍HOLD。噪声14/27部分完成，13NOT_RUN，无偏完整噪声比较未形成，禁止汇总部分均值。原方法、事后诊断、弱幅度控制与取消批次都应完整公开。

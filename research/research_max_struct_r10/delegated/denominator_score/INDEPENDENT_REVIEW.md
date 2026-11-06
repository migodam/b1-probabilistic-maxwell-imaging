# 低秩分母实现与理论的独立复核

**源码层面未发现分母材料导数遗漏、Woodbury伴随次序错误、Galerkin误实现为GP或错误沿用四颜色独立更新。** 本审查只阅读实际源码和理论，未重跑机器恒等式或新增反演；最终数学资格与新颖性判断归root。六条历史DEV评分已独立完成，不能将其称为新确认或blind结果。

## 已核对的数学对应

| 项目 | 实现位置 | 复核结论与条件 |
|---|---|---|
| 实际K | `src/lowrank_denominator.py:38–50` | Galerkin分支先取GU，再取U(U*GU)，故BU*=PGP；right分支才是GP。U由公开G/S固定右子空间构造，材料和GT不参与选基。 |
| Woodbury正解 | `:52–69` | A=D_a(I−FU*)，F=D_(m/a)B；代码先计算b/a，再加FT⁻¹U*(b/a)，与理论一致。没有交换D_m与P。 |
| 伴随 | `:65–67` | A⁻*=D_a⁻¹(I+UT⁻*F*)，代码先F*rhs、T伴随求解，再右式整体除a；a在当前实m、实β合同下为实数，次序正确。不能把此代码无条件推广到复m/复β。 |
| D_m导数 | `:71–79`；理论式(3) | J、Q、β、U、Ω在一次材料读出内固定。de=−D_hQ，dA=D_h(βI−K)，故dw=−A⁻¹D_h[Q+(βI−K)w]。代码的`h=Q+beta*w-self.action(w)`保留了分母导数。不存在把A冻结求错目标梯度的问题。 |
| 源metric与先验 | `:74–79` | 在Ω为固定Hermitian源metric的既有合同下，`w@Os`和A伴随解给出所写实梯度；`gprior`另加，未漏先验导数。若将来Ω随m变化，则本公式不含其导数，须另立合同。 |
| 非局部材料优化 | `:86–110` | 使用完整576维有界L-BFGS-B，最多8次迭代；普通CIE种子、旧点和优化结果按同一完整候选辅助目标比较。未调用原A+的四颜色条件独立更新。继承的OneHop对象仍可用于旧诊断，不能据此误判候选更新仍是FourColor。 |
| 分母近似差 | `THEORY_DERIVATIONS.md:76–82` | `R_K−R_G=−A_K⁻¹D_m(G−K)R_G`符号与次序正确，要求两分母均可逆；并未声称仅靠小谱尾即可控制GT误差。 |
| 非凸反例 | 理论`:88–94` | 所列一维二阶导数与结论相符，且文中明确不将其当具体Maxwell成像反例。 |

## 零集、优化和适用范围

理论的零集结论应始终读为：**β≠0且A_K(m)可逆时**，`βA_K⁻¹e=0 iff e=0`。当前冻结β为正、m为实合法材料，a=1+βm为正；但这本身没有证明小矩阵T在整个材料盒内处处可逆。`factor()`目前直接LU，没有给出全盒T的逆界/条件界；一次有限运行成功不能替代这一资格。理论原文已经保留“只要A_K可逆”的条件，因此这是适用范围，未发现其被偷换为无条件定理。若未来真的遇到奇异T，SciPy的LU warning/非有限值也需作为失败事实保留，不能靠末端回退推导一个全局证书。

`accepted`只要求优化返回有限值且不大于旧点/普通种子中较小的候选损失；即使L-BFGS-B未达到收敛，也可以返回一次已改善读出。`material_inner_success`保留了优化器自己的success，故应分别报告“接受了改善”与“优化收敛”。这种接受不保证原Maxwell Φ下降，更不保证GT改善；物理选优和corrector合同仍不可省。

这里改变的是固定J下的材料metric和完整材料读出规则，没有改变完整物理Green、FFT/current算法或材料自由度。静态投影降到19/28维的是K的秩，不是将材料反演限制为19/28个未知数，也不增加物理信息。

## history标签已由实际adapter补齐

历史DEV入口`experiments/denominator_run.py`在`DENOMINATOR_SETUP.json`写出旧`state_term_material`仅为legacy onehop诊断。fresh入口的`_prepare()`自身未重复此字段，但已核对实际`src/twod_adapter.py:305–312`：custom分支为每个iteration写入

`state_term_material_semantics=LEGACY_ONEHOP_DIAGNOSTIC_NOT_CUSTOM_MATERIAL_OBJECTIVE`

并将`custom_material_objective_field`指定为`material_final_block_loss`，execution receipt也说明相同语义。因此现有adapter已经补标签，不构成当前fresh运行遗漏；本审查未热改root代码。后续汇总应读取`material_final_block_loss`，不能将继承的旧onehop `state_term_material`画成新分母目标。

## 成本和比较口径

每次`value_gradient()`分解一次r×r的T；正解与伴随各一个ns列小系统，计数为一小LU、2ns个小三角RHS。例：历史b20d Twofold记录33772个小LU与810528个小RHS，恰为33772×24；它们与1491个原物理576维LU、35736个物理复RHS分列，不能互换或合并成“同一LU”。

`prepare_bases()`确实为每条候选都构造G/S谱与投影正交化，连G19也保留这笔setup；setup有CPU与SVD计数。在线低秩乘法统计不是所有BLAS乘法/字节的逐项完整账，但进程CPU/wall包含实际操作，不应拿局部计数直接推出复杂度实测。auxiliary execution、外层auxiliary与总进程计时互相嵌套，不能相加重复收费。

六条新结果与本表全部六条R9 cie_dual/global/A+参考均记录Darwin arm64、host=`migodamdeMac-mini.local`、Python3.13.13，因此可以作**同记录主机的历史描述性cost比较**。它们不是本轮同时重跑，也没有控制两次执行的后台负载、硬件频率等条件，不能升级为严格同负载速度优势。两条R6 CIE-single参考未记录完整execution_platform，保持unknown，不将其混入上述同主机分组。fresh同队列、同物理合同的比较仍是本轮主要成本证据。R6 CIE-single是成熟旧实现参考，R9同code fixed .001的历史DEV无CIE-single；此缺口不能改名掩盖。

## 历史结果对主张的限制

六条候选在b20d与40ac最终均通过原材料门，但G19、G28和Twofold28接近，**没有证据把收益称为Twofold接收信息特有增益**。b20d原R9 A+本已通过，新的NRMSE与IoU之间也有取舍。40ac的确改善了GT材料，并非仅残差变小，但Twofold/G19在corrector前NRMSE仍略高于0.15，最终通过包含继承物理corrector贡献。当前合理表述是公开历史DEV上的材料metric候选结果；fresh结果必须独立保留PASS/FAIL，不按DEV或freshGT再选rank/阈值/臂。

# 固定核费用修订：rowG2代数与共同stored I/O

新增src/p3r4/decision_refresh_fast.py及tools/run_r4_interval_chain_fast.py。本修订只优化确定性代数、缓存和证据持久化实现，不提出新定理、不更换原统计target/核；原两对已封账链和全部结果文件不改。尚未启动任何付费物理。

以Tight界为原样基础：预存row_g2[i]=∑j|Gij|²，利用严格代数恒等式||diag(dχ)G||F=sqrt(dot(|dχ|²,row_g2))。每次只处理576行标量，避免临时576×576扰动矩阵。数学界、左逆缺陷项、两端势区间决策、beta和max4Richardson均不变；float64不同归约次序有约1e−16差异，不把其称作定向区间证书。

Fast对象拥有readonly的χ0/j0/E0/G/S/incident/R0副本，2项LRU cache按χ dtype及完整字节标识，cached constants为不可改MappingProxy。每次归约、key/向量copy bytes、命中、矩阵动作、setup copy/row预计算及wall均记账。R0仍须真正通过付费solve_anchor(I576)获得；R0矩阵应用只记MVP和wall/内存，不冒充新的求解器RHS。

经父调度器明确批准，core.save增加进程局部R4_ARRAY_CODEC=stored选项，对应np.savez lossless ZIP_STORED；缺省保持np.savez_compressed，未知值报错。仅新fast pair子进程显式启用stored，其他进程与旧实验默认不变。Book SOURCE_REVISION/FROZEN_ARM/COMMAND记录codec，stored array写次数/字节/时间有进程统计。没有改Target求解或后端数学。

固定parents R4DEV_00040…00047，每parent1024步，full66 pCN beta=.002，seed20261025131。Full与interval两个独立子进程使用相同初值与SeedSequence([seed,parent_index])，每步normal66后uniform1；运行中不读另一arm。两arm共同采用stored模式及下述相同物理端点数据合同，Full是优化后的实际执行对照，不从旧压缩实验推断收益。

每个提案在求解前立即保存x/q/χ/logU尝试，包括后来失败的提案。每个物理材料端点的最终current/field/prediction/nativeq/χ/backward-residual必须立即可达。Full直接引用Target已经写好的native endpoint，runner不再重复保存相同大数组；精确Full状态也不重复拷贝大field arrays。interval仅保存自己的最终近似端点大数组，公共数组字段与单位一致。每步拒绝/重复状态保存endpoint引用、完全相同的float64 x及判决/区间；最终CHAIN保留完整x路径与全部accept标志，不把拒绝丢掉。

同一物理端点的中间Richardson精化只立即追加残差范数/白化误差/Φ包络等proof scalars，不再写每次中间大数组。原R0/base及背景身份保存，可由原χ与累计refinement次数确定性重建。当一个旧accepted approximate CURRENT在之后步骤被进一步精化时，允许原位置更新**这个相同材料**的最新近似数组；ENDPOINT_VERSIONS逐次追加canonical array版本和accuracy iteration，旧标量proof不删除。所有overwrite次数、写入bytes及时间公开。这是一个canonical文件而非“永远只写一次”的虚假断言；paidFull文件及任何closed job绝不覆盖。

如果旧近似端点后来真正付费Full，canonical索引转到新的不可改paid endpoint；旧近似文件可保留为superseded证据，所有数组身份一并封存。新物理端点从不因状态重复而漏记，未完成的物理尝试保留attempt与错误；失效时停止而非筛样本。

硬上限仍为pair200000复RHS：Full arm65536，interval arm134464(包括8×576冷R0和所有真实fallback)，每arm1200秒与3GiB。BookII/physics.lock管理。没有GT、概率reference或posterior资格声明。source input identities与所有canonical/latest/attempt/CHAIN/reconstruction数组身份在关闭后核对，不做SHA。

父调度器后续release才允许运行：

```json
{"released":true,"experiment":"interval_chain_fast_v1","campaign_id":"P3LR4-20261025"}
```

入口`<python> tools/run_r4_interval_chain_fast.py --release RELEASE_PATH`只协调三个独立子进程(interval/full/compare)，并为两个arm显式设置stored。显式单arm调试要求调用者也设置R4_ARRAY_CODEC=stored；正式全费用判决需要pair协调器记录ARM_PROCESS_WALL.json。外部进程计时覆盖cold import、R0/setup/copy、全部大数组与标量proof/账本、最终gate写盘及退出，不用“计算小计”替代整个作业。

两个arm Book都CLOSED COMPLETE、PID不同、输入一致及数组身份未变后，事后检查全部accept flags、x路径、最终当前/候选势区间，以及**每一个中间proof**对相同χ的独立Full Φ包含关系。只有0flag/path/enclosure miss且interval外部整进程wall相对共同优化Full节约≥20%，才PASS_FIXED_KERNEL_GE20_PERCENT_WALL_SAVING；任一不满足即NO_GO，不扩时长、阈值或换数据救结果。该PASS仍只是固定核实现/费用结论，不说明混合、MCSE或FULL后验合格。

FAST_EQUIVALENCE_CHECK.json记录0新RHS/状态/训练的纯检查：30个norm case最大扰动范数差5.8e−17、η差2.9e−16；首4个已付tight候选的保存数组重放包络一致且包含其独立Full；2项cache及readonly检查、默认压缩/无损stored/非法codec拒绝、同端点版本替换与Full不重复大数组均通过。只读旧证据，不重跑旧物理、不以这个测试宣称20%加速。

import json,sys
from pathlib import Path
OUT=Path(__file__).resolve().parent
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=json.loads((OUT/'RESULTS.json').read_text());cost=json.loads((OUT/'OPTIMIZED_COST_BENCHMARK.json').read_text());ledger=json.loads((OUT/'COST_LEDGER.json').read_text());rows={r['family']:r for r in R['summary']};cm={(r['arm'],r['batch']):r for r in cost['rows']};gates={r['family']:r for r in R['final_cost_gates']};selected=['strong_CNN','original_K32_frozen','deterministic_frozen','deterministic_matched','direct_moments_matched','distilled_w64_n4096','distilled_w64_n32768','distilled_w128_n32768'];names=['Strong CNN','K32 frozen','Det. frozen head','Det. matched head','Direct GT moments','Distill64 / 4096','Distill64 / 32768','Distill128 / 32768']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'figure.dpi':150})
fig,axs=plt.subplots(1,3,figsize=(15,5.2),gridspec_kw={'width_ratios':[1.25,1,1]});ax=axs[0]
for i,(family,name) in enumerate(zip(selected,names)):
 vals=np.array(rows[family]['seed_MSE'])*1e3;ax.scatter(vals,[i]*3,s=17,color='#6b7280');ax.plot(rows[family]['seed_mean_parent_MSE']*1e3,i,'s',color='#0d7490',markersize=6)
ax.set_yticks(range(len(names)),names);ax.invert_yaxis();ax.set_xlabel('Full-grid complex MSE × 1000');ax.set_title('Development VAL256, 3 seed losses');ax.grid(axis='x',alpha=.2)
fams=['distilled_w64_n4096','distilled_w64_n32768','distilled_w128_n32768'];vals=[gates[f]['VAL_retained_fraction']*100 for f in fams];axs[1].bar(range(3),vals,color=['#86b6c6','#418aab','#145b7d']);axs[1].axhline(80,color='#b34439',ls='--',label='80% frozen gate');axs[1].set_xticks(range(3),['64 / 4096','64 / 32768','128 / 32768']);axs[1].set_xlabel('Width / teacher TRAIN count');axs[1].set_ylabel('Retained incremental risk gain (%)');axs[1].set_title('Capacity + data ablation');axs[1].legend(frameon=False)
for i,v in enumerate(vals):axs[1].text(i,v+1,f'{v:.1f}%',ha='center')
labels=['Optimized K32','Det. Q16','Distill64','Distill128'];keys=['optimized_original_K32','deterministic_q16','distilled_w64_n32768','distilled_w128_n32768'];x=np.arange(4);provider=[cm[(k,64)]['provider_warm_per_parent_seconds']*1e6 for k in keys];full=[cm[(k,64)]['full_warm_per_parent_seconds']*1e6 for k in keys];axs[2].bar(x-.17,provider,width=.34,label='Provider only',color='#418aab');axs[2].bar(x+.17,full,width=.34,label='Provider + features + head',color='#b39452');axs[2].set_yscale('log');axs[2].set_xticks(x,labels,rotation=24,ha='right');axs[2].set_ylabel('CPU µs / parent, batch64 (log)');axs[2].set_title('Common paid anchor acquisition excluded');axs[2].legend(frameon=False,fontsize=8)
fig.suptitle('Track B: useful moment features compress; locked confirmation remains resource-limited',fontsize=13);fig.tight_layout();fig.savefig(OUT/'TRACK_B_SUMMARY.png',bbox_inches='tight');fig.savefig(OUT/'TRACK_B_SUMMARY.svg',bbox_inches='tight');plt.close(fig)
phase='\n'.join(f"| {p['region']} | {100*p['ROI_cell_fraction']:.2f}% | {p['risks']['strong_CNN']:.7f} | {p['risks']['original_K32_frozen']:.7f} | {p['risks']['deterministic_matched']:.7f} | {p['risks']['distilled_w128_n32768']:.7f} |" for p in R['phase'])
table='\n'.join(f"| {name} | {rows[family]['seed_mean_parent_MSE']:.8f} | {rows[family]['ensemble_image_MSE']:.8f} |" for family,name in zip(selected,names))
costtable='\n'.join(f"| {label} | {cm[(k,64)]['provider_warm_per_parent_seconds']*1e6:.2f} | {cm[(k,64)]['full_warm_per_parent_seconds']*1e6:.2f} |" for k,label in zip(keys,labels))
nested=[]
for arm in ['raw_y','decoder_only','physical_residual','same_spectrum_random']:
 rr=[r for r in R['nested_RRR']['records'] if r['arm']==arm];nested.append(f"| {arm} | {rows['nested_RRR_'+arm]['seed_mean_parent_MSE']:.8f} | {np.mean([r['selected']['nested_baseline_OOF_MSE'] for r in rr]):.8f} | {np.mean([r['selected']['nested_corrected_OOF_MSE'] for r in rr]):.8f} |")
report=f'''# Track B：材料矩与风险表示，已完成 development，独立确认受预算限制

本轮没有增加物理调用。最有信息量的结果是：Gaussian 模型输出的四个材料矩可以用一个 **214,272 参数**、仅消费原始测量的模型蒸馏；其三 seed 验证风险保留原 Gaussian fusion 相对强 CNN 新增收益的 **{100*gates['distilled_w128_n32768']['VAL_retained_fraction']:.2f}%**，额外 provider 的 batch64 推断费用约为优化 K32 的 **{100*gates['distilled_w128_n32768']['provider_cost_ratios']['64']:.2f}%**。这是已有 VAL256 上的压缩候选，**不是新的独立 TEST 确认，也不是 posterior 概率合格证据**。完整 provider + 特征整理 + risk head 的费用仍约为 K32 的 **{100*gates['distilled_w128_n32768']['fullpipeline_cost_ratios']['64']:.2f}%**；共同物理 anchor 成本尚未删除。

## 1. 恢复资产与实际定义

只读恢复 R4 真正的 Gaussian checkpoint、三份原 fusion head、三份 strong CNN、TRAIN32768 / VAL2048、付费的 TRAIN1024 / VAL256 强 anchor，以及原 C2 decoder。原 x 的变换是 `x = (x_normalized * x_scale + x_mean) @ frame_inverse`。实际 checkpoint 恢复后的场–属性 covariance 为零，因此生产路径可以将属性与场的期望相乘；这一结论来自实际 frame/covariance 核验，不能从网络坐标 diagonal 字样推断。

Gaussian encoder 保留 **682,943** 参数。原测量为 source-major Re32/Im32，噪声 sigma = 0.00119211798514519，白化为 sqrt(2)/sigma。所有训练与推断均保持原模型、温度 0.5、C2 quintic 与两个 probit 属性区间。没有替换 sigmoid、重标噪声或引入新 prior。

`material_moments.py` 实现任意联合 Gaussian 下的一维场积分：先在给定场值时计算属性 probit 的一阶、二阶和实虚交叉 CDF，再对 C2 饱和函数及平方积分。保留场–属性 covariance 和属性之间 covariance；退化场方差取极限。`gaussian_image` 返回完整四通道，可选逐像素 Re–Im covariance。它没有声称提供完整跨像素 posterior sampler。

原 R4 K32 使用 empirical variance `ddof=0`，其期望是 population variance 的 31/32；确定性路径计算 population variance，是大 K 极限，**没有乘 31/32**。因此冻结 head 的替换是明确的输入分布改变诊断；匹配重训才测试改变后的表示链。

## 2. 数值资格与失败留痕

- 相关/强相关/退化/尾部/平台边缘案例：确定性材料均值及边缘方差对独立 adaptive quadrature 的绝对误差小于 4e-15；每例 1,000,000 MC 的均值误差小于 1.54 个 MC 标准误，联合 covariance 也检验。
- 真正恢复的 Gaussian provider：8 个 VAL parents 各 131,072 MC 作为积分参考，全部使用 4096 分块；材料四通道 RMSE 约 2e-4 至 7e-4，符合有限 MC 误差量级。
- 实际 VAL256 上 Q16、Q32、Q64 均与 Q128 的 float32 四通道结果一致。生产训练使用稳定 Q64；Q16 是本 checkpoint/该验证范围的数值快路径，不能据此宣称任意极窄 Gaussian 都可固定低阶积分。
- 宽方差截断多项式 recurrence 曾发生 2.23e-5 消减误差，已改为条件不良时稳定积分；高相关属性 CDF 固定节点曾误差 1.03e-5，已用角度变换并在 |rho| > .95 时 adaptive fallback。初次失败文件全部保留。
- **一次初始 high-K 全像素 materialization 的 RSS 达到 5.77 GB，违反 3 GiB 限制**，已作为失败保留；随后改为 streaming，重跑 peak 1.50 GB。此后训练主任务观察到的峰值低于 3 GiB，但没有声称证明连续 host + MPS 联合分配总量的硬上界。

## 3. 四个必要比较与直接压缩修订

风险是每 parent 的 full-grid complex MSE，再对三个训练 seed 的 loss 取平均。ensemble 图像风险另列，不混用。新 head seed 固定为 2026102601 / 02 / 03；冻结历史 head seed 为 2026102501 / 02 / 03。

| 表示/处理 | 三 seed 平均 loss | 三 seed ensemble 图像 loss |
|---|---:|---:|
{table}

确定性矩 + 原 head 从 0.00714550 降至 0.00704106；匹配重训为 0.00703415。这说明此处去除 K32 扰动有验证收益，**不证明任意非线性 head 去 MC 噪声必然改善**。独立训练 seed 与历史 head seed 不相同，因而纯输入替换与新训练增量应分开解释。

原监督材料矩 provider 实际已经很小：**120,192 参数**，并非 0.68M。它只保留 43.61% 的 fusion 增量，因此保留这一负面证据。随后在同一第四表示臂内预声明 objective/capacity 修订：用 Gaussian 确定性四通道作 teacher，raw-y → 64 → 64 → 4ROI，**103,552 参数**。teacher 只读 joint TRAIN 的 y,parent，不 materialize GT。4096 teacher 数据仅保留 44.47%；按冻结的严格正 VAL 均值差 gate，扩大既有 TRAIN 到 32768 后保留 67.60%。最后一组预声明 width128/128（实际 214,272 参数）达到 92.96%。各 provider 与匹配 head 均固定 4000 updates、三 seed；没有继续搜索容量或更新步数。width128 协议的参数估算 214144 漏计128个 bias，实际架构未变，实际值在 checkpoint/结果中为准。

这里“删除 encoder”指在线蒸馏模型替代原 Gaussian 682,943 参数 encoder。离线 teacher 仍依赖那个已训练模型、其 latent-supervised TRAIN32768，以及历史数据和物理支出，不能称为免费或无监督生成真实 posterior。

## 4. RRR：训练谱收益没有转化成数据外收益

使用输入 covariance 白化后的输出–输入 cross-covariance SVD，固定 rank {{2,4,8,16}} 和 ridge {{.001,.01,.1,1}}；**没有对 GT covariance 做 PCA 然后称其可预测**。比较 raw-y、decoder-only、physical residual、与物理输入同谱的随机 measurement projection，feature normalizer/PCA 只使用相应 TRAIN。

原 correction-only OOF 以全 TRAIN 训练的冻结 CNN 为 m0，并条件于全 TRAIN 的无监督输入 PCA；这一窄结果保留，不能称 end-to-end crossfit。随后实际完成 18 个 baseline CNN fit 的诚实嵌套二折：每 outer heldout512 的 baseline 仅用另512训练；RRR 用于 outerTRAIN 的 residual 由两个 inner256 baseline 交叉产生。所有 baseline 固定最后4000更新，不通过 outer heldout 或 global VAL 选 checkpoint。C 独立读取全部18训练/heldout parent 集，确认没有 supervised overlap。

| RRR 表示 | 最终 VAL loss | nested baseline OOF loss | nested corrected OOF loss |
|---|---:|---:|---:|
{chr(10).join(nested)}

全部仍由 VAL 选择 rank2/ridge1；nested OOF correction 对 baseline 变差，没有约2%的物理专属增量。结论仅限冻结的线性 rank/feature/数据制度，不能升级成所有非线性可压缩表示均不存在。

剩余边界：原几何 frame 使用历史全 TRAIN online J，是无监督 conditional/transductive geometry，不是完全 fold-inductive 新几何；最终 VAL 的 m0 是冻结 fullTRAIN CNN，训练 residual 来自 halfTRAIN OOF baseline，数据制度存在已声明差别。VAL 选定 rank/ridge 后的 nested OOF 分数仍是次要、经选择的诊断。三个确定性 RRR 行是配对三 baseline seed，而非三次独立随机纠正学习；随机控制与 fold 划分才有单独 RNG。

## 5. 区域与属性诊断

统一用原 y-only anchor 的 field 阈值 ±1 划分 ROI；不是用 GT 或 posterior 的相区。

| anchor 区域 | ROI占比 | Strong CNN 每单元 MSE | K32 | Deterministic matched | Distill128 |
|---|---:|---:|---:|---:|---:|
{phase}

`RESULTS.json` 另给逐单元改善与 pure-inclusion 属性误差。属性诊断是 ensemble 材料在 anchor-pure-inclusion 像素的平均，含 anchor 相区污染，不能当成纯 posterior 属性估计；pure-host 材料本身不能识别两种 inclusion 属性。

## 6. 实测成本与门槛

为防止旧 API 重复执行 encoder 让对照虚慢，优化 K32 只按 parent/batch 执行一次 Gaussian encoder，再按原 parent seed 和 chunk8 RNG 产生32 draw。与存储的原四通道最大差 1.01e-6（预声明2e-6数值容差），不改变原 Gaussian law/frame/decoder。原 API、vectorized 参考和优化 receipt 均保留。

以下统一采用最终 `OPTIMIZED_COST_BENCHMARK.json` 的 **CPU 单线程、batch64、每 parent 微秒**；30 次原始 provider/full timing、p10/p90、环境、冷加载与首次推断均已保存，不能混用不同 receipt 的有利分子和分母。

| 路径 | provider only µs | provider+特征+head µs |
|---|---:|---:|
{costtable}

Distill128 在现有 VAL 上通过原先的“保留至少80%新增风险收益、额外 provider 费用≤1/3”两项点估计：**92.96% 与 {100*gates['distilled_w128_n32768']['provider_cost_ratios']['64']:.2f}%**。如果另要求完整 provider+特征+head 也≤1/3，batch64 不通过：比值 **{100*gates['distilled_w128_n32768']['fullpipeline_cost_ratios']['64']:.2f}%**。二者必须分开报告，不能删除共同 head/feature 开销，也不能把 provider 加速比例当成含 Maxwell anchor 的系统加速。

新训练 fit GPU wall 合计 **{ledger['GPU_training_fit_seconds']:.1f}s**；whole-job 保守 GPU 预留 wall **{ledger['GPU_reserved_whole_job_seconds_conservative']:.1f}s**；从冻结协议到关闭的实际受控间隔 **{ledger['controlled_wall_since_frozen_protocol_seconds']:.1f}s**，均在24h范围。所有新增物理 RHS = **0**。完整 ledger 含失败与所有 objective/data/capacity/fold 修订，未覆盖或重写 R4 ledger。

## 7. 最终状态与停止边界

- `B0_NUMERICAL`: **PASS**，指近似 Gaussian law 的 decoder 积分，不指真实 posterior。
- `B2_COMPRESSION_DEVELOPMENT`: **PROVISIONAL_PASS**（原“额外 provider”定义）；完整 pipeline≤1/3 的更强定义 **FAIL**。
- `PHYSICS_SPECIFIC_RRR_INCREMENT`: **NOT_ESTABLISHED**；冻结线性控制没有收益。
- `FRESH_256_PLUS_256_CONFIRMATION`: **RESOURCE_LIMITED / NOT_RUN**。相同 R4 native96+probit64 anchor 的历史512 parent成本约1,934,840 RHS，远超独立 B cap32768；没有用廉价降级 anchor 冒充公平确认，没有以旧已暴露 R4 TEST 重命名盲测。
- 两个 prospective finalists **deterministic_matched**、**distilled_w128_n32768** 与全部候选 checkpoint 已在 `SELECTION_FROZEN.json` 冻结，新 truth 未生成。停止本分支架构/容量/步数扩张；后续独立确认必须遵守既定预算/协议。

以上是可审查的局部算法结果，不是新的基本 Gaussian 积分理论、完整 posterior 资格、物理任务表示原创性、或最终发表接受结论。最终科学判断由主 Codex 进行。

![Track B summary](TRACK_B_SUMMARY.png)

主要机器证据：`RESULTS.json`、`METRICS.csv`、`COST_LEDGER.json`、`SELECTION_FROZEN.json`、`MOMENT_VALIDATION.json`、`JOINT_CDF_VALIDATION.json`、`OPTIMIZED_COST_BENCHMARK.json`、`NESTED_RRR_SUMMARY.json`、各 seed checkpoint 与逐100 update loss JSONL。所有细节只写本 B 子目录。
'''
(OUT/'REPORT_ZH.md').write_text(report)
print('report written; summary figures rendered')

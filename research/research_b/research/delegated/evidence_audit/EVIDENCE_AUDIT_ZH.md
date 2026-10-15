# Research B：历史证据定向审计

核对日期：2026-10-15（Asia/Singapore）。这是只读证据审计，不是 E0 action-ranking 实验，没有新增 Maxwell、训练、SHA 检查或历史结果修改。主会话负责新理论、最终新颖性与科学判断。

## 可直接用于本轮的结论

1. Joint P2 已有单频补偿 witness 与第二频率降低部分事件歧义的信号，**尚未实现主动 action optimization**。这足够形成待检验假设，不能冒称 Research B 的实验结果。
2. A17 的后续 closure 已运行有限 nonlinear reconstruction；早期“nonlinear 未运行”的说法过时。它仍没有建立相对 unchanged receiver 的加速。A23/A24 的失败结果不应被相对 Maxwell 恒等式或单个代数检查覆盖。
3. A30 明确分离了 candidate capacity 与 measurement-conditioned selection 的失败。其 bounded compensation 没有达到原协议的清晰→模糊要求，不能作为本轮补偿普遍存在的先验证明。
4. P2 的 accepted fresh reference 与 full-Mie residual correction 是不同资产。前者是连续插值目标的数值积分；后者只有两个固定 observation、有限已访问 full 样本。二者均不是全支持 exact full posterior。

## 正负证据表

| 来源 | 实际结果 | 适用范围与不得外推的结论 |
|---|---|---|
| P2 非线性补偿 | r=.25R 中心空腔：单频白化距离 raw 1.529652，ε 补偿 1.295221，ε+gain .900170；双频共享 gain 下 2.467944 | 合法 best-found full-Mie endpoints；不是补偿优化全局最小，也不是结构唯一性证明 |
| P2 小空腔/环形负控制 | r=.08R 双频补偿后仍仅 .10750；ring h=.08R 单频补偿后约4.45，双频约8.02 | 增加频率不保证消除弱可辨识；ring 并非任意容易补偿 |
| P2 局部机制 | raw/projected 半径幂次斜率 2.9778116 / 2.9747005 | 不支持普遍 r³ 消去后 r⁵ 的主导机制，原主张已停止 |
| P2 第二频率 | large-cavity 的 solid 事件 proxy 质量 .162591→.000431；large-central .677377→.966140；solid case 的 solid .407769→.833822 | 同支持下 C→D 加信息的数值信号；没有 action 排名与新策略优势，D 未有全支持 full 后验认证 |
| P2 独立几何 context | 128观测、64parent；weighted-vs-single-particle paired regret −.247748，95%区间[−.302452,−.193408] | B1 reference-derived particles 的 simulated context handoff；不是 C1-only 成功，更不是主动 Maxwell MAE 收益 |
| A17 closure | 原36单元35改善/1相同；额外预定99单元95改善/1相同/1失败/2未运行；10个完整额外对象 GN-gap 中位下降50.7% | 冻结规则下历史暴露资产、固定实际 current rank；不是99独立物体、strict blind 或同费用比较 |
| A17 nonlinear | 三完整配对材料误差相对改善10.84%/0.22%/7.16%；child-wall 比2.68/1.74/1.46；第四配对未完成 | 有限真实重建后果，付费更慢；GN-gap 改善本身不能替代最终材料指标 |
| A17 efficiency | top1/最多3交换/geometry cache：原 exchange 时间降低26.89%/9.64%，仍为 receiver 的1.958/1.318倍；Jd 仅占原净增量1.19%/3.82% | 主要瓶颈是 bank/anchor/workspace 及 trajectory；verification ROM NOT_RUN，不能预设 ROM 能提速本轮 design |
| A23 最终 | H1/H2/H3全部 FAIL、STOP_NO_EXPANSION；204输出79材料可行/125拒绝；full-validation 1 RUN、19 REJECTED、184 NOT_RUN | 四个历史三维 vector-DDA 场景；8个 nominal A23 条件0合格，未建立跨家族 ≥10% 材料改善或同质量加速 |
| A23 成本与局部正面 | cached Schur 精确重用成立，完整 cold 准备降幅仅约0.35–0.86%；2003的 A23 raw 改善7.824%/5.817%但材料拒绝 | 代数/局部 raw headroom 不能覆盖完整失败；A23_EXECUTIVE_DECISION 是早期理论建议，实际结果以正式 report/gate/最终一致性审计为准 |
| A24 最终 | CURRENT_COMPRESSED_FEEDBACK_NO_GO；G1 PASS，G2 FAIL（6个有效任务条件超过5%），G3 HOLD，G5未达质量区；224行136OK/4FAILED/84初始化失败 | 24²二维 scalar TM、四个暴露parent；weak 子组 full-feedback改善3.98%仅描述性，3D/NN未运行 |
| A30 capacity | 12复杂parent full-hull offline oracle 改善13.12%，L-weighted mean恶化3.32%；固定支持sparse64容量改善35.69% | truth 只用于 offline 拟合/评分；没有相应在线重建器，不能当生成/成像算法收益 |
| A30 compensation | 29/32合法条件，最小full-F白化距离8.84，复杂条件最小9.18；符合nominal≥3且compensated≤1的parent为0（要求6/12）；E2未运行 | 只否定当前局部见证的继续资格；非全局无补偿证明。q≥1仅为保守充分界失效，不代表物理奇异 |

## P2 reference 的实际数值误差

- 原 coefficient-uniform prior（包括零几何原子），R=1，同心径向四维 latent；kR=1.6/2.4、8角度×2偏振、每频32实通道，Re全部在前、Im随后。
- proper truncated real amplitude gain N(1,.15²) 限于[.7,1.3]，同一个材料世界跨频共享；每实通道方差σ²。P2 在有限频带内固定复ε，没有现实 dispersion 认证。
- `COMMON_RESPONSE_BANK.npz` 是每频625个实际 full-Mie 节点的连续 Chebyshev 插值资产；节点不是离散 prior。112独立prior probes×2freq 的最大响应相对误差2.986e−4，最大白化差.014871；有限点而非全支持界。
- accepted fresh：全部128观测 event/mean/SD 做13/17核对；12/39/48升级17并对照21后无未解决项。CDF是positive interval partition积分，全128网格单调缺陷0；CDF13/17细化仅固定前8加升级3，直接反演quantiles仍只前8。
- 已披露不同积分partition与旧13阶的最大log-evidence差.00998733、event质量差.00550557；前8的partition CDF13/17差.00097375、conditional CDF差.00414377。这些不是 accepted 全域误差上界。
- large-cavity/solid full-minus-proxy correction 分别384/192 paid samples；最大event改变量.00015128/.00010187，最大residual MC SE4.583e−5/6.114e−5。bootstrap未覆盖模型外mass、全域solver、cubature与漏解。

## 缓存实际可得性与重要冲突

本工作树 `Gaussian/B1/b1-research-b` 的 HEAD 为 `7ed4faef9bc3dd214d4de689c2203cc383f7495f`。

| 文件 | 实际位置/状态 | 解释 |
|---|---|---|
| COMMON_RESPONSE_BANK.npz | 当前工作树 `research/joint_p2/results/`，存在；fields_1.6/2.4均(5,5,5,5,32) | 两频 paid full 节点及 proxy coefficients，不是任意材料的缓存full响应 |
| FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz | **当前工作树缺失**；实际保存在 `Gaussian/B1/b1-joint-p2/research/joint_p2/results/` | 可只读引用/复制已有资产，不能降级旧13阶后继续称accepted |
| FRESH_FINAL_OBSERVATIONS.npz | 当前工作树存在，y形状(128,64)，F_parent(64,64) | 两频旧观测已在历史中使用时，不得当作新独立测量；只读检查header，未展开truth数值 |
| solid/large_cavity_FULL_CORRECTION.npz | 当前工作树存在；F_full分别(192,32)/(384,32) | **仅原单频32维full预测**，不能假称已缓存双频full predictions |
| packages/solid/large_cavity_samples.npz | 当前工作树存在，latent(8256,4)/(8288,4) | 多数来自proxy条件采样，corrected event mass重赋权；不是全部full-Mie draws，也不是已验证的完整候选池 |

历史预算P2为3459/4000 attempts、受控wall2191.3570846秒，余541**不属于Research B**。所有历史ledger、Gate和配置保持原样。

## 实际来源身份

完整逐文件路径、大小、mtime、revision归属见 `SOURCE_RECORDS.json`；不进行SHA校验。P2以本工作树HEAD为Git身份；A23当前HEAD `a47cda424dfb55792e69d952a25f54ce9a2a5613`，实际pilot源码在报告中标为 `8e7d8cf209c6e54f55567b940658e1742aa4cc55`；A24当前HEAD `ca80ca2a29867e1b03589061dd5c17688c9acbe0`；A30 `717d3d26a290248bd2b9aeb55f3d03694b407b84`。

A17本地顶层不是Git checkout；读取的是local closure/efficiency最终证据。既有发布receipt声明exchange commit `56097653e6a814399581c3759fc3e3e8f0c2284c`、efficiency commit `df61a1b533682a883cdaf65259b4528c354bd832`（历史baseline `f5bd7d0aa0e229110c425f1becb6384e09663196`）。这只是receipt提供的来源身份，**本次没有验证当前文件字节与这些提交相同**。

最近邻核对见 `FOCUSED_LITERATURE_AUDIT.md`。本审计不作全球novelty裁决。

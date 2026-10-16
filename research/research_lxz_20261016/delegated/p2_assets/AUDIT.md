# LXZ：已付费 Joint P2 资产与目标只读审计

审计日期2026-10-16。只读检查源代码、JSON、NPZ headers和CPU checkpoint metadata；没有调用forward、proxy prediction、posterior inference、采样或训练，没有访问网页、检查SHA或启动campaign。全部30份NPZ资产和15份C1权重均可用，逐字段shape/dtype、source IDs及绝对路径见 [assets.json](assets.json)。旧文件未改变。

## 源版本与可直接复用入口

- B1 LXZ根目录：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-research-lxz`，当前commit `0b9024b11b1bbb338ae334a873e731fb46e4bb66`。
- `src/b1_joint`和`research/joint_p2`相对P2冻结commit `7ed4faef9bc3dd214d4de689c2203cc383f7495f`的只读git diff为空；没有把后来A/B/G/T研究误当P2 reference。
- 原P2工作目录：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-joint-p2`。三个较大资产不在新worktree，但仍在这里：`C1_TRAINING_DATA_REF13.npz`、`C1_FINAL_DATA.npz`、`FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz`。无需重新生成、复制或重付physics。
- C1根目录：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/C1/c1-joint-p2`，commit `83056881e89329212b16f5be47a4d0dbc4465943`。权重实际存在于`p2/results/models/`，不是只有metadata或Release说明。
- 当前LXZ冻结范围来自`research/research_lxz_20261016/configs/PROTOCOL.json`。它选定`solid/large_cavity/annular_cavity/combined`；本审计沿用这四个暴露观测，不更改已启动任务。最初只读建议中的`prior_boundary`保留为备选，`gain_negative_control`为额外nuisance负控制。

## 唯一共同概率target

主target为`joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1`。原四维材料坐标`x=(z_c,z_s,ε′,ε″)`的proper prior是独立uniform `[-3,1]²×[2.2,2.8]×[.04,.16]`，相对于`dx`。Decoder在`src/b1_joint/material.py`：两个不相交compact radial bumps，`u<0`为严格真空，其余为共享复ε；外半径1。无topology输入。

Coefficient≥−1时相应空腔半径恰为0，inactive coefficient仍有完整先验体积。实心几何prior atom=.25，但因为ε连续，不能给一张精确完整材料图非零概率。改用物理半径积分必须保留Jacobians与zero atoms，不能另用radius-uniform prior。

Response bank的物理坐标是`(r_c,h_s,ε′,ε″)`，两频1.6/2.4各625已付费nodes。每轴五个Chebyshev nodes，意味着每轴最高degree4的tensor interpolation；不是离散材料prior，也不是旧A32的固定H0/H1/H2专家。`Response.predict`只评价保存系数，不导入物理solver。Full target `full-vector-Mie-matching-n14:C:bounded-shared-gain`另列，不能将proxy likelihood叫作新full-Mie评估。

C target只用kR=1.6的32个real channels：八角度20/40/…/160度、S1/S2共16复数，先全部Re后全部Im。每个实通道独立Gaussian噪声，sigma来自固定实心ε=2.5+.1i response RMS的.12倍。当前观测JSON保存64个双频channels，C应读取`y[:32]`和`sigma[:32]`；D才使用64个channels且同一个gain跨频共享。

Likelihood位于`src/b1_joint/likelihood.py`，对实gain proper truncated Normal（均值1、untruncated SD .15、范围[.7,1.3]）精确边缘化。保留rank-one Gaussian logdet和posterior/prior截断区间概率比；不是profile minimum，也不是每频独立gain。所有四个推荐观测均用相同C target。

特别区分三个坐标：原材料`x`；prior-normal坐标`u=Φ⁻¹((x−lower)/width)`；C1 MDN的logit坐标`t=log(x−lower)−log(upper−x)`。`u`在原prior下为标准Gaussian；`t`不服从这个standard Gaussian prior。二者不能互换而仍称同一prior。

## 已付费数组：shape与用途

完整绝对路径在assets.json，不以当前worktree缺少大文件为由重跑。

| 资产 | 保存shape | 正确用途与限制 |
| --- | --- | --- |
| `COMMON_RESPONSE_BANK.npz` | 两频`fields/coeff=[5,5,5,5,32]`，coarse coeff `[3,3,3,3,32]` | 已付费响应与连续proxy；mesh axes是物理坐标 |
| 七个`*_REFERENCE_PARTICLES.npz` | `latent[8192,4]`, `weights[8192]`, `reference_log_density[8192]`, `y/sigma[32]`, event mass[12] | order17 stratified quadrature-node resamples，inactive coordinates条件uniform；不是iid posterior参考、连续CDF或真实采样ESS |
| 七个`packages/*_samples.npz` | 约8192–8320行latent/weights/event IDs | 条件proxy粒子；solid/large-cavity事件权重可带full-residual修正，不能说整条joint sample law已full-corrected |
| `TRAIN_OBSERVATIONS.npz` | y[256,64], truth[256,4], parent[256], F_parent[128,64] | 128暴露parents×两噪声，已有physics |
| `VALIDATION_OBSERVATIONS.npz` | y[64,64], truth[64,4] | 32暴露parents×两噪声 |
| `FRESH_FINAL_OBSERVATIONS.npz` | y[128,64], truth[128,4] | 原P2 fresh，现已公开暴露；不能再叫independent new validation |
| `C1_TRAINING_DATA_REF13.npz` | train_y[256,34], teacher[256,512,4], validation_y[64,34], validation_teacher[64,512,4] | 额外posterior监督，teacher/setup费用非免费 |
| `C1_FINAL_DATA.npz` | reference_samples[128,4096,4], log density[128,4096], event mass[128,12] | 旧13阶resampled-node评价，CDF仅node诊断 |
| `FRESH_CONTINUOUS_REFERENCE_ACCEPTED.npz` | grid[4,65], CDF[128,4,65], conditional_CDF[128,12,4,65], mean/SD[128,4] | 最终positive-partition continuous reference；另外保存每row实际order |
| `solid/large_cavity_FULL_CORRECTION.npz` | full/proxy F[192/384,32]，event残差[12,16/32] | 付费full-minus-proxy residual reference；conditional samples/quantiles要分别标target |
| 三个`CORRECTED_PROPOSAL_*.npz` | latent[64,4], full F[64,32], known logq/weights, bootstrap[2000,12] | 一次暴露large-cavity C，64 calls/方法，不是独立普遍排名 |
| 两个`SMC_*.npz` | latent[64,4], weights[64] | 历史generic SMC，不能冒充精确参考 |

## C1权重的适用范围

15个checkpoint来自五方法×三seeds。可先用`full_distill-26101321.pt`，另保留同seed diagonal/single-Gaussian作为真正对照。Full配置为D=4、condition_dim34、8 components、两个96 hidden layers，head weight `[120,96]`；checkpoint含train-only mean/scale、normalization_fitted、FP64语义bounds和目标identity。代码不接topology/state/memory；旧P1的二维权重不适用。

Condition必须是`[y32,log(sigma),.15]`，而且原模型只支持32个同sigma实通道。不能直接传64通道双频、不同gain prior、不同material decoder或新prior。这样的变动须另建target并适配模型；proposal即使可用也需正确density correction。

`c1p2.density`使用logit Gaussian mixture及精确Jacobian得到`dx`归一化density，Cholesky正对角确保数学SPD；raw logits和machine endpoint标记保留。`b1_joint.proposal_correction.load_conditioned_proposal`还校正FP32 softmax export约1e−8级全局categorical roundoff，使sampler与logq相同，不能省略qratio，也不能在评估网格上重归一化。

旧C1 full/diagonal event MAE最终约.013215/.013028，区间跨零；latent KL .216/.275仍是13阶node expectation，不是材料field KL或accepted physical KL。Accepted physical grid-CDF约.099/.102、重要conditional约.142/.162，说明normalized proposal仍远非完整posterior。这些权重可以提出跨区域candidate，不可将其raw frequency自动当物理posterior概率。

## 四个现有暴露观测

都取`DEVELOPMENT_CASES.json`（seed261013201）中的C-slice，精确JSON pointer、y/sigma和参考在assets.json。当前LXZ protocol已冻结以下四项。

| ID / JSON index | 主要order17 proxy事件质量 | 科学用途 |
| --- | --- | --- |
| solid /0 | event0=.407769，6=.116304，9=.410802 | 实心atom与连续中心空腔族竞争；有付费full residual reference |
| large_cavity /2 | event0=.162591，6=.071867，9=.677377 | 非唯一但偏向较大中心空腔；有付费full residual reference与64-call proposal对照 |
| annular_cavity /3 | event2=.519907、8=.086533、11=.117430等 | 环形/中心连续补偿，更多重要区域；不自动叫多分离峰 |
| combined /4 | event2=.501669、8=.086331、10=.079982、11=.124422等 | 组合truth但后验支持其他几何区域；检验generator是否只追truth标签 |

备选`prior_boundary`（index5）有event2=.519023、8=.109706、11=.280949，适合下一项有授权的prior边界诊断；本审计没有因此改动当前四观测。`gain_negative_control` truth gain1.2在C prior内；A/B假定gain1则是故意误设，不能拿来作公平known-calibration优势比较。

## 参考误差、ridge与mode结论

七例C的13/17 event最大差1.35e−5、CDF约1.17e−4；相同prior的物理坐标5/9 adaptive integration收敛，较tensor17 event差约4.80e−6。Coefficient坐标adaptive最初失败仍保留。112 prior probes的两频full/proxy比较最大相对field差2.986e−4、白化差.014871；这是有限点检查，不是uniform bound。

Fresh128的event/mean/SD均13/17检查，12/39/48先失败后升级17、经21检查通过，无unresolved；CDF细化仍只固定8加3例，fixed65-grid不是连续supremum。原truncated CDF非单调.00265042、新positive partition的旧新差.01476339与event差.00550557均保留。不可将节点CDF或旧近似升级为精确参考。

目前证据建立的是**结构不同、响应可补偿、多个互斥材料区域有显著posterior mass**，没有建立这些区域各自是分离density modes。有限bins、MDN components、不同truth类别、低residual endpoints都不是mode证据。Inactive coefficient有精确flat fibers，中心/环形半径可形成连续族；P2主报告也明确允许ridge。没有发现保存的density-valley、路径barrier或连通性证书。因此LXZ可以比较已声明事件与CDF，但不能仅据现有event mass宣称E3的“两个分离mode”条件满足。

最后，P2为已封口历史阶段：3459/4000 attempts、controlled jobs2191.357s，剩541不是本任务新增授权；旧B1 384/384保持不动。本审计只定位已有资产，新增physics/proxy evaluations/samples均0。NPZ header读取最初遇到NumPy private helper不可用，改用public版本化header reader后完成；没有改旧文件或补跑物理。


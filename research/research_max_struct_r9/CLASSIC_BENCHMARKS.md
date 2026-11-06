# 原始Austria二维采集合同比较（R9 phase2）

运行状态：**RUNNING，16个任务已冻结；尚未读取本轮评分GT**。这是已公开论文对象的开发比较，不是新盲测确认集。原R9第一阶段113次反演的结论与成本账本保持冻结；本扩展单列 `RESULTS/classic_extension_xinan/PHASE2_COST_LEDGER.json`，保留其中主Windows与部分Mac的分项。

硬件迁移：用户在运行中明确要求优先远端i9 P核。现按原科学合同完整重跑**全部16项Windows主cohort**，使用独立 `RESULTS/classic_extension_xinan/`；固定四个物理P核逻辑编号0/2/4/6，各一个单BLAS线程worker，进入NumPy前设置亲和性。本任务CPUonly、不使用GPU、不改其他进程。Mac派发器已停，两个当时在跑的子任务继续完成；最多9个Mac终点作为**单独的部分平台诊断及已付成本**，余7项为用户迁移后的NOT_RUN，不标算法失败。主费用比较不把Windows与Mac混在同一方法表中。两例Windows16项及遗留Mac两子任务全部terminal前，仍不读GT。迁移配置为 `configs/classic_extension_xinan.json`。

## 固定问题与对照

复用R8原始150²合成观测，未重生成、未换噪声。两个400 MHz Austria条件是εr=3.5原论文正例，以及εr=4已知接近方法极限的条件。完整2 m域、20个整圆平面波、40个半径3.75 m整圆接收点；64²反演，4096个独立实材料χ=εr−1∈[0,20]，从空气χ=0和零APCS开始。DPCS秩15；原mf/β：ε3.5用5/7/10与6/2/.5，ε4用5/8/10与8/2/.5。原频率键系数精确零填充跨stage传输。

六个匹配臂：sourcewise CIE single、sourcewise Dual、ordinary GlobalDual、A global、A+ global、direct physical。前三种与A/A+同smooth=1e−4、TV=1e−4先验。A和A+的current步骤与ordinary GlobalDual完全相同；A只换材料目标，A+再换材料读出，不增加优化自由度。另有独立paper NIE、paper LSIE参考，二者均不加材料prior，**不同合同，单列，不称作者代码**。

每stage最多500次并使用论文δ2D=0.01/0.001/0.001（内部平方阈值0.0001/1e−6/1e−6）。sourcewise用每源相对系数变化平方均值，global用整体Frobenius相对变化平方。所有stage材料、current、系数和物理Phi留存；候选只依初值及已付stage端点的原物理Phi选取。随后同一连续L-BFGS-B最多300 accepted iterations；direct也是300。保留连续history中的80步检查点，不以重启代替连续迭代。CAP与优化器收敛分别报告。

物理模型继承R8等面积圆盘脉冲、解析自项与padded FFT，GMRES每RHS最多1000次，检查相对线性残差≤2e−9，失败后同一离散算子exact dense LU。新adapter只缓存不可变dense Green，避免每个材料点重新计算Hankel矩阵；所有臂一致，失败迭代、组装、LU和全部RHS计费。不能把工程加速解读为科学改进。

## 限定与评分口径

原论文来源为[2016 NIE](https://doi.org/10.1109/TAP.2016.2535492)，详见 `delegated/classic_extension/AUDIT.md`。εr=4不是原文承诺稳定成功的条件；原文连无噪初跑也有环内伪影。所有论文对象都是公开开发对象，不能作为未见过几何的泛化证据。

原150²→64²仍有离散差异：R8已测ε3.5、400 MHz的64→150 clean场相对误差约9.96%，与10%噪声同量级。固定范数噪声只定义数据扰动，既不自动给出高斯似然覆盖率，也不包含模型误差上界。物理残差达到噪声量级不等于材料恢复；χ接近20的幅度偏差不能被IoU掩盖。

所有16项完成／显式失败、预测封存后，才单独离线构造公开几何真值并写 `sealed_truth/classic_extension/`。评分并列：论文逐像素相对ε的RMS、contrast relativeL2、R7兼容amplitudeNRMSE，以及Otsu／固定半对比度／绝对0.3阈值IoU。不新增总体二元成功阈值。包含各stage、80步和最终300上限内端点；共用εr=1–21全幅色轴，另有注明饱和的形态图。

在GT评分前另冻结 `configs/classic_extension_paid_budgets.json`：六个匹配臂的300/600/900 **inclusive CPU秒**前缀，只使用实际保存的stage、auxiliary、80步和final状态（以及真实使用过的规定零初值）。先按已付公开Phi选可用incumbent并保存选择记录，再读GT；另列最新可用保存点及覆盖缺口。不插值材料、不按GT选择迭代。checkpoint内部CPU加上 `wait4总CPU−最终内部CPU` 的保守上移，明确标记包含输出尾部；final使用完整wait4 CPU。300/600/900间没有保存材料的区间，真实优化器incumbent未知，不能把粗检查点结果当作精确等CPU截断实验。实际费用仍支付该前缀内全部工作，不只支付所选较早材料点。

## 复现入口

- 配置：`configs/classic_extension_main.json`；16项冻结于启动前。
- 公开数据：`data_public/classic_extension/`，从R8逐字节文件复制，运行只依仓库相对路径；原来源路径／size／mtime用于来源记录，不用于默认运行。未做SHA检查。
- 主程序：`experiments/classic_extension.py`，adapter：`src/classic_extension.py`；import复用仓库内R8/R7/R6代码，不依赖原工作树的绝对数组路径。
- 新任务入口：`python experiments/classic_extension.py solve --task-id <manifest task_id>`，从仓库根使用study前缀；已有terminal输出会拒绝覆盖。
- 一次有限队列：`python research/research_max_struct_r9/experiments/classic_extension.py queue --lanes 2`，最多两个进程、各1个BLAS线程，worker日志和wait4 inclusive CPU见 `RESULTS/classic_extension/queue/`。
- 若从原R8首次打包：`freeze --r8-repo <R8仓库根>`；已冻结配置拒绝覆盖。已经随包提供Y时，无需再次打包。
- 队列全部terminal之后：依次执行 `score`、`figures`、`ledger`。score内部先验证全部16项封存，再进入GT；只读既有预测，拒绝覆盖已有评分。
- Windows主cohort用 `delegated/classic_extension/remote_transport.py` 的prepare/run/status/fetch；SSH全程保持连接，远端有限FIFO为 `remote_fifo.py`。只传必要代码与公开Y/config，首包113,578bytes，无GT。抓回后加 `--output-root research/research_max_struct_r9/RESULTS/classic_extension_xinan` 执行score/figures/ledger。Windows完整child CPU来自GetProcessTimes user+kernel，Mac完整child CPU来自wait4，遗留两子任务如未被旧父进程回收则明确只报内部计时下界及未计尾部，不伪造外部精确CPU。

三维另属原3λ完整矢量合同，不以本二维结果代替。2019公式纯文本抽取曾丢掉根号，R9审计已按原MathML纠正；`data_public/classic_extension/PAPER_CONTRACT_CORRECTIONS.json` 留存更正，R8冻结证据不回写。

## 当前结果

待固定队列全部terminal后填入。没有基于本轮GT选预算、候选、stage或方法。

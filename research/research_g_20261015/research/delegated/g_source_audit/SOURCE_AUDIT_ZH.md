# Research G：只读源资产与表示空间核对

核对时间：2026-10-15 10:09 UTC。新增物理、训练、优化、SHA检查均为0。只读查看历史报告、协议、实现和Git元数据；没有打开locked/test truth数组、A30 `OFFLINE_TARGETS.npz`、P2真值或修改T/A/C1目录。下述不是T当前实验成功的声明。

## 可复用资产与版本

|资产|实际本地版本/位置|可复用内容|不能据此声称|
|---|---|---|---|
|A30|`Gaussian/A30/a30-maxwell-material-manifold`，HEAD `717d3d26a290248bd2b9aeb55f3d03694b407b84`，工作树干净|已冻结24×24标量TM/VIE后端；TRAIN-only Haar/PCA表示；已保存容量与补偿结果|一般3D vector Maxwell、多解在线生成成功、truth-fit系数可供在线搜索|
|C1公开源码|`Gaussian/C1/c1-physics-state-generative-imaging/publication/repository`，HEAD `83056881e89329212b16f5be47a4d0dbc4465943`，工作树干净|P1/P2报告、accepted continuous comparison、normalized bounded MDN源码、训练冻结JSON|P2模型已训练为A30图像生成器，或C1已恢复准确完整posterior|
|C1实际P2权重|`Gaussian/C1/c1-joint-p2`，同HEAD `83056881e89329212b16f5be47a4d0dbc4465943`，工作树干净|15个真实`.pt`，包括`p2/results/models/full_distill-26101321.pt`，各有JSON训练回执|published源码目录有全部权重；本次读取没有执行模型/重训|
|Research A上游说明|本G工作树继承的`research/research_a_20261015/delegated/history/SUMMARY.md`|A17实际源快照位置、冻结GN与非线性成本限制|为A17编造独立Git revision|
|Research T|`Gaussian/B1/b1-research-t`，HEAD `f1b674eb4893e39126701ed48f3d3dcef3098c3c`，工作树干净；`research/research_t_20261015`为空|当前只能定位仓库，未发现Git提交的T physical handoff|接口已交付、corrector已冻结、合法基础解已具备|

外层`c1-physics-state-generative-imaging`本身不是Git仓库；不能把其路径误报为已验证Git revision。所有revision都是Git溯源标识，没有计算或核对任何SHA256。

## A30的合法原forward identity

实际协议为`results/a30/FROZEN_PROTOCOL.json`，原source A29=`5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4`、A28=`c7819fa142bb70f9920835e8e91ba145f15ed8e3`。`src/a30/physics.py`只从未修改vendored A29的`a24.physics`/`a24.contracts`导入，加载已付费known-background与N3缓存；本审计没有实例化adapter。

冻结对象是二维24×24、8源/32接收器、单标量TM，物理网格side=1、发射半径.75、接收半径.8、全孔径。材料为χ，Reχ∈[-.5,4]、Imχ∈[0,.5]，每格mass=1/576，complex128/float64。材料质量坐标、current/field按sqrt(cell-volume)归一化；`L=I-diag(χ)G`，`current=χ*field`，receiver是`S@current`。源/接收器顺序由保存descriptor精确核对。

`a24.contracts.pack`与`a30.compensation.pack`均为**逐source `[Re q, Im q]`**，不是P2的全Re后全Im。Re/Im噪声各σ²/2，白化sqrt(2)/σ。gain固定1，不能偷换为P2径向Mie的bounded gain。Research A E2继承该目标时相容阈值为sqrt(χ²_512,.95)=23.785449004057032；只有target确实相同时才可继承，不能从一个数字反向证明identity匹配。

## GT-only capacity与在线实验必须分开

`A30_E0_E1_DECISION.md`、`README.md`、实际`CAPACITY_SUMMARY.json`一致：12复杂parent先平均weak/ring，再parent配对统计。Sparse64达到相对N3改善.3568962378，探索性95%区间[.2950424211,.4281295418]；TRAIN-PCA64为.2190742349，[.1525004051,.2804729484]。full convex hull约13.12%，实际L-weighted mean反而约-3.32%。

Sparse/PCA离散支持、基底和起点是TRAIN/N3先冻结，但**拟合连续系数和挑选最佳起点的目标使用truth**。Sparse64不是在线重建器，也不保证这些拟合材料通过观测full-forward阈值。`offline_best_s_tail`还另使用truth挑支持，且不强制actual material box；其容量类与实际fixed-support box-constrained QP不同。不能把前者系数、最优support、truth-fit输出当生成器训练标签/初始化。

实际A30 E1只有29/32合法局部补偿对；最小full-F白化距离8.84，复杂条件9.18，没有满足原目标的parent。E1未获进入E2资格，旧结论不改变。本次Research G是新授权，不继承A30剩余额度或扩展其失败主张。

## 连续系数可改变结构，但表示仍受限制

完整材料field是权威对象，低维编码不是世界身份。连续field系数跨越foreground阈值可使一个连通区域开洞/闭合，不必输入孔数；这只是表示可能性，不是Maxwell相容性。

实际A30 Sparse使用24→12→6→3 Haar和3×3粗层DCT，预算16/32/64指独立实坐标。每个start的top-s支持根据N3或TRAIN anchor系数选定，之后不扩展；因此未被该支持表达的细节可能根本不可达。TRAIN-PCA是64个独立parent的未中心化SVD，其inactive位置被精确置零；TRAIN不支持的位置和span外形变不能由连续系数自动补全。两者可有结构变化，但不能冒称完整1152实维材料域。

旧atlas primary还有明确TRAIN数字拓扑保持规则，不能用作本轮无类别生成主架构。level-set允许拓扑变化，但依赖固定TRAIN相位partition、.5格软化、有限DCT位移与区域contrast范围，不能保证任意新孔洞或边界可达。其当前负容量结果是有限非凸优化结果，不是全局表示不可能定理。

C1 P1 P0仅13个vector dipoles与两个材料系数，P2仅四维`(z_c,z_s,Reε,Imε)`。P2 decoder规定同心球中中心及固定半径.65的环形空腔，inactive系数造成非单射和零半径原子；不输入topology label仍不能生成偏心洞、任意多物体或A30全图。P2 density的8个mixture components无拓扑身份，但这不能撤销decoder先验支持。

## C1最新结果的实际含义

已完整读取START_HERE、P1/P2报告、accepted continuous comparison和bounded MDN实现。P2 full/diag-distill accepted event MAE约.013215/.013028，full-minus-diag差+.000187，95%区间[-.000580,.000880]，未建立full事件优势。full重要事件内grid-CDF约.141761、边际约.098906；旧latent KL≈.216仍是13阶估计，不是accepted physical-field KL。

`p2/c1p2/density.py`给出proper bounded logistic Gaussian mixture与Jacobian；`TRAINING_SUMMARY.json`冻结15模型、三个seed、每个3000updates；full_distill seed26101321选择2500 checkpoint。这是measurement-only、34维condition、四维radial proxy target，不是A30 condition或2D full-field生成器。P1缓存physics-state/memory未有跨seed稳定增益，不据此追加GRU/Transformer训练。

冻结proposal的合法full-Mie入口在B1 `src/b1_joint/proposal_correction.py`：明确加载seed26101321指定模型，`q=.9 MDN+.1 proper uniform prior`并计算exact dx logq；没有把residual softmax当posterior。其一次64-draw×3模型development修正不是跨对象优越性。G可保留这类proposal接口设计，但不能直接将其权重迁移到A30。

## A17只按上游说明保留

Research A上游审计记录A17 closure_final没有独立Git；源快照为`Gaussian/A17/CLOSURE_R1/remote_snapshots/closure_final_closed_v1/source`，Exchange baseline记录56097653。十个额外对象冻结GN gap中位下降50.7%不是材料MAE下降。三组非线性pair材料误差改善10.84%、.22%、7.16%，child wall比2.68/1.74/1.46，均更慢。G不能将frozen-model-change公式变成跨材料finite-world certificate或加速保证。

## 对正式新物理比较的交接事实

当前时间快照未找到T提交的目标/config、corrector、validator、基础解及独立cost identity；所以这份源审计**没有建立E1所需接口齐备的证据**。历史A30后端存在，不等于G已经冻结相同target/corrector并取得4个合格parent。历史A E2合法base只3/8，所有ring失败；唯一额外成员与base改变少于6格，不能作为新结构成功。

G现在可以做零物理E0、连续材料proposal和archive开发；能否用完全相同历史冻结后端启动正式E1，由主Codex核对完整physical handoff与合法base后决定。本代理没有发消息给T，也没有依据未提交文件启动实验。

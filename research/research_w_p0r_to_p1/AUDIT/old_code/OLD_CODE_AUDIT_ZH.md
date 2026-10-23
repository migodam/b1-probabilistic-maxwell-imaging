# Research W 旧 FEM / 搜索程序只读审计

审计范围：旧包 `Gaussian/Research_W/research-w-continuation` 的程序、已保存 JSON/NPZ 与账本。对照 B1 固定发布提交 `483d260c` 的 `research/research_w_continuation/`；`search.py`、`inverse_support.py`、`fem.py`、`experiments.py` 四文件均逐字节相同。当前工作树 HEAD 为 `33471215`，本审计不把 HEAD 当作原始实验的提交。所有旧文件保持原样。

**本审计新增 Maxwell / FEM / JVP / VJP / LU / RHS 动作均为 0。** 只运行本目录 `analyze_saved_only.py`，导入 NumPy、SciPy 的连通分量函数、标准库；没有导入 `wq`、没有实例化任何正演对象。未运行旧实验或测试套件。

## 1. 主结论与证据边界

旧搜索没有进行至收敛的普通反演。每方法只有一个初始提案、一个 Full-J 正确器方向，以及最多四个有限半步。其结果为 `d_cal=19.9123–49.8925`，证明这些短程候选不满足原 `d_cal<=1`，不能据此给出不可辨识性、唯一性、普通 FWI 已失败或 Task-aware 有优势的结论。

能确定的缺陷包括：原错误版本正确器锚定被拒绝提案；修正版仍用不连续 top-384 支持投影处理无约束全域方向；原最小特征 0.1 没有被实现；没有记录投影后的预测下降与支持变化；所谓多启动只是一条 BP 邻域扰动。能确定的工程问题包括每个新候选销毁正演后端及 LU、反复组装所有频率，以及原 preflight 的 DOF 字段误标。

**不能确定“大残差就是某一个导数 bug 造成”。** 已保存固定离散 FD 与伴随证据反而支持导数实现。选定三个校准端点的细化变化仅 0.000586–0.004319，远小于其 19.9–47.1 的拟合差；至少对这三个端点，普通场离散误差不支持成为几十量级残差的主要解释。未对所有候选重做求解，不将此推广成全候选数值资格。

以下源代码行号均相对旧包根；结果字段保留历史值，不据本审计修改。

## 2. 历史事实：搜索实际做了什么

|项目|源证据|结果与限制|
|---|---|---|
|共同初值|`src/wq/search.py:71–80,114–118`|背景 Born 小矩阵求解后，归一化无穷范数强制不超过 0.35，再做非负/面积投影；不是完整非线性初始化求解。|
|普通步|`search.py:21–25,118,122`|全域无约束 ridge 步，归一化每分量最多 0.2；之后才做非凸约束投影。|
|方法长度|`search.py:124–145`|无外层迭代、无收敛判据、无多频延拓或分支继续。一次 proposal 和一次 corrector，接受后立即结束该方法。|
|有限回溯|`search.py:138–143`|仅 α=1,0.5,0.25,0.125；首次改善即退出；最小 α 仍失败则 no-op，不继续缩小、不换方向。|
|多启动标签|`search.py:121–126`|`MULTISTART_FWI=random_smooth(...,.35)` 仍从共同 `z` 出发，再统一把整个 proposal 缩到最多 0.2；不是多个独立初始化各自运行至收敛。|
|先做任务方向|`search.py:117–118`|所有方法开始前已支付任务态 J。共同 `shared_rhs=160` 含校准 80 与任务 80；不能把单方法增量 RHS / 时间称为端到端公平成本。|
|未完成 context|`search.py:111–113`|只检查剩余搜索相位额度至少 800；不是原子整 context 预留。历史余 480 时第四 context 直接 NOT_RUN。搜索相位全局账本仍保留首次错误版本费用。|

结果定位：`results/p0b/search-corrected/RESULTS.json`，以及 `P0_B_GATE_DECISION.md:3–18`。本审计重新用保存的候选 data 与 `results/wall-pml/ARRAYS.npz` 算距离，12 个值最大误差为 `7.105427357601002e-15`。不是转录或白化标量计算错误。

|seed / band|初始 BP d_cal|FULL_J|RANDOM_MATCHED|MULTISTART_FWI|TASK_AWARE|
|---|---:|---:|---:|---:|---:|
|2026102201 / wide|45.354386|45.354386|45.354386|45.354386|45.354386|
|2026102201 / narrow|56.984289|29.942024|19.912259|32.065442|26.713549|
|2026102202 / wide|51.208123|49.892531|48.699200|49.892531|47.129469|

第一行四个候选数组与共同 BP **逐元素完全相同**，并非四个不同的近等价材料。对应所有 α 的正确器轨迹也相同（例如 `RESULTS.json:11–49` 与 `66–104`）；最低 trial 仍为 46.344392，高于 BP。该行不能作为四算法独立到达同一解的证据。

## 3. 支持投影：确定的实现缺口与可检验机制

### 3.1 top-K 改写了实际方向，局部模型没有跟踪这种变化

`inverse_support.py:108–141` 先把 contrast/sigma 裁剪至 `[0,8]` / `[0,.15]`，将强制空气内整像素置零，再按 `contrast²+(8/.15*sigma)²` 排序，只留 `floor(3/(12/1536))=384` 个像素。这个规则只表达幅度、被迫空气和保守面积；没有任何最小宽度、曲率、孔洞、连通形态或尺度约束。

而 `search.py:21–25,119–120` 使用未投影的 `r+J*a` 来造方向、挑 task 方向正负。实际进入 FEM 的是 `legal(z+a)` (`search.py:19,126,140`)。投影可能同时删除原活跃坐标并引入新坐标，实际位移是 `normed(legal(z+a))-z`，不必等于 `a`。因此未投影模型下降不保证投影后下降。这是程序控制流确定的失配机制；本审计没有保存各首次 proposal 的 J/r/raw proposal，不能反推出每次失败的精确因果份额。

保存数组给出的客观支持变化如下（全部见 `saved_candidate_diagnostics.csv`）：

* 12/12 个最终候选正好保留 384 个活跃像素。
* 最佳窄带 RANDOM_MATCHED 相对 BP 新增 246、删除 246 个像素，即最终支持中约 64.1% 来自支持交换。FULL_J 为 244/244；TASK_AWARE 为 224/224。这些有限更新不应被解释为沿单一平滑材料支流的微小切向移动。
* 12 个候选有 20–47 个四邻域连通分量，每个有 1–13 个单像素分量。四邻域分量是格点描述；它本身不是用户所定义所有几何最小特征规则的严格认证。
* 格宽 dx=4/48=0.083333、dy=3/32=0.09375，均小于冻结的 0.1。单像素结构没有被禁止。Gaussian 提案平滑 (`search.py:32–34`) 不会让最终 top-K 投影保留这一尺度约束。这个缺口已经在历史 `P0_B_GATE_DECISION.md:14` 与 `docs/IMPLEMENTATION_REVISIONS.md:12` 公开，不能事后放宽阈值来恢复见证资格。

### 3.2 面积规则偏保守而且易被误读

投影按整像素收费 (`inverse_support.py:132–140`)，实际物理面积按扣除 query 的精确交叠面积计算 (`95–105,144–148,166–168`)。所以 384 个像素最终常显示小于 3 的面积（如 2.965833），并非预算计错或偷偷多占区域；只是这套投影没有使用所有允许的部分像素面积。clip+按裁剪后能量保留也仅是声明过的启发式，不宜宣称是任意越界输入的精确最小距离投影。

### 3.3 新实现的明确修正方向

1. 先把 `minimum_feature=.1` 写成可执行、独立于 truth 的几何/表示约束；同时定义其在边界、窄颈、孔洞与 query 剪切处的判定，输出合法性证据。仅把 nx 改小不足以自动保证所有特征条件。
2. 将固定支持内的连续优化与支持交换分开；在固定合法支持内求带边界/信赖域的步，支持交换必须产生合法的新几何再接受。不把无约束 ridge 后 top-K 当作可微延拓。
3. 对每个 trial 记录真实投影后位移 `delta`、`||r+J*delta||`、实际下降、支持新增/删除、触界数和接受锚点。在正式物理预算内用这些量区分局部模型失真与非线性失真。
4. 只有校准已达标、合法材料条件与端点数值资格通过后，才将 task 方向用于沿近等价集的分支探索。大 residual 初值处的广义比不是纤维证据。

以上是新实现建议，未在旧包执行，也没有据此获得任何改进结果。

## 4. 正确器控制流：旧错误已修、短搜索的局限仍在

历史错误路径仍可用 `main('search')` 调到：若首次 proposal 被拒绝，`search.py:127–129` 虽保存 BP 为 best，`133–136` 只在 mode=`search-corrected` 回退锚点；旧 mode 的 else 继续在被拒绝 theta 上建正确器。其后 trial 与 best 比较是保守的，**不会自动把差 proposal 保存为胜者**，但浪费了纠正方向和预算，也污染了方法比较。历史日志 `results/jobs/p0b-search-001/OUTPUT.log` 与 `docs/IMPLEMENTATION_REVISIONS.md:9` 已将此轮排除，保留 760 search RHS 与 133.673416 process CPU 秒。

修正版 `search.py:133–137` 在拒绝时恢复共同 theta/data/z/r/shared_j，逻辑一致；接受时使用该 proposal 的缓存 J/data，也一致。没有发现修正版混用别的 theta 的 J 的证据。

第一 seed 宽带 FULL_J 的首次 proposal 就是共同 ridge 步。被拒后回退 BP，再构造同一 corrector，因此 α=1 是同一个已拒候选。由于精确 theta 缓存命中，该重复 trial 不付新 RHS，但仍重复一条已知不改善的决策；其 method RHS 为 140，而其他从不同提案回退者为 160。该 20 RHS 差是缓存/重复候选导致，不能解释为更高效的算法证明。

建议新入口显式禁用旧 `search` 错误路径，维护唯一 accepted state（theta/data/J/metadata），并把首次 proposal、回退与全部 trial 写入完整轨迹。不要悄悄覆盖旧结果。若开展普通反演，明确设置多次 accepted 迭代、约束的一阶停止指标及实际预算，避免继续沿用“多启动 FWI”标签暗示充分优化。

## 5. 导数与 FEM：支持正确性的证据、未覆盖的部分

### 5.1 没有找到可支持“J 符号/尺度/共轭出错”的证据

* `search.py:17–20,40,65`：残差 Re/Im 与 J 统一除 .01；J 列乘 `[8,.15]`，与归一化 z 坐标匹配。
* `inverse_support.py:231–288`：volume 像素 J 包含 free-space 直接接收场和 transpose 解形成的散射接收场；sigma 块是 `i/f` 倍。这里用 `trans='T'` 是复双线性 J 构造，不应机械改为 Hermitian。
* `fem.py:677–691`：volume JVP 显式包含 `delta_epsilon*state.field` 和 `chi*delta_field` 两项，没有漏掉对接收体积分显式材料项的导数。
* `fem.py:720–743`：VJP 使用共轭接收 Green 和 `trans='H'`，匹配实内积伴随；与上述复 J 的 transpose 用途不同。
* `results/fem-derivatives/RESULTS.json`：实伴随缺陷 `3.150179205e-14`，最佳中心差分相对误差 `1.488888079e-9`，volume reciprocity `3.904804731e-15`。
* `results/p0b/preflight/RESULTS.json:3–8`：像素 FD `3.754391893e-9`。本审计直接由保存 FD/Jdirection 数组复算得到同值。

这些证据针对固定离散材料映射。preflight 的方向由 `v=theta*random*.1` 生成 (`search.py:97–102`)，本审计证实其 768 个非零参数以外方向全为零。因此它**没有**检查从空气出生新材料、top-K 支持交换、投影导数、几何/形状导数，也没有证明所有物理频率/新候选均准确。

`fem.py:641` 的私有 `_state_action` 将所有非 `'N'` 情形按 Hermitian 计算 residual。目前仅被 N/JVP 和 H/VJP 调用，故不是观测到的错误；新实现若让它接受 `'T'`，必须新增真正 `matrix.T` 分支。当前像素 J 独立在 `inverse_support.py:260–270` 对 `'T'` 使用正确的 matrix.T residual。`volume_pixel_jacobian` 也没有主动拒绝 probe state，但历史 caller 明确使用 volume (`search.py:53`)，不是本次结果的错误来源。

### 5.2 线性解残差不是拟合残差

`fem.py:569–576` 检查 `||A u-b||/||b||`，默认 tolerance=1e-9 (`297`)；`search.py:85,123,145` 的 `d_cal` 是 `||F(theta)-y||/.01`。两者不是同一种残差，不能把“20–50”读成 LU 失效，也不能用机器精度线性残差证明材料拟合成功。

参考世界使用连续 wall 系数 `experiments.py:13–21` 与 interface-fitted 网格；搜索使用完整像素系数 `inverse_support.py:151–181` 与 pixel-fitted 网格。这是有意的表示差异，不是自动构成代码 bug。表示下是否存在低 residual 的合法材料没有在旧短搜索中证明；对“像素表达误差导致下限”的说法仍须新的、有预算的证据。

### 5.3 端点误差足以排除一个过强解释，但不覆盖 task 的所有情况

`results/p0b/finer-endpoints/RESULTS.json` 的三个选定 cal 端点变化为 0.0043188、0.0005858、0.0036904，d_cal 仍为 19.9118、45.3544、47.1294。它们不支持把几十的失配主要归因于场网格误差。原 unfitted task 端点变化为 0.198–0.270，超过 .1；其后 fitted 圆周 256→512 检查为 0.000644–0.000923 (`results/p0b/fitted-task-endpoints/RESULTS.json`)。这一修正没有令校准材料变成见证，更没有给所有 task/所有候选建立连续误差界。

## 6. 缓存、资源与复现元数据

### 6.1 数值缓存未发现错配，保留策略却很昂贵

`Forward.evaluate` 只在 task=None 且 theta **完全相等** 时返回 data/J (`search.py:44–45`)；否则 `clear()` 清空 backend、theta/data/J，再逐频构造 FEM (`48–53`)。每个频率结束后 backend 被 pop、state 被释放 (`66`)；`cached=False` 分支和 `states` 列表实际上不形成跨频或跨材料复用。task 分支结束又 clear (`69`)。这降低驻留内存，却意味着每个新 trial 重建全部 basis、背景稀疏矩阵、Green、LU。

底层 FEM 的材料缓存键为完整 quadrature epsilon 精确相等 (`fem.py:456–458`)，失败前主动清空旧 LU (`459–464`)；state 持有自己的 factor/matrix (`609–619`)，导数校验 backend token。这些是正面的 fail-closed 设计。瓶颈在上层主动丢弃后端，不是发现了“错误材料复用旧 LU”。

首次 proposal 在知道是否接受前已经付完整 J (`search.py:126`)。五频每个 F+J 付 `(4+12)*5=80 RHS`，只做 F 是20 RHS。若拒绝并回退共同锚点，这 60 transpose RHS 不再用于正确器。新实现可以在明确预算规则下先 F 判断接受，接受后才在确切同一 state 求 J；需要显式管理内存和一次性 state，不能靠增加隐藏求解来声称节省。

### 6.2 历史 preflight 字段误标

`results/p0b/preflight/RESULTS.json:8` 写 `ndofs:3795`。同作业账本 `results/p0b/PHYSICS_LEDGER.jsonl:1` 明确写 P3 `ndofs:33415`、nelements7344；3795 是该 mesh 的 vertex 数。发布 source `search.py:103` 已改成 `mesh_vertices=fw.mesh.nvertices`，但旧保存行未改。它是元数据字段名错误，不是 FD 用了3795 DOF 的证据。新报告需同时给出 mesh vertices 与 FEM DOFs，保留这项历史差异。

### 6.3 预算与结果轨迹的建议

`SearchBook` (`src/wq/accounting.py:45–47`) 保持四 phase 独立上限；每次数值动作前 reserve，不能把旧失败成本退回。`search.py:111–113` 的整 context 可用性检查并非事务式预留；若未来允许并发物理作业，必须在共同账本上添加真正 context/method 预留，不能只在前端 sum 后检查。历史未见并发绕过，因此此处仅标将来控制流风险。

`search.py:94,141,145` 保存最终候选与 corrector trial，但没有首次 proposal 的 objective、拒绝原因、原始方向、投影后预测、支持差、accepted anchor 标识与梯度停止量。靠现有 NPZ 无法精确重建每次拒绝到底源自支持跳变、非线性、幅度限制还是候选数不足。新 campaign 若要解释原因，必须把这些小型诊断随每次真实动作记录，而不是在失败后补推因果。

## 7. 纠错优先次序与允许的结论

|优先级|修正|旧证据能说什么|
|---|---|---|
|P0|落实全部合法材料条件，尤其 minimum_feature；停用错误 legacy corrector 入口|原候选没有完整合法性资格；旧错误路径确实存在且已排除比较。|
|P0|把普通拟合做到有记录的多次接受迭代，设置真实预算与约束停止指标|旧残差只说明短程候选失败，不说明普通 FWI 已收敛或物理模型无解。|
|P0|基于投影后位移记录预测/实际下降，明确 accepted anchor|支持大交换已从保存数组确认；对某次拒绝的因果解释尚缺轨迹。|
|P1|合法支持内连续步与合法支持交换分离；明确多启动数量与覆盖范围|当前四标签不构成充分分支覆盖或统计优势。|
|P1|正确 cache/state 生命周期、F-only 首次接受筛选、精确成本归因|可减少冗余设置与被拒 proposal 的未用 J；实际节省尚未测量。|
|P1|只对需要的真实新接口补方向/伴随校验，保留现有通过证据|没有依据将旧几乎机器精度的固定支持 FD 判成失败；它没有覆盖投影与形状。|

**本审计最终状态：READ_ONLY_COMPLETE；NEW_PHYSICS_NOT_RUN；修正后的收敛性/可辨识性结果 NOT_ESTABLISHED。** 所有建议均待主线程以新冻结协议和实测结果判断；本审计不改变旧 P0-Q/P0-B 判决、不启动 P1。

## 本目录产物

* `OLD_CODE_AUDIT_ZH.md`：本报告与准确源行定位。
* `saved_candidate_diagnostics.csv`：12 个保存候选的 residual 复算、支持差、格点分量与幅度。
* `saved_only_audit.json`：源逐字节对照、保存 preflight 复算、原账本分组摘要；早期 reserve 记录缺 pid，不能把 pid 分组解释成完整每作业成本。
* `analyze_saved_only.py`：可复查的纯保存数据分析脚本；没有 PDE 调用。

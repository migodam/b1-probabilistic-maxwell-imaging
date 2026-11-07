# R11 vector完成复现与历史边界

本文件描述已经完成的生产者与只读复算入口。本轮科学已结束；不在当前目录重新启动任何Maxwell反演、benchmark、corrector扩展或新设计。建议在归档副本复算，保留原始预测与public seals。

## 真实入口与冻结输入

工作树：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r1/r11_worktree`。实际Mac Python：`/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python`；Windows：`D:/python/python.exe`，历史remote repo `D:/AI/MAX_STRUCT_R10_VECTOR_20261107`。

真实依赖为现有R10 `vector_grid/vector_cuda_model/vector_physical/vector_cie/vector_cuda_dense/vector_resources`及R9资源job helper；R10不改写。物理是12³ dyadic CM+RR，3N interleaved电流与每cell一个实χ。新代码只在R11 `src/vector_denominator.py`、`src/vector_cie.py`材料回调、薄backend/storage/resource包装；forward未改变。

正式合同：`configs/vector_confirmation_contract.json`、`vector_confirmation_material_design.json`、`vector_confirmation_queue.json`。源基在Y生成前由公开Einc与G Einc构造，relative SVD阈值1e−10、实际rank120；不是全G谱截断。源基与传播遗漏诊断、唯一必要full/lagged/z新梯度接线检查均在本目录独立receipt保留，不需重复。

## 已完成科学次序与源码证据

1. ordinary DEV6（default库）→一次公开3状态×default/cusolver=6LU完整value/gradient/原方程残差→new DEV6（选cusolver）→fresh50。probe嵌套在new DEV6，库差不是纯算法速度。DEV12先封存后仅旧曝光GT评分。
2. 新观测生成：8个12³ dense-LU/all60与预选c13/c18的18³ GMRES/all60；3% iid CN无样本重归一化。10个public observations只含Y与公开noise metadata；private generation仅写GT，反演从未读GT。
3. host最初50GPU计划，后25/25，再基于已付冷首臂38.6min改35GPU/15Mac；Mac15完成后按实际10个spool臂中位21.11min和条件授权将完全未启c18 matched+cross10转回Mac。最终25/25，每对象五方法与两grid同host。所有旧manifest/config版本与时刻/预测保留；最终实际成本替代早期ETA。
4. Mac c11五臂旧RAM版自然结束；仅空闲controller受控drain，c14/c15与随后c18两组用单临时memmap spool保存历史J/Q视图，数学/FP64/cache-hit/预算不变。20个spool均在原数组/前缀封存后释放并删临时文件，原输出完整。
5. GPU original25完成15臂后，c17 Dual在315 accepted/317完成requests之后、下一FW之前nvidia-smi10s查询超时。原attempt0额外320×5184LU/19200primal/19020adj RHS全计；80/160/240快照保留，315态与拟牛顿历史未持久化。root授权一次同初值同预算独立重跑。
6. retry1 c17 Dual自然320封存；下一c17Global路径预建明确sentinel，既有exclusive mkdir在其初始FW前受控结束。finally全部LU/primal/adj计数等于唯一已完成arm，证明下一臂精确0FW。raw FileExistsError不是算法失败。
7. 后9原未启任务在新进程：仅查询超时杀自建nvidia-smi helper并有限重试；成功external=[]快照缓存30s，过期/外部出现先暂停并重查，最大新增外部探测延迟30s。资源等待改善不是科学算法提速。新增seal+归档后stop-request接口，本轮未使用该停止请求。
8. 全50预测/400prefix公共联合封存后，root独立逐值审查并于UTC2026-11-07T07:53:45.607075写`audit/VECTOR50_ROOT_POSTSEAL_AUTHORIZATION.json`，之后才运行唯一fresh scorer。原public seals的GTfalse保持历史记录，不回写。

源码/输入快照与成本证据在`remaining_transport/SOURCE.zip`及其各overlay、`c18_append_transport/SOURCE.zip`、`gpu_recovery_attempt1_transport/SOURCE.zip`、`gpu_recovery_attempt2_transport/SOURCE.zip`、各PROCESS/FROZEN_CONFIG/SOURCE_EVIDENCE。source快照/size/mtime与registry保留，不新增hash检查。

## 只读复算与离线产物

原评分数据已在`confirmation_postseal_scores/`：50主表、400冻结前缀、已付轨迹、固定accepted80/160/240/320、CN99与物理残差、全30图。fresh scorer真实入口为`experiments/vector_postseal_score.py`；外部wait4入口`experiments/vector_fresh_postseal_external.py`，必须有后来root独立授权与all50公共barrier。当前已完成，不在原输出目录再次跑；脚本exclusive输出保护会拒绝覆盖。

只从已评分CSV复算机器汇总/本报告与去重账本（0新FW/0GT数组访问），在归档副本的worktree中：

```sh
VECTOR_PYTHON="/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/.venv_nn/bin/python"
"$VECTOR_PYTHON" research/research_max_struct_r11/delegated/vector/compile_completed_analysis.py
"$VECTOR_PYTHON" research/research_max_struct_r11/delegated/vector/finalize_completed_costs.py
```

`compile_completed_analysis.py`只读已授权评分CSV/公开合同，输出completed_analysis与REPORT；不按GT重选状态，不新增成功门。`finalize_completed_costs.py`调用原metadata costcompiler并纳入本次报告生产者，所有重复元数据编译run的CPU历史保留。总账scope含本worker全部vector DEV/setup/diagnostic/generation/fresh/失败/transport/helper/score/plot/summary；root环境清单与2D/理论/发布单列排除。原failed wholephase替代15个arm而非叠加；Macwait4含子孙/RSShelper；含null gaps的平台已知CPU下界不称精确总费，Mac/Windows秒不相加称算法速度。5184 full-LU与32/120 small-LU、17496 iterative工作量分开。

科学历史runner的实际接口是`experiments/vector_run.py --config ... --queue ... --offset ... --max-tasks ... --output ...`。原50、Mac分组、recovery1/2的真实命令已完整保存在各controller PROCESS_RECEIPTS/ACTIVE、transport manifests与FROZEN_CONFIG；历史transport脚本固定exclusive输出，不能在已完成archive里直接重跑。若未来明确要求独立科学重现，须复制成隔离输出并保持相同公式、初值、数据、rank120、预算/FP64；本轮没有授权追加这种重跑。

## 路径与统计限制

`ATTEMPT_PATH_REGISTRY.json`将原c17 Dual历史canonical路径解析到archived attempt0，将recovery1 sealed路径映射到当前canonical；历史receipt路径原样保留。recovery2是后9个原未启动任务的执行代际，不是九次科学重跑。共50完成预测+1个基础设施失败attempt，0FW sentinel/取消waiting controller不算科学attempt。

8matched材料是规定设计有限集，不是iid总体；cross只重复其中2个材料，不扩N。主8表不删c17；同成本敏感性整例排c17，因唯一retry prefix非端到端。全部连续指标、分箱bias与图保留；CN99资格、IoU连续值、主轮廓可辨认与幅值准确分别解释。没有新科学调用、bootstrap或新3D成功阈值；本轮分析与交付后结束。

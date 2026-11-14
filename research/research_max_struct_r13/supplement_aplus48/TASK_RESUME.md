# ACTIVE CONTINUATION: finish publication; do not restart experiments

All184 new runs and all requested analysis are COMPLETE. Read RESEARCH_REPORT_ZH.md, CHAT_HANDOFF.md, THEORY2_MATERIALS.md. Arc andMM final48 nearly identical; Arc faster; NN effects solver-dependent; costs/dynamics audited. New112 sealed/accepted/authorized/scored; original params and selections unchanged. New derived summary-label correction only, no prediction change. Independent table audit148checks0issues.

Artifacts prepared by src/package_supplement.py:6assets total1215403300B (see audit/RELEASE_PACKAGE_RECEIPT.json), all on external disk, not yet uploaded. Website built locally (94actualfigures/8galleries), no downloads until actualmanifest. Remaining: selective source/compacttables/figures Git commit+branch push; create new same-repo supplementRelease, upload all6assets, verify public headers/range read; createDOWNLOADS+ARTIFACTS; rebuildsite; add top-levelR13 addendum/banner preserving oldbody/results; normal pushmain; read back publicsite/images/report/downloads (CUA ifused must completion_summary). No new algorithm/data/network training/CUDAport/benchmarks.

Current branch codex/r13-aplus48-20261112, origin/main still initial2e32f0626 at last fetch,3localfreezecommitsahead. Huge olduntrackedR13 results unrelated mustNOTbecommitted; stage only supplement allowlist and explicit mainaddenda. Exclude private_operational, archives asGitblobs, D:alias trees, gradientslargearrays fromGit; allarrays goRelease. Figures/compactCSV+JSON inGit.

# 本次任务与续接

用户请求：把原A+补到相同48材料，与CIE、Dual-CIE、G19、DR同成本比较，参数/选根不变；比较A+与DR擅长对象和失败动力学，寻找算法/理论启发。

必须完成：未插桩原A+缺失32份；原16复用清单；全48配对raw/corrected/primary及成本曲线；失败/成功机制诊断及合理边界；原始数据、费用、中文报告和网页/GitHub发布。不得看到结果后调参或称新盲测；不得把辅助状态插桩的额外开销藏入原主比较。

保存位置：R13/supplement_aplus48；隔离分支codex/r13-aplus48-20261112。保持旧R13 seal/RESULTS/configs及R4–R12不变。所有大型文件只在外盘或任务专用远端目录。

## Completed science / publication in progress — 2026-11-14

User's active task: finish unchanged A+ on the same48, reconcile historical algorithms/NN/corrector with native dynamics and Theory2 material, use remote resources reasonably, diagnose low busy-P frequency, and publish the complete supplement. Do not stop at an interim summary.

- Main32 extension COMPLETE, downloaded, public-only acceptance32/0issues, sealed and authorized. Original16 reused explicitly; no frozen R13 predictions changed. comparison48 analysis complete: primary A+-MM NRMSE .1521001711; equal-actual-CPU paid prefixes .1574596474 versus CIE .1946844599, Dual .1797826225, G19 .1822155553, DR .1767194422. Exposed existing48 follow-up, not fresh blind.
- 40 native diagnostics COMPLETE and separately charged. Downloaded to dynamics_run, acceptance40/0issues, separate authorized seal. Native analysis complete561snapshots/844LU/18924complexRHS/0issues; parent interpretation and evidence summaries in progress. Observer replays never replace primary results.
- Historical ledger found R13 default A+ uses quintic-MM without R7's final arc/guard options. Existing main comparison remains unchanged and explicitly named R13-MM.
- Registered fixed-stack112 under configs/aplus_stack112.json and STACK_COMPARISON_CONTRACT.md before execution, commit0a8c39de1. Remote D:/AI/MAX_STRUCT_R13_STACK_20261114 running original P8/BLAS1/120CPU/72+48:48uniform R7-readout Arc controls +64NN interfaces (CIE/Dual/MM/Arc on fixed16). Reuses frozen R8 predictions, no training or tuning. Live foreground session86544. Remote status42finished/0errors at latest check; do not launch duplicate work.
- After112 terminal: private_operational/download_stack.py -> stack_run; public acceptance/seal; parent authorization authorized_stack_scoring with four seals; src/compare_stack.py. Keep newArc/MM columns separate, compare NN same16 and raw/corrected, count historical training402.558889MacCPU separately from no-new-training.
- CPU read-only audit confirms loaded P8 low frequency; AC already BestPerformance, no explicit cap/99%/disabledboost/EcoQoS. Temperature/packageW/hardwarelimitreason/OEMmode unavailable. No settings changed. New E/HT authorization should not change frozen comparison hardware midway.
- Remaining: finish112, analyze actual native evidence and parent-review figures, unified Chinese report/Theory2/version+NN+corrector/cost/failures/reproduce, new raw-data Release plus supplement webpage, GitHub main normal fast-forward (no reset/force), live readback images/downloads. Exclude private_operational, private SSH logs, local D: manifest-alias trees. Big files external/remote only. No new NN training/new random objects/new algorithm cycle.

Prior progress log below is historical and superseded by this status.

## Historical progress log
（追加）

- 原C2源码/运行时复制到 `D:/AI/MAX_STRUCT_R13_APLUS48_20261112`，原legacy依赖 `D:/AI/MAX_STRUCT_R11_DEV_20261107`。32缺失A+已全部COMPLETE，0失败；120CPU/72+48/8原P核/原选根不变。
- Git冻结补32配置提交 `2e6916aac5e6984292a65773b40ec2150087055d`，未改旧R13任何seal/RESULTS/config。新本地分支codex/r13-aplus48-20261112。
- 初次detached启动被SSH关闭终止，发生在queue创建前；改为foreground有限队列后完成。费用/故障记录在audit。执行session12851已返回0，564.826s wall。
- 完整主结果ZIP远端125854107B，下载脚本 private_operational/download_primary.py 正在执行session22478，落外盘SUP/archives，之后抽取到SUP/RESULTS、SUP/queues，拷贝原公开输入到SUP/data_public。
- 第一次Windows默认文本编码导致注释误显示，explicitUTF8重读后8个实际继承模块源码与本地完全一致；没有SHA验证。
- 没有打开新32评分。下一步：下载完成后parent/src/audit_predictions.py --root SUP --manifest SUP/configs/aplus48_extension32.json --output SUP/audit/PUBLIC_ACCEPTANCE；parent/src/seal_results.py到SUP/audit/aplus48_extension32_seal.json。父验收后才将新seal authorized_scoring=true，写SUP/audit/COMPARISON_AUTHORIZATION.json (authorized_comparison_scoring=true,scope FOLLOWUP_EXISTING48,root parentR13,supplement_root SUP)。旧core48与mechanism16 seal保持不动。
- 比较脚本 SUP/src/compare48.py：原16复用＋新32，5方法48×5view/5prefix；支持原CPUprefix及同actualCPU/commonphysicalLU公开minPhi规则（configs/comparison_rules.json已登记，非GT选根）。输出RUN_INDEX供后续物理诊断。运行目标NEW SUP/analysis/comparison48，尚未执行。
- r13_full_current完成 SUP/src/diagnostic_instrumentation.py：从原C2目录或ZIP创建NEW诊断runtime，仅复制src/runtime，插入copy-only hooks、nativehistory、最多20状态，固定10检查点＋每stage首末＋最多2个公开raw/primary迭代目标。参数/预算不变但记录费用改变CPU停止，所以不能回填主结果。实现说明 delegated/dynamics_audit/IMPLEMENTATION_HANDOFF.md。
- 比较评分后按ANALYSIS_CONTRACT.md选8个事后机制对象：A+减DR误差差值两端各2，共同困难2、共同好2（去重完整ID顺延）。用独立诊断root跑5方法（40份，单独费用），公有raw/proposal迭代传diagnostic_target_iterations便于验证重放一致。尚未部署/运行。
- r13_physics在写 SUP/src/analyze_physical_prefixes.py（等RUN_INDEX+父授权）：所有5×48 raw/corrected/primary+5prefix的原Phi/grad/box mapping，扣重计LU/RHS；不是新优化或独立物理验证。尚未运行。
- 数学诊断已在ANALYSIS_CONTRACT.md：u=J-mQ，v=Q-E-GJ，J-Jphys=A_m^-1(u+m v)，端口差S(...)；不重验证恒等式，用于实际缺陷/端口放大、A+与普通CIE材料梯度/真实物理方向分析。需要完整native状态后进一步写native诊断分析。
- 已知公开cohort差异：旧16有8sixsource、4half12、4full12及3crossgrid；补32有8sixsource、12half12、12full12及9crossgrid，不能把旧16均值直接推广全48。原材料家族/variant分配经过固定随机permutation，并非indexmod6；当前不读新GT。
- 要继续到：完整48结果/配对/失败动力学、公式解释、raw/corrector与cost、所有原始证据与图、中文报告和同GitHub/Pages发表并读回。不要以这条状态更新结束。

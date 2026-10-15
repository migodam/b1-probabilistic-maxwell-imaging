# 复现与数据布局

本仓库保持历史A/B/P2，G位于research/research_g_20261015。现有Gaussian/.venv_nn Python/NumPy/SciPy/Matplotlib，BLAS threads1，nice10，无安装GPU/新框架。完整源包、raw/corrected/independently-validated χ/current/field/data、observations、所有失败及cost随同仓库Release ZIP；解压到仓库根保持相对路径。

code/e0_trajectory_audit.py、e0_archive_repair.py与e0_history_signal.py只读旧保存资产；cache replay不能用于产生新物理速度/反事实搜索主张。原资产revision见SOURCE_MANIFEST。历史A30源和输入需要同级Gaussian/A30/a30-maxwell-material-manifold；不得因找不到资产偷偷新跑solver。

code/run_g_e1.py是唯一G付费入口，复用原A30与A corrector，必须使用冻结协议；已完成RESULTS存在时拒绝重跑。不要重置COST_LEDGER或覆盖旧协议来获得新额度。

code/offline_g_score.py必须在独立进程、online冻结后运行，只有它读已暴露OFFLINE_TARGETS；不由online模块导入。online数组与manifest只读。首次reader路径错误保留在results/E1_OFFLINE_path_bridge_failure_v1，后来的评分没有重算物理。

code/export_g_episodes.py对真实保存数组生成GSearchEpisode侧车和合法候选索引，不改A/T native schema。软件35项合同测试通过只证明预算/archive/接口语义，不是成像成功；check_assets只查路径/array key，不是新物理验证。

从头复现实验须新独立输出/授权，所有states/RHS/failed trials计费；普通报告/数组阅读不触发求解。不能将本阶段未用余额自动借给新3D或RL。

# 冻结边界和补充分析

- Stage1的初始gauge/切向/梯度审计通过后开始576单元开发。生产helper、ridge正规方程及A17共同anchor illumination derivative的补充审计随后完成。最终16项PASS不能倒填为全部在首次Stage1启动前完成。补充审计没有推翻生产梯度；误把current tangent-span projector称为illumination derivative的标签被明确纠正，两个检查分开保存。
- Stage2使用带LU/RHS预算计数的版本；Stage1执行的是45秒wall-only版本。旧10任务stage1_plan.json是未执行占位，不在64个实际任务的分母中。
- 新鲜Stage4的六方法、数据生成/评分合同和参数于生成前冻结。其主推理使用与Stage2相同的core；后加的Schur、coordinate-only和joint-current探针是独立CLI，只作用于已披露开发对象，不改变新鲜六方法。
- 数学证明、离线诊断、聚合脚本、图和网页可在冻结后整理；不得改变已封存的预测、参数或目标。附加方法如未执行，保留NOT_RUN；执行也不冒充本轮新鲜预注册方法。
- Stage2固定总预算是名义45秒分配。不可抢占的物理/Full-J/SVD及返回incumbent导致实际超支；K4最多47.78秒、K8最多50.24秒。实际LU/RHS未超过上限。不得写成严格等45秒墙钟胜出。
- `calibration_adjoint_RHS`是`adjoint_RHS`的子标签，不能重复相加。辅助current normal solve与原Maxwell LU/伴随也必须分列。

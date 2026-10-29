# 失败、修复与未成立条件

## 数学与算法失败

- 原E2固定weak幅度加固定材料半径使strong更新可行空间过窄。原结果保留；另开REV1解除幅度约束，B/C达到8/8有限资格，但没有20%费用收益。
- receiver-dark current的选择性lift不能满足多源共同材料；8/8记录MATERIAL_LIFT_FAIL，不作为合法候选，也不外推所有SOM失败。
- Fourier/Wavelet映回原66支持有显著投影损失；不是单纯同prior可逆换基。
- 新E5在绝对数据门下生成大量合法候选，但固定精细GT覆盖仅6/32。19/32的实际候选仿射容量下界仍超过epsilon，另7个有仿射容量但有限点未命中。GT最优仿射点可能非法，不能当算法。
- 固定云uniform与likelihood加权图像均未胜共同anchor。完整posterior不准不是对此唯一或必要的解释；方向任务失配与有限点放置已有直接证据。
- 旧relative竞争门和prior增量门会排除本批全部GT，不能当GT集合置信门。旧门结果不修改，新目标另行冻结绝对噪声门。
- 正确低维prior积分数值通过；任务包络均值误差与费用未胜普通强求积。严格外包未建立。

## 执行错误与保存

1. E5首次启动出现重复role日志参数，在任何physics/GT生成前失败。v1代码、冻结文件、日志与窄修复diff存于audit/possibility_runtime_fix1；新增RHS0。
2. 编排器解析venv symlink误选base Python，缺numpy，生成前失败；改用venv原始绝对路径。日志、0RHS与实际启动wall保留。
3. E4 underflow/事件报告修复、数值资格v1—v3与费用保留于其独立audit；以canonical v3为最终数值资格，不删除早期输出。
4. 文档/图库助手两次读取公开结果schema错误保留于reproducibility记录，未改变模型、预测或冻结评分。

E5实际物理job无solver action异常；bracket未找到目标半径18次/固定臂、18次/动态臂是明确预算内搜索失败，不等于数学不可行。随机臂0次bracket失败却有409/512未通过最终数据资格。所有拒绝与失败费用在完整账本中。

## 未执行不是成功

没有新NN训练、noise/layout扩展、W迁移、RTO/RML主要物理对照、Full66可信posterior、全域纤维覆盖或连续PDE认证。新32对象揭晓后没有根据GT扩大半径、改threshold、重选checkpoint或追加候选。

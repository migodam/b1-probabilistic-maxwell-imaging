# Research G 本轮迭代记录

1. **假设与失败诊断，E0，新增物理0。** 旧first-member archive有654个duplicate改善某个保留代表的full residual；raw新颖大量在corrector后返回旧邻域。区分archive质量丢失与真实结构发现不足。
2. **算法修改。** 只把first-member换成pairwise-safe、truth-free的低full-residual替换；success ledger不删，必须保持与其它active代表的分离。没有使用truth、孔数标签或改变相容阈值。
3. **同池回放，新增物理0。** 369次质量替换，原结构邻域不丢；没有新增物理世界。固定缓存回放不能模拟重走adaptive policy后的反事实轨迹。
4. **付费假设检验。** 预先提交355990a，三合格历史对象、四可用方法，各三proposal；first/quality两轮保持同seed、动作库、corrector、状态与RHS上限。各轮144新states，公共base另3states；两轮合计291states/6936复RHS。
5. **结果与理论修订。** 小pilot中质量替换未触发，两轮相同；一个近结构门槛的边界world重复命中，去重只有1。全部Best-of-K仍是base。archive质量修复成立不意味着搜索效率增加，单次bandit命中不支持历史学习机制。
6. **停止。** 不用余额补有利案例或反复调阈值。T新物理机制未作为G臂运行；发布时T已有公开交接提交，留给独立兼容接入和下一轮共同corrector比较。

离线reader第一次路径错误记录保留；只修reader后完成评分，没有重跑物理。两轮paired重放不是独立种子或独立对象复现。

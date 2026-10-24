# 第一64个L1：预GT探索性优化补充诊断

两项choices均先封存并封账；算法字节与DEV as-run相同，所有模型、45risk/26density主输出、原320case/callback flags不变。此处只重放Full-J梯度/ratio、KKT与目标进展，不作GT风险评分。

|lambda|新RHS|Full状态|接受更新|拒绝试探|Phi变化中位|目标变化中位|GN decrement原/末中位|
|---|---:|---:|---:|---:|---:|---:|---|
|0|15624|469|371|98|-1.1499183|-1.1499183|5.52331/4.89802|
|1|15128|679|303|376|-4.21001295e-06|-5.43398571e-07|0.00107859/3.83092e-05|

lambda1的decrement中位很小，说明多数初始点在GN曲率尺度上已近MAP；逐父对象的最大值/梯度仍必须保留，不能由中位宣告全部全局收敛。lambda0的Proper TR确能下降数据Phi，原unscaled SLSQP+callback-only返回停滞不能当强优化基线。所有TR仍受相同16状态/256RHS上限，停止不是全局最优证书。

实际新TR 30752/32768RHS，其中Full8=9184、acceptedJ32=21568。受控22.424秒，CPU21.036秒。cached初始J来自原96+64上下文，64对象实际context拥有240304RHS，不作为免费冷启动。best-query复用原80640RHS查询池，新增Full0。

best-query ordinary 64对象中64个Phi比原返回更低，mu256case中0个找到更低可行mu目标。它叫BEST_VALIDATED_QUERY_CANDIDATE，不叫优化器接受步或收敛解。原outputs和成本保留；没有mu缩放、额外迭代或参数搜索。

旧callback的Python触发条件是major iteration增加，不是下降/可行证书。Full-J TR的精确梯度、GN decrement与trust ratio另行逐对象保存；GN不是完整非线性Hessian，不据此宣称全局最优。任何后续GT附录只能读取封存choices，64个父对象为统计单位，不能把两lambda臂当128个独立对象。补充发生于方法冻结之后、GT之前，仅探索性，不追认confirmatory primary结论。

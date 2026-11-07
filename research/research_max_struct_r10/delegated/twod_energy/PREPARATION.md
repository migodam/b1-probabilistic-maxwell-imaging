# Fresh energy脚本已准备，未运行

入口为`experiments/analyze_twod_energy.py`，复用`src/energy_diagnostics.py::Geometry/physical_current`。只做语法编译；没有读取fresh GT或计算fresh结果。

程序须先通过与material scorer相同的root显式授权/168任务全终态门，再验证material scoring已COMPLETE和其公开seal仍一致；其后才打开GT。不能绕过“先完成评分，再energy诊断”的顺序。

预算与复用：

- 原24个GT若已有合法j则直接复用；否则每case最多一次既有576维全物理forward，新增上限为**24 LU、288 source RHS、0 adjoint RHS**，独立标为OFFLINE_SCORING，绝不进入inverse或在线选优。
- GT current建立带source大小/mtime和公开operator来源的本地offline缓存，防止后续失败恢复重复支付同一GT forward；缓存不作为inverse输入，授权门之前不读取。
- 168个final physical j从保存prediction读取，缺失明确记录，不加forward。辅助J只取最后实际joint快照；对应stage_physical必须与joint的m逐元素相同，才能计算physical误差、gap与有符号cross。不同m或缺j即保留缺口，不补正演。
- 固定公开G/S只构造一个Geometry实例，复右基、实材料控制、两序投影、数值交集/remainder与传播cross共同复用。
- `Geometry.affines`按`(rank,shape,mf)`缓存系数range SVD，跨case和七arm共享。oracle按`(case,rank,mf)`去重，JD按`(case,rank)`缓存，投影后的GT current也缓存，不为每arm重复计算oracle或SVD。

三个历史诊断tau .03/.1/.3全部输出，主G强/Twofold及paper图固定使用已冻结的.1；不按fresh恢复表现挑门槛。paper图优先全24对象汇总的绝对带能量＋百分比，并标有效n；同case actual FFT/APCS capacity与oracle floor另图。floor只约束辅助affine current，不约束材料导出的final physical current或GT恢复。

预定输出：逐状态能量与aux/physical/gap/cross、全tau单G与Twofold表、同case/rank/mf oracle容量表、按regime/method的pooled能量、明确缺口和独立offline-scoring成本。原variance/solver预算、候选、当前算法和FFT空间均不改。

根授权之前不得执行。命令形式与material scorer一致，并须显式提供实际postseal授权说明；默认输出`delegated/twod_energy/scored/`，paper图追加到`FIGURES/twod_comparison/`的新文件名，既有图不覆盖。

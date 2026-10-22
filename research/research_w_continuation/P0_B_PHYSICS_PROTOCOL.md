# P0-B 开始前冻结协议

2026-10-23，读取P0_Q_QUALIFICATION.json为PASS后冻结；首次搜索前写入。原物理、噪声、两个seed、宽/窄带、24任务probes、目标总材料质心/质量门槛全部不变。P0-Q与P0-B使用两个独立账本。

## 数据、未知量与共同初值

固定reference仅提供高精度cal均值和含probe8均值；算法不得读reference材料、mask、孔洞数、wall_values或mesh轮廓。高精度真实墙体网格仅资格/终点审核使用，不进入inverse模块。

逆变量为全域48×32像素两实场：contrast∈[0,8]和conductivity∈[0,.15]，跨频共享。cal合同在查询T实际为空气，网格包含T边界。精确空气允许拓扑生灭，无已知wallmask；最多占据面积3。所有方法采用同一top-area材料投影、同一正则化、同一trust限制与corrector。不将像素边界表达限制解释为连续材料全局不存在。最小特征受固定像素/平滑proposal先验影响，必须报告实际细小结构，而非假称完全任意拓扑。

初值仅由cal测量的background Born/adjoint与已声明先验给出。所有比较含该共同no-op初值。绝不以reference truth挑初值或接受步。

## 四方法和计算矩阵

四个context=两个seed×两套原带宽。每个context搜索预算800 RHS，其中共享当前cal full F/J和固定probe8任务F/J最多160 RHS；各方法最多160 RHS。

- FULL_J：普通全材料GN拟合提案，不加角度限制。
- RANDOM_MATCHED：同像素、同prior尺度的固定随机平滑材料提案。
- MULTISTART_FWI：第二个合法测量条件/先验初始化，接相同FWI正确器。
- TASK_AWARE：最大化当前完整任务切线/可见切线的正则广义比；只作提案。没有nuisance投影，主仪器已校准。

每臂使用相同Full-J trust-region corrector，最多4次标准半步line trials；在线以完整cal数据失配+同一弱Gaussian/ridge先验接受，绝不使用任务真值。最大归一化材料步幅0.2；corrector通过全材料Jacobian构造，先验ridge0.02，null/task比ridge0.01。坐标按contrast/8和sigma/.15归一化。无结果驱动超参数调整。全部拒绝/失败计费；当预算不足，余行写NOT_RUN。

只有cal的数值距离<=1且终点精度可判定时才允许GO-MECHANISM。对最优合法端点先支付固定probe8任务响应；这只是全24任务类sup的下界，不冒充已评估全部任务。值得验证的近等价端点再评估全部24任务及独立细化。若没有cal相容端点，保存best-found目标响应和cal残差，不自动用更多搜索额度救结果。

cal合格且任务响应显著不同的端点必须做未知目标joint补偿；用T内总材料质心差>=.1λ、质量±20%判定，仪器先固定。校准可分解不当作joint证据。若没有cal候选，joint为NOT_RUN_NO_CALIBRATION_WITNESS。

## 费用和数值

硬上限8000 complex RHS、2000材料频率LU、7200 process CPU秒；search/endpoint/validation/repair分别40/30/20/10%。搜索辅助导数实际预约12接收转置RHS/频率，多源前向4RHS/频率；纯矩阵小解CPU计费但不冒充Maxwell LU。

搜索采用独立高阶FEM全网格材料接口，体积分接收输出；这一新逆表示须先做真实pixel-J差分检查。最终GO端点必须重新高精度审核；不将原两墙qualification自动迁移至每个新候选。数值误差不能靠调大sigma掩盖。

计算结果是受限表示/预算的见证发现实验；没有找到只记NOT_FOUND_WITHIN_BUDGET，不证明唯一性。P0-C/P1/P2/P3本轮不执行。

## 运行后实施审阅（不修改以上冻结合同）

最终审阅发现原配置`materials.minimum_feature=0.1`未被48×32像素、Gaussian平滑提案和top-area投影严格执行。像素边长0.08333/0.09375，可能产生小于限制的孤立结构。box、passivity、面积和跨频律的检查不能替代此几何先验条件。记为PROTOCOL_IMPLEMENTATION_GAP，所有当前候选不具备正向见证资格；不事后改变原材料类。全部cal距离本已大于1，故没有任何GO被撤换或改门槛，也不能据失败排除原合法材料类中的物理多解。新一轮应先修正这一实施缺口，而非利用剩余验证费用重做搜索。

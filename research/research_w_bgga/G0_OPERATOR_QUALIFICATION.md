# G0：算子实现已就绪，Maxwell资格尚未执行

## 已经实际运行

固定公共PSLG网格：wall ROI为[-1.5,1.5]×[.5,1.25]，target积分ROI为[-1.2,1.2]×[1.4,2.8]。网格key不含未知目标位置、半径或材料。P3 scalar H1场；复128/实64。原被动材料、PML、单位线源与白化保持。

圆盘按384边形与每个固定三角形相交，积分P3质量及解析incident forcing。不是把目标边界刻入网格，也不是把随机目标换成真值mask。几何测试包括目标域边缘位置、不同半径、常量/线性/四次incident；四次例明确能分辨直接积分与错误P3节点插值。

测得网格28,903 DOF/28,420 free DOF/6,369 triangles；公共ports 2,257；示例目标79个结构活跃自由度。目标多边形面积0.04523691560091064，分割积分面积差约−6.94e−18。该面积一致性不认证电磁场误差。

散射场增量实现保留`b_X`。小复矩阵测试核对普通full、resolvent、独立外域Schur、共享port Schur、背景预条件GMRES、同空间goal-Galerkin；删掉b_X的负控制确实不一致。这是合成代数验证，**不是新Maxwell实验**。

## 尚未执行

`tools/g0_qualify.py`在任何新物理预约前检查AUTHORIZATION。已冻结12个几何/contrast资格条件；每个先做公共固定mesh完整solve，再用原interface-fitted P3独立核验。第一条白化差超过.1即停止后续cache主张，至多另登记一次限定修复；不会因失败提高噪声。低损耗高contrast是事前stress，不预称near resonance。

四类各3个信息/传播stress在`G0_CASES.json`列出公开构造，类别是待检验假设。正式Jacobian/数值秩、absolute Schur谱、joint covariance、有限端点和全部task传播尚未执行；该完整stress驱动尚未交付，不把资格脚本冒充完成整个G0。

## 成本障碍

完整目标ROI的Green images估计1,026,303,040字节，核心81,504,784字节。准备至少2,257 full RHS；每个全目标仅2源RHS。普通Schur可复用同一已付费port。所有接口都有相同数值门，不给BGGA更松误差。

端口矩阵/低秩core的数值条件、真实forward误差、独立网格、伴随及跨目标全波一致性均为NOT_RUN。不得把这里的同离散数学恒等式或小的积分误差写成continuum certificate。

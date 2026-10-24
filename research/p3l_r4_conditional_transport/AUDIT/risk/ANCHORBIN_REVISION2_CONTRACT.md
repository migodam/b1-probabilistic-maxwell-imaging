# 最终相矩修订2：在线anchor质量倾斜初始化

新增risk_anchorbin.py、train_risk_anchorbin.py。没有训练或物理启动；旧训练器和模型文件不改。evaluate_r4_risk.py仅追加新format加载路由：loader=anchorbin可显式指定；原r4入口也通过format-aware dispatcher识别新format，旧检查点仍调用原RiskActionPredictor。

目标只检验初始化/优化差异：标签继续是同一joint-parent的精确等h质量箱联合δ事件标签，loss、AdamW(lr5e−4,weight_decay1e−3)、batch64、有放回父样本抽样、梯度clip5、验证间隔100更新、update0可选及max3000更新均保持。三个固定seed为2026102501/2/3；K8/16仍6114/6386参数，不增加分支/自由参数，不独立化属性和相位的条件相关性。

raw native在线q0给出δ0，raw complexχ0给出两个分量的相比例。两者可识别时Re/Im比例须在2e−5浮点容差内一致且落在原[0,1]相支持内；否则报错。共同比例用δ加权最小二乘形式计算以减少小虚部舍入误差。绝不从GT、标签x或标准化image通道恢复属性/相位。显式q0/χ0缺失即拒绝。

基准箱矩Mj0=δ0·clip(Kh0−j,0,1)。对每个合法最大分量B，构造p0=1−M0/B，pj=(M(j−1)−Mj)/B，pK=M(K−1)/B；它们非负、和1。冻结eps=1e−4，pbase=(1−eps)p+eps/(K+1)，输出质量为softmax(log(pbase)+CNN_correction)。CNN最后1x1头仍零初始化，故update0相矩及材料图接近原anchor。完整K+1质量倾斜保持与旧单调头相同的自由度，输出继续单调、有界、按1/K积分，是材料风险动作而非posterior样本。

精确算术下update0材料图=(1−eps)χ0+eps(χhost+B/2)。因此初始化并非逐位等于anchor，有可声明的eps平滑偏差。δ→0时未平滑联合矩→0；固定eps使最终初始化仍保留eps·B/2的小宿主偏移，不能声称零增量时输出严格为host。原P1属性下界为正，零δ只是数学退化测试。没有靠改eps、学习率或额外更新救结果。

命令：`PYTHONPATH=src <python> -m p3r4.train_risk_anchorbin --train JOINT_TRAIN --validation JOINT_VALIDATION --anchor-train RAW_CONTEXT_TRAIN --anchor-validation RAW_CONTEXT_VALIDATION --kind phase8 --seed 2026102501 --updates 3000 --out NEW_RUN --device cpu`。支持相同spectral/decoder_only/random context模式。base-mass始终由未归一化raw q0/χ0独立传入，几何输入置零不删除合法anchor基准。

新checkpoint明确format=p3r4_anchorbin_risk_v2、anchor_epsilon、moment_rule、raw online contract和最终相矩family停止规则。独立AnchorBinPredictor输出原全576complex raw图；既有独立评价器按原矩形另存clipped，旧format加载路径保持。父调度器决定实际训练与最终科学解释；若validation无增益，本次是最终实质phase revision，应停止家族而非继续换规则。

ANCHORBIN_REVISION2_CHECK.json：32个源MaterialPrior decoder anchor、K8/16的质量和误差≤3.58e−7；update0相对声明平滑公式误差≤5.96e−7(按float32 eps×B量级核对)，相对原anchor的全网格RMS约4.29e−5；joint label差≤4.89e−6，包含float32 A/比例计算差异。单调、有界、有限梯度、零/近零δ极限、raw参数必需、Re/Im不一致拒绝及独立predictor精确重放均通过。实际训练/validation仍NOT_RUN，无新RHS/状态。

# 固定Spline模型矩作CNN风险上下文

新增risk_fusion.py、train_risk_fusion.py、tools/r4_fusion_moments.py及evaluate_r4_risk的独立loader=fusion路由。未启动特征生成、训练或物理。目标是检验较粗的模型分布矩能否提供CNN上下文；不是把风险较差的Spline样本均值称作posterior estimator或校准结果。

特征生成：`PYTHONPATH=src <python> tools/r4_fusion_moments.py --checkpoint FROZEN_RAWY_SPLINE --observations ONLINE_DIRECTORY_OR_PACKED_Y --parents PARENTS_JSON --out NEW_MOMENTS_DIRECTORY --device cpu`。manifest为parents固定列表。只接受rawy512/full66 Spline，K32固定，每parent seed由固定sample-seed20261025121与parent身份生成，sampleindex0…31不挑选。每次网络最多8draws，逐parent经原decoder计算全24×24的meanRe/meanIm/varianceRe/varianceIm，保存moment_image、原生成x_samples、K/seed/checkpoint身份及MOMENTS_COMPLETE gate。packed joint分支只读parent/y，绝不打开True x/material/s标签。固定模型矩是模型pushforward的经验矩，不是真后验矩。

训练：`PYTHONPATH=src <python> -m p3r4.train_risk_fusion --train JOINT_TRAIN --validation JOINT_VALIDATION --anchor-train STRONG_PACKED_TRAIN --anchor-validation STRONG_PACKED_VALIDATION --moments-train FROZEN_TRAIN_MOMENTS --moments-validation FROZEN_VAL_MOMENTS --control moments --seed 2026102501 --updates 4000 --out NEW_RUN --device cpu`。最终1024/256强context由父调度器提供；不存在时不生成或替代context。控制为moments/permuted，各3seed2026102501/2/3，共6fits，max4000更新。TRAIN-only拟合新增4通道各自的均值/尺度；validation不参与。唯一GT读取为材料MSE监督，不读GT latent。

原8个在线anchor/geometry属性通道加4个模型矩通道，卷积12→16→16→16、1x1双分量zero-head，6418参数，输出是anchor+材料残差。permuted控制仅在TRAIN标准化后对同parent的4个新增通道用同一个固定空间置乱，seed20261025122+parent，保留通道幅度与同parent信息；原8通道不动，不借外来parent。初始输出逐位等于float32 anchor。AdamW lr5e−4、wd1e−3、batch64、clip5，验证含update0且每100步；保存raw/clipped risk动作。

FusionPredictor必须显式接收moment数组、parent身份和原允许的online anchor，不能由True x或GT构造特征。评价spec使用loader=fusion、checkpoint及moments_directory；冻结前验证对应每parent封存K32矩及同一个Spline身份，线上读取moment_image并通过原允许anchor接口预测。缺失矩走明确error/fallback协议，不偷偷换别的样本。旧field/anchorbin加载路径保持。

FUSION_MATH_CHECK.json：6418参数、zero-head精确anchor、置乱幅度逐parent保留且确定、有限梯度、独立加载器精确重放、缺失moment拒绝，均PASS。新物理/训练更新为0。没有实际validation/development增益证据；若development无增益则停止该融合机制，不按locked结果调整。

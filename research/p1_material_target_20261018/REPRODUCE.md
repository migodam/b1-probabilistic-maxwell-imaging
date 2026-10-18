# 代码、数组与复现

核心模型与protocol在1674a7f、测量语义在4f71c0a826f9db6217e5081bd56300716152db24固定，均早于材料输出。执行脚本/线上隔离接口及后续报告由本阶段最终commit保存。原数据生成seeds和任务尺度没有按结果修改。

src/p1_material/prior.py为独立先验及解析Dg；proposals.py为prior可逆提案和candidate-only retraction；target.py为proper likelihood及受控A30接口。无需NN、GPU或新框架。沿用已有NumPy/SciPy/Matplotlib环境Gaussian/.venv_nn。

已经执行：`PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python tools/prepare_p1_material_target.py`。原目录MATERIAL_STAGE_CLOSED阻止覆盖本轮样本；分析复现应复制配置/结果到新目录。该脚本不导入物理backend。

复现材料样本时，先将Release解压到仓库根目录。以下独立调用直接使用冻结prior/seed重建128份calibration材料，输出在新的复现目录，不触碰已封账campaign，不调用物理：

```python
import json
from pathlib import Path
import numpy as np
from p1_material.prior import MaterialPrior
p = Path('research/p1_material_target_20261018')
cfg = json.loads((p / 'configs/PROTOCOL.json').read_text())
points = json.loads((p / 'inputs/geometry.json').read_text())['points']
prior = MaterialPrior(cfg, points)
q = prior.draw(np.random.default_rng(cfg['calibration_seed']), 128)
chi = np.array([prior.decode(a)[0] for a in q])
out = Path('publication_local/p1_reproduction')
out.mkdir(parents=True, exist_ok=False)
np.savez_compressed(out / 'calibration.npz', latent=q, chi=chi)
```

用`PYTHONPATH=src python`运行；development和material-check分别使用protocol中对应seed和32个样本。只重绘已保存数组可运行`python tools/publish_p1_material.py`；该入口不运行材料优化或Maxwell。它重新生成网页/图/本地证据包，不上传、不增加物理额度。

Release包含calibration/dev/material_check三个192样本数组、noise draws、24个材料运动数组、128步prior-only proposal、原已付费背景以及所有记录。offline_materials中q与材料只属于数据生成/离线分析；线上反演只读configs/ONLINE_TARGET和inputs/online。当前online observation目录为空，不能从空目录推定已做过成像。

旧F3只读取sigma metadata与已支付背景，没有重算旧P0、读取旧locked truth、执行SHA256或额外smoke。首次记录字段冲突的失败保留；恢复复用原材料/noise数组，未重新抽取有利parent。

目标、测度和原始定义见MODEL_AND_PRIOR；真实阶段命令和费用见PHYSICS_NEXT_RUN。P1的材料law可以独立复用，但若改prior/noise/支持或采集，要新建target/protocol，不混合概率和实验结果。

未来真实阶段另需现有A30后端，固定版本见SOURCE_INVENTORY；跨机器通过`B1_P1_A30_ROOT=/path/to/a30-maxwell-material-manifold`指定其根目录。包含复用的已付费background不等于从零部署没有setup成本。该外部backend未在本轮重新开发或重新支付。

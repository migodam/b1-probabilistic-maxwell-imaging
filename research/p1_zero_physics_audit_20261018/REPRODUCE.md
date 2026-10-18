# 复现与信息隔离

源码：tools/p1_zero_online.py、p1_zero_offline.py和publish_p1_zero_audit.py。依赖NumPy/SciPy/Matplotlib，复用现有CPU环境，不安装新框架。

先取得本轮证据Release与原[P1真实物理证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p1-physics-qualification-20261018)，解压到同一B1 checkout根目录。本轮原始数组包含材料metric、方向与decoder试验、隔离offline pair；原物理Release提供被复用的24个receiver J和26个已付费state。没有重新生成物理或暗中把它算免费。

只重绘本轮已有结果：`PYTHONPATH=src python tools/publish_p1_zero_audit.py`。它会重建图、网页与本地证据zip，不上传、不调用Maxwell。在线和离线分析原目录已封账，直接重跑会拒绝覆盖。

若要独立复算全部零物理分析，在仓库根目录使用`PYTHONPATH=src:tools python`执行：

```python
from pathlib import Path
import shutil
import p1_zero_online as online
import p1_zero_offline as offline
replay = Path('publication_local/p1_zero_audit_replay').resolve()
replay.mkdir(parents=True, exist_ok=False)
shutil.copytree('research/p1_zero_physics_audit_20261018/configs', replay / 'configs')
online.OUT = replay
online.main()
offline.OUT = replay
offline.main()
```

这只复算cached linear algebra、decoder与已有prediction，创建独立的零物理分析记录；它不删除旧账本、不刷新任何物理额度。online先封账，再运行offline。online没有访问GT端点、offline_materials、旧LU、backend或新solver的入口。输入读取清单见ONLINE_INPUT_MANIFEST，GT的单独读取位置在offline源码。

本轮parent A/B和原数据均已暴露，不能称独立blind validation。所列新物理实验只有冻结协议，没有获授权的execution overlay或已运行结果。本轮不提供自动支付的future runner，也不修改旧physics ledger。未来若获授权，应另建独立stage/ledger，并保留本轮和原26/976证据。

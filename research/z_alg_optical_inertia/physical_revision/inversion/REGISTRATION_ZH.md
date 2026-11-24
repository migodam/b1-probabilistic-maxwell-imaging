# 本轮新增数据与调用入口

| 数据族 | 独立材料数 | 观测 ID | 使用边界 |
|---|---:|---|---|
| 原始 optfresh 全环复算 | 6 | `optfresh_{01..06}__pw360_10x20__pair1pct`、`__pair2pct`、`__n10` | 96²，400 MHz，2 m 方域，10 发 20 收。原 2% 数据和独立 10% 复用，1% 为此前冻结的配对缩放。这六件材料已有完整暴露史 |
| 电小片过渡算例 `compact_tm` | 1 | `compact_tm_g045`、`compact_tm_g060`、`compact_tm_g075`、`compact_tm_g095` | 12²，400 MHz，0.125 m 方域，10 发 20 收。Codex 按本轮设计新定义材料及测量。训练与留出阵列的角度相差半步。四档复场是预定配对噪声而非四件物体 |

原始观测和算子由数据加载函数从 `synthetic/public` 和 `submission/curve/observations` 读取。电小片公开数据只有算子、复场和半径，放在 `compact_synthetic/public`。生成材料、留出端口与生成回执在 `compact_synthetic/offline_truth`。逐级放宽回放在主比较评分后开始，使用相同四份公开数据。 [原始协议](PROTOCOL.md) · [电小片布置](SYNTHETIC_DESIGN.md) · [曝光后回放](COMPACT_FOLLOWUP.md)

共享全波和伴随求解与前期最小峰值研究采用同一 96² FFT 入口。电小片用相同标量 TM 圆盘基函数的稠密矩阵，以及点解区间残差反校验。原生更新空间的区别和实际费用在逐例回执中可查。 [调用命令与数据](RESULTS_EN.md)

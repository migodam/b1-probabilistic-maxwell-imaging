# 响应校准交付索引

主报告：`RESEARCH_REPORT_ZH.md`。18行：`OBSERVATION_SUMMARY.csv/json`。

- 解析19组：`math/A1_MATH_RESULTS.json`，7 checks PASS，另3个积分 unit tests PASS。
- 三空间：S0=(0,1)，S1=(0,64)，S2=(12,27)；6 independent parents ×3 paired noise。posterior 是显式截断2D新模型。
- S1 regular P1：V3真实1强+1弱坐标、det/prior/normal residual；σ/4 cross mean=.0003388 SD、Z=.0006413。
- S1 host P0：rank0 prior mass=.22086048；直接 paid-node readback posterior host=.92652/.94128/.95314。前两档 mean convergence仍不达1%SD。
- S2 normal P0：rank2 chart3→5；σ/4 cross mean=1.59e−6 SD、Z=4.82e−6。rank2没有弱坐标。
- 归一化连续 Gaussian 的数值 TV/KL 见 `v3_strata/*GAUSSIAN_BASELINES.json`、`remaining_spaces/*RESULT.json`。没有 chart TV=0声明。
- 全域覆盖、漏质量、分支 certificate、全局 Δ、同精度成本冠军均未建立，扩大 HOLD。

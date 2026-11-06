# v04 追加原物理预算诊断

分类：POSTHOC_BUDGET_DIAGNOSTIC，已曝光development，不能充作独立确认。三臂从各自原physical_accepted80公开材料出发，一致重置L-BFGS历史；原18预测和canonical不变，未重新造Y/噪声或重跑aux。工人未读GT。

目标保持Phi3=mean_s ||S J_s−Y_s||²/||Y_s||²，无1/2、R0、χbox[0,8]；同一G/E/S、60源。它不是全Y归一化，也不是iid同σ白化ML。

|arm|m80 Phi3|追加Phi3|PG∞|追加accepted/requests|global相对数据残差|
|---|---:|---:|---:|---:|---:|
|sourceCIE|0.00098868388|0.000953504147|6.88e-06|80/81|0.0302585|
|GlobalDual|0.000987525063|0.000955191814|5.44e-06|80/81|0.0302882|
|AplusGlobal|0.000980661414|0.000951046733|9.09e-06|80/81|0.0302182|

三臂都以新增80步预算上限终止，PG仍约5.4–9.1e−6，不能称收敛/驻点或唯一性。公开Phi继续下降只说明原80步不是最终物理最小值，不能据此推断材料恢复；材料评分由root在三预测共同seal之后独立做。选择追加v04的依据是公开PG/fit，不是GT排名。

新增已测CPU下界：Mac 1.359581s，XINAN 502.750000s。一整GPU科学子进程含coldsetup/三臂试探/保存，计243LU+14580正向/14580伴随RHS；adjoint复用每个材料的原LU。query子CPU和transport另计，nested三臂/设置/stream-event不重加；原18prefix为非additive引用。手工实现与未单计progress/query等仍null。

图：FIGURES/vector/V04_CONTINUATION_PUBLIC_PHI.png；费用：RESULTS/vector_continuation_v04/ACCOUNTING_FINAL.json；公开表：同目录PUBLIC_ENDPOINTS.csv。

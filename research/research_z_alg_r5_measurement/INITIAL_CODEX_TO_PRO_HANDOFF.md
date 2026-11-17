# R5 → GPT Pro：严格理论/原创性复核交接

2026-11-18。首先读 RESEARCH_REPORT_ZH.md、EXPERIMENT_PROTOCOL.md、THEORY_DERIVATIONS.md，再查原始 results。不要把本交接当作扩大实验的授权。

## 当前决定

NO-GO for Gate2/complex host/unknown support/3D/PEC/576 expansion. Gate1六新目标全部运行，未替换；没有全材料box或噪声可行集覆盖。论文主张未成立。不是数学不可能性证明；保留有用局部证书与反例。

## 最重要的新证据

- 历史双根对偶预测与完整VIE误差1.29e−15；旧源组合gap3.47e−18；可行外部线源gap最大4.42e−3。这是历史调试，不是新盲测。
- T0两合法实根在clean/noisy-central方程均有严格局部Krawczyk框。新源local-box gap下界0.0040682/0.0056466。NB=0的旧源不可消歧对整个严格数据纤维成立。
- 同一T0的reciprocal64找全相同两根并复制主动源恢复。没有商代数独占成功集合。
- T3无噪Q初始化确有数值优势，但根认证失败、带噪零合法点根；不能把它作为稳健测量优势。
- 带噪+固定0.1%半径失配：一源max-min、固定源、Fisher、reciprocal-active各3/6 RMS≤0.1；直接补基+联合TRF64为4/6。
- 预定模型误差2e−5被4/6目标超过；所有max-min点根硬筛都丢空。局部数值tube保留不等于真材料准确：T5保留一候选但RMS0.204。
- T4的原始补基readout误差16.06，联合校正0.05793；强基线本身也有失败，必须保留。

## 请逐项审查，而非再次建议扫参

1. 在当前rank-one独立cell、单接收、k=2s结构下，能否给整个明确材料box的有限解覆盖条件？必须处理节点、pole、投影无穷远、重根、系数误差。D=20浮点端点齐全不是证明。
2. 对含旧观测误差的可行材料集合，如何做候选邻域补集排除？w0=tB†必须随旧数据移动。不能只把中央根扩大成未经验证的椭球。
3. 除标准线性代数/有限超平面分离外，剩下哪些真正非平凡、可发表且可证明的Maxwell结论？不得重新命名Schur、Woodbury、Fisher或Gram匹配。
4. 若ordinary reciprocal/homotopy同样枚举根，Q的价值是否只剩较低小模型发现成本？这与物理测量优势不同，应分别归因。
5. 本轮理想线源辐射Gram是否足够定义研究场景？真实loaded antenna power/hardware mismatch没有测试，不应外推。
6. 若没有新的全域定理，是否应只保留认证歧义/不可辨识性的技术记录，而停止主动消歧方法论文？请给明确判断，不自动建议继续。

## 审计边界

原来源/目标/预算冻结不变。全部GT已在封存后评分，目前附件包含GT便于独立复核，后续重放必须叫exposed replay。局部证书针对binary64 G/E/S/y；正式安全框位于qualification.json，其K严格包含检查处理了报告端点舍入。没有Hankel积分/连续PDE误差证书。

成熟homotopy因工具不可用NOT_RUN；不能声称已战胜它。引用核查涉及metadata/abstract级证据和HTTP403/429记录，详见inputs。Gate2及D924 probe按停止条件NOT_RUN，不属于遗漏待自动执行。

## 建议优先看的具体反例

T0/noisy/Q.json与reciprocal64.json、两active_trace.json：独占性失败。
T2/noisy：Q找到中央实根距GT RMS7.33，新增一源联合后仍1.91；残差合法不意味着材料正确。
T3/noisy/direct.json对maxmin_trace.json：补基校正0.0869对一源12.53。
T4/noisy/direct.json：满秩但病态，raw失效。
T5/noisy/maxmin_trace.json：局部相容却不准确。

最终科学裁决由Codex主线程作出；worker只负责提取、普通代码与独立审查。此轮未发布远程仓库/网站，未修改冻结历史结果。

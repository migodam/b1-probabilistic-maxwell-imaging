# B1 最终研究裁决

**GO_LIMITED_TOPK_PROTOTYPE**：已实现从散射数据到3D vector-Mie模板材料事件概率、Top-K代表图、第三分支余量与不确定性。竞争主例50.64%单空腔、40.43%双空腔、8.93%实心。它是有明确proper prior的条件材料posterior，不是多图美化。

**HOLD_GENERAL_GEOMETRY_AND_CERTIFICATION**：几何固定、模板支持有限，位置没有被推断。原生全网格4个reference失败；第二轮经验插值核验通过，但无全域余项证书。24parent的校准只是开发检查，部分分箱存在偏差、CDF细化未认证，无posterior-rank SBC或仪器校准。

**选择积分，暂停当前importance扩张**：相同720数据上，名义cold75的Gauss5通过720/720完整参考门槛；来自单竞争y0的冻结importance27+48只通过124/720。这个结论限定在冻结proposal/已暴露低维panel，不是一般高维方法的NO-GO。

**NO_GO_STRONG_TAP_CLAIM_AT_THIS_STAGE**：没有新Maxwell-specific认证余项/可变量几何理论或同概率精度部署加速。A17只作generic affine签名诊断，没有native交换，rank2/材料2维不压缩。现有组合来自已知Mie/Bayes/importance/插值工具。大型NN、beam confidence当posterior、无密度修正的频率、以residual赢者替代evidence均停止。

**预算收口**：新物理384/384，完整历史不删除，原包二维74另外保留，合计458解析材料–频率评估。GPU/NN/locked0。后续新物理或扩展support需要用户另行授权；本阶段代码/结果/报告/失败/成本完成本地交付。下一阶段优先增加一个可积分的几何变量或proper其他分支，重建reference，再考虑likelihood-informed/tempered/RJ方案。

详细PASS/FAIL/HOLD/NOT_RUN与证据见[GATE_DECISION.json](GATE_DECISION.json)、[中文报告](RESEARCH_REPORT_ZH.md)、[独立native review](research/delegated/probability_review/REVIEW_NATIVE.md)。公开发布未执行。HTML现场渲染被本地file访问策略拒绝，记NOT_RUN；8张静态科学图已视觉核对。

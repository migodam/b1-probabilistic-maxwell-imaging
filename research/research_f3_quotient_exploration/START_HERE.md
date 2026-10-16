# Research F3：商坐标、有限噪声探索与P0停止裁决

**已运行并封账。原P0的局部Explorer研究到此收束：商先验正确、两例新驻点已获得、weak背景发现任务尺度竞争端点；未发现新拓扑或高阶反馈增量，不训练NN。**

[中文报告](RESEARCH_REPORT_ZH.md) · [网页版](report.html) · [实际材料与路径图](gallery.html) · [全部数值](METRICS.json) · [账本](RESOURCE_LEDGER.jsonl) · [原始数组Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-f3-quotient-20261017)

- 独立推导并核验非均匀商先验；40万次原权重pushforward抽样符合解析CDF。
- 仅4个新state使两个原F1端点通过新商坐标驻点标准；旧Gate2仍FAIL。
- weak的普通native候选RMS=.00527556、ΔΦ=1.14060；ring无合格竞争端点。D高阶方法没有优于普通native。
- 3个数值不同的合格端点彼此只差RMS .000639–.001100，无拓扑变化，不能称3个posterior modes。
- 追加9个既有端点的中间点：均在原noise球，但全部ΔΦ>2；未建立连续likelihood-competitive路径。
- 合计55/64 states、712/2048 complex RHS、17.641秒受控wall；55个端点立即保存。GPU/Mie/3D/NN均0。

下一步选择：另立有物理理由的P1及可核验ambiguity目标；本轮不启动P1、不再追加P0优化。完整五项决策见GO_NO_GO_DECISION。

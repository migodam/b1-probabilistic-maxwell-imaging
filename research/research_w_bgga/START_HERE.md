# Research W-BGGA：实现与零PDE preflight

**当前状态：代码/合同/几何benchmark完成，新全波G0等待独立资源授权。没有宣布G0/G1/G2通过。**

- [中文完整报告](RESEARCH_W_BGGA_REPORT_ZH.md)
- [实际公共网格与端口图](figures/PUBLIC_MESH_AND_PORTS.png)
- [待授权的G0方案与费用预测](RESOURCE_REQUEST.md)
- [可复现命令](REPRODUCE.md)
- [分层裁决](SCIENTIFIC_DECISION.md)、[机器可读状态](GATE_DECISION.json)

实际完成：32页附件提取/核心公式核读，指定两个Git版本来源审计，2950原子及原始数组定位，固定公共目标网格、cut-polygon P3装配、显式材料相关RHS的增量式、普通full/Schur/resolvent/Krylov/goal-ROM接口，以及有限prior likelihood区间、条件lift和唯一factor版本逻辑。

实际网格28,903 DOF、公共端口2,257，示例目标触及79个DOF。公共Green images约1.03GB，准备需2,257 RHS；这不是实测全波加速。`results/jobs/`保留全部成功和失败记录。新FE=0、新物理RHS=0。测试数量/最终结果见 [代码检查](CODE_CHECK.json) 与最终测试receipt。

旧W明确封账、旧W-SLAM资源分scope；用户本轮第0/9节明确预算仅为建议，不能冒充既有授权。因此未借旧E4余额运行新的PDE。申请只限G0最多4,000 RHS、500 FE、2 CPU小时，不捆绑G1/G2或新锁定测试。

网页版发布不等于科学通过。旧 oracle joint witness、有限posterior成功和一般A0失败均保留原身份。

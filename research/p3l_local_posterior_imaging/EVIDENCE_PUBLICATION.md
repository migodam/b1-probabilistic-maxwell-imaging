# P3-L：原始证据发布回执

**已公开，科学结论仍为 REDIRECT。** 本回执记录科学封账后的上传，不修改封账账本、原始数组或研究裁决。

- [完整中文报告与图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/p3l_local_posterior_imaging/)
- [本轮证据 Release](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3l-local-posterior-20261022)
- [轻量关键结果包，约17MB](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/p3l-local-posterior-20261022/P3L_key_results.zip)：128对象风险摘要、协议、预先固定的前4个对象材料/图像及必要数组。
- [机器可读公开读取回执](WEB_READBACK.json)

完整现存原始证据为61卷，压缩后23,969,420,103字节，逻辑文件28,313,674,687字节，共1,767,734个文件。链、重复拒绝、端点、模型权重、失败记录和完整追加账本均保留。早期内存中断造成的3个端点缺口仍按原记录披露，没有重跑补齐或重置费用。

## 版本与科学含义

证据来源固定于 `08f56796a0b786a839b4880992244d1a87460d33`。`CAMPAIGN_CLOSURE.json` 是上传前的科学关闭记录；本文件与WEB_READBACK是后续发布回执，版本时序不同不改变科学结论。原始材料图和近似均值不能视为已校准posterior；11组四链Full参考全部未合格。

## 下载与恢复

普通阅读只需网页和关键结果包。完整复核时，在Release下载所有 `P3L_evidence_part*.tar.gz`、`P3L_EVIDENCE_MANIFEST.json`、`P3L_EVIDENCE_INVENTORY.jsonl.gz` 及 `package_p3l_evidence.py`。把清单命名为MANIFEST.json、解压目录表为INVENTORY.jsonl，再恢复到一个不存在的新目录：

```sh
cp P3L_EVIDENCE_MANIFEST.json MANIFEST.json
gzip -dc P3L_EVIDENCE_INVENTORY.jsonl.gz > INVENTORY.jsonl
python3 package_p3l_evidence.py --bundle-dir /absolute/downloaded-bundle --restore /absolute/new-restored-directory
```

恢复按路径、文件大小和分片偏移重组，不运行Maxwell、不覆盖原实验。物理依赖及原P1 geometry见[复现说明](REPRODUCTION.md)。读取证据不消耗新物理预算；已封账余额不能作为下一轮授权。

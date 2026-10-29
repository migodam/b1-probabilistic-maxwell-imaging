# R3 GitHub发布记录

已沿用既有 `migodam/b1-probabilistic-maxwell-imaging` 仓库及其GitHub Pages，没有创建新仓库。报告、源码、原始输出、图像、成本与失败证据位于 `research/research_z_alg_r3/`。同期R10的新提交与导航保留，未覆盖冻结R2/R9研究。

- [中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r3/)
- [latent与characteristic subspace推导](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r3/theory/LATENT_AND_CHARACTERISTIC_SUBSPACES.html)
- [匿名原始Markdown](https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_z_alg_r3/RESEARCH_REPORT_ZH.md)
- [源码与全部证据](https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_z_alg_r3)
- [Release阅读包、完整复现包及清单](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r3-20261030)

科学内容提交 `9d19da42564adf9d6b22f81a781c89c70d3e959c`。Pages部署运行40933218273成功。19项源码/Release资源和8项网页/图像通过匿名完整GET及逐字节比较；范围见 `audit/PUBLIC_RAW_AND_RELEASE_READBACK.json` 与 `audit/PUBLIC_SITE_READBACK.json`。不是对全仓库逐文件验证，也不构成科学主张的额外证明。

完整ZIP 8,846,094字节，257成员；阅读ZIP 2,011,542字节，82成员。完整包含最小不可变R2物理依赖、冻结输出、公开观测和评分后仅含m的GT；不包含私有生成seed。后续发布回读与封账收据另外归档，不递归改写已冻结压缩包。

科学结论仍为有限实验PASS_WITH_LIMITS，首创性及一般优势HOLD；latent扩展只完成推导，未实现Full66。最终新鲜四例的3/4与2/4以事前数据门衡量，top1由原完整物理残差选择，不使用best@K替换。所有限制见中文报告与独立科学审查。

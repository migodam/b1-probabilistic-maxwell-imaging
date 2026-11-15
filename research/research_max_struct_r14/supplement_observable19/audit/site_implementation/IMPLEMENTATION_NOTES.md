# 补充研究静态页面实现

仅新增本补充的 src/build_site.py 与本审计目录。借用父 R14 的 Renderer、CSS、MathJax 配置、已核验 DOWNLOADS 筛选和链接审计函数；父文件只读，不生成其缓存文件。无需新依赖或站点框架。

实际入口（从 r14_worktree 工作目录）：

```sh
python3 research/research_max_struct_r14/supplement_observable19/src/build_site.py --check-only
```

正文和分析完成后，主线程才运行正式构建：

```sh
python3 research/research_max_struct_r14/supplement_observable19/src/build_site.py
```

也可通过 `--root` 指定本补充研究的绝对目录。脚本默认根目录由自身位置确定，不依赖当前工作目录。需要保留父 `research/research_max_struct_r14/src/build_site.py`，因为渲染能力直接复用它。

正式构建写七个约定文档 HTML、galleries.html、evidence.html 和本目录 BUILD_RECEIPT.json。缺少正文时明确标记缺失；不生成 benchmark 模板或补写结论。Markdown 相对链接限制在本补充目录，链接到已有文档时改用对应 HTML；根外或私密链接抑制且留回执。父 R14 的科学引用请作者采用已核验的绝对网页 URL，不用 ../ 链接。

图和表只检查 analysis 及其直接分析子目录、这些目录中的实际 figures/FIGURES，以及补充根 FIGURES。不递归扫描数值目录，不读取大型数值 JSON、图像数组、RESULTS 或真值。所有发现图片都进入图集，不按结果挑选；同名 PDF 仅供下载。保留原图中的 case 标签，图注显示原文件名和完整来源相对路径。CSV 仅列下载链接，不再读取或计算指标；analyze_results 原 index.html 作为独立分析入口。

只列 DOWNLOADS.json 明确标记已核验的真实 https release 附件；没有文件或核验标记时显示未发布，不猜 URL，也不联网核验。源码及冻结配置有本地入口，原始大档通过实际 release 附件提供。

静态验收见 STATIC_SELF_CHECK.json。small_fixture 是本目录下的小型合成呈现测试，不是研究数据或正式网页。测试包含公式原 LaTex/下划线、缺件状态、相对链接、根外/私密链接抑制、未核下载筛除、图片/PDF去重。第一次测试发现 Mac 大小写目录别名导致重复卡片；已改为使用实际目录项名称并通过复测。

截至 INVENTORY_CHECK.json，真实根目录仅有理论澄清和冻结实验合同，尚无报告、复现、费用、失败、交接正文或分析图表。本任务未构建真实网站、未运行科学程序、未联网或发布。待主线程完成正文、最终分析与实际 DOWNLOADS 后构建，检查 BUILD_RECEIPT 中本地链接零错误，并视觉审阅公式和全部24对象图，再按既有 Pages 流程发布。

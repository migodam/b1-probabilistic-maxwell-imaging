# R14 静态研究网页构建说明

交付源文件：`src/build_site.py`。本次只写源文件和说明，没有执行网页构建、浏览器验收或发布，没有读密封GT、运行优化器/评分器或新科学实验。

## 父线程在报告完成后的使用

在已有兼容Python环境运行：

```text
python research/research_max_struct_r14/src/build_site.py --root research/research_max_struct_r14
```

程序只依赖Python标准库。视觉排版延续R13补充页的浅底、蓝链接、中文宽正文、可横向滚动表格和完整原尺寸图；没有新框架、服务、仪表盘或外部JS。输入Markdown中的公式以可复制的TeX原文显示，避免离线外部渲染依赖；父线程可以据实际需要后续改善公式排版。

主报告`RESEARCH_REPORT_ZH.md`完整渲染到index.html；其余实际存在的理论、轨迹、benchmark、费用、失败、复现、交接、论文与合同各有可读HTML及原Markdown下载。缺失报告不会被声称完成，缺失主报告时首页明确说明尚未生成，SITE_BUILD_RECEIPT记录缺失项目。

图集只遍历真实study中的analysis/*/FIGURES、audit图文件和根FIGURES，禁止跟随符号链接，不扫描D盘别名或整个工作树。以路径分组为暴露48、独立24、controls、3D、fullstate/root和其他，保留完整文件名、可识别caseID、原图及存在的配对PDF链接。路径不能判断科学类型时保留在“其他”，不臆造方法/胜负。PDF-only图以下载链接呈现，不创建假缩略图。未来新增campaign使用标准目录名会被自动发现。

`evidence.html`列出现有报告、src代码、configs、analysis/audit的公开Markdown/CSV。仅对已有analysis摘要CSV读取前60行作原值预览，不重算指标、不挑例、不生成科学结论。原始状态/输入NPZ、进程日志、任意JSON不会被自动暴露；完整数据、状态、代码包由父线程审定data closure后，通过实际发布附件提供。

DOWNLOADS.json只展示HTTPS `/releases/download/`附件，要求顶层`status=PUBLISHED_AND_ANONYMOUSLY_READ_BACK`或`verified=true`，或单asset明确`verified=true`/`verification_status=PASS|VERIFIED|ANONYMOUSLY_READ_BACK`。name/url使用父线程实际填写值；没有占位URL。网页不代替匿名网络readback，也不把本地构建等同发布。无清单或无核验资产时明确说明。

## 本地链接与隐私范围

Markdown原HTML一律转义；链接/图片仅由renderer生成。只允许HTTP(S)外链与解析后仍在实际study内的真实文件，拒绝越界路径、file/javascript等scheme、私有目录及别名。自动生成导航只连现存文档与本次生成页；相邻已有Pages通过仓库既有路由接入，由父线程发布时更新原入口，不生成新站点地址。

永久禁止原始链接：`sealed_truth`、`private_operational`、`sources/CHAT_READBACK`及含credentials/secrets的路径。程序不读取这些树；不因data closure通过而自动解禁。MD中指向这些文件的链接显示为不可公开文本说明。不会直接链接原始DOWNLOADS.json以免间接暴露未核验的占位URL。

构建后只对本次生成HTML解析href/src，验证本地目标真实存在且未越界/私有，并记录全部外链供父线程网络核验。结果写`SITE_BUILD_RECEIPT.json`；broken或private/outside链接使退出状态非零。该检查不扫描虚拟D盘、不检查SHA，不证明视觉质量，也不验证外站内容。

父线程需实际打开最终网页和相关原图，检查长公式、表格、所有case页面、跨页导航、下载附件及匿名可访问性后发布。正文科学判断、结果完整性与发布验收始终由父线程负责。若重新构建，当前生成页会更新；旧不再收录的图库页未自动删除，不再由新导航指向，父线程发布清单应以最新SITE_BUILD_RECEIPT.generated_pages为准。

## 公式与主图修订（最新行为）

当前改用固定MathJax3.2.2 CDN `https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js`；这是静态页的一段可选公式渲染JS，无新框架。`\(...\)`置于span保留delimiter，`\[...\]`与`$$...$$`置于可横向滚动div保留delimiter；fenced code仍是pre/code并被MathJax跳过。CDN不可用或禁用JS时，数学原始TeX仍在正文可读，不丢公式。支持$...$作为行内delimiter；作者正文若含货币$需转义。前文“公式只原文、无外部JS”的说明已被此修订取代。[官方配置文档](https://docs.mathjax.org/en/v3.2/web/configuration.html)

导航新增实际存在时的BASELINE_AUDIT、THEORY_STATUS、PRIOR_ART_AND_NOVELTY、METHODS_AND_ABLATIONS。主暴露48图集优先`FIGURES/exposed48`的8页可读色条版本。`analysis/exposed48_v2/FIGURES/materials_*`和all48_materials PDF只进单独archive页面，从evidence页链接，不再次混入主画廊；其其他汇总曲线仍按实际文件保留。bugged `analysis/exposed48_v1`整体排除主画廊和自动evidence列表/预览，仅允许作者FAILURES.md显式提供追溯链接。SITE_BUILD_RECEIPT记录归档旧图与排除v1文件。

本轮仅在/tmp生成一个很小的公式HTML自检（检查delimiter、div/span、code隔离、CDN与新增导航名单），未构建真实站点、未浏览验收大图、未发布。主报告完成后由父线程执行构建，并实际检查公式CDN加载、数学布局和最终完整材料图。包扫描器的父线程修正及473.875CPU失败回执未由本轮改动，package算法保持父线程版本。

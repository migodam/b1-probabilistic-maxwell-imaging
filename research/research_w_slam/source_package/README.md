# Research W-SLAM 研究交付包

版本：2026-10-24 UTC。

## 先读什么

1. `Research_W_SLAM_Report_ZH.pdf`：完整中文报告，24个正文章节、理论推导、反例、原创性裁决、E0–E7与A0实验门，以及固定研究证据和外部文献。
2. `Research_W_SLAM_Report_ZH.md`：与PDF内容对应的可编辑原稿。
3. `Literature_Audit.md`：27行原创性矩阵、逐篇模型与访问范围、未能取得全文的限制。
4. `Codex_Experiment_Prompts.md`：可独立移交的七阶段任务，包含默认物理参数、输出结构、预算、强基线和停止条件。
5. `experiment_contract.example.json`：机器可读的候选实验合同；状态是草案。`null` 表示尚待预注册，不表示可以自行忽略。真实墙体形状和孔洞数没有预填给算法。

`References_Main.md`与主报告的S编号一致；文献审计附录使用独立L编号。图件、图注及可再生成源代码位于`figures/`与`code/`。

## 本次做了与没有做的工作

已完成：固定仓库与相关旧理论审计、原始论文定位和方法范围核查、条件数学推导、反例、两个小型NumPy代数核验、报告与执行协议。

没有运行新的Maxwell PDE反演、宽材料模式搜索、measurement-only在线地图、主动采集、训练、速度基准或三维实验。附图中的恢复和冻结曲线是明确标注的有限维解析例子。D标签表示在声明假设下的条件推导，不表示十四项独立文献原创。

研究建议是以动态补偿与有限可辨识性为主线，先做最小任务证书应用。旧measurement-only锚点仍未取得，source-progressive和A17/ROM的独立净收益尚无证据。若普通smoother/direct posterior在相同精度下更便宜，应采用它们并停止相应加速主张。

## 复现代数核验

需要Python 3和NumPy；画图另需Matplotlib。可在本目录执行：

```bash
python code/observability_toy.py
python code/freezing_toy.py
python code/make_research_figures.py
```

脚本按自身路径定位输出目录，无需任何真实墙体、外部网络或PDE求解器。核验结果写入`results/`，随机种子为20261025；该整数只是固定种子，不是报告日期。浮点结果在不同BLAS/平台上可有最后若干位差异，断言采用明确容差。

关键已核验结果：共享状态例子单次信息0、0而联合为2；奇异B反例中错误伪逆给4/3，正确联合信息为0（浮点舍入范围）；100组SPD及100组resolvent/端口代数检查通过。它们不证明任意Maxwell材料类中的物理可实现性、后验校准或计算节约。

## 重建PDF

主PDF使用Pandoc、XeLaTeX、标准TeX包和Noto Sans CJK SC字体。提供一个本地Unicode CJK字体文件，例如`NotoSansCJKsc-Regular.otf`：

```bash
python code/build_report_pdf.py --font-path /absolute/path/NotoSansCJKsc-Regular.otf
```

构建器不下载字体或安装依赖。需安装DejaVu Sans/Mono作为拉丁扩展与代码字体。PDF已嵌入展示所需字体，阅读PDF不需要这些构建依赖。

## 下一轮从哪里开始

只先执行`Codex_Experiment_Prompts.md`中的统一约束和Prompt 0。E0和E1共用预算账本，不重复重置；完成几何、源/噪声尺度、前向资格和预注册后，再按gate进入下一阶段。未定外域/PML、源归一化、边界规则、事件和rank接受规则必须在正式数据生成前冻结。

保存所有失败、原始数据ID和成本；旧因子重算后同步修正粒子权重或信息系统。不要把相同测量的重复数值运算当成新的观测，不用GT支撑或旧两墙先验替代未知支持实验。

## 证据与文件完整性

所有外部论文和固定Git证据都通过报告中的原始URL定位，未复制外部论文全文。`results/environment_versions.json`记录本次运行环境；`results/delivery_validation.json`记录交付检查范围。文件校验值见`SHA256SUMS.txt`，该清单不包含自身。

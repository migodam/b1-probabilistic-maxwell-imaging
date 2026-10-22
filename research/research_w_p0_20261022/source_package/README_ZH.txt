Research W — 未知复杂墙体的 Maxwell 传播可辨识性与概率成像
研究日期：2026-10-22

一、建议先读

Research_W_Report_ZH.pdf 是完整中文报告。
第1节给出中心判断与证据分级；第3节是可辨识性及已执行反例；
第4节是四块传播签名、联合后验压缩和风险；第8–9节是既有项目及文献审计；
第10节是P0/P1；第12节回答最终立项、首选路线和停止条件。

本报告推荐：先研究有限单侧观测下的任务传播可辨识性，把概率压缩作为建设性后续。
不预设按深度递进、A17、Flow或Diffusion必须成为最终方案。
本文给出的标准恒等式、可证推论、解析反例、待检验机制和未解决问题均作了区分。
正确的数学推导不等于已建立论文首创性。

二、实际做过什么

1. 核验五个指定GitHub仓库的核心原始决策和报告，并固定commit。
   原报告的读取范围和未取得的L/Z/X、P3-L等资料在PDF中明确列出。
   没有重新运行全部历史项目实验。
2. 定向核验48条项目或文献来源，审计DBIM、layer stripping、blind Green、
   T-matrix、Marchenko/redatuming、Bayesian目标降维和已有A17-X/F理论重叠。
3. 独立推导任务纤维条件、条件性有限延拓、四块传播误差界、
   保持联合先验的后验误差界，以及错误固定和重复likelihood的问题。
4. 实际运行非磁、无损、无限平面介质板的一维精确Maxwell实验：
   完整内部反射、固定参考面双程传播、PEC目标级联、目标位置补偿和带宽扫描。
   传输矩阵与独立Fresnel公式核验通过；能量/幺正性缺陷在机器精度量级。

没有执行新的二维/三维成像、复杂墙体posterior采样、大网络训练或加速benchmark。
P0-B与P1是下一轮预先定义的实验建议，不能引用为已经成功的结果。

三、解析实验复现

依赖：Python 3、NumPy、SciPy、Matplotlib。
从本目录运行：

    python analytic_witness.py

会重新生成analytic_results中的JSON、CSV、PDF与PNG。
实验长度以中心频率的真空波长lambda0归一化。
墙1：n=2、d=1/4；墙2：n=3、d=1/6；mu=mu0；相同前表面。
每个复数据的实部与虚部噪声标准差分别为0.01；每档带宽9个频点。
这是已知两均值的检验距离计算，不是校准过的材料posterior。

关键预期结果：
中心反射约1.63e-16；同参考面双程响应差模为1；相位差pi/3。
PEC位置相差lambda0/12后联合反射差约2.45e-16。
0.5%带宽的wall-only距离0.887，而该联合补偿对距离已为3.038。
因此必须保留“完整目标数据可以打破特定墙体歧义”的负面限制。
本例不证明有限横向支撑、不规则墙、现实宽带或任意目标的普遍非唯一性。

四、报告源文件与重新编译

report_parts包含全部LaTeX章节和参考文献，main.tex为入口。
fonts/README.txt给出使用的Noto Serif CJK SC字体及下载入口。
已交付PDF嵌入了所需字形；阅读PDF不需要另装字体。
重新编译需要XeLaTeX、常见TeX Live包及两份字体文件：

    bash build_report.sh

图由analytic_results/exact_maxwell_witness.pdf引入。
SOURCE_MANIFEST.json保留所有参考文献键及链接；项目链接固定到本轮commit。
REPRODUCIBILITY.json记录本次运行的软件版本与输出元数据。

五、下一轮执行

P0_P1_PROTOCOL_ZH.txt是便于执行的提要；完整定义和限定以PDF第10节为准。
先冻结材料空间、任务、采集、噪声、预算和门槛，再搜索有限全波见证。
禁止用truth-selected candidates、oracle权重、未知真值fallback或新增未来观测
作为原固定数据的在线性能。

若只得到candidate ensemble，就按candidate ensemble交付。
只有明确先验、正确联合likelihood、合法权重和独立概率验证，才能使用posterior表述。


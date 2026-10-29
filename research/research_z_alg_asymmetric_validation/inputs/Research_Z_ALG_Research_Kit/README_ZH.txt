Research Z-ALG 研究复现包
日期：2026-10-29

一、主结论和阅读顺序

完整中文报告为单独交付的 Research_Z_ALG_Report_ZH.pdf，共38页。
本包用于审计推导、复现微型实验及启动下一轮 Codex 工作。

建议先读报告的研究判断、共享材料条件定理、两cell物理见证、带噪E4对照和原创性判断，再读 CODEX_PROTOCOL_ZH.txt。
结论为 CONDITIONAL GO：建立研究线，但当前尚未证明非对称系统的全解覆盖或一般计算加速。

二、文件结构

report_source/                  完整LaTeX源及已生成科学图
zalgo_maxwell/two_cell_witness.py
zalgo_maxwell/two_cell_witness.json
                                冻结两cell TM/VIE、商变量反演、暗电流不可实现反例
zalgo_probability/e4_two_cell_posterior.py
zalgo_probability/e4_two_cell_posterior.json
                                原材料与quotient的有限噪声参考积分
zalgo_probability/analytic_probability_checks.py
zalgo_probability/analytic_checks.json
                                QD、三次根权重、TV和单位方形解析概率检查
zalg_algebra_e0.py
zalg_algebra_e0.json             三次映射两条复闭路的根置换及解析结论汇总
CODEX_PROTOCOL_ZH.txt            完整下一轮执行协议
SOURCE_MANIFEST.json             R6/R7冻结来源、完整commit和文件路径
INPUT_AUDIT.json                 本次读取材料及SHA256、未独立核验范围
CLAIM_LEDGER.json                本轮结果资格与尚未完成事项
ENVIRONMENT.json                本轮Python与主要库版本
SHA256SUMS.json                  本包文件校验值（不包括其自身）
requirements.txt                Python依赖

三、复现命令

在解压后的 Research_Z_ALG_Research_Kit 目录内，使用具备 numpy、scipy、matplotlib 的Python：

python zalgo_maxwell/two_cell_witness.py
python zalgo_probability/e4_two_cell_posterior.py
python zalgo_probability/analytic_probability_checks.py
python zalg_algebra_e0.py
python report_source/figures/make_two_cell_structure.py

前三个脚本会更新同目录的结果JSON；若要保留交付时原始结果，请先复制工作目录。
E0脚本把两条monodromy闭路的数值结果输出到标准输出；对应汇总JSON还包含报告中证明的解析根与作用恒等式。
图形脚本只读取现有结果并绘图，不重新执行Maxwell求解。
不同平台的末位浮点值和时间可能不同；不要用逐字节结果相等代替数值与条件检查。

四、编译报告

需要 XeLaTeX（fontspec、amsmath、amsthm、mathtools、graphicx、booktabs、longtable、tabularx、enumitem、hyperref、fancyhdr、needspace）以及 DejaVu Sans Mono。
在研究包根目录建立 tmp/pdf_fonts/，放入官方字体 NotoSerifCJKsc-Regular.otf。
官方下载地址：
https://raw.githubusercontent.com/notofonts/noto-cjk/main/Serif/OTF/SimplifiedChinese/NotoSerifCJKsc-Regular.otf

字体未随本包重复分发。来源为Noto CJK项目，使用时遵循其OFL许可。
然后切换到 report_source/，连续执行两遍：

xelatex -interaction=nonstopmode -halt-on-error report.tex

若修改了分页或目录，按编译提示再执行一遍。最终报告已完成38页视觉检查及公式/表格复核；源码中使用的相对字体路径与上述目录一致。

五、结果资格

1. 两cell实验是真实Hankel Green的圆盘pulse-basis中心配点离散TM/VIE，不是Born替代模型，也不是连续Maxwell全解认证。
2. 材料交换的观测保持由解析恒等式证明。完整两根结论是条件定理：商线性系统可逆、两根合法、resolvent可逆。具体Hankel数值没有严格区间包围。
3. 冻结two_cell_witness.json中的symbolic/certified措辞须按上述条件定理理解；数值行列式非零诊断并非独立区间证书。报告正文明确区分解析证明和浮点见证。
4. 暗电流反例使用另一个源合同，不能与两轴向源的交换对称实验混为一体。
5. E4使用[0.5,4]^2上的rho=x1/27.5625、8个实化噪声分量各sigma=0.001。它与原Full66 prior不同，也与analytic_checks中的单位方形toy不同。
6. E4的32/64阶收敛和两个积分协议一致是数值参考，不是严格区间积分误差界。
7. 几何覆盖、材料可实现性、正确概率和计算收益分别声明。当前没有非对称全分支覆盖、Full66连续精确纤维或一般加速结论。
8. R6/R7、原Research Z和TSOM的核查范围在报告与INPUT_AUDIT中列明；旧大型实验未重跑。

六、费用

各实验JSON包含独立费用账本。E4总计15360次解析likelihood计算，并有5次2×2数值分解、10个复RHS用于独立校验。
解析公式仍有算术费用。不要把微型示例的毫秒级时间外推为大规模Maxwell加速。
每次重跑是新的一次实验费用；后续Codex协议要求统一记录构建、失败、认证和概率计算。

七、近期先例

Everett 2026（DOI:10.1093/gji/ggag432）是电磁材料多项式多解路径搜索的直接近邻。
报告给出独立的往返路径定理：沿gamma*sin(2*pi*t)的参数线段往返，在全程正则唯一提升下必返回原解。
奇点或Newton跳支可能产生另一个合法端点，但不自动建立循环群、正则monodromy、全覆盖或后验质量。

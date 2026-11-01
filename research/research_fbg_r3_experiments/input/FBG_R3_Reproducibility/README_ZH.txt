Research FBG-R3-T — 复算与研究附件
2026-10-31

这份附件包含本轮自行编写的数学核验、实际运行结果、图、研究备忘录和报告源文件。
所有程序都是小型有限维代数核验。没有新的24×24/576参数成像实验，没有新算法性能基准。
PDF中的历史成像数字来自原报告；本轮代码不重新生成这些历史图像。

一、目录

deliverables/
  Research_FBG_R3_T_Report_ZH.pdf    完整中文研究报告
  CODEX_FBG_R3_HANDOFF_ZH.txt        可直接交给Codex Sol Max的执行任务
  FBG_R3_Algorithm_Comparison.csv   18算法×14列完整比较矩阵，UTF-8 BOM

checks/
  root_math_checks.py              反馈梯度、源协变、有界Gram、同jet反例
  check_a17_math.py                固定B导数、profile材料步、0/2/4、有限响应
  branches_checks.py               法向导数、有限端口、两单元分支与弦线
  supplement_checks.py             dressing等价、软释放、改基、normal trace
  exact_two_cell.py                Python Fraction精确两单元根与无极点证书
  branch_corridor.py               由既有两单元数据生成示意图
  *.json / root_checks/            本轮真实运行结果
  branch_corridor.png              本轮标准绘图库生成的图

memos/
  agent_a17.md / agent_branches.md / agent_cie.md
  algorithm_release_final.md / algorithm_branch_final.md
  以上是独立研究线的中间推导和详细算法。保留它们用于追溯；若早期措辞与
  最终报告不同，以PDF和final算法稿为准。某条研究线未取得某资料，不代表
  主研究者最终没有取得；最终来源审计在PDF中。

report_sources/
  report_body.md / report_complete.md / matrix_appendix.md
  report_template.tex / report_complete.tex / build_report.py
  自行生成的图、封面和完整算法CSV。
  fonts/只包含开放字体许可和下载说明，不附第三方论文全文。

ENVIRONMENT.json                   运行依赖版本与研究来源固定版本
MANIFEST_SHA256.json                附件内部文件的SHA-256

二、复算顺序

需要Python 3.10+、NumPy、Matplotlib；精确两单元脚本只用标准库。
进入checks目录后运行：

  python root_math_checks.py
  python check_a17_math.py
  python branches_checks.py
  python supplement_checks.py
  python exact_two_cell.py
  python branch_corridor.py

脚本在自身附近写结果，重跑会更新对应JSON或图。请先复制附件或保留原JSON用于对照。
部分脚本的完整数值结果打印到stdout；浮点尾数可随BLAS/平台略变。
不要为极小差异反复改随机种子。exact_two_cell.py的有理检查应精确通过。

关键预期：
- 反馈梯度恒等式约1e-15，中心差分约1e-10；
- CIE理想读出差：一般方向约二阶、物理切向直线约四阶、精确物理路径为零；
- profile H/g/step公式与直接消元到双精度一致；
- raw/dressed同span给相同材料步；
- 运输后的ridge不受改基影响，裸系数ridge会改变问题；
- p=2等权完整复源模式组的normal trace占比为0.75，不能选模；
- 两单元二次式严格为(m1-2)(764347498747*m1-3470613075400)=0，
  两根合法且精确响应相同；
- 图中的电流弦线数据误差与材料重算后的物理误差必须分开。

三、报告重建（可选）

数学正文使用XeLaTeX，封面使用ReportLab和开源Noto轮廓；图由Matplotlib生成。
需要pandoc、xelatex及相应LaTeX包，另需要Python reportlab、fonttools、pypdf。
先按report_sources/fonts/FONT_SOURCES.txt取得四个开放字体，再在report_sources中运行：

  python build_report.py

输出写入相邻output目录。正文公式保持可搜索文本；封面中文为精确字体轮廓。
附件已包含最终PDF，单纯阅读和数学复算不需要重建PDF或下载字体。

四、范围与使用边界

1. 同jet反例没有证明bounded反馈优于CIE；它证明局部jet不足以决定有限搜索。
2. finite-chord公式假设两个端点已知；它不是未知另一端点的求解器。
3. Q6/Q20是已存在Z-ALG条件代数核，新的成像算法实现留给下一轮Handoff。
4. 本轮没有实现或运行报告中完整的受控释放/空间块分支成像算法。
5. 不把这些候选集合解释为概率、后验或全局解覆盖。
6. 没有打包下载的外部论文正文；可访问来源和SOURCE UNAVAILABLE边界在报告中。

五、固定原研究版本

主仓库：migodam/b1-probabilistic-maxwell-imaging
SHA：6659aa49dff2535996d40962791e9a5a6e974adc
A17：migodam/a17-maxwell-certified-exchange
SHA：2ac81b989a2f13d62df9d5a7abd690a33c2189c7

最终研究建议：优先检验材料可实现性法向与任务相关电流释放；小支持精确多根
提案作为独立竞争路线；任意exp/tanh/Kfield保留为机制对照，不先组合为大算法。

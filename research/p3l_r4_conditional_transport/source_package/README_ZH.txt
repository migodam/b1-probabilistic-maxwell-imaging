Research P3L-R3：Maxwell 弱几何输运、材料后验体积与风险导向计算
报告日期：2026-10-23

文件说明
1. Research_P3LR3_ZH.pdf：65 页完整中文报告，含主报告、历史对照、三个数学推导附录和原始文献。
2. CODEX_EXECUTION_PROMPT_ZH.txt：可独立交给 Codex 的完整下一阶段任务书。
3. source/Research_P3LR3_ZH.tex：可编辑、可编译的完整报告源码。
4. source/fonts/：中文字体和字体许可。Noto Serif CJK SC 采用 SIL Open Font License；详见 OFL.txt。
5. SOURCES.json：三个推导附录的来源元数据，正文的完整文献在 PDF / TeX 中。
6. MANIFEST.json：包内各文件的 SHA-256、报告范围与固定版本。

研究范围与证据
固定项目版本：6eadb4690c962aaa95821fc4056b677f0f924a6b
项目：https://github.com/migodam/b1-probabilistic-maxwell-imaging
主要研究目录：research/p3l_r2_risk_posterior_exploration
输入交接 ZIP 的 SHA-256：9486199125f4c72e12f17af1240c924603d6f7044a69dd980f29454d1adb9ab3
本报告使用所提供证据包、固定历史报告、原 P1 配置/源码及已检索原始文献。没有声称重新运行或独立重现历史 Maxwell 数值实验。
本轮为理论与算法设计：没有新增 Maxwell 实验、训练、采样或 GPU 作业。
所有新算法性能、阶段预算与放大门槛均为待验证的研究设计。完整后验资格尚未由已有实验建立。

阅读路线
先读“研究判断与阅读路线”及第1—4节，再读算法第6节、验证与 campaign 第7—8节、研究决策第9节。
第11—12节给出来源、新颖性边界和历史负结果。附录 A/B/C 展开几何、概率输运、材料任务与计算误差推导。
任务书明确区分未授权时的 Stage 0 静态准备和获得新授权后的实验 campaign，不以本报告代替实验授权。

重新编译
环境：XeLaTeX / TeX Live，fontspec、unicode-math、fvextra、titlesec、fancyhdr、longtable、booktabs 等标准宏包。
在 source 目录内先后执行两次以下命令：
  xelatex -interaction=nonstopmode -halt-on-error Research_P3LR3_ZH.tex
中文字体通过相对路径 fonts/ 读取；数学字体为 Latin Modern Math。少数作者名的扩展拉丁字符采用 DejaVu Serif。
PDF 已检查全页缩略图及公式/表格密集页；最终编译无 Overfull、缺失字形或未定义引用。

包内不重复包含用户提供的原始交接 ZIP，也不包含未授权生成的新物理数据。

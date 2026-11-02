Research MAX-STRUCT — 理论研究源包
日期：2026-11-02

内容与证据边界
Research_MAX_STRUCT_ZH.tex 是完整中文报告的独立 LaTeX 源，不依赖未打包的 Markdown。
template.tex 是排版模板参考。SOURCE_MANIFEST.json 含冻结版本、36项来源和未运行任务。
CODEX_SOL_MAX_HANDOFF_ZH.txt 是完整实验执行计划；新实验尚未执行。
本包不包含用户仓库副本、外部论文全文或任何凭据。

生成 PDF
需要 Python 3、XeLaTeX，以及 tex 文件声明的常见 TeX 包。
1. 按 assets/FONT_SOURCE.txt 放置字体。
2. 在本目录运行：python3 build_pdf.py
已经提供 figures/tensor_excitation_spectrum.png，无需重新绘图。
不同 TeX 版本可改变换页，但不改变数学内容。

理论示意图
需要 numpy 和 matplotlib。运行：python3 make_figure.py
图来自正文的精确局部激发谱公式和条件界，不是 Maxwell 成像结果。

唯一新运行的小模型
需要 numpy。运行：python3 dynamics_toy.py
脚本在同目录写 dynamics_toy.json。随包 JSON 是本轮实际结果。
这是四材料、三源、三接收器的复对称随机代数 VIE，不能冒充 Hankel–Green
或三维 Maxwell 成像实验。它只检查给定条件下的局部交替映射推论与源基协变性反例。
不要由这个小模型外推成熟 CIE 的全局性能、真实物理混沌或材料恢复成功率。

后续实验
严格按交接文件先做资格审查和最小决定性实验。禁止在线读取测试 GT 选解，
不先启动大规模扫参或 NN 训练。否证同样是有效产出。

完整性
FILE_HASHES.json 给出包内主要文件的 SHA-256。

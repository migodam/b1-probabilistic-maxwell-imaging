"""Render the pilot pages using the existing research site's typography and navigation."""
from pathlib import Path
import sys

PILOT = Path(__file__).resolve().parents[1]
PAPER = PILOT.parent / "z_alg_optical_inertia"
sys.path.insert(0, str(PAPER))
from build_site import html_page, parser

ARCHIVE = "https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-neural-dual-pilot-v1/Z_ALG_Neural_Dual_Pilot.zip"
PUBLISHED_PAPER = "https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/z_alg_optical_inertia/ip_final/"
MAIN = [
    ("RESEARCH_REPORT_ZH.md", "index.html", "研究结果"),
    ("FOLLOWUP_RESULTS_ZH.md", "interpolation.html", "连续权重与接收角度"),
    ("MATHEMATICAL_NOTES_ZH.md", "theory.html", "数学依据"),
    ("REPRODUCE.md", "reproduce.html", "复现"),
    ("NEXT_PRO_HANDOFF.md", "handoff.html", "给下一轮研究的交接"),
]
DETAILS = [
    ("EXPERIMENTAL_RESULTS.md", "tables.html", "逐组数值与费用"),
    ("PILOT_DESIGN.md", "design.html", "盲测设计"),
    ("SUPPLEMENT_DESIGN.md", "supplement.html", "配对噪声"),
    ("METHODS_TRAINING_ZH.md", "training.html", "训练"),
    ("COARSE_METHOD.md", "gradients.html", "谱梯度与固定对照"),
    ("NATIVE_VERIFICATION.md", "native.html", "原生区间核验"),
    ("ADVISOR_ZH.md", "advisor.html", "机制与创意讨论"),
    ("GEOMETRY_INERTIA_ZH.md", "state.html", "状态负方向"),
    ("followup/THEORY_ZH.md", "port-theory.html", "端口矩阵与接收孔径"),
    ("followup/NESTED_SEGMENT_COVER.md", "segment.html", "权重区间覆盖"),
    ("followup/MASKED_PORTS_REPLAY.md", "receivers.html", "接收角度"),
    ("followup/GALERKIN_COVER_ZH.md", "currents.html", "端点电流空间"),
    ("followup/MASKED_SEGMENT.md", "arcs.html", "连续接收弧"),
    ("followup/TRANSFERRED_WEIGHT_ZH.md", "transfer.html", "权重跨材料检验"),
    ("followup/ONLINE_TRANSFER_ZH.md", "online.html", "逐件调整权重"),
]


def main():
    routes = {source:name for source,name,label in MAIN+DETAILS}
    nav = [(page,label) for _,page,label in MAIN]
    details = [(page,label) for _,page,label in DETAILS]
    for source, page, _ in MAIN+DETAILS:
        md = (PILOT/source).read_text()
        title = next(line[2:].strip() for line in md.splitlines() if line.startswith("# "))
        content = parser().render(md, dict(source=source, root=PILOT, routes=routes))
        html = html_page(title, content, page, nav, details, source, ARCHIVE,
                         stylesheet="../z_alg_optical_inertia/research.css",
                         extra_links=[(PUBLISHED_PAPER, "已封存的英文投稿稿")])
        # The shared template's project name in the browser tab belongs to Paper 1.
        html = html.replace(" · Z-ALG 光学惯性", " · Maxwell 材料峰值认证")
        (PILOT/page).write_text(html)


if __name__ == "__main__":
    main()

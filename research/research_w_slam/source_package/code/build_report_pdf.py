#!/usr/bin/env python3
"""Build the Chinese report with Pandoc and XeLaTeX.

Requires pandoc, xelatex, a Unicode CJK OTF/TTF font, and normal TeX packages.
Example: python build_report_pdf.py --font-path /path/NotoSansCJKsc-Regular.otf
No network downloads are performed.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("--font-path", type=Path, required=True)
parser.add_argument("--output", type=Path)
parser.add_argument("--keep-tex", type=Path)
args = parser.parse_args()
root = Path(__file__).resolve().parents[1]
font = args.font_path.resolve()
if not font.is_file():
    raise SystemExit(f"Font does not exist: {font}")
for exe in ["pandoc", "xelatex"]:
    if not shutil.which(exe):
        raise SystemExit(f"Missing executable: {exe}")
output = (args.output or root / "Research_W_SLAM_Report_ZH.pdf").resolve()
header = r"""
\usepackage{fontspec}
\setmainfont[Path=FONTDIR/,BoldFont=FONTFILE,BoldFeatures={FakeBold=1.8},ItalicFont=FONTFILE,ItalicFeatures={FakeSlant=0.12}]{FONTFILE}
\setsansfont[Path=FONTDIR/,BoldFont=FONTFILE,BoldFeatures={FakeBold=1.8}]{FONTFILE}
\setmonofont[Scale=0.84]{DejaVu Sans Mono}
\newfontfamily\LatinFallback{DejaVu Sans}
\XeTeXlinebreaklocale "zh"
\XeTeXlinebreakskip=0pt plus 1pt
\usepackage{microtype}
\usepackage{titlesec}
\usepackage{fancyhdr}
\usepackage{etoolbox}
\usepackage{fvextra}
\usepackage{xurl}
\definecolor{InkBlue}{HTML}{163C50}
\definecolor{SoftTeal}{HTML}{245F6B}
\titleformat{\section}{\Large\bfseries\color{InkBlue}}{}{0em}{}
\titleformat{\subsection}{\large\bfseries\color{InkBlue}}{}{0em}{}
\titleformat{\subsubsection}{\normalsize\bfseries\color{SoftTeal}}{}{0em}{}
\titlespacing*{\section}{0pt}{18pt}{9pt}
\titlespacing*{\subsection}{0pt}{13pt}{6pt}
\titlespacing*{\subsubsection}{0pt}{10pt}{4pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\setlength{\emergencystretch}{4em}
\setlength{\headheight}{16pt}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\small\sffamily\color{InkBlue}RESEARCH W-SLAM}
\fancyhead[R]{\small\sffamily 理论审查与实验协议}
\fancyfoot[L]{\scriptsize\sffamily 2026-10-24 UTC}
\fancyfoot[R]{\small\thepage}
\renewcommand{\headrulewidth}{0.3pt}
\renewcommand{\footrulewidth}{0pt}
\AtBeginEnvironment{longtable}{\footnotesize\setlength{\tabcolsep}{3pt}\renewcommand{\arraystretch}{1.25}}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\},breaklines,breakanywhere,fontsize=\small}
\RecustomVerbatimEnvironment{verbatim}{Verbatim}{breaklines,breakanywhere,fontsize=\small}
\AtBeginDocument{\hypersetup{colorlinks=true,linkcolor=InkBlue,urlcolor=SoftTeal,citecolor=InkBlue}}
\let\originaltableofcontents\tableofcontents
\renewcommand{\tableofcontents}{\originaltableofcontents\clearpage}
""".replace("FONTDIR", font.parent.as_posix()).replace("FONTFILE", font.name)
with tempfile.TemporaryDirectory(prefix="wslam_pdf_") as temp:
    tmp = Path(temp)
    (tmp / "header.tex").write_text(header)
    (tmp / "inline_code.lua").write_text('''
function Code(el)
  return pandoc.RawInline("latex", "\\\\nolinkurl{" .. el.text .. "}")
end
''')
    tex = tmp / "report.tex"
    subprocess.run([
        "pandoc", str(root / "Research_W_SLAM_Report_ZH.md"), "-o", str(tex),
        "--standalone", "--from=markdown+tex_math_dollars", "--to=latex",
        "--toc", "--toc-depth=2", "--resource-path", str(root),
        "--include-in-header", str(tmp / "header.tex"),
        "--lua-filter", str(tmp / "inline_code.lua"),
        "--variable=fontsize:11pt", "--variable=papersize:a4",
        "--variable=geometry:margin=19mm", "--variable=classoption:titlepage",
        "--variable=linestretch:1.12", "--variable=toc-title:目录",
        "--variable=colorlinks:true", "--pdf-engine=xelatex",
    ], check=True, cwd=root)
    cover = r"""
\begin{titlepage}
\thispagestyle{empty}
\vspace*{18mm}
{\sffamily\large\color{SoftTeal}Research W / 理论与算法研究\par}
\vspace{16mm}
{\fontsize{30}{36}\selectfont\bfseries\color{InkBlue}Research W-SLAM\par}
\vspace{10mm}
{\fontsize{21}{32}\selectfont\bfseries 从子空间可观测性到\\可回访的概率 Maxwell 地图\par}
\vspace{12mm}
{\large 未知墙体与目标联合反演\\动态补偿、任务风险与有限非线性反馈\par}
\vspace{13mm}
{\color{SoftTeal}\hrule height 1pt}
\vspace{8mm}
{\normalsize 全文24章 · 条件推导与反例 · 原始文献审计\\可直接移交的分阶段实验协议\par}
\vfill
{\small 本次完成文献、条件推导与代数核验。\\新增全波成像、成本基准与三维实验尚未运行。\par}
\vspace{7mm}
{\normalsize 2026年10月24日（UTC）\par}
\end{titlepage}
"""
    tex_text = tex.read_text().replace("č", r"{\LatinFallback č}")
    if "\n\\maketitle\n" not in tex_text:
        raise SystemExit("Pandoc title command was not found")
    tex.write_text(tex_text.replace("\n\\maketitle\n", "\n" + cover + "\n", 1))
    # TeX runs from the report root so relative figure paths remain valid.
    for run in range(3):
        proc = subprocess.run([
            "xelatex", "-interaction=nonstopmode", "-halt-on-error",
            f"-output-directory={tmp}", str(tex),
        ], cwd=root, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        (tmp / f"xelatex_{run}.log.txt").write_text(proc.stdout)
        if proc.returncode:
            if args.keep_tex:
                args.keep_tex.mkdir(parents=True, exist_ok=True)
                for f in tmp.iterdir():
                    if f.is_file(): shutil.copy2(f, args.keep_tex / f.name)
            print(proc.stdout[-12000:])
            raise SystemExit(proc.returncode)
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(tmp / "report.pdf", output)
    if args.keep_tex:
        args.keep_tex.mkdir(parents=True, exist_ok=True)
        for f in tmp.iterdir():
            if f.is_file(): shutil.copy2(f, args.keep_tex / f.name)
    print(output)

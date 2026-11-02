"""Compile the complete reviewed LaTeX manuscript with XeLaTeX."""
from pathlib import Path
import subprocess

root = Path(__file__).resolve().parent
font = root / "assets" / "NotoSerifCJKsc-Regular.otf"
if not font.is_file():
    raise SystemExit("请先按 assets/FONT_SOURCE.txt 下载字体，保存到指定位置。")
for _ in range(3):
    subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error",
                    "Research_MAX_STRUCT_ZH.tex"], cwd=root, check=True)
print(root / "Research_MAX_STRUCT_ZH.pdf")

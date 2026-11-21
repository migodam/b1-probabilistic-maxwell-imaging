"""Render the tracked handbook for the existing B1 Pages site (Python-Markdown)."""
import csv
from html import escape
from pathlib import Path
import re
from urllib.parse import urlsplit

import markdown


HERE = Path(__file__).resolve().parent
PAGES = {
    "README.md": ("index.html", "实验手册"),
    "DATASETS_ZH.md": ("datasets.html", "数据设计"),
    "PROJECT_FAMILIES_ZH.md": ("families.html", "研究系列"),
    "OPERATORS_ZH.md": ("operators.html", "算子与设备"),
    "DESIGN_AND_CONTEXT_ZH.md": ("design.html", "运行与接续"),
    "BASELINES_AND_ERRATA_ZH.md": ("errata.html", "基线错题本"),
    "WALKTHROUGHS_ZH.md": ("walkthroughs.html", "实例与复算"),
    "STORAGE_ZH.md": ("storage.html", "磁盘管理"),
    "CHAT_HANDOFF.md": ("handoff.html", "交接"),
}
REPO = "https://github.com/migodam/b1-probabilistic-maxwell-imaging"
PUBLIC_DIR = "research/experiment_handbook"


def internal_link(match):
    url = match.group(1)
    path, separator, anchor = url.partition("#")
    if path in PAGES:
        url = PAGES[path][0] + separator + anchor
    return f'href="{url}"'


def write_page(path, title, content, source):
    nav = "".join(f'<a href="{file}"' + (' aria-current="page"' if file == path else '')
                  + f'>{name}</a>' for file, name in [*PAGES.values(), ("catalog.html", "查目录")])
    html = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)} · Maxwell 逆散射实验手册</title><link rel="stylesheet" href="handbook.css"></head>
<body><header><a href="../">B1 Maxwell 研究</a><nav aria-label="手册章节">{nav}</nav></header>
<main>{content}</main><footer><a href="{REPO}/blob/main/{PUBLIC_DIR}/{source}">查看本页来源</a> ·
<a href="{REPO}/tree/main/{PUBLIC_DIR}">代码与原始登记</a></footer></body></html>
"""
    (HERE / path).write_text(html, encoding="utf8")


def registry():
    with (HERE / "COHORTS.csv").open(encoding="utf8", newline="") as stream:
        cohorts = list(csv.DictReader(stream))
    with (HERE / "OPERATORS.csv").open(encoding="utf8", newline="") as stream:
        operators = list(csv.DictReader(stream))
    parts = ['<h1>数据与算子目录</h1><p>查找名称、来源、物理类别或研究轮次。每个数据条目列出可用于求解的输入和已知用途。</p>',
             '<label for="term">查找条目</label> <input type="search" id="term" placeholder="例如 Fresnel、R13 或三维"><p id="count" aria-live="polite"></p>',
             '<h2>数据集与观测版本</h2>']
    for r in cohorts:
        source = r["source"]
        link = (f'<a href="{escape(source, quote=True)}">原始来源</a>'
                if urlsplit(source).scheme in ("http", "https") else escape(source))
        parts.append(f'<article data-entry><h3><code>{escape(r["dataset"])}</code> · {escape(r["population"])}</h3>'
                     f'<p>{escape(r["origin"])} · {escape(r["kind"])}</p><p>{escape(r["contract"])}</p>'
                     f'<p>在线输入　{escape(r["online_fields"])}<br>实验用途　{escape(r["role"])}</p>{link}</article>')
    parts.append('<h2>算子与实现入口</h2>')
    for r in operators:
        link = (f'<a href="{REPO}/blob/main/{escape(r["repo_file"], quote=True)}">源码</a>'
                if r["repo_file"] else escape(r["workspace_file"]))
        parts.append(f'<article data-entry><h3><code>{escape(r["operator"])}</code> · {escape(r["device"])}</h3>'
                     f'<p><code>{escape(r["entry"])}</code></p><p>{escape(r["workload"])} {escape(r["accounting"])}</p>{link}</article>')
    parts.append("""<script>const rows=Array.from(document.querySelectorAll('[data-entry]'));
const term=document.getElementById('term'), count=document.getElementById('count');
function update(){let found=0; for(const row of rows){row.hidden=!row.textContent.toLowerCase().includes(term.value.trim().toLowerCase());if(!row.hidden)found++;}
count.textContent='显示 '+found+' / '+rows.length+' 个条目';}term.addEventListener('input',update);update();</script>""")
    write_page("catalog.html", "数据与算子目录", "\n".join(parts), "COHORTS.csv")


def main():
    for name, (file, title) in PAGES.items():
        body = markdown.markdown((HERE / name).read_text(encoding="utf8"),
                                 extensions=["tables", "fenced_code", "toc"])
        body = re.sub(r'href="([^"]*)"', internal_link, body)
        body = body.replace("<table>", '<div class="wide"><table>').replace("</table>", "</table></div>")
        write_page(file, title, body, name)
    registry()


if __name__ == "__main__":
    main()

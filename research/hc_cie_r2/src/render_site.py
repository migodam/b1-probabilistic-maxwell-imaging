"""Render the study's Markdown reports and preserve their mathematical notation."""
import html
import re
from pathlib import Path
from urllib.parse import urlsplit

import markdown


ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "https://github.com/migodam/b1-probabilistic-maxwell-imaging"
RELEASE = f"{REPOSITORY}/releases/tag/hc-cie-r2-20261128"
SPECIAL = {
    "RESEARCH_REPORT_ZH.md": "index.html",
    "CHAT_HANDOFF.md": "handoff.html",
    "README.md": "contents.html",
}
HANDBOOK = {"README.md": "index.html", "OPERATORS_ZH.md": "operators.html"}
MATH = re.compile(r"\$\$(.+?)\$\$|\\\[(.+?)\\\]|\\\((.+?)\\\)|(?<![\\$])\$([^\n$]+?)\$(?!\$)", re.S)
# Keep shell variables in the original code examples.
FENCES = re.compile(r"(^[~\x60]{3}.*?^[~\x60]{3}[ \t]*$)", re.M | re.S)


def page_name(source):
    return SPECIAL.get(source, Path(source).stem.lower() + ".html")


def body(source):
    math = []

    def preserve(match):
        expression = next(group for group in match.groups() if group is not None)
        display = match.group(1) is not None or match.group(2) is not None
        token = f"HCXMATHPLACE{len(math)}END"
        delimiters = (r"\[", r"\]") if display else (r"\(", r"\)")
        math.append((token, delimiters[0] + html.escape(expression) + delimiters[1]))
        return token

    sections = FENCES.split(source)
    for i in range(0, len(sections), 2):
        sections[i] = MATH.sub(preserve, sections[i])
    result = markdown.markdown("".join(sections), extensions=["fenced_code", "tables", "sane_lists"])
    for token, expression in math:
        result = result.replace(token, expression)

    def local_link(match):
        link = html.unescape(match.group(1))
        target = urlsplit(link)
        if target.scheme or target.netloc:
            return match.group(0)
        path = target.path
        if path.startswith("../experiment_handbook/") and Path(path).name in HANDBOOK:
            path = "../experiment_handbook/" + HANDBOOK[Path(path).name]
        elif path.endswith(".md") and "/" not in path:
            path = page_name(path)
        if target.fragment:
            path += "#" + target.fragment
        return 'href="' + html.escape(path, quote=True) + '"'

    result = re.sub(r'href="([^"]+)"', local_link, result)
    result = result.replace("<table>", '<div class="scroll"><table>').replace("</table>", "</table></div>")
    return result


def main():
    pages = sorted(ROOT.glob("*.md"))
    nav = [("index.html", "研究结论"), ("theory.html", "证明"), ("results.html", "逐案证据"),
           ("reproduce.html", "复现"), ("handoff.html", "交接")]
    for source in pages:
        content = source.read_text()
        title = html.escape(content.splitlines()[0].removeprefix("# "))
        filename = page_name(source.name)
        menu = " ".join(f'<a href="{path}"' + (' aria-current="page"' if path == filename else '')
                        + f'>{name}</a>' for path, name in nav)
        lang = "en" if source.name == "PAPER_DRAFT.md" else "zh-CN"
        markup = f"""<!doctype html>
<html lang="{lang}"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title} · HC-CIE-R2</title><link rel="stylesheet" href="study.css">
<script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script></head>
<body><header><a href="../">Maxwell 逆散射研究</a><span>HC-CIE-R2</span><nav>{menu}</nav></header>
<main>{body(content)}</main><footer><a href="contents.html">完整目录</a>
<a href="{REPOSITORY}/blob/main/research/hc_cie_r2/{source.name}">Markdown 原文</a>
<a href="{RELEASE}">数值附件</a></footer></body></html>"""
        (ROOT / filename).write_text(markup + "\n")


if __name__ == "__main__":
    main()

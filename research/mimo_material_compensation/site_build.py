"""Build the Chinese research pages from the scientific records. Requires Markdown."""
from pathlib import Path
from urllib.parse import urlsplit
import html
import os
import re
import markdown

ROOT = Path(__file__).resolve().parent
ASSETS = "https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/mimo-material-compensation-20261128/"
ARCHIVE = ASSETS + "mimo-material-compensation-20261128.zip"
REPO = "https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/mimo_material_compensation"
PAGES = {"RESEARCH_REPORT_ZH.md":"index.html", "NEXT_PRO_HANDOFF.md":"handoff.html",
         "REPRODUCE.md":"reproduce.html", "theory/THEORY_ZH.md":"theory/index.html",
         "literature/PRECEDENTS_ZH.md":"literature/index.html", "reviews/RESPONSE_ZH.md":"reviews/index.html"}
MENU = [("index.html","研究报告"),("theory/index.html","数学推导"),("reproduce.html","原始数据与复算"),("handoff.html","研究交接")]


def allowed(doc):
    p = doc.relative_to(ROOT)
    if "inputs" in p.parts or doc.name == "TASK_RESUME.md":
        return False
    return "reviews" not in p.parts or doc.name in ("RESPONSE_ZH.md", "public_response.md")


def dest(doc):
    rel = doc.relative_to(ROOT).as_posix()
    return ROOT/PAGES.get(rel, str(Path(rel).with_suffix(".html")))


def rel(path, folder):
    return os.path.relpath(path, folder).replace(os.sep, "/")


def emit(doc):
    out = dest(doc)
    raw = doc.read_text()
    title_match = re.search(r"^#\s+(.+)$", raw, re.M)
    title = html.escape(title_match.group(1).replace(chr(96), "") if title_match else doc.stem)
    formulas = []
    def reserve(match):
        formulas.append(html.escape(match.group(0), quote=False))
        return "XMATHPLACEHOLDER" + str(len(formulas)-1) + "ENDX"
    # Markdown escape processing would otherwise drop the TeX delimiters.
    raw = re.sub(r"\\\[.*?\\\]|\\\(.*?\\\)|\$\$.*?\$\$|(?<![\\\w])\$(?![\s$])[^\n$]+?(?<!\s)\$(?!\w)", reserve, raw, flags=re.S)
    md = markdown.Markdown(extensions=["tables","fenced_code","sane_lists","toc"],
                           extension_configs={"toc":{"slugify":lambda value,sep: re.sub(r"[^\w\-]+","-",value.lower()).strip("-")}})
    body = md.convert(raw)
    body = re.sub(r"XMATHPLACEHOLDER(\d+)ENDX", lambda m: formulas[int(m.group(1))], body)
    def link(m):
        label,raw_url = m.group(1),html.unescape(m.group(2))
        url = urlsplit(raw_url)
        if url.scheme or url.netloc or not url.path:
            return m.group(0)
        target=(doc.parent/url.path).resolve()
        if target.suffix.lower()==".md" and target.is_file() and allowed(target):
            target=dest(target)
        if target.suffix in {".zip",".npz",".npy"}:
            special={"prior-core.zip":"prior-core.zip","matched_amplitude_energy_archive.zip":"matched_amplitude_energy_archive.zip"}
            href = ASSETS+special[target.name] if target.name in special else ARCHIVE
        elif not str(target).startswith(str(ROOT)+"/"):
            return m.group(0)
        else:
            href = rel(target,out.parent)
            if url.fragment: href += "#"+url.fragment
        return label+'="'+html.escape(href,quote=True)+'"'
    body = re.sub(r'(href|src)="([^"]+)"',link,body)
    body=body.replace("<table>",'<div class="table-wrap"><table>').replace("</table>","</table></div>")
    subheadings = re.findall(r'<h2 id="([^"]+)">([^<]+)</h2>',body)
    toc = '<details class="toc"><summary>本页内容</summary><ol>'+''.join('<li><a href="#'+key+'">'+name+'</a></li>' for key,name in subheadings)+'</ol></details>' if len(subheadings)>3 else ''
    pos = body.find("</h1>")
    if toc and pos>=0: body=body[:pos+5]+toc+body[pos+5:]
    nav = ''.join('<a'+(' aria-current="page"' if out==ROOT/p else '')+' href="'+rel(ROOT/p,out.parent)+'">'+s+'</a>' for p,s in MENU)
    math = r"""<script>window.MathJax={tex:{inlineMath:[['\\(','\\)'],['$','$']],displayMath:[['\\[','\\]'],['$$','$$']]},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']},svg:{fontCache:'global'}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-svg.js"></script>"""
    page = f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="color-scheme" content="light"><title>{title} · 微波材料反演</title><link rel="stylesheet" href="{rel(ROOT/'site.css',out.parent)}">{math}</head>
<body><header><div class="topline"><a class="brand" href="{rel(ROOT/'index.html',out.parent)}">微波材料反演</a><nav aria-label="项目导航">{nav}</nav></div></header>
<main><article>{body}</article></main><footer><div><a href="{ARCHIVE}">完整研究包</a><a href="{REPO}">公开源码</a><a href="https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/">研究目录</a></div></footer></body></html>"""
    out.write_text(page)


if __name__=="__main__":
    docs = sorted(p for p in ROOT.rglob("*.md") if allowed(p))
    for doc in docs: emit(doc)
    print(f"Rendered {len(docs)} research pages")

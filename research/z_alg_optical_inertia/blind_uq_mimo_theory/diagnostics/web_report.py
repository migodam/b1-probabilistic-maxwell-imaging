# /// script
# dependencies = ["markdown==3.11.4"]
# ///
"""Render the research Markdown into the static B1 project site.

Run this file with the study requirements installed, or with `uv run`.
No publication, Markdown edits, or changes to the study results.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from html import escape, unescape
from html.parser import HTMLParser
import os
from pathlib import Path
import re
import unicodedata
from urllib.parse import quote, unquote, urlsplit, urlunsplit

import markdown


ROOT = Path(__file__).resolve().parents[1]
PROJECT = "https://github.com/migodam/b1-probabilistic-maxwell-imaging"
WEBSITE = "https://migodam.github.io/b1-probabilistic-maxwell-imaging/"
TAG = "z-alg-blind-uq-mimo-20261127"
RELEASE = f"{PROJECT}/releases/tag/{TAG}"
DOWNLOADS = f"{PROJECT}/releases/download/{TAG}"
ARCHIVES = {
    "study": f"{DOWNLOADS}/blind-uq-mimo-study.zip",
    "operators": f"{DOWNLOADS}/blind-uq-mimo-operators.zip",
    "scoring": f"{DOWNLOADS}/blind-uq-mimo-scoring-truths.zip",
}

CSS = """
:root{color-scheme:light;--paper:#fbfaf6;--ink:#172f39;--muted:#465e66;--accent:#006e74;--deep:#075967;--rule:#ccd7da;--tint:#e6efee;--code:#edf1ef;--white:#fff;--max:77rem;--read:50rem}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;color:var(--ink);background:var(--paper);font:1rem/1.82 system-ui,-apple-system,"PingFang SC","Hiragino Sans GB","Microsoft YaHei",sans-serif;overflow-wrap:break-word}
a{color:var(--accent);text-decoration-thickness:.075em;text-underline-offset:.19em}
a:hover{color:var(--deep);text-decoration-thickness:.12em}
:focus-visible{outline:3px solid var(--accent);outline-offset:4px;border-radius:3px}
.skip{position:fixed;left:1rem;top:-8rem;background:var(--white);padding:.55rem .8rem;z-index:4}
.skip:focus{top:1rem}
.site-header{border-bottom:1px solid var(--rule);background:var(--paper)}
.top,.hero,.layout,.footer-inner{width:min(100%,var(--max));padding-inline:clamp(1rem,3.3vw,2.7rem);margin-inline:auto}
.top{display:flex;align-items:center;flex-wrap:wrap;gap:.4rem 1.8rem;padding-block:.7rem}
.brand{display:inline-flex;align-items:center;gap:.65rem;font-weight:700;white-space:nowrap}
.brand::before{width:1rem;height:.33rem;background:var(--accent);content:""}
.site-nav{display:flex;align-items:center;gap:.12rem .3rem;flex-wrap:wrap}
.site-nav a{display:block;padding:.39rem .61rem;color:var(--muted);font-size:.92rem;white-space:nowrap;text-decoration:none;border-radius:4px}
.site-nav a:hover{background:var(--code);color:var(--deep)}
.site-nav a[aria-current=page]{background:var(--tint);color:var(--deep);font-weight:650}
.hero{padding-block:clamp(1.9rem,5vw,4.5rem) clamp(1.35rem,3vw,2.7rem)}
.hero .rule{width:3.2rem;height:.25rem;background:var(--accent);margin:0 0 1rem}
.hero h1{max-width:var(--read);margin:.24em 0 0;font-size:clamp(2rem,5.1vw,3.7rem);line-height:1.18;letter-spacing:-.04em;text-wrap:balance}
.kicker{color:var(--muted);font-size:.95rem;margin:0;line-height:1.45}
.hero-links{display:flex;flex-wrap:wrap;gap:.32rem 1.35rem;margin-top:1.12rem;font-size:.93rem}
.layout{display:grid;grid-template-columns:minmax(12.2rem,16.2rem) minmax(0,1fr);gap:clamp(1.45rem,4vw,3.6rem);align-items:start;padding-bottom:4rem}
.layout.no-toc{grid-template-columns:minmax(0,var(--read));justify-content:center}
.sidebar{position:sticky;top:1.2rem;border-top:3px solid var(--rule);padding-top:.55rem;max-height:calc(100vh - 3rem);overflow:auto}
.sidebar details>summary{font-weight:660;color:var(--deep);cursor:pointer}
.sidebar ol{list-style:none;margin:.5rem 0 0;padding-inline-start:0}
.sidebar ol ol{padding-inline-start:.75rem;margin-top:.15rem}
.sidebar li{margin:.25rem 0}
.sidebar a{display:block;color:var(--muted);font-size:.88rem;line-height:1.5;padding:.23rem .42rem;text-decoration:none;border-radius:3px}
.sidebar a:hover{background:var(--tint);color:var(--deep)}
.article{min-width:0;max-width:var(--read);border-top:1px solid var(--rule);padding-top:1rem}
.article>p:first-child{font-size:1.09rem}
.article p{margin:.85em 0}
.article h2,.article h3,.article h4{line-height:1.38;letter-spacing:-.025em;scroll-margin-top:1.25rem}
.article h2{font-size:clamp(1.43rem,3.2vw,1.95rem);margin:1.9em 0 .6em;border-top:1px solid var(--rule);padding-top:.8em}
.article h2:first-child{margin-top:.1em;border-top:0;padding-top:0}
.article h3{font-size:1.19rem;margin:1.53em 0 .52em}
.article h4{font-size:1.05rem;margin:1.4em 0 .4em}
.article strong{color:#173e47}
.article ul,.article ol{padding-inline-start:1.45rem}
.article li{margin:.4em 0}
.article blockquote{margin:1.15rem 0;border-inline-start:.22rem solid var(--accent);background:var(--tint);padding:.3rem 1rem;color:var(--muted)}
.article hr{border:0;border-top:1px solid var(--rule);margin:2.3rem 0}
.article figure{margin:1.55rem 0}
.article img{height:auto;max-width:100%;vertical-align:middle}
.article figure img{display:block;background:var(--white);border:1px solid var(--rule)}
.article figcaption{color:var(--muted);font-size:.91rem;padding-top:.5rem;line-height:1.58}
.table-scroll{max-width:100%;overflow-x:auto;margin:1.1rem 0 1.5rem;border:1px solid var(--rule);border-radius:4px;scrollbar-color:var(--muted) var(--paper)}
.table-scroll table{border-collapse:collapse;width:100%;min-width:35rem;font-size:.92rem;line-height:1.56;font-variant-numeric:tabular-nums}
.table-scroll th,.table-scroll td{padding:.56rem .72rem;text-align:left;vertical-align:top;border-bottom:1px solid var(--rule)}
.table-scroll th{background:var(--tint);font-weight:700}
.table-scroll tr:last-child td{border-bottom:0}
.table-scroll tbody tr:nth-child(even) td{background:#f1f4f0}
code{background:var(--code);padding:.17em .29em;border-radius:2px;font: .88em ui-monospace,"SF Mono",Menlo,Consolas,monospace}
pre{margin:1rem 0;overflow:auto;max-width:100%;background:var(--code);border:1px solid var(--rule);padding:1rem;border-radius:4px;line-height:1.5}
pre code{padding:0}
.formula{overflow-x:auto;max-width:100%;margin:1rem 0;padding:.25rem .08rem}
.formula mjx-container[jax=CHTML][display=true]{margin:.55rem 0!important}
.math-inline{white-space:normal}
.archive-link small{color:var(--muted);font-size:.79em;margin-left:.25em}
.page-end{display:flex;gap:.4rem 1rem;flex-wrap:wrap;border-top:1px solid var(--rule);margin-top:2.3rem;padding-top:1rem;font-size:.92rem}
.footer{border-top:1px solid var(--rule)}
.footer-inner{display:flex;gap:.4rem 1.5rem;flex-wrap:wrap;justify-content:space-between;padding-block:1.15rem 1.8rem;color:var(--muted);font-size:.91rem}
@media(max-width:800px){.layout{grid-template-columns:minmax(0,1fr);gap:.92rem}.sidebar{position:static;max-height:initial}.sidebar ol:first-child{columns:2}.sidebar ol ol{columns:1}.layout.no-toc{justify-content:normal}.article{max-width:none}}
@media(max-width:460px){.table-scroll table{min-width:31rem}.sidebar ol:first-child{columns:1}.top{gap:.2rem}.site-nav a{padding:.34rem .45rem}.hero{padding-top:1.45rem}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}
@media print{.site-header,.sidebar,.hero-links,.page-end,.footer{display:none}.layout{display:block}.article,.hero{max-width:none}.table-scroll{overflow:visible}.table-scroll table{min-width:0}a{color:inherit}}
"""

MATHJAX = r"""
<script>
window.MathJax = {
  tex: {inlineMath: [['\(', '\)']], displayMath: [['\[', '\]']], processEscapes: true},
  options: {skipHtmlTags: ['script','noscript','style','textarea','pre','code']}
};
</script>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js" defer></script>
"""


@dataclass(frozen=True)
class Page:
    source: str
    target: str
    label: str
    group: str


PAGES = [
    Page("BENCHMARK_REPORT_ZH.md", "index.html", "实验报告", "研究报告"),
    Page("THEORY_REFINEMENT.md", "theory.html", "理论", "理论与证明"),
    Page("REPRODUCE.md", "reproduce.html", "复现", "实验复现"),
    Page("NEXT_PAPER_BLUEPRINT.md", "handoff.html", "研究交接", "研究交接"),
    Page("OBSERVATIONS.md", "observations.html", "观测与材料", "数据台账"),
    Page("bounds/RESULTS_ZH.md", "bounds/index.html", "区域界与暗电流", "量化结果"),
    Page("bounds/roi_mimo_RESULTS_ZH.md", "bounds/roi_mimo.html", "MIMO 区域峰值", "量化结果"),
    Page("imaging/RESULTS_ZH.md", "imaging/index.html", "盲成像结果", "量化结果"),
    Page("imaging/README_ZH.md", "imaging/method.html", "盲成像方法", "方法与输入"),
    Page("ports/PORTS_ZH.md", "ports/index.html", "互耦阵列与实测", "量化结果"),
    Page("ports/ROI_COUPLING_ZH.md", "ports/roi_coupling.html", "互耦下的区域峰值", "量化结果"),
    Page("ports/CLASS_CERTIFICATE_ZH.md", "ports/class.html", "整类材料分离", "理论与复核"),
    Page("results/aperture_control/REPORT_ZH.md", "results/aperture_control/index.html", "孔径配对控制", "实验补充"),
    Page("results/aperture_joint/REPORT_ZH.md", "results/aperture_joint/index.html", "同读数孔径控制", "实验补充"),
    Page("diagnostics/PROFILE_DESIGN.md", "diagnostics/profile.html", "物理可行范围设计", "实验设计"),
]

REVIEW_PAGE = Page("reviews/RESPONSE_ZH.md", "review.html", "评审", "独立评审")
PUBLIC_REVIEW_FILES = {"RESPONSE_ZH.md", "neon_02.md", "neon_03.md", "neon_04.md"}
PAGES += [
    Page("reviews/neon_02.md", "reviews/review-02.html", "第二次公开评审", "独立评审"),
    Page("reviews/neon_03.md", "reviews/review-03.html", "第三次公开评审", "独立评审"),
    Page("reviews/neon_04.md", "reviews/review-04.html", "最终数值核查", "独立评审"),
]

# Files below the report tree that are deployed as static assets. Binary arrays
# and source files retain their names as a ZIP member hint.
STATIC_SUFFIXES = {".png", ".jpg", ".jpeg", ".webp", ".gif", ".svg", ".pdf", ".html"}
PUBLIC_SOURCE_SUFFIXES = {".csv", ".json", ".md", ".py", ".tex", ".tsv", ".txt"}
PUBLIC_ARRAY_NAMES = {"mechanism_witness.npz", "mechanism_class_collision.npz"}
BLOCKED_TAGS = {"script", "style", "iframe", "object", "form"}
MATH_EXPR = re.compile(
    r"(?P<code>```[^\n]*\n[\s\S]*?^```[^\n]*$|~~~[^\n]*\n[\s\S]*?^~~~[^\n]*$|(?P<tick>`+)[^\n]*?(?P=tick))"
    r"|(?P<block>\\\[[\s\S]*?\\\]|\$\$[\s\S]*?\$\$)"
    r"|(?P<inline>\\\([^\n]*?\\\)|(?<!\\)\$[^\n$]+?(?<!\\)\$)",
    re.MULTILINE,
)


def href(origin, destination):
    return quote(os.path.relpath(destination, Path(origin).parent), safe="/-_.~")


def slug(value, separator):
    text = unicodedata.normalize("NFKC", unescape(value)).strip().casefold()
    text = re.sub(r"[^\w\s-]", "", text)
    return re.sub(r"[-\s]+", separator, text).strip(separator)


def protect_math(source):
    entries = {}

    def capture(match):
        if match["code"] is not None:
            return match["code"]
        kind = "block" if match["block"] is not None else "inline"
        formula = match[kind]
        tex = formula[2:-2] if formula.startswith(("\\", "$$")) else formula[1:-1]
        token = f"mJxScience{len(entries):06d}Token"
        entries[token] = (kind, tex)
        if kind == "block":
            line = source.rfind("\n", 0, match.start()) + 1
            prefix = source[line:match.start()]
            if not prefix.strip():
                return f"\n\n{prefix}{token}\n\n"
        return token

    return MATH_EXPR.sub(capture, source), entries


def restore_math(html, entries):
    for token, (kind, tex) in entries.items():
        tex = escape(tex, quote=False)
        if kind == "block":
            display = f'<div class="formula">\\[{tex}\\]</div>'
            if f"<p>{token}</p>" in html:
                html = html.replace(f"<p>{token}</p>", display)
            else:
                html = html.replace(token, f'<span class="math-inline">\\({tex}\\)</span>')
        else:
            html = html.replace(token, f'<span class="math-inline">\\({tex}\\)</span>')
    return html


def archive_for(path):
    name = path.as_posix()
    if name.startswith("ports/public/operators/"):
        return "operators"
    if name.startswith("ports/private/"):
        return "scoring"
    return "study"


class Rewrite(HTMLParser):
    def __init__(self, page, mapping, raw_relative):
        super().__init__(convert_charrefs=False)
        self.page = page
        self.mapping = mapping
        self.raw_relative = raw_relative
        self.output = []
        self.archive_labels = []
        self.suppress = 0

    def location(self, value, image=False):
        split = urlsplit(value)
        if split.scheme:
            return (value, None) if split.scheme in {"http", "https", "mailto"} else ("", None)
        if not split.path:
            return value, None
        if split.path.startswith("/"):
            if re.match(r"/(Volumes|Users|tmp|private)/",split.path):
                return "", None
            return value, None
        candidate = (ROOT / self.page.source).parent.joinpath(unquote(split.path)).resolve()
        if not candidate.is_relative_to(ROOT):
            if not candidate.is_relative_to(ROOT.parents[2]):
                return "", None
            url = href(self.page.target,os.path.relpath(candidate,ROOT))
            return urlunsplit(("", "", url, split.query, split.fragment)), None
        path = candidate.relative_to(ROOT)
        addressed = self.mapping.get(path.as_posix())
        if addressed:
            url = href(self.page.target, addressed)
            return urlunsplit(("", "", url, split.query, split.fragment)), None
        suffix = path.suffix.lower()
        if path.parts[:1] == ("reviews",):
            allowed = PUBLIC_REVIEW_FILES | {"review-02.html", "review-03.html", "review-04.html"}
            if path.parent != Path("reviews") or path.name not in allowed:
                return "", None
        eligible = ("private" not in path.parts and
                    path.name!="ORIGINAL_REQUEST.txt" and not any(part.startswith(".") for part in path.parts))
        # Long per-iteration traces are distributed intact with the raw observations.
        is_trace = (path.parts[:2] == ("imaging", "results") and suffix == ".json"
                    and (ROOT/path).is_file() and (ROOT/path).stat().st_size > 200_000)
        published_source = eligible and not is_trace and (ROOT/path).is_file() and (
            suffix in PUBLIC_SOURCE_SUFFIXES or path.name in PUBLIC_ARRAY_NAMES or self.raw_relative)
        if image or suffix in STATIC_SUFFIXES or published_source:
            url = href(self.page.target,path)
            return urlunsplit(("", "", url, split.query, split.fragment)), None
        return ARCHIVES[archive_for(path)], path.as_posix()

    def render_tag(self, tag, attrs, closed=False):
        attributes = dict(attrs)
        if tag == "a" and "href" in attributes:
            target, member = self.location(attributes["href"])
            if target:
                attributes["href"] = target
            else:
                del attributes["href"]
            self.archive_labels.append(bool(member))
            if member:
                attributes["class"] = "archive-link"
                attributes["title"] = "下载归档，文件 " + member
        elif tag == "img" and "src" in attributes:
            attributes["src"], _ = self.location(attributes["src"], image=True)
            attributes["loading"] = "lazy"
            attributes["decoding"] = "async"
        if tag == "table":
            self.output.append('<div class="table-scroll" role="region" tabindex="0" aria-label="可横向滚动的数据表">')
        rendered = "".join(
            f' {k}="{escape(v or "", quote=True)}"' for k, v in attributes.items()
            if not k.lower().startswith("on")
        )
        self.output.append(f"<{tag}{rendered}{' /' if closed else ''}>")

    def handle_starttag(self, tag, attrs):
        if tag in BLOCKED_TAGS:
            self.suppress += 1
        elif not self.suppress:
            self.render_tag(tag, attrs)

    def handle_startendtag(self, tag, attrs):
        if not self.suppress and tag not in BLOCKED_TAGS:
            self.render_tag(tag, attrs, closed=True)

    def handle_endtag(self, tag):
        if tag in BLOCKED_TAGS:
            self.suppress = max(0,self.suppress-1)
            return
        if self.suppress:
            return
        if tag == "a" and self.archive_labels and self.archive_labels.pop():
            self.output.append("<small>（归档内）</small>")
        self.output.append(f"</{tag}>")
        if tag == "table":
            self.output.append("</div>")

    def handle_data(self, value):
        if not self.suppress:
            self.output.append(value)

    def handle_entityref(self, name):
        if not self.suppress:
            self.output.append(f"&{name};")

    def handle_charref(self, name):
        if not self.suppress:
            self.output.append(f"&#{name};")


def local_link(origin, target, label, active=False):
    current = ' aria-current="page"' if active else ""
    return f'<a href="{href(origin,target)}"{current}>{escape(label)}</a>'


def toc_markup(tokens, entries):
    bits = []
    for item in tokens:
        label = escape(unescape(item["name"]))
        for key, (_, tex) in entries.items():
            label = label.replace(key, escape(tex,quote=False))
        child = toc_markup(item["children"], entries)
        anchor = quote(item["id"],safe="-_.~")
        bits.append(f'<li><a href="#{anchor}">{label}</a>{child}</li>')
    return "<ol>"+"".join(bits)+"</ol>" if bits else ""


def add_figures(html):
    def wrap(match):
        tag = match[1]
        description = re.search(r'\balt="([^"]*)"', tag)
        caption = f"<figcaption>{description[1]}</figcaption>" if description and description[1] else ""
        return f"<figure>{tag}{caption}</figure>"

    return re.sub(r"<p>\s*(<img [^>]+>)\s*</p>",wrap,html)


def render_markdown(source, page, mapping, raw_relative):
    protected, entries = protect_math(source)
    converter = markdown.Markdown(
        extensions=["tables", "fenced_code", "sane_lists", "toc"],
        extension_configs={"toc":{"slugify":slug,"toc_depth":"2-3"}},
        output_format="html",
    )
    fragment = converter.convert(protected)
    heading = re.search(r'<h1(?: id="([^"]*)")?>(.*?)</h1>', fragment, flags=re.DOTALL)
    if heading:
        title = restore_math(heading[2], entries)
        heading_id = heading[1] or ""
        fragment = fragment[:heading.start()] + fragment[heading.end():]
    else:
        title, heading_id = escape(page.label), ""
    rewriter = Rewrite(page, mapping, raw_relative)
    rewriter.feed(fragment)
    rewriter.close()
    article = restore_math(add_figures("".join(rewriter.output)), entries)
    toc = toc_markup(converter.toc_tokens, entries)
    return title, heading_id, article, toc


def htmldoc(page, title, heading_id, article, toc):
    title_plain = re.sub(r"<[^>]+>", "", title)
    label = escape(page.group)
    sep = "" if page.target=="index.html" else "相容区间与互耦阵列 · "
    tabs = [
        ("index.html","实验报告"), ("theory.html","理论"), ("reproduce.html","复现"),
        ("review.html","评审"), ("handoff.html","研究交接"),
    ]
    nav = "".join(local_link(page.target,dest,name,page.target==dest) for dest,name in tabs)
    downloads = " ".join(f'<a href="{ARCHIVES[key]}">{name}</a>' for key,name in (
        ("study","代码与实验数据 ZIP"), ("operators","原生算子 ZIP"), ("scoring","离线评分真值 ZIP")
    ))
    related = [
        ("bounds/index.html","区域界与暗电流"), ("bounds/roi_mimo.html","MIMO 峰值区间"),
        ("imaging/index.html","盲成像"), ("ports/index.html","互耦阵列与实测"),
        ("observations.html","观测与材料"),
    ]
    context_links = " ".join(local_link(page.target,dest,name) for dest,name in related)
    if page.target in ("index.html","reproduce.html","handoff.html"):
        branch = quote(TAG,safe="-_.~")
        context_links += f' <a href="{PROJECT}/tree/{branch}/research/z_alg_optical_inertia/blind_uq_mimo_theory">GitHub 源文件</a>' 
    sidebar = ''
    if toc:
        sidebar = f'<aside class="sidebar"><details open><summary>本页目录</summary><nav aria-label="本页目录">{toc}</nav></details></aside>'
    layout = "layout" if sidebar else "layout no-toc"
    return f'''<!doctype html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="theme-color" content="#fbfaf6">
<title>{sep}{title_plain}</title>
<style>{CSS}</style>
{MATHJAX if 'math-inline' in article or 'class="formula"' in article else ''}
</head>
<body>
<a class="skip" href="#article">跳到正文</a>
<header class="site-header"><div class="top">
<a class="brand" href="{WEBSITE}">Maxwell 逆散射研究</a>
<nav class="site-nav" aria-label="研究导航">{nav}<a href="{RELEASE}">GitHub 发布</a></nav>
</div></header>
<main>
<div class="hero"><div class="rule" aria-hidden="true"></div>
<p class="kicker">{label} · 稳态盲成像与有限端口</p>
<h1{f' id="{escape(heading_id,quote=True)}"' if heading_id else ''}>{title}</h1>
<div class="hero-links">{downloads}</div></div>
<div class="{layout}">{sidebar}
<article class="article" id="article">
{article}
<nav class="page-end" aria-label="配套报告">{context_links}</nav>
</article></div>
</main>
<footer class="footer"><div class="footer-inner">
<span>固定 Green 离散模型 · 2026 年 11 月</span><a href="{RELEASE}">版本与完整归档</a>
</div></footer>
</body></html>
'''


def review():
    source = ROOT/REVIEW_PAGE.source
    if source.is_file():
        return source.read_text()
    notes = ["# 独立评审与修改"]
    for number, name in (("02","第二次评审"),("03","第三次评审")):
        path = f"neon_{number}.md"
        if (source.parent/path).is_file():
            notes.append(f"[{name}]({path})")
    return "\n\n".join(notes)


def main():
    args = argparse.ArgumentParser(description=__doc__)
    args.add_argument("--raw-relative",action="store_true",
                      help="use direct relative links to existing source and data files")
    options = args.parse_args()
    sources = {page.source:(ROOT/page.source).read_text() for page in PAGES
               if (ROOT/page.source).is_file()}
    sources[REVIEW_PAGE.source] = review()
    if PAGES[0].source not in sources:
        sources[PAGES[0].source] = (
            "# 盲区域的相容区间与互耦阵列\n\n"
            "[理论与证明](THEORY_REFINEMENT.md) · [区域界](bounds/RESULTS_ZH.md) · "
            "[MIMO 区域峰值](bounds/roi_mimo_RESULTS_ZH.md) · [盲成像](imaging/RESULTS_ZH.md) · "
            "[互耦阵列与实测](ports/PORTS_ZH.md)\n"
        )
    if "REPRODUCE.md" not in sources:
        sources["REPRODUCE.md"] = (
            "# 实验复现\n\n[观测清单](OBSERVATIONS.md) · [区域界](bounds/RESULTS_ZH.md#复现) · "
            "[MIMO 区域峰值](bounds/roi_mimo_RESULTS_ZH.md#复现) · "
            "[盲成像实现](imaging/README_ZH.md) · [依赖列表](requirements.txt)\n"
        )
    pages = [p for p in [*PAGES,REVIEW_PAGE] if p.source in sources and sources[p.source].strip()]
    mapping = {p.source:p.target for p in pages}
    for page in pages:
        title, heading_id, content, toc = render_markdown(
            sources[page.source],page,mapping,options.raw_relative)
        target = ROOT/page.target
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_text(htmldoc(page,title,heading_id,content,toc),encoding="utf-8")
        print(page.target)


if __name__ == "__main__":
    main()

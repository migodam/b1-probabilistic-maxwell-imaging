"""Render the Chinese research documents as standalone, linked HTML pages."""

import argparse
from html import escape
from pathlib import Path, PurePosixPath
import posixpath
import re
from urllib.parse import quote, unquote, urlsplit, urlunsplit

try:
    from markdown_it import MarkdownIt
except ImportError as exc:
    raise SystemExit("Install markdown-it-py or run: uv run --with markdown-it-py==4.2.0 python build_site.py") from exc


ROOT = Path(__file__).resolve().parent
PROJECT = "https://migodam.github.io/b1-probabilistic-maxwell-imaging/"
RESEARCH_ZIP = (
    "https://github.com/migodam/b1-probabilistic-maxwell-imaging/"
    "releases/download/z-alg-secants-20261122/Z_ALG_Secants_2026-11-22.zip"
)
PAGES = (
    ("RESEARCH_REPORT_ZH.md", "index.html", "研究报告"),
    ("MATHEMATICAL_THEORY.md", "theory.html", "数学论证"),
    ("EXPERIMENTAL_RESULTS.md", "experiments.html", "实验结果"),
    ("PRIOR_ART.md", "literature.html", "相关文献"),
    ("reviews/RESPONSE_ZH.md", "review.html", "Neon 的建议"),
    ("NEXT_PRO_HANDOFF.md", "handoff.html", "研究交接"),
)
DETAILS = (
    ("highdim/reflection/README_ZH.md", "symmetry.html", "镜像与停测"),
    ("highdim/reflection/PORT_GAINS_ZH.md", "calibration.html", "端口标定"),
    ("highdim/sector/README_ZH.md", "contrast.html", "材料峰值"),
    ("highdim/sector_geometry/README_ZH.md", "hidden-currents.html", "接收不到的电流"),
    ("highdim/sector_geometry/DESIGN_ZH.md", "receiver-design.html", "增加接收器"),
    ("highdim/sector_inertia/README_ZH.md", "multiple-illuminations.html", "多源共同约束"),
    ("multiplicity/aperture/README.md", "aperture.html", "接收角与频率"),
    ("REPRODUCE.md", "reproduce.html", "复现命令"),
    ("TASK_RESUME.md", "scope.html", "原任务与分工"),
)


def line(state, number):
    start = state.bMarks[number] + state.tShift[number]
    return state.src[start:state.eMarks[number]]


def block_math(state, start, limit, silent):
    opening = line(state, start).rstrip()
    delimiters = [(r"\[", r"\]"), ("$$", "$$")]
    for left, right in delimiters:
        if opening.startswith(left):
            break
    else:
        return False

    first = opening[len(left):]
    if right in first and not first.endswith(right):
        return False
    if first.endswith(right):
        pieces, finish = [first[:-len(right)]], start
    else:
        pieces = [first]
        for finish in range(start + 1, limit):
            candidate = line(state, finish).rstrip()
            if candidate.endswith(right):
                pieces.append(candidate[:-len(right)])
                break
            pieces.append(candidate)
        else:
            return False
    if silent:
        return True
    token = state.push("math_block", "div", 0)
    token.block = True
    token.content = "\n".join(pieces).strip("\n")
    state.line = finish + 1
    return True


def inline_math(state, silent):
    start = state.pos
    sample = state.src[start:state.posMax]
    for left, right in ((r"\(", r"\)"), (r"\[", r"\]"), ("$$", "$$"), ("$", "$")):
        if sample.startswith(left):
            break
    else:
        return False

    body_start = start + len(left)
    end = state.src.find(right, body_start, state.posMax)
    while end >= 0 and right in ("$", "$$"):
        nearby = state.src[body_start:end]
        if (len(nearby) - len(nearby.rstrip("\\"))) % 2 == 0:
            break
        end = state.src.find(right, end + len(right), state.posMax)
    if end <= body_start or "\n" in state.src[body_start:end]:
        return False
    if left == "$" and (state.src[body_start].isspace() or state.src[end - 1].isspace()):
        return False
    if not silent:
        token = state.push("math_inline", "span", 0)
        token.content = state.src[body_start:end]
        token.meta["display"] = left in (r"\[", "$$")
    state.pos = end + len(right)
    return True


def math_display(tokens, i, options, env):
    return f'<div class="math-display" tabindex="0">\\[{escape(tokens[i].content)}\\]</div>\n'


def math_inline(tokens, i, options, env):
    token = tokens[i]
    start, end = (r"\[", r"\]") if token.meta["display"] else (r"\(", r"\)")
    return f'<span class="math-inline">{start}{escape(token.content)}{end}</span>'


def rewrite(url, env):
    parts = urlsplit(url)
    if parts.scheme or parts.netloc or not parts.path:
        return url
    path = unquote(parts.path)
    if path.startswith("/"):
        try:
            path = Path(path).relative_to(env["root"]).as_posix()
        except ValueError:
            return url
    else:
        parent = PurePosixPath(env["source"]).parent
        path = posixpath.normpath(posixpath.join(parent.as_posix(), path))
    target = env["routes"].get(path, path)
    if target == path and parts.path.endswith("/") and not target.endswith("/"):
        target += "/"
    return urlunsplit(parts._replace(path=quote(target, safe="/-._~!$&'()*+,;=:@")))


def heading_anchors(state):
    counts = {}
    for token, inline in zip(state.tokens, state.tokens[1:]):
        if token.type != "heading_open":
            continue
        title = inline.content.casefold()
        base = re.sub(r"[^\w\s-]", "", title, flags=re.UNICODE)
        base = re.sub(r"\s+", "-", base).strip("-") or "section"
        n = counts.get(base, 0)
        token.attrSet("id", base if n == 0 else f"{base}-{n}")
        counts[base] = n


def linked_image(renderer, tokens, i, options, env):
    token = tokens[i]
    target = rewrite(token.attrGet("src"), env)
    token.attrSet("src", target)
    token.attrSet("loading", "lazy")
    token.attrSet("decoding", "async")
    markup = renderer.image(tokens, i, options, env)
    inside_link = sum(t.nesting for t in tokens[:i] if t.tag == "a") > 0
    if inside_link or not target:
        return markup
    return f'<a class="image-original" href="{escape(target, quote=True)}">{markup}</a>'


def parser():
    md = MarkdownIt("commonmark", {"html": False}).enable("table")
    md.block.ruler.before("fence", "math_block", block_math, {"alt": ["paragraph", "blockquote"]})
    md.inline.ruler.before("escape", "math_inline", inline_math)
    md.core.ruler.after("inline", "heading_anchors", heading_anchors)
    md.renderer.rules["math_block"] = math_display
    md.renderer.rules["math_inline"] = math_inline
    original_render = md.renderer.renderToken

    def link(tokens, i, options, env):
        token = tokens[i]
        token.attrSet("href", rewrite(token.attrGet("href"), env))
        return original_render(tokens, i, options, env)

    def table_open(tokens, i, options, env):
        return ('<div class="table-scroll" tabindex="0" role="region" aria-label="表格，可横向滚动">\n'
                + original_render(tokens, i, options, env))

    def table_close(tokens, i, options, env):
        return original_render(tokens, i, options, env) + "</div>\n"

    md.renderer.rules["link_open"] = link
    md.renderer.rules["table_open"] = table_open
    md.renderer.rules["table_close"] = table_close
    md.add_render_rule("image", linked_image)
    return md


def html_page(title, content, this_page, nav, source, zip_link, zip_label):
    navigation = []
    for href, label in nav:
        active = ' aria-current="page"' if href == this_page else ""
        navigation.append(f'<a href="{href}"{active}>{label}</a>')
    nav_html = "".join(navigation)
    return f"""<!doctype html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(title)} · Z-ALG</title>
  <link rel="stylesheet" href="research.css">
  <script>
    window.MathJax = {{
      tex: {{inlineMath: [['\\\\(', '\\\\)']], displayMath: [['\\\\[', '\\\\]']]}},
      options: {{skipHtmlTags: ['script','noscript','style','textarea','pre','code']}}
    }};
  </script>
  <script async src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js"></script>
</head>
<body>
  <a class="skip-link" href="#content">跳到正文</a>
  <header class="site-header">
    <div class="site-shell">
      <div class="project-links">
        <a href="{PROJECT}">B1 总项目</a>
        <a href="{PROJECT}research/">研究主页</a>
        <a class="download" href="{escape(zip_link, quote=True)}">{escape(zip_label)}</a>
      </div>
      <nav aria-label="本轮研究页面">{nav_html}</nav>
    </div>
  </header>
  <main class="site-shell" id="content">
    <article>{content}</article>
    <footer><a href="{escape(source, quote=True)}">原始 Markdown</a></footer>
  </main>
</body>
</html>
"""


def render_documents(root, handoff_zip=None, require_all=False):
    pages = PAGES + DETAILS
    present = [row for row in pages if (root / row[0]).is_file()]
    missing = [row[0] for row in pages if row not in present]
    if require_all and missing:
        raise SystemExit("缺少文稿：" + "、".join(missing))
    if not present:
        print("目录中暂无研究文稿。")
        return []
    if missing:
        print("跳过尚未写入的文稿：" + "、".join(missing))
    local_archives = sorted(root.glob("*.zip"))
    if handoff_zip:
        zip_link, zip_label = handoff_zip, "本轮交接 ZIP"
    elif len(local_archives) == 1:
        zip_link, zip_label = quote(local_archives[0].name), "本轮交接 ZIP"
    else:
        zip_link, zip_label = RESEARCH_ZIP, "本轮交接 ZIP"
    routes = {source: out for source, out, _ in present}
    nav = [(out, label) for source, out, label in PAGES if source in routes]
    renderer = parser()
    files = []
    for source, out, nav_label in present:
        text = (root / source).read_text(encoding="utf-8")
        heading = re.search(r"^#\s+(.+?)\s*$", text, flags=re.MULTILINE)
        title = heading.group(1) if heading else nav_label
        content = renderer.render(text, {"source": source, "root": root, "routes": routes})
        destination = root / out
        destination.write_text(html_page(title, content, out, nav, source, zip_link, zip_label), encoding="utf-8")
        files.append(destination)
    print("已生成 " + "、".join(file.name for file in files))
    return files


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--handoff-zip", metavar="URL", help="Use this release URL or page-relative ZIP path.")
    cli.add_argument("--require-all", action="store_true", help="Report a missing document instead of skipping it.")
    args = cli.parse_args()
    render_documents(ROOT, args.handoff_zip, args.require_all)

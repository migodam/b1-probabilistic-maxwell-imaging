#!/usr/bin/env python3
"""Render supplement prose and existing figures only; no science or publication."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import time
from datetime import datetime, timezone
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
# Reuse the frozen parent's renderer/style/math/download checks; never modify it.
PARENT_RENDERER = ROOT.parent / 'src/build_site.py'
spec = importlib.util.spec_from_file_location('r14_public_renderer', PARENT_RENDERER)
base = importlib.util.module_from_spec(spec)
exec(compile(PARENT_RENDERER.read_text(encoding='utf8'), str(PARENT_RENDERER), 'exec'), base.__dict__)
DOCS = [
    ('RESEARCH_REPORT_ZH.md', 'index.html', '补充研究报告'),
    ('THEORY_CLARIFICATION_ZH.md', 'theory_clarification_zh.html', '公式与理论澄清'),
    ('EXPERIMENT_CONTRACT.md', 'experiment_contract.html', '冻结实验合同'),
    ('REPRODUCE.md', 'reproduce.html', '复现说明'),
    ('COST_LEDGER.md', 'cost_ledger.html', '费用总账'),
    ('FAILURES.md', 'failures.html', '失败与范围'),
    ('CHAT_HANDOFF.md', 'chat_handoff.html', '研究交接'),
]
BLOCKED = {'sealed_truth', 'private_operational', 'results', 'predictions',
           'chat_readback', 'secrets', 'credentials'}


def allowed(path, root):
    try:
        rel = path.relative_to(root)
        return not any(part.lower() in BLOCKED for part in rel.parts) and base.safe_file(path, root)
    except ValueError:
        return False


class Renderer(base.Renderer):
    def target(self, url, source_base, output_base, image=False):
        u = urlsplit(url)
        if not u.scheme and u.path:
            p = (source_base / unquote(u.path)).resolve()
            if not allowed(p, self.root):
                self.suppressed.append(dict(url=url, reason='OUTSIDE_OR_NONPUBLIC_SUPPLEMENT'))
                return None
        return super().target(url, source_base, output_base, image)


def direct_files(folder, root):
    """Inspect only explicit small public directories, with no recursive walk."""
    if not folder.is_dir() or folder.is_symlink():
        return []
    return sorted(p for p in folder.iterdir() if allowed(p, root))


def inventory(root):
    # analyze_results writes analysis/<run>/figures and CSVs/index.html.
    # Also accept analysis/figures or uppercase FIGURES, without traversing data.
    analysis = root / 'analysis'
    locations = [analysis]
    if analysis.is_dir() and not analysis.is_symlink():
        locations += sorted(p for p in analysis.iterdir()
                            if p.is_dir() and not p.is_symlink()
                            and p.name.lower() not in BLOCKED)
    images = set()
    tables = set()
    source_indexes = set()
    pdfs = set()
    figure_dirs = [p for parent in [root] + locations
                   if parent.is_dir() and not parent.is_symlink()
                   for p in parent.iterdir()
                   if p.name.lower() == 'figures' and p.is_dir() and not p.is_symlink()]
    # Actual diagnostic figures are also direct children of analysis/<run>.
    figure_dirs += locations
    for folder in figure_dirs:
        for p in direct_files(folder, root):
            if p.suffix.lower() in base.IMAGE_EXT:
                images.add(p)
            elif p.suffix.lower() == '.pdf':
                pdfs.add(p)
    for folder in locations:
        for p in direct_files(folder, root):
            if p.suffix.lower() == '.csv' and p.name != 'PAID_HISTORY.csv':
                tables.add(p)
            elif p.name == 'index.html':
                # Original analyzer indexes link the full paid-history table,
                # archived in release rather than duplicated into static Git.
                pass
    # Same-stem PDF is a download on the image, never a second gallery card.
    standalone_pdfs = pdfs - {p.with_suffix('.pdf') for p in images}
    return dict(images=sorted(images), tables=sorted(tables),
                analysis_indexes=sorted(source_indexes), standalone_pdfs=sorted(standalone_pdfs))


def page(title, body, docs):
    links = [('index.html', '研究报告'), ('galleries.html', '全部实际图'), ('evidence.html', '数据与下载')]
    links += [(target.name, label) for source, label, target in docs if target.name != 'index.html']
    return ('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1"><title>'
            + base.esc(title) + '</title><style>' + base.CSS + '</style>' + base.MATHJAX_HEAD
            + '</head><body><main><header><h1>' + base.esc(title) + '</h1><nav>'
            + ''.join('<a href="' + href + '">' + base.esc(label) + '</a>' for href, label in links)
            + '</nav></header>' + body
            + '<footer>R14 · supplement_observable19。页面只呈现已有文字、图与文件；建站不新增科学结果。'
              '公式加载MathJax CDN；无JavaScript时保留LaTex原文。</footer></main></body></html>\n')


def audit(pages, root):
    result = base.audit_links(pages, root)
    for p in pages:
        parser = base.Links(); parser.feed(p.read_text(encoding='utf8'))
        for url in parser.links:
            u = urlsplit(url)
            if not u.scheme and u.path and not allowed((p.parent / unquote(u.path)).resolve(), root):
                result['forbidden_private_or_outside_links'].append(dict(page=p.name, url=url))
    return result


def build(root, check_only=False):
    root = root.resolve()
    cpu, wall = time.process_time(), time.perf_counter()
    docs = [(root / src, label, root / dest) for src, dest, label in DOCS]
    present = [src.name for src, _, _ in docs if allowed(src, root)]
    missing = [src.name for src, _, _ in docs if not allowed(src, root)]
    files = inventory(root)
    receipt = dict(snapshot_utc=datetime.now(timezone.utc).isoformat(), root=str(root),
                   parent_renderer=str(PARENT_RENDERER), present_documents=present, missing_documents=missing,
                   inventory={k:[str(p.relative_to(root)) for p in v] for k,v in files.items()},
                   science_executed=False, network_accessed=False, publication_performed=False)
    if check_only:
        receipt['status'] = 'INVENTORY_ONLY_NO_SITE_WRITTEN'
    else:
        renderer = Renderer(root, {src.resolve():dest for src, _, dest in docs if allowed(src, root)})
        pages = []
        for src, label, dest in docs:
            if allowed(src, root):
                body = '<article>' + renderer.markdown(src.read_text(encoding='utf8'), src.parent, dest.parent) + '</article>'
                body += '<p><a href="' + base.relurl(src, root) + '">原始Markdown</a></p>'
            else:
                body = '<p class="empty">' + base.esc(src.name) + ' 尚未提供；此页不表示报告或实验已完成。</p>'
            if dest.name == 'index.html':
                body += '<p><a href="galleries.html">查看全部实际图</a> · <a href="evidence.html">查看CSV与已核验下载</a></p>'
            dest.write_text(page(label, body, docs), encoding='utf8'); pages.append(dest)
        gallery = '<p class="note">下列图只来自本补充目录已有文件，保留原图标签与完整对象页面；不按恢复表现挑选。相同名称PDF只作下载。</p><div class="gallery">'
        for p in files['images']:
            href = base.relurl(p, root)
            gallery += '<figure><a href="' + href + '"><img loading="lazy" src="' + href + '" alt="' + base.esc(p.stem) + '"></a><figcaption>' + base.esc(p.stem) + '<span class="path">' + base.esc(p.relative_to(root)) + '</span>'
            pdf = p.with_suffix('.pdf')
            if allowed(pdf, root): gallery += '<a href="' + base.relurl(pdf, root) + '">PDF下载</a>'
            gallery += '</figcaption></figure>'
        gallery += '</div>'
        if not files['images']: gallery += '<p class="empty">实际分析图尚未生成。</p>'
        for p in files['standalone_pdfs']:
            gallery += '<p><a href="' + base.relurl(p, root) + '">' + base.esc(p.name) + '（PDF）</a></p>'
        target = root / 'galleries.html'; target.write_text(page('全部实际图', gallery, docs), encoding='utf8'); pages.append(target)
        download_body, links, skipped = base.downloads(root)
        body = '<h2>已核验公开附件</h2>' + download_body
        for key, label in [('tables','实际分析CSV'),('analysis_indexes','原始分析页面')]:
            body += '<h2>' + label + '</h2>'
            body += ''.join('<p><a href="' + base.relurl(p,root) + '">' + base.esc(p.relative_to(root)) + '</a></p>' for p in files[key]) or '<p class="empty">尚无文件。</p>'
        body += '<h2>本补充研究代码与冻结配置</h2>'
        for folder, suffixes in [(root/'src',{'.py'}),(root/'configs',{'.json'})]:
            for p in direct_files(folder,root):
                if p.suffix in suffixes: body += '<p><a href="' + base.relurl(p,root) + '">' + base.esc(p.relative_to(root)) + '</a></p>'
        target = root/'evidence.html'; target.write_text(page('数据与下载',body,docs),encoding='utf8'); pages.append(target)
        receipt.update(status='BUILT_NOT_PUBLISHED', generated_pages=[p.name for p in pages],
                       verified_downloads=links, skipped_downloads=skipped,
                       suppressed_prose_links=renderer.suppressed, link_audit=audit(pages,root))
    receipt.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    dest=root/'audit/site_implementation'; dest.mkdir(parents=True,exist_ok=True)
    (dest/('INVENTORY_CHECK.json' if check_only else 'BUILD_RECEIPT.json')).write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    if not check_only and (receipt['link_audit']['broken_local_links'] or receipt['link_audit']['forbidden_private_or_outside_links']):
        raise RuntimeError('Site link audit failed; see BUILD_RECEIPT.json')
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=ROOT,help='Supplement study directory only')
    p.add_argument('--check-only',action='store_true',help='Inventory prose/figure filenames; do not write site pages')
    args=p.parse_args()
    print(json.dumps(build(args.root,args.check_only),ensure_ascii=False,indent=2))


if __name__=='__main__':main()

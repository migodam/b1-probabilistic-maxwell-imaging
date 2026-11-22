#!/usr/bin/env python3
"""Render the submission report, companion pages, and links from the B1 index."""
import argparse
import csv
import json
from pathlib import Path
import re
import runpy
import sys
from urllib.parse import quote

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
ROOT = STUDY.parents[1]
sys.path.insert(0, str(STUDY))
from build_site import DETAILS as LEGACY_DETAILS, PAGES as LEGACY_PAGES, RESEARCH_ZIP, html_page, parser

DEFAULT_TAG = 'z-alg-optical-inertia-submission-20261122'
DATA_NAME = 'Z_ALG_Optical_Inertia_Submission_2026-11-22.zip'
LATEX_NAME = 'Z_ALG_Optical_Inertia_LaTeX_2026-11-22.zip'
MANUSCRIPT = 'manuscript_ip.pdf'
SITE_PATH = 'research/z_alg_optical_inertia/submission/'
DEFAULT_CHECKOUT = ROOT / 'research/delegated/z_alg_secant_publication'
PAGES = (
    ('RESEARCH_REPORT_ZH.md', 'index.html', '定稿报告'),
    ('single/REPORT_ZH.md', 'coherent_comparison.html', '单源与联合对照'),
    ('curve/README_ZH.md', 'noise_curve.html', '原生全环噪声曲线'),
    ('REPRODUCE.md', 'reproduce.html', '复现导航'),
    ('NEXT_PRO_HANDOFF.md', 'handoff.html', '交接'),
)
DETAILS = (
    ('diagnostic/RESULTS_EN.md', 'diagnostic.html', 'Gauss–Newton 对照'),
    ('highcaps/REPORT.md', 'highcaps.html', '高阈值证明'),
    ('typeset/README.md', 'typeset.html', '排版与投稿要求'),
)
MARKDOWN_NOTICE = (
    '[英文投稿稿与定稿报告](submission/index.html) · '
    '[本次复现入口](submission/reproduce.html)\n'
)
LEGACY_END = '<!-- /Z_ALG_OPTICAL_INERTIA_PUBLICATION -->'
RELEASE_ROOT = 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/'
COHORT = 'z_alg_optical_fullring_paired30'
OPERATOR = 'z_alg_optical_fullring_native'


def release_url(tag, file):
    return f'{RELEASE_ROOT}{quote(tag, safe="")}/{file}'


def render_source(root, source, out, nav, detail_nav, routes, asset, release_tag=None, **options):
    text = (root / source).read_text(encoding='utf-8')
    if release_tag and release_tag != DEFAULT_TAG:
        text = text.replace(release_url(DEFAULT_TAG, ''), release_url(release_tag, ''))
    match = re.search(r'^#\s+(.+?)\s*$', text, flags=re.MULTILINE)
    title = match.group(1) if match else nav[0][1]
    content = parser().render(text, {'source': source, 'root': root, 'routes': routes})
    (root / out).write_text(
        html_page(title, content, out, nav, detail_nav, source, asset, **options),
        encoding='utf-8',
    )


def render_submission(tag, *, preview=False):
    if not preview:
        required = [HERE / file for file, _, _ in PAGES + DETAILS]
        required.append(HERE / MANUSCRIPT)
        missing = [path.name for path in required if not path.is_file()]
        if missing:
            raise FileNotFoundError('定稿网页缺少 ' + '、'.join(missing))
    pages = [row for row in PAGES if (HERE / row[0]).is_file()]
    details = [row for row in DETAILS if (HERE / row[0]).is_file()]
    routes = {source: out for source, out, _ in pages + details}
    routes.update({f'../{source}': f'../{out}' for source, out, _ in LEGACY_PAGES + LEGACY_DETAILS})
    nav = [(out, label) for _, out, label in pages]
    extra = [('../index.html', '上一轮光学惯性报告')]
    if not preview:
        extra.extend(((MANUSCRIPT, '英文论文 PDF'), (release_url(tag, LATEX_NAME), 'LaTeX 源包')))
    for source, out, _ in pages + details:
        render_source(HERE, source, out, nav, [(p, l) for _, p, l in details], routes,
                      release_url(tag, DATA_NAME), release_tag=tag,
                      stylesheet='../research.css', extra_links=extra,
                      english=out in ('diagnostic.html', 'highcaps.html', 'typeset.html'))


def render_legacy_entry():
    source, out, _ = LEGACY_PAGES[0]
    markdown = STUDY / source
    text = markdown.read_text(encoding='utf-8')
    if MARKDOWN_NOTICE not in text:
        start, rest = text.split('\n\n', 1)
        markdown.write_text(start + '\n\n' + MARKDOWN_NOTICE + '\n' + rest, encoding='utf-8')
    routes = {md: dest for md, dest, _ in LEGACY_PAGES + LEGACY_DETAILS}
    nav = [(dest, label) for _, dest, label in LEGACY_PAGES]
    nav.append(('submission/index.html', '英文定稿与复现'))
    details = [(dest, label) for _, dest, label in LEGACY_DETAILS]
    render_source(STUDY, source, out, nav, details, routes, RESEARCH_ZIP)


def insert_once(original, start, end, anchor, block, gap=""):
    if start in original:
        before, remainder = original.split(start, 1)
        _, after = remainder.split(end, 1)
        return before + block + after
    if original.count(anchor) != 1:
        raise ValueError(f'Cannot locate unique B1 publication marker {anchor}.')
    left, right = original.split(anchor, 1)
    return left + anchor + gap + block + right


def stage_indexes(checkout):
    start = '<!-- Z_ALG_OPTICAL_INERTIA_SUBMISSION -->'
    end = '<!-- /Z_ALG_OPTICAL_INERTIA_SUBMISSION -->'
    title = 'Z-ALG · 有限孔径散射数据的介电峰值认证'
    desc = '英文投稿稿、原生 96² 分块对照、配对噪声曲线与严格证书，保留定位和实测的失败范围。'
    for file, prefix in (('index.html', SITE_PATH), ('research/index.html', SITE_PATH.removeprefix('research/'))):
        block = (
            f'{start}<section class="card"><h2>{title}</h2><p>{desc}</p><p>'
            f'<a href="{prefix}">中文定稿报告</a> · <a href="{prefix}{MANUSCRIPT}">英文论文 PDF</a>'
            f' · <a href="{prefix}reproduce.html">复现与下载</a> · <a href="{prefix}handoff.html">研究交接</a>'
            f'</p></section>{end}'
        )
        target = checkout / file
        target.write_text(insert_once(target.read_text(encoding='utf-8'), start, end, LEGACY_END, block),
                          encoding='utf-8')
    for file, prefix in (('README.md', SITE_PATH), ('research/README.md', SITE_PATH.removeprefix('research/'))):
        block = (
            f'{start}\n## {title}\n\n{desc}\n\n'
            f'[中文定稿报告]({prefix}index.html) · [英文论文 PDF]({prefix}{MANUSCRIPT}) · '
            f'[复现与下载]({prefix}reproduce.html) · [研究交接]({prefix}handoff.html)\n{end}'
        )
        target = checkout / file
        text = target.read_text(encoding='utf-8')
        anchor = LEGACY_END if LEGACY_END in text else '<!-- /Z_ALG_SECANT_TASKS_PUBLICATION -->'
        target.write_text(insert_once(text, start, end, anchor, block, gap='\n\n'), encoding='utf-8')


def registries(directory, filename, keys, rows):
    path = directory / filename
    with path.open(encoding='utf-8', newline='') as stream:
        existing = {tuple(row[key] for key in keys) for row in csv.DictReader(stream)}
    fresh = [row for row in rows if tuple(row[:len(keys)]) not in existing]
    if not fresh:
        return
    with path.open('a', encoding='utf-8', newline='') as stream:
        csv.writer(stream, lineterminator='\n').writerows(fresh)


def stage_handbook(checkout):
    directory = checkout / 'research/experiment_handbook'
    study = SITE_PATH.rstrip('/')
    website = f'https://migodam.github.io/b1-probabilistic-maxwell-imaging/{SITE_PATH}'
    registries(directory, 'COHORTS.csv', ('dataset',), [(
        COHORT, '项目自建；沿用已公开的 Z-ALG 光学惯性六件材料', '二维仿真',
        '6 件独立材料；30 份同噪声方向配对全环观测',
        '96², 400 MHz, 10 平面波 × 20 接收；0.5%、1%、2%、5%、10%；2% 原样复用，10% 异于独立旧抽样。',
        'B1 网页公开观测或定稿 Release 的完整科研 ZIP',
        '公开 Y、原生 Green FFT 频谱、B、S 和各观测扩大 Frobenius 半径',
        '既有材料已曝光；配对噪声曲线、联合与单源松弛对照',
        study, study, 'curve/README_ZH.md', 'curve/observation_audit.json', website+'noise_curve.html',
    )])
    manifest = HERE / 'curve/observations/manifest.json'
    observations = json.loads(manifest.read_text(encoding='utf-8'))['observations']
    registries(directory, 'OBSERVATIONS.csv', ('dataset', 'object_id', 'observation_id'), [
        (COHORT, row['object_id'], row['id'],
         f"native96_PW360_10x20_pair_{row['noise_pct']:g}pct",
         'curve/observations/'+row['file'],
         f"{row['object_id']}__pw360_10x20__n02；2% 逐位复用，其余沿该噪声方向重标度")
        for row in observations
    ])
    verify = study + '/curve/verify.py'
    registries(directory, 'OPERATORS.csv', ('operator',), [(
        OPERATOR, 'curve/verify.py --cap --arm --out', 'CPU',
        '原生 96² 标量 TM 的共享照明不定惯性残差界；支持全局权或四象限空间权和同球配对比较。',
        '逐例保留原生 FFT、区间几何、源谱回放和未越线原因；高阈值单独提供不读 Y 的接收核障碍。',
        verify, verify,
    )])
    source = directory / 'README.md'
    original = source.read_text(encoding='utf-8')
    entry = (
        f'|全环噪声变动时还能排除什么峰值？| `{COHORT}` 的 6 件材料和 30 份配对复观测'
        '|逐件比对标量、空间乘子、噪声球和高阈值限制。'
        '[配对曲线](../z_alg_optical_inertia/submission/noise_curve.html)|'
    )
    if entry not in original:
        anchor = '|能否用较少发射排除镜像材料？|'
        original = original.replace(anchor, entry + '\n' + anchor, 1)
    for key, label in (('COHORTS.csv', '个数据集、子集及派生视图'), ('OPERATORS.csv', '个可调用的算子入口')):
        with (directory / key).open(encoding='utf-8', newline='') as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        original = re.sub(r'\d+(?= ' + label + ')', str(count), original, count=1)
    source.write_text(original, encoding='utf-8')
    renderer = runpy.run_path(str(directory / 'render_site.py'))
    body = renderer['markdown'].markdown(original, extensions=['tables', 'fenced_code', 'toc'])
    body = re.sub(r'href="([^"]*)"', renderer['internal_link'], body)
    body = body.replace('<table>', '<div class="wide"><table>').replace('</table>', '</table></div>')
    renderer['write_page']('index.html', '实验手册', body, source.name)
    renderer['registry']()


if __name__ == '__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--release-tag', default=DEFAULT_TAG)
    cli.add_argument('--preview', action='store_true', help='Render pages currently available without altering site navigation.')
    args = cli.parse_args()
    render_submission(args.release_tag, preview=args.preview)
    if not args.preview:
        render_legacy_entry()

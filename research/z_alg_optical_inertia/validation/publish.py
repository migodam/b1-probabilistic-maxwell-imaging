#!/usr/bin/env python3
"""Assemble the revised manuscript, validation pages, and B1 navigation locally."""
import argparse
import csv
import json
from pathlib import Path, PurePosixPath
import re
import runpy
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
SUBMISSION = STUDY / 'submission'
ROOT = STUDY.parents[1]
sys.path.insert(0, str(SUBMISSION))
from build_release import EXCLUDE, EXCLUDE_SUFFIX, SOURCE_FILE, flat_source
from publish_pages import (
    DETAILS as OLD_DETAILS, PAGES as OLD_PAGES, LEGACY_DETAILS, LEGACY_PAGES,
    DEFAULT_CHECKOUT, insert_once, registries, release_url, render_source,
)

TAG = 'z-alg-optical-inertia-validation-v1'
DATA = 'Z_ALG_Optical_Inertia_Validation.zip'
LATEX = 'Z_ALG_Optical_Inertia_LaTeX.zip'
PDF = 'manuscript_ip.pdf'
SITE = 'research/z_alg_optical_inertia/validation/'
BASE = ROOT / '.artifacts/z_alg_optical_submission_release/Z_ALG_Optical_Inertia_Submission_2026-11-22.zip'
OUTPUT = ROOT / '.artifacts/z_alg_optical_validation_release'
EIGOR = 'https://home.cc.umanitoba.ca/~lovetrij/EMILab/research/EIGOR/Tomographic_Dataset_1.zip'
SRAP = 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/srap-cie-array-study-20261121/SRAP_CIE_arrays_20261121.tar.gz'
PAGES = (
    ('RESEARCH_REPORT_ZH.md', 'index.html', '独立验证'),
    ('srap/README_ZH.md', 'srap.html', 'SRAP 22 件'),
    ('max_struct/REPORT_ZH.md', 'maxstruct.html', 'MAX-STRUCT 7 件'),
    ('manitoba/REPORT_ZH.md', 'manitoba.html', 'Manitoba 实测'),
    ('reconstructions/README_ZH.md', 'endpoints.html', '历史重建端点'),
    ('transfer/README_ZH.md', 'transfer.html', '校准转移'),
    ('REPRODUCE.md', 'reproduce.html', '复现与下载'),
    ('NEXT_PRO_HANDOFF.md', 'handoff.html', '研究交接'),
)
OLD_FIGURES = (
    STUDY/'figures/fine_peaks.pdf', SUBMISSION/'curve/full_ring_noise.pdf',
    SUBMISSION/'single/partition_comparison.pdf', SUBMISSION/'diagnostic/figure_gn_diagnostic.pdf',
)
PUBLIC_NUMPY = {Path('srap/public_observations.npz'), Path('max_struct/public_inputs.npz')}
PUBLIC_SUFFIX = {'.md', '.html', '.py', '.tex', '.bib', '.cls', '.json', '.csv', '.png', '.pdf', '.svg'}
PRIVATE_PARTS = {'sealed_truth', 'private', 'raw', 'raw_eigor', '.artifacts', *EXCLUDE}
MANITOBA_COHORT = 'manitoba_eigor2010_dataset1'
MANITOBA_AUDIT = HERE / 'manitoba/physical_audit.json'
VALIDATION_DEPENDENCIES = f'''
## 独立验证修订版

`{SITE}` 保留三批验证的公开复算代码、数值表、图、论文和交接。
本包沿用首版研究附件的已有依赖。

- SRAP 的证书复算读取本包的 `srap/public_observations.npz` 和首版几何、乘子。
  独立评分或重建输入需要原 SRAP 回放与已公开参考，从 [原 SRAP Release]({SRAP}) 取得。
  若某历史依赖未随首版完整包提供，按原 Release 的相对路径放回工作区。
- MAX-STRUCT 的证书复算读取本包的 `max_struct/public_inputs.npz`。历史封存材料沿用原权限，
  本包提供原合成全波可行性和历史重建的逐例标量。
- Manitoba 重算从 [采集方原 ZIP]({EIGOR}) 开始。原测量、原始文件与制备相干数组
  由采集方档案提供或在本地按脚本生成。本包提供派生标量、图及从合法原始输入重算的程序。
'''


def missing_files():
    required = [HERE / src for src, _, _ in PAGES]
    required += [HERE / SOURCE_FILE, MANITOBA_AUDIT, *OLD_FIGURES, SUBMISSION/'typeset/iopjournal.cls']
    return [str(path.relative_to(STUDY)) for path in required if not path.is_file()]


def render_pages(preview=False):
    available = [(src, page, label) for src, page, label in PAGES if (HERE / src).is_file()]
    routes = {src: page for src, page, _ in available}
    routes.update({f'../{src}': f'../{page}' for src, page, _ in LEGACY_PAGES + LEGACY_DETAILS})
    routes.update({f'../submission/{src}': f'../submission/{page}' for src, page, _ in OLD_PAGES + OLD_DETAILS})
    nav = [(page, label) for _, page, label in available]
    extras = [('../submission/index.html', '首版英文投稿稿')]
    if not preview:
        extras.extend(((PDF, '修订稿 PDF'), (release_url(TAG, LATEX), 'LaTeX 源包')))
    for src, page, _ in available:
        render_source(HERE, src, page, nav, (), routes, release_url(TAG, DATA),
                      stylesheet='../research.css', extra_links=extras)


def is_public(path):
    relative = path.relative_to(HERE)
    return (not path.is_symlink() and path.is_file() and not PRIVATE_PARTS.intersection(relative.parts)
            and not path.name.endswith(EXCLUDE_SUFFIX)
            and (relative in PUBLIC_NUMPY or path.suffix in PUBLIC_SUFFIX))


def make_latex(output):
    filename = output / LATEX
    pending = filename.with_suffix('.partial')
    with TemporaryDirectory(prefix='flat-validation-', dir=output) as folder:
        stage = Path(folder)
        flat_source(HERE / SOURCE_FILE, stage)
        kept = []
        for fig in OLD_FIGURES:
            dest = stage / fig.name
            if not dest.exists():
                shutil.copy2(fig, dest)
                kept.append(fig.name)
        if kept:
            with (stage/'README.txt').open('a', encoding='utf-8') as stream:
                stream.write('\nFirst-version figures retained in this package\n' + '\n'.join(kept) + '\n')
        source_files = sorted(stage.iterdir())
        subprocess.run(['tectonic', '-X', 'compile', SOURCE_FILE], cwd=stage, check=True)
        with ZipFile(pending, 'w', ZIP_DEFLATED) as bundle:
            for asset in source_files:
                bundle.write(asset, asset.name)
        # The revised PDF comes from the same flattened sources as the source ZIP.
        shutil.copy2(stage / PDF, HERE / PDF)
    pending.replace(filename)
    shutil.copy2(HERE / PDF, output / PDF)


def old_page_notices(checkout):
    start = '<!-- Z_ALG_VALIDATION_LINK -->'
    end = '<!-- /Z_ALG_VALIDATION_LINK -->'
    old_relative = STUDY.relative_to(ROOT)
    notices = {}
    for old_dir, file, target in (
        ('', 'index.html', 'validation/index.html'),
        ('', 'paper.html', 'validation/index.html'),
        ('submission', 'index.html', '../validation/index.html'),
        ('submission', 'handoff.html', '../validation/index.html'),
    ):
        relative = old_relative / old_dir / file
        html = (checkout / relative).read_text(encoding='utf-8')
        opening = re.search(r'<article(?:\s+[^>]*)?>', html)
        if opening is None:
            raise ValueError(f'Missing article container in {relative}')
        banner = (f'{start}<p><a href="{target}">查看加入独立验证的修订版</a>'
                  '。本页保留首版的研究记录。</p>{end}')
        notices[relative.as_posix()] = insert_once(html, start, end, opening.group(0), banner)
    return notices


def make_data(base, output, artifacts, notices):
    filename = output / DATA
    pending = filename.with_suffix('.partial')
    with ZipFile(base) as previous, ZipFile(pending, 'w', ZIP_DEFLATED) as package:
        for name in previous.namelist():
            if (not name.endswith('/') and name not in notices and name != 'DEPENDENCIES.md'
                    and not name.startswith(SITE) and 'sealed_truth' not in PurePosixPath(name).parts
                    and not name.endswith('/r10_twod_0001_offline_reference.npz')):
                package.writestr(name, previous.read(name))
        inherited = previous.read('DEPENDENCIES.md').decode('utf-8')
        package.writestr('DEPENDENCIES.md', inherited.rstrip() + '\n' + VALIDATION_DEPENDENCIES)
        for name, html in notices.items():
            package.writestr(name, html)
        for path in artifacts:
            package.write(path, SITE + path.relative_to(HERE).as_posix())
    pending.replace(filename)


def update_indexes(checkout):
    start = '<!-- Z_ALG_OPTICAL_INERTIA_VALIDATION -->'
    end = '<!-- /Z_ALG_OPTICAL_INERTIA_VALIDATION -->'
    anchor = '<!-- /Z_ALG_OPTICAL_INERTIA_SUBMISSION -->'
    title = 'Z-ALG 光学惯性 · 独立验证修订版'
    description = '既有 SRAP、MAX-STRUCT 材料迁移，Manitoba 实测误差、校准转移及各自适用范围。'
    for directory, href in ((checkout, SITE), (checkout/'research', SITE.removeprefix('research/'))):
        pages = (
            ('index.html',
             f'{start}<section class="card"><h2>{title}</h2><p>{description}</p><p>'
             f'<a href="{href}">中文验证报告</a> · <a href="{href}{PDF}">修订英文稿 PDF</a> · '
             f'<a href="{href}reproduce.html">复现与下载</a> · <a href="{href}handoff.html">研究交接</a>'
             f'</p></section>{end}', ''),
            ('README.md',
             f'{start}\n## {title}\n\n{description}\n\n'
             f'[中文验证报告]({href}index.html) · [修订英文稿 PDF]({href}{PDF}) · '
             f'[复现与下载]({href}reproduce.html) · [研究交接]({href}handoff.html)\n{end}', '\n\n'),
        )
        for file, block, gap in pages:
            path = directory / file
            original = path.read_text(encoding='utf-8')
            path.write_text(insert_once(original, start, end, anchor, block, gap=gap), encoding='utf-8')


def update_handbook(checkout):
    directory = checkout / 'research/experiment_handbook'
    public_study = SITE.rstrip('/')
    raw_directory = '.artifacts/z_alg_validation/manitoba'
    from_raw = '../../../' + public_study
    report = json.loads(MANITOBA_AUDIT.read_text(encoding='utf-8'))
    cases = report['cases']
    registries(directory, 'COHORTS.csv', ('dataset',), [(
        MANITOBA_COHORT, 'University of Manitoba，EIGOR Tomographic Dataset 1', '二维等效实测',
        '4 件介质实物、每件 7 频点，共 28 份观测',
        '3–6 GHz 步进 0.5 GHz；24×23 原复端口；固定奇偶源组保留 288/552 端口，40² 标量 TM 参考。',
        '从 University of Manitoba 官方 ZIP 取得原数据；B1 发布派生标量、验算和代码。',
        '原始复场、训练空场和铝柱、留出空场和钢柱、等效 Green / B / S 和固定端口组。',
        '已曝光实物；校准误差敏感性与可转移材料上限，误差条件须单独报告。',
        raw_directory, '', from_raw+'/manitoba/FREEZE.md', from_raw+'/manitoba/physical_audit.json', EIGOR,
    )])
    raw_name = 'Tomographic_Dataset_1.zip'
    registries(directory, 'OBSERVATIONS.csv', ('dataset', 'object_id', 'observation_id'), [
        (MANITOBA_COHORT, row['target'], row['id'], f"EIGOR1-{row['frequency_GHz']:g}GHz-24x23",
         raw_name, report['members'][row['target']]['member'] +
         f"；{row['frequency_GHz']:g} GHz；原 552 通道，证书冻结取 288 通道")
        for row in cases
    ])
    physical = public_study + '/manitoba/physical.py'
    registries(directory, 'OPERATORS.csv', ('operator',), [(
        'manitoba_eigor_native_calibration', 'physical.py / scan.py', 'CPU',
        'EIGOR Data 1 等效二维复增益标定、固定接收行、40² 严格阈值证书。',
        '单列空场和金属留出误差、有限模型网格差与区间惯性。', physical, physical,
    )])
    markdown_path = directory / 'README.md'
    original = markdown_path.read_text(encoding='utf-8')
    entry = (
        '|独立实验室的材料峰值结论如何受校准约束？|'
        f'`{MANITOBA_COHORT}` 的 4 件实物与 28 份频率观测|'
        '[校准误差分解与复算](../z_alg_optical_inertia/validation/manitoba.html)|'
    )
    if entry not in original:
        anchor = '|能否用较少发射排除镜像材料？|'
        original = original.replace(anchor, entry+'\n'+anchor, 1)
    for key, caption in (('COHORTS.csv', '个数据集、子集及派生视图'), ('OPERATORS.csv', '个可调用的算子入口')):
        with (directory / key).open(encoding='utf-8', newline='') as stream:
            count = sum(1 for _ in csv.DictReader(stream))
        original = re.sub(r'\d+(?= ' + caption + ')', str(count), original, count=1)
    markdown_path.write_text(original, encoding='utf-8')
    renderer = runpy.run_path(str(directory / 'render_site.py'))
    html = renderer['markdown'].markdown(original, extensions=['tables', 'fenced_code', 'toc'])
    html = re.sub(r'href="([^"]*)"', renderer['internal_link'], html)
    html = html.replace('<table>', '<div class="wide"><table>').replace('</table>', '</table></div>')
    renderer['write_page']('index.html', '实验手册', html, markdown_path.name)
    renderer['registry']()


def stage(checkout, artifacts, notices):
    for path in artifacts:
        dest = checkout / SITE / path.relative_to(HERE)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    for relative, content in notices.items():
        (checkout / relative).write_text(content, encoding='utf-8')
    update_indexes(checkout)
    update_handbook(checkout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-zip', type=Path, default=BASE)
    parser.add_argument('--assets-out', type=Path, default=OUTPUT)
    parser.add_argument('--checkout', type=Path, default=DEFAULT_CHECKOUT)
    parser.add_argument('--preview', action='store_true', help='Generate local pages only from available sources.')
    parser.add_argument('--list-missing', action='store_true', help='Print unfinished source paths without building.')
    args = parser.parse_args()
    missing = missing_files()
    if not args.base_zip.is_file():
        missing.append(str(args.base_zip))
    if args.list_missing:
        print('\n'.join(missing) if missing else '定稿输入齐全。')
        return
    if args.preview:
        render_pages(preview=True)
        return
    if missing:
        raise SystemExit('仍缺少定稿输入：\n'+'\n'.join(missing))
    args.assets_out.mkdir(exist_ok=True, parents=True)
    make_latex(args.assets_out)
    render_pages()
    notices = old_page_notices(args.checkout)
    artifacts = [p for p in sorted(HERE.rglob('*')) if is_public(p)]
    make_data(args.base_zip, args.assets_out, artifacts, notices)
    stage(args.checkout, artifacts, notices)
    for file in (PDF, DATA, LATEX):
        path = args.assets_out / file
        print(f'{path.resolve()} ({path.stat().st_size:,} bytes)')
    print(f'已同步 B1 工作树：{args.checkout / SITE}')


if __name__ == '__main__':
    main()

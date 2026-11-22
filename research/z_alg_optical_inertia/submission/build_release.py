#!/usr/bin/env python3
"""Build the two release ZIPs and stage the manuscript site in the existing B1 checkout."""
import argparse
import re
import shutil
import subprocess
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from zipfile import ZIP_DEFLATED, ZipFile

from publish_pages import (
    DATA_NAME, DEFAULT_CHECKOUT, DEFAULT_TAG, HERE, LATEX_NAME, MANUSCRIPT, ROOT, SITE_PATH, STUDY,
    render_legacy_entry, render_submission, stage_handbook, stage_indexes,
)

OLD_RELEASE = ROOT / '.artifacts/z_alg_optical_release'
BASE_ZIP = OLD_RELEASE / 'Z_ALG_Optical_Inertia_2026-11-22.zip'
OUTPUT = ROOT / '.artifacts/z_alg_optical_submission_release'
R10 = Path('research/delegated/max_struct_r1/r10_worktree/research/research_max_struct_r10')
R10_PUBLIC = R10 / 'data_public/twod/r10_twod_0001.npz'
R10_PRIVATE = R10 / 'sealed_truth/twod/r10_twod_0001.npz'
R10_LOCAL_REFERENCE = Path('diagnostic/r10_twod_0001_offline_reference.npz')
SOURCE_FILE = 'manuscript_ip.tex'
TEX_FILE = re.compile(r'\\(includegraphics\*?|input|include|bibliography|addbibresource)(?:\s*\[[^\]]*\])?\s*\{([^{}]+)\}')
GRAPHIC_EXTS = ('.pdf', '.png', '.jpg', '.jpeg', '.eps')
EXCLUDE = {'.DS_Store', 'Thumbs.db', '__pycache__', '.git'}
EXCLUDE_SUFFIX = ('.pyc', '.pyo', '.aux', '.log', '.toc', '.out', '.partial', '.synctex.gz')
RELEASE_FILES = {DATA_NAME, LATEX_NAME}
RELEASE_DEPENDENCIES = '''

## 英文定稿新增依赖

定稿从原光学惯性完整附件继承所有冻结观测与源码。新增计算、区间回执、网页和论文位于
`research/z_alg_optical_inertia/submission/`。旧 24² 体检额外带入 MAX-STRUCT R10 的同名
公开合成观测，位于其原项目的 `data_public/twod/`。原 R10 封存材料及其本地复制均沿用之前
的发布边界。完整端口重算、区间验证和峰值排除只需公开输入。图和离线标量评分沿用已记录的结果。
先执行不读真值的 `diagnostic/run.py`，再用 `diagnostic/audit.py` 独立复核端口。

打开 `research/z_alg_optical_inertia/submission/index.html` 浏览定稿，或阅读同目录
`REPRODUCE.md` 独立复算。图、已冻结提案和双精度到区间的证书都按原相对路径提供。
'''


def files(directory):
    return [path for path in sorted(directory.rglob('*'))
            if path.is_file() and not EXCLUDE.intersection(path.relative_to(directory).parts)
            and not path.name.endswith(EXCLUDE_SUFFIX) and path.name not in RELEASE_FILES
            and path != HERE / R10_LOCAL_REFERENCE]


def resolve_tex_file(request, extensions, directories):
    name = Path(request.strip())
    if name.is_absolute():
        raise ValueError(f'TeX input must be relative: {request}')
    choices = (name,) if name.suffix else tuple(name.with_suffix(ext) for ext in extensions)
    for folder in directories:
        for option in choices:
            path = folder / option
            if path.is_file():
                return path.resolve()
    raise FileNotFoundError(f'TeX dependency not found: {request}')


def flat_source(source, stage):
    copied = {}
    root_dirs = (HERE, HERE/'figures', STUDY/'figures', HERE/'curve', HERE/'single', HERE/'diagnostic', HERE/'typeset')

    def put(path, target):
        if target in copied and copied[target] != path:
            raise ValueError(f'TeX source package name collision: {copied[target]} and {path}')
        if target in copied:
            return target
        if not re.fullmatch(r'[A-Za-z0-9_]+\.[A-Za-z0-9]+', target):
            raise ValueError(f'Rename the TeX input to an IOP-compatible basename: {target}')
        copied[target] = path
        if path.suffix == '.tex':
            raw = path.read_text(encoding='utf-8')
            directories = (path.parent,) + root_dirs

            def dependency(match):
                kind, arg = match.groups()
                graphics = kind.startswith('includegraphics')
                bibtex = kind == 'bibliography'
                exts = GRAPHIC_EXTS if graphics else ('.bib',) if bibtex or kind == 'addbibresource' else ('.tex',)
                names = []
                for requested in (arg.split(',') if bibtex else (arg,)):
                    path_ref = resolve_tex_file(requested, exts, directories)
                    new_name = put(path_ref, path_ref.name)
                    names.append(Path(new_name).stem if bibtex or kind == 'include' else new_name)
                prefix = match.group(0)[:match.group(0).rfind('{')]
                return prefix + '{' + ','.join(names) + '}'

            # Keep commented manuscript code unchanged.
            lines = []
            for line in raw.splitlines(keepends=True):
                front, tail = re.split(r'(?<!\\)(?=%)', line, maxsplit=1) if re.search(r'(?<!\\)%', line) else (line, '')
                lines.append(TEX_FILE.sub(dependency, front) + tail)
            (stage / target).write_text(''.join(lines), encoding='utf-8')
        else:
            shutil.copy2(path, stage / target)
        return target

    put(source.resolve(), SOURCE_FILE)
    for name in ('iopjournal.cls', 'orcid.pdf'):
        put((HERE/'typeset'/name).resolve(), name)
    other = sorted(p.name for p in stage.iterdir())
    (stage/'README.txt').write_text(
        'Inverse Problems LaTeX sources\n\n'
        'Extract all files to one directory and run\n'
        '  tectonic -X compile manuscript_ip.tex\n'
        'or\n'
        '  latexmk -pdf manuscript_ip.tex\n\n'
        'The LaTeX class and optional ORCID icon are distributed by IOP Publishing.\n'
        'Publisher archive\nhttps://publishingsupport.iopscience.iop.org/wp-content/uploads/2025/07/ioplatextemplate.zip\n\n'
        + '\n'.join(other) + '\n', encoding='utf-8',
    )


def latex_zip(output):
    pending = output.with_suffix('.partial')
    with TemporaryDirectory(prefix='flat-latex-', dir=output.parent) as folder:
        stage = Path(folder)
        flat_source(HERE / SOURCE_FILE, stage)
        archive_members = files(stage)
        subprocess.run(['tectonic', '-X', 'compile', SOURCE_FILE], cwd=stage, check=True)
        with ZipFile(pending, 'w', compression=ZIP_DEFLATED) as archive:
            for path in archive_members:
                archive.write(path, path.name)
    pending.replace(output)


def research_zip(base, output):
    overlays = {path.relative_to(ROOT).as_posix(): path for path in files(HERE)}
    for path in (STUDY/'index.html', STUDY/'build_site.py', STUDY/'RESEARCH_REPORT_ZH.md', ROOT/R10_PUBLIC):
        overlays[path.relative_to(ROOT).as_posix()] = path
    pending = output.with_suffix('.partial')
    with ZipFile(base) as previous, ZipFile(pending, 'w', compression=ZIP_DEFLATED) as archive:
        for name in previous.namelist():
            if (name not in overlays and name not in ('DEPENDENCIES.md', R10_PRIVATE.as_posix())
                    and not name.endswith('/')):
                archive.writestr(name, previous.read(name))
        inherited = previous.read('DEPENDENCIES.md').decode('utf-8')
        archive.writestr('DEPENDENCIES.md', inherited.rstrip() + '\n' + RELEASE_DEPENDENCIES)
        for name, path in sorted(overlays.items()):
            archive.write(path, name)
    pending.replace(output)


def stage_site(checkout):
    for path in files(HERE):
        destination = checkout / PurePosixPath(SITE_PATH) / path.relative_to(HERE)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
    published = checkout / STUDY.relative_to(ROOT)
    for name in ('index.html', 'build_site.py', 'RESEARCH_REPORT_ZH.md'):
        shutil.copy2(STUDY / name, published / name)
    stage_indexes(checkout)
    stage_handbook(checkout)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--base-zip', type=Path, default=BASE_ZIP)
    cli.add_argument('--assets-out', type=Path, default=OUTPUT)
    cli.add_argument('--checkout', type=Path, default=DEFAULT_CHECKOUT)
    cli.add_argument('--release-tag', default=DEFAULT_TAG)
    args = cli.parse_args()
    args.assets_out.mkdir(parents=True, exist_ok=True)
    render_submission(args.release_tag)
    render_legacy_entry()
    latex_zip(args.assets_out/LATEX_NAME)
    research_zip(args.base_zip, args.assets_out/DATA_NAME)
    shutil.copy2(HERE/MANUSCRIPT, args.assets_out/MANUSCRIPT)
    stage_site(args.checkout)
    for asset in (DATA_NAME, LATEX_NAME, MANUSCRIPT):
        path = args.assets_out/asset
        print(f'{path.resolve()} ({path.stat().st_size:,} bytes)')
    print(f'已写入 B1 网页及导航：{args.checkout/SITE_PATH}')


if __name__ == '__main__':
    main()

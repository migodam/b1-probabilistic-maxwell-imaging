"""Publish this bounded study inside the existing B1 public repository."""
from pathlib import Path
import sys,subprocess,shutil,json,zipfile,re
root=Path(__file__).resolve().parents[1]
repo=Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/Gaussian/B1/b1-probabilistic-maxwell-imaging")
work=root.parent/'publish-w-p0'
relative=Path('research/research_w_p0_20261022')
tag='research-w-p0-20261022'
def run(args,cwd=None):
    return subprocess.check_output(args,cwd=cwd or root,text=True,stderr=subprocess.STDOUT).strip()

def main():
    archive=root/'publication/Research_W_P0_ARRAYS.zip';archive.parent.mkdir(exist_ok=True)
    arrays=sorted((root/'results').glob('qualification*/WALL*.npz'))
    arrays += [root/'figures/raw/02_geometry.csv',root/'figures/raw/03_first_source_field.csv']
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in arrays:z.write(p,p.relative_to(root))
    if not work.exists():
        run(['git','fetch','origin'],repo)
        run(['git','worktree','add','-b','research-w-p0-20261022',str(work),'origin/main'],repo)
    dest=work/relative
    if dest.exists():raise RuntimeError('Publication destination already exists; inspect instead of overwriting')
    def ignored(directory,names):
        banned={'__pycache__','.DS_Store','.git','.venv','publication','.pytest_cache'}
        for name in names:
            p=Path(directory)/name
            if name.endswith('.pyc') or (name.startswith('WALL') and name.endswith('.npz')):banned.add(name)
            if p in arrays:banned.add(name)
        return banned
    shutil.copytree(root,dest,ignore=ignored)
    # Public root deliverable is a packaging snapshot of the authoritative
    # append-only results ledger, not a new accounting or hash operation.
    shutil.copy2(root/'results/PHYSICS_LEDGER.jsonl',dest/'PHYSICS_LEDGER.jsonl')
    (dest/'publication').mkdir()
    (dest/'publication/README.md').write_text('''# Publication

Scientific state: P0-A PASS; P0-B UNRESOLVED / UNQUALIFIED_SOLVER. No witness search or P1.

Large wall fields and grid CSVs: [Research_W_P0_ARRAYS.zip](https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-w-p0-20261022/Research_W_P0_ARRAYS.zip). Extract relative paths in this project root. Source, figures, small diagnostic arrays, decisions and packaging-time ledgers are in Git.

The Release also carries final publication/readback receipts and the final append-only ledger, including post-commit publication costs. No file hash checks were run. The root PHYSICS_LEDGER.jsonl is the packaging snapshot; results/PHYSICS_LEDGER.jsonl is the live authoritative path in the working project.
''')
    patterns=[re.compile(r'-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----'),re.compile(r'gh[pousr]_[A-Za-z0-9]{25,}'),re.compile(r'sk-[A-Za-z0-9]{35,}')]
    matches=[]
    for p in dest.rglob('*'):
        if p.is_file() and p.suffix in {'.py','.md','.txt','.json','.jsonl','.html','.sh','.tex','.log'}:
            value=p.read_text(errors='ignore')
            if any(pat.search(value) for pat in patterns):matches.append(str(p.relative_to(dest)))
    if matches:raise RuntimeError('Sensitive-pattern review needed: '+','.join(matches))
    nav=work/'research/README.md'
    text=nav.read_text()
    marker='research_w_p0_20261022'
    if marker not in text:
        nav.write_text(text+'\n\n## Research W P0 — 2026-10-22\n\n[执行入口](research_w_p0_20261022/START_HERE.md) · [中文报告](research_w_p0_20261022/RESEARCH_REPORT_ZH.md) · [图文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_w_p0_20261022/)\n\n一维解析回归PASS；二维前向网格资格未过，P0-B为UNRESOLVED。没有合格墙体歧义见证、搜索优势或P1结论。\n')
    run(['git','add',str(relative),'research/README.md'],work)
    status=run(['git','diff','--cached','--stat'],work)
    print(status)
    run(['git','commit','-m','Publish Research W P0 analytic reproduction and unresolved 2D qualification'],work)
    sha=run(['git','rev-parse','HEAD'],work)
    print(run(['git','push','origin','HEAD:main'],work))
    record=dict(repository='migodam/b1-probabilistic-maxwell-imaging',commit=sha,subdirectory=str(relative),release_tag=tag,
                arrays_asset_bytes=archive.stat().st_size,arrays_files=len(arrays),sensitive_pattern_matches=matches,
                science_status='P0_B_UNRESOLVED_NUMERICAL_QUALIFICATION')
    (root/'publication/PUBLISH.json').write_text(json.dumps(record,indent=2))
    notes=root/'publication/RELEASE_NOTES.md'
    notes.write_text(f'''Research W P0, source commit `{sha}`. P0-A independent 1D Maxwell PASS. P0-B remains UNRESOLVED / UNQUALIFIED_SOLVER: finest-grid full-data whitened differences 0.951–1.229 > 0.1. No candidate search, joint ambiguity witness, posterior, NN or 3D results.

258 actual discrete states, 266 physical LU/core factorizations, 1184 complex RHS. Failures and repairs retained. The array ZIP contains original/corrected wall fields and large per-grid CSVs; extract relative paths inside `{relative}`. Other source/evidence is committed in that subdirectory. Final publication costs and anonymous verification receipts are separate assets.
''')
    print(run(['gh','release','create',tag,'--repo',record['repository'],'--target',sha,'--title','Research W P0: analytic PASS, 2D qualification UNRESOLVED','--notes-file',str(notes),str(archive)],work))
    print(json.dumps(record))

if __name__=='__main__':main()

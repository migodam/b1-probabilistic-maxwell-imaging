"""Copy only this campaign's public text/code/figures into existing B1 tree.

No network, no Git mutations, no hashes, no new physical calls. Raw arrays are
separate release artifacts. Other campaign directories are never altered; this
campaign's public stage and its research-index entries may be updated.
"""
from pathlib import Path
import argparse,shutil,json,time,sys,os,stat
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from wslam.accounting import BudgetBook
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repository',type=Path,required=True);ap.add_argument('--partial',action='store_true');args=ap.parse_args()
    before=time.process_time();started=time.monotonic();repo=args.repository.resolve()
    if not (repo/'.git').exists():raise RuntimeError('Existing repository/worktree required')
    dest=repo/'research/research_w_slam';dest.mkdir(parents=True,exist_ok=True)
    allowed={'.md','.json','.jsonl','.csv','.py','.txt','.png','.pdf','.html'}
    selected=[]
    for path in ROOT.iterdir():
        if path.is_file() and path.suffix.lower() in allowed:selected.append(path)
    for dirname in ['src','tools','tests','docs','source_package','legacy_inputs','inputs','figures','research/delegated']:
        for path in (ROOT/dirname).rglob('*'):
            if path.is_file() and path.suffix.lower() in allowed and '__pycache__' not in path.parts:
                selected.append(path)
    # Top-level summaries and the original ledger are evidence, not raw arrays.
    for path in (ROOT/'results').glob('*'):
        if path.is_file() and path.suffix.lower() in {'.json','.csv','.md'}:selected.append(path)
    for dirname in ['e0','e1','a0','pilot_summary_tables']+([] if args.partial else ['e2e3/inference','e2e3/inference_attempt_001_failed','e2e3/finite_bank','e2e3/bank_validation','e2e3/posterior_endpoint_audit','e2e3/negative_controls','e4','e5']):
        for path in (ROOT/'results'/dirname).rglob('*'):
            if path.is_file() and path.suffix.lower() in {'.json','.csv','.md'}:selected.append(path)
    # Preserve actual job provenance, including failures. Private executable paths
    # inside these scientific records remain historical evidence, not runtime deps.
    for path in (ROOT/'results/jobs').rglob('*'):
        if path.is_file() and path.name in {'RECEIPT.json','OUTPUT.log'}:selected.append(path)
    selected=[p for p in sorted(set(selected)) if not p.is_symlink() and p.resolve().is_relative_to(ROOT.resolve())
              and p.relative_to(ROOT).as_posix() not in {'README.md','PUBLICATION_FILES.json','RELEASE_LAYOUT.json'}]
    adaptations=[]
    for path in selected:
        target=dest/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        if target.is_symlink():raise RuntimeError('Refuse publication destination symlink: '+str(target))
        if target.exists() and not target.stat().st_mode & stat.S_IWUSR:
            # copy2 preserves read-only legacy input modes. Only the generated
            # public copy is made writable for replacement; source stays frozen.
            target.chmod(target.stat().st_mode|stat.S_IWUSR)
        shutil.copy2(path,target)
        if path.relative_to(ROOT).as_posix()=='tools/run_job.py':
            # Change only the public copy. The active campaign wrapper and its
            # current finite-bank child keep their exact original source/path.
            text=target.read_text()
            old="(ROOT/'src',ROOT/'vendor',ROOT.parent/'research-w-continuation'/'.deps')"
            new="(ROOT/'src',ROOT/'vendor')"
            if old in text:
                target.write_text(text.replace(old,new))
                adaptations.append(dict(path='tools/run_job.py',change='public copy removes old sibling .deps preference; installed requirements supply dependencies',active_source_unchanged=True))
    shutil.copy2(ROOT/'START_HERE.md',dest/'README.md')
    note='实验进行中：E0/E1/A0已执行，概率参考尚在计算' if args.partial else '已完成的物理、概率与成本结果；限制及负结果保留'
    index=repo/'research/index.html'
    marker='<!-- WSLAM_START -->';end='<!-- WSLAM_END -->'
    card=marker+'<section><h2><a href="research_w_slam/">Research W-SLAM</a></h2><p>'+note+'</p></section>'+end
    if index.exists():
        text=index.read_text()
        if marker in text:
            a=text.index(marker);b=text.index(end,a)+len(end);text=text[:a]+card+text[b:]
        else:text=text.replace('</main>',card+'</main>') if '</main>' in text else text+card
        index.write_text(text)
    readme=repo/'research/README.md'
    text=readme.read_text() if readme.exists() else '# Research\n'
    entry='- [Research W-SLAM](research_w_slam/)：'+note+'。'
    if 'research_w_slam/' not in text:text+='\n'+entry+'\n'
    else:text='\n'.join(entry if line.startswith('- [Research W-SLAM](research_w_slam/)') else line for line in text.splitlines())+'\n'
    readme.write_text(text)
    manifest=[dict(path=str(p.relative_to(ROOT)),bytes=(dest/p.relative_to(ROOT)).stat().st_size,
                   source_bytes=p.stat().st_size,copied_verbatim=not any(a['path']==p.relative_to(ROOT).as_posix() for a in adaptations)) for p in selected]
    manifest.append(dict(path='README.md',source_path='START_HERE.md',bytes=(dest/'README.md').stat().st_size,copied_verbatim=True))
    # Metadata only: no release ZIP is created or uploaded, and arrays are never
    # opened. Archive members must preserve these stage-relative paths on unpack.
    release=[]
    copied={p.relative_to(ROOT).as_posix() for p in selected}
    for dirname in ['inputs','results','figures','source_package','legacy_inputs']:
        for path in (ROOT/dirname).rglob('*'):
            if path.is_file() and not path.is_symlink() and path.suffix.lower() in {'.npz','.npy','.png'} and path.resolve().is_relative_to(ROOT.resolve()):
                rel=path.relative_to(ROOT).as_posix()
                if rel not in copied:release.append(dict(path=rel,bytes=path.stat().st_size))
    (dest/'RELEASE_LAYOUT.json').write_text(json.dumps(dict(status='ARRAY_LAYOUT_MANIFEST',
        stage_root='research/research_w_slam',archive_members='stage-relative paths with no extra directory prefix',
        unpack_relative_to='research/research_w_slam/',files=sorted(release,key=lambda r:r['path']),
        snapshot_only=True,snapshot_of_terminal_research=not args.partial,release_verification_recorded_separately=True,contains_hash_checks=False,
        warning='Public text/code alone lacks these arrays; do not trigger a new paid solve to replace a missing saved input.'),indent=2))
    (dest/'PUBLICATION_FILES.json').write_text(json.dumps(dict(status='PARTIAL_RUNNING' if args.partial else 'COMPLETE_RESEARCH_SOURCE_SNAPSHOT',files=manifest,
        generated_outputs=['PUBLICATION_FILES.json','RELEASE_LAYOUT.json'],public_copy_adaptations=adaptations,
        raw_arrays='Separate release layout in RELEASE_LAYOUT.json; archive not created by this helper; no hashes generated or checked'),indent=2))
    cpu=time.process_time()-before
    accounting=dict(kind='nonphysical_publication_preparation',wall_seconds=time.monotonic()-started,partial=args.partial,files=len(manifest),new_physical_calls=0)
    if os.environ.get('WSLAM_STAGE'):
        accounting.update(internal_cpu_seconds_informational=cpu,CPU_accounting='enclosing run_job receipt only; no nested CPU charge')
    else:
        accounting.update(process_cpu_seconds=cpu,CPU_accounting='standalone helper measured process CPU')
    event=BudgetBook().event('E0E1','prepare_publication',**accounting)
    print(json.dumps(dict(destination=str(dest),files=len(manifest),event_id=event['id'])))
if __name__=='__main__':main()

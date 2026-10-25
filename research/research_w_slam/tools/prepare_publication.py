"""Copy only this campaign's public text/code/figures into existing B1 tree.

No network, no Git mutations, no hashes, no new physical calls. Raw arrays are
separate release artifacts. Existing research directories are never altered.
"""
from pathlib import Path
import argparse,shutil,json,time,resource,sys
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
    for dirname in ['src','tools','tests','docs','source_package','legacy_inputs','figures','research/delegated']:
        for path in (ROOT/dirname).rglob('*'):
            if path.is_file() and path.suffix.lower() in allowed and '__pycache__' not in path.parts:
                selected.append(path)
    for dirname in ['e0','e1','a0']+([] if args.partial else ['e2e3/inference','e2e3/bank_validation','e2e3/negative_controls','e4','e5']):
        for path in (ROOT/'results'/dirname).rglob('*'):
            if path.is_file() and path.suffix.lower() in {'.json','.csv','.md'}:selected.append(path)
    for path in sorted(set(selected)):
        target=dest/path.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
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
    if 'research_w_slam/' not in text:text+='\n- [Research W-SLAM](research_w_slam/)：'+note+'。\n'
    readme.write_text(text)
    manifest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size) for p in sorted(set(selected))]
    (dest/'PUBLICATION_FILES.json').write_text(json.dumps(dict(status='PARTIAL_RUNNING' if args.partial else 'FINAL_CANDIDATE',files=manifest,raw_arrays='Release; no hashes generated or checked'),indent=2))
    event=BudgetBook().event('E0E1','prepare_publication',kind='nonphysical_publication_preparation',process_cpu_seconds=time.process_time()-before,wall_seconds=time.monotonic()-started,partial=args.partial,files=len(manifest),new_physical_calls=0)
    print(json.dumps(dict(destination=str(dest),files=len(manifest),event_id=event['id'])))
if __name__=='__main__':main()

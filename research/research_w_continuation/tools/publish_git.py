"""Publish the already-reviewed new directory without modifying old P0."""
from pathlib import Path
import argparse, subprocess, json
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--repository-root',required=True);p.add_argument('--asset',required=True);a=p.parse_args()
repo=Path(a.repository_root).resolve();asset=Path(a.asset).resolve()
def run(args,capture=False):
    if capture:return subprocess.check_output(args,cwd=repo,text=True).strip()
    subprocess.run(args,cwd=repo,check=True)
allowed=('research/research_w_continuation/','research/README.md','research/index.html')
if run(['git','remote','get-url','origin'],True)!='https://github.com/migodam/b1-probabilistic-maxwell-imaging.git':raise RuntimeError('Unexpected repository')
run(['git','add','--','research/research_w_continuation','research/README.md','research/index.html'])
staged=run(['git','diff','--cached','--name-only'],True).splitlines()
if not staged or any(not (s.startswith(allowed[0]) or s in allowed[1:]) for s in staged):raise RuntimeError('Staged scope exceeds new continuation publication')
run(['git','commit','-m','Research W continuation: qualify 2D full-wave and retain bounded search limitations'])
commit=run(['git','rev-parse','HEAD'],True)
run(['git','push','origin','HEAD:main'])
run(['gh','release','create','research-w-continuation-20261023',str(asset),'--repo','migodam/b1-probabilistic-maxwell-imaging','--target',commit,'--title','Research W continuation: qualified forward model, no ambiguity witness','--notes-file',str(root/'RELEASE_NOTES.md')])
row=dict(commit=commit,repository='https://github.com/migodam/b1-probabilistic-maxwell-imaging',
         release='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-continuation-20261023',
         scope='New continuation, research reading entries only; old P0 unchanged')
(root/'PUBLICATION.json').write_text(json.dumps(row,indent=2));print(json.dumps(row),flush=True)

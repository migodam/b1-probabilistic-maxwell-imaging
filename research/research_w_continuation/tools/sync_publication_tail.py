"""Refresh public text/artifacts and disclose costs later than report snapshot."""
from pathlib import Path
import argparse, json, shutil
root=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser();p.add_argument('--destination',required=True);a=p.parse_args();dest=Path(a.destination).resolve()
old=[json.loads(x) for x in (root/'PHYSICS_LEDGER.jsonl').read_text().splitlines()]
old_cost={(r['stage'],r.get('label')) for r in old if r.get('event')=='cost'}
allrows=[];tail=[];charged={}
for stage,path in [('P0-Q','results/PHYSICS_LEDGER.jsonl'),('P0-B','results/p0b/PHYSICS_LEDGER.jsonl')]:
    rows=[json.loads(s) for s in (root/path).read_text().splitlines() if s.strip()]
    charged[stage]=sum(r.get('charged_cpu_seconds',0) for r in rows)
    for row in rows:
        r=dict(stage=stage,**row);allrows.append(r)
        if r.get('event')=='cost' and (stage,r.get('label')) not in old_cost:tail.append(r)
(root/'PHYSICS_LEDGER.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in allrows))
(root/'PUBLICATION_TAIL_COST.json').write_text(json.dumps(dict(scope='Additional measured jobs after narrative/COST_SUMMARY snapshot; publication push/readback receipt may follow as Release evidence.',additional_rows=tail,additional_charged_cpu_seconds=sum(r.get('charged_cpu_seconds',0) for r in tail),current_stage_charged_cpu_seconds=charged,physics_actions_added_by_this_script=0),indent=2))
extensions={'.py','.md','.json','.jsonl','.txt','.csv','.png','.svg','.html','.log'}
files=[]
for src in sorted(root.rglob('*')):
    if not src.is_file() or src.suffix not in extensions:continue
    rel=src.relative_to(root)
    if any(s in ('.deps','.venv','__pycache__','.git','mpl','.DS_Store') for s in rel.parts):continue
    target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,target)
    files.append(dict(path=rel.as_posix(),bytes=src.stat().st_size))
(dest/'PUBLICATION_FILES.json').write_text(json.dumps(dict(files=files,scope='Public text, code, plots and metrics. NPZ files in separate Release asset; no dependency installs.'),indent=2))
print(json.dumps(dict(public_files=len(files),tail_cpu_seconds=sum(r.get('charged_cpu_seconds',0) for r in tail),stage_cpu_seconds=charged)))

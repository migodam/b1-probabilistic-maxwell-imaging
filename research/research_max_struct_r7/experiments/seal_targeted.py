"""Public-only targeted80 sealing. No truth access or hash operations."""
from pathlib import Path
from datetime import datetime,timezone
import json
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
def stamp(p):
    p=Path(p);s=p.stat();return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)
def main():
    manifest=ROOT/'configs/targeted_final_manifest.json';tasks=json.loads(manifest.read_text())['tasks']
    expected={(f'r7_target_{i:02d}',arm) for i in range(20) for arm in ('a_plus','cie_dual','global_dual','cie_single')}
    assert len(tasks)==80 and {(t['case_id'],t['arm']) for t in tasks}==expected
    entries=[]
    for task in tasks:
        directory=REPO/task['output_dir'];rec=json.loads((directory/'RECEIPT.json').read_text())
        assert rec['status'] in ('SEALED_PUBLIC_PREDICTION','FAILED','WORKER_FAILED','CANCELLED','HARD_TIMEOUT'),task['task_id']
        if rec['status']=='SEALED_PUBLIC_PREDICTION':
            assert all((directory/name).is_file() for name in ('prediction.npz','prediction.json','BUDGET_STATES.npz'))
        entries.append(dict(task=task,status=rec['status'],files=[stamp(p) for p in sorted(directory.iterdir()) if p.is_file()],public_files=[stamp(REPO/p) for p in task['public_files']]))
    evidence=[stamp(manifest),stamp(ROOT/'data_public/targeted/cases.json'),stamp(ROOT/'configs/targeted_frozen.json'),stamp(ROOT/'TARGETED_PROTOCOL.md'),stamp(ROOT/'configs/targeted_single_amendment.json')]
    evidence.extend(stamp(REPO/p) for p in sorted({t['config'] for t in tasks}))
    seal=dict(status='ALL_TARGETED80_PUBLIC_OUTPUTS_SEALED',created_utc=datetime.now(timezone.utc).isoformat(),manifest=str(manifest),evidence=evidence,tasks=entries,truth_access='NONE',hash_checks='NONE',objects=20,tasks_count=80)
    path=ROOT/'audit/targeted/PUBLIC_SEAL.json'
    with path.open('x') as f:json.dump(seal,f,indent=2)
    print(json.dumps(dict(status=seal['status'],objects=20,tasks=80,path=str(path))))
if __name__=='__main__':main()

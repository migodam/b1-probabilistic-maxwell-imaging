"""Object-level hardware assignment; no truth or prediction values are read."""
from pathlib import Path
import json,time,zipfile
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]

def main():
    all_tasks=json.loads((ROOT/'configs/warm_all_manifest.json').read_text())
    remote_ids={f'r8_test_{i:04d}' for i in range(6,42) if i%4!=0}
    local=[];remote=[];already=[]
    for t in all_tasks['tasks']:
        if t['case_id'] in remote_ids:
            if (REPO/t['output_dir']).exists():raise RuntimeError('Remote object already started locally: '+t['task_id'])
            remote.append(t)
        elif (REPO/t['output_dir']/'RECEIPT.json').exists():
            rec=json.loads((REPO/t['output_dir']/'RECEIPT.json').read_text())
            if rec['status']!='SEALED_PUBLIC_PREDICTION':raise RuntimeError('Local task still running/failed '+t['task_id'])
            already.append(t['task_id'])
        else:local.append(t)
    for name,tasks in [('warm_local_remaining',local),('warm_xinan',remote)]:
        (ROOT/f'configs/{name}_manifest.json').write_text(json.dumps(dict(config=all_tasks['config'],tasks=tasks),indent=2))
    paths=set(REPO.glob('research/**/*.py'))
    paths={p for p in paths if not {'__pycache__','audit','delegated'}.intersection(p.relative_to(REPO).parts)}
    for t in remote:
        paths.update(REPO/p for p in t['public_files'])
        if t.get('initialization_file'):paths.add(REPO/t['initialization_file'])
    paths.update([ROOT/'configs/solver.json',ROOT/'configs/warm_xinan_manifest.json'])
    target=ROOT/'audit/transfer_warm_xinan.zip'
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=3) as z:
        for p in sorted(paths):z.write(p,p.relative_to(REPO))
    out=dict(status='OBJECT_LEVEL_HARDWARE_ASSIGNMENT_BEFORE_GT_SCORING',rule='Indices >=6 and not divisible by4 to XINAN; other objects stay Mac. Whole objects, all arms and extra initialization ablations share hardware.',remote_objects=sorted(remote_ids),remote_tasks=len(remote),local_remaining_tasks=len(local),already_sealed_local=already,truth_access='NONE',numerical_change='NONE; run_case added Windows peak-memory metadata only',resource_policy='XINAN read-only inventory found no Python process, GPU idle, CPU2.9%; native CPU8 P-core lanes, affinity0/2/.../14, GPU not used',archive=str(target.relative_to(REPO)),bytes=target.stat().st_size,archive_files=len(paths),cost_comparison='Pair within object/hardware; stratify absolute seconds by hardware, no mixed-CPU speed claim',old_queue_status='Parent-only SIGINT stops dispatch; child solves finish normally and individual receipts retained')
    (ROOT/'audit/HARDWARE_ASSIGNMENT.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))

if __name__=='__main__':main()

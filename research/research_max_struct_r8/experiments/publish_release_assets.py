"""Upload the already reviewed R8 archive assets; no scientific computation.

Uses an existing gh-authenticated draft release. Resume skips only same-name,
same-byte-count assets; this is transfer bookkeeping, not a content certificate.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path
import resource
import subprocess
import time

REPO = 'migodam/b1-probabilistic-maxwell-imaging'
ROOT = Path(__file__).resolve().parents[1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--archive-dir', type=Path, required=True)
    ap.add_argument('--tag', default='research-max-struct-r8-20261106')
    ap.add_argument('--workers', type=int, default=4)
    args = ap.parse_args()
    index = json.loads((ROOT/'delegated/cost_package/ARCHIVE_INDEX.json').read_text())
    assets = index['built_archives']
    existing = json.loads(subprocess.check_output([
        'gh','release','view',args.tag,'--repo',REPO,'--json','assets']))['assets']
    existing = {v['name']:v['size'] for v in existing}
    started = datetime.now(timezone.utc).isoformat()
    wall = time.perf_counter(); cpu = time.process_time()
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    child0 = child.ru_utime + child.ru_stime
    receipt = {'status':'UPLOADING','started_utc':started,'tag':args.tag,
               'repository':REPO,'workers':args.workers,'assets':[],
               'scientific_runs':0,'content_hash_checks':'NONE',
               'timing_scope':'this upload process and gh children; pre-list/setup and human/tool review excluded'}
    dest = ROOT/'audit/PUBLICATION_UPLOAD_RECEIPT.json'
    if dest.exists():
        with (ROOT/'audit/PUBLICATION_UPLOAD_HISTORY.jsonl').open('a') as f:
            f.write(json.dumps(json.loads(dest.read_text()))+'\n')
    def save():
        dest.write_text(json.dumps(receipt,indent=2)+'\n')
    def upload(a):
        path = args.archive_dir/a['name']
        if path.stat().st_size != a['bytes']:
            return {'name':a['name'],'status':'LOCAL_SIZE_CHANGED','bytes':path.stat().st_size}
        if a['name'] in existing:
            return {'name':a['name'],'bytes':a['bytes'],
                    'status':'EXISTING_SIZE_MATCH' if existing[a['name']]==a['bytes'] else 'REMOTE_SIZE_MISMATCH'}
        w = time.perf_counter()
        proc = subprocess.run(['gh','release','upload',args.tag,str(path),'--repo',REPO],
                              capture_output=True,text=True)
        return {'name':a['name'],'bytes':a['bytes'],
                'status':'UPLOADED' if proc.returncode==0 else 'UPLOAD_FAILED',
                'returncode':proc.returncode,'wall_seconds':time.perf_counter()-w,
                'output':(proc.stdout+proc.stderr)[-3000:]}
    save()
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending = {pool.submit(upload,a):a for a in assets}
        for future in as_completed(pending):
            receipt['assets'].append(future.result()); save()
            print(len(receipt['assets']), '/', len(assets), receipt['assets'][-1]['status'],
                  receipt['assets'][-1]['name'], flush=True)
    ok = all(a['status'] in ('UPLOADED','EXISTING_SIZE_MATCH') for a in receipt['assets'])
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt.update(status='ALL_ASSETS_UPLOADED' if ok else 'UPLOAD_INCOMPLETE',
                   completed_utc=datetime.now(timezone.utc).isoformat(),
                   wall_seconds=time.perf_counter()-wall,
                   cpu_seconds=time.process_time()-cpu+child.ru_utime+child.ru_stime-child0,
                   uploaded_bytes=sum(v['bytes'] for v in receipt['assets'] if v['status']=='UPLOADED'))
    save()
    if not ok: raise SystemExit(1)

if __name__ == '__main__': main()

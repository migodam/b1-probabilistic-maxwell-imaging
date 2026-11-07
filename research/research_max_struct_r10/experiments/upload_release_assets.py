"""Finite GitHub release upload with whole-process publication accounting.

Run only after the existing repository release has been created. No hashes,
scientific computation, background service, overwrite or automatic retry.
"""
from pathlib import Path
import concurrent.futures
import json
import resource
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = 'migodam/b1-probabilistic-maxwell-imaging'
TAG = 'research-max-struct-r10-20261107'


def main():
    started, cpu = time.perf_counter(), time.process_time()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    out = ROOT / 'audit' / 'publication_upload'
    out.mkdir(exist_ok=False)
    files = sorted((ROOT / 'release_assets').glob('*.zip'))
    files += [ROOT / 'audit' / 'PUBLICATION_MANIFEST.json', ROOT / 'audit' / 'PACKAGING_RECEIPT.json']

    def upload(path):
        t = time.perf_counter()
        command = ['gh', 'release', 'upload', TAG, str(path), '--repo', REPO]
        with (out / (path.name + '.log')).open('w') as log:
            p = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        return {'name': path.name, 'bytes': path.stat().st_size,
                'returncode': p.returncode, 'wall_seconds': time.perf_counter()-t}

    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(upload, files))
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    receipt = {'status': 'COMPLETE' if all(x['returncode'] == 0 for x in rows) else 'UPLOAD_FAILURES_RETAINED',
               'repo': REPO, 'tag': TAG, 'files': rows,
               'CPU_seconds': time.process_time()-cpu + after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
               'wall_seconds': time.perf_counter()-started, 'execution_platform': 'Mac',
               'CPU_scope': 'This uploader and all gh child processes once; no scientific or packaging CPU included',
               'new_forward_calls': 0, 'monetary_cost': None,
               'monetary_cost_reason': 'No hardware/network/tool billing rate available'}
    (ROOT / 'audit' / 'PUBLICATION_UPLOAD.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k != 'files'}))


if __name__ == '__main__':
    main()

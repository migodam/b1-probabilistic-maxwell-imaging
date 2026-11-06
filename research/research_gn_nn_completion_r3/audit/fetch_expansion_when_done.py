"""Collect the one finite remote expansion batch after its terminal receipt."""
from pathlib import Path
import json, subprocess, sys, time
S=Path(__file__).resolve().parents[1]
log=S/'audit/EXPAND1024_REMOTE.log'
while True:
    lines=log.read_text(errors='replace').splitlines() if log.exists() else[]
    final=None
    for line in reversed(lines[-10:]):
        try:r=json.loads(line)
        except ValueError:continue
        if r.get('selected')==896 and 'status'in r:final=r;break
    if final is not None:break
    time.sleep(30)
with (S/'audit/EXPAND1024_FETCH.log').open('w')as f:
    p=subprocess.run(['python3','-B',str(S/'src/remote_jobs.py'),'fetch'],cwd=S,stdout=f,stderr=subprocess.STDOUT)
(S/'audit/EXPAND1024_REMOTE_DONE.json').write_text(json.dumps(dict(batch=final,fetch_exit_code=p.returncode,unix=time.time()),indent=2))

"""Sequential official physical campaign; no labels scored, no model selection."""
from pathlib import Path
import sys,subprocess,json,time,hashlib
OUT=Path(__file__).resolve().parents[1]
def main():
    freeze=json.loads((OUT/'configs/FREEZE.json').read_text())
    for entry in freeze['files']:
        p=Path(entry['path'])
        if hashlib.blake2b(p.read_bytes()).hexdigest()!=entry['blake2b']:
            raise ValueError('frozen source/config changed: '+str(p))
    replay=json.loads((OUT/'CHECKPOINT_REPRODUCTION.json').read_text())
    for group in ('A','B'):
        for script,role in [('generate_observations.py','generate'),('run_anchors.py','anchors')]:
            cmd=[sys.executable,'-u','-B',str(OUT/'src'/script),group]
            dest=OUT/'raw/jobs'/f'{group}_{role}.log'
            print(json.dumps(dict(event='START',group=group,role=role,log=str(dest))),flush=True)
            # Append only: retries never erase an interrupted or failed attempt.
            with dest.open('a') as log:
                result=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
            if result.returncode:raise RuntimeError(f'{group} {role} stopped ({result.returncode}); inspect preserved log')
            print(json.dumps(dict(event='COMPLETE',group=group,role=role)),flush=True)
    print(json.dumps(dict(status='ALL_512_FULL_FIDELITY_PHYSICS_COMPLETE',scoring=False,time=time.time())),flush=True)
if __name__=='__main__':main()

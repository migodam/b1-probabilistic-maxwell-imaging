"""Serial guarded queue; no truth or scores are opened here."""
from pathlib import Path
import argparse,subprocess,sys,time,json,datetime,os
ROOT=Path(__file__).resolve().parents[1]
METHODS=['FULL_GN','FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART']
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--noise',type=int,choices=[0,1],default=0);a=ap.parse_args();clock=json.loads((ROOT/'RUN_CLOCK.json').read_text());cut=datetime.datetime.fromisoformat(clock['no_new_heavy_after_utc'].replace('Z','+00:00')).timestamp();log=ROOT/'results'/f'CORE_QUEUE_r{a.noise}.jsonl'
    if not (ROOT/'METHODS_FROZEN.json').exists():raise RuntimeError('Definitive methods must be frozen before blindqueue')
    for casebase in ['blind_normal','blind_cavity','blind_strong']:
        case=f'{casebase}_r{a.noise}';inp=ROOT/'inputs/online'/case
        for method in METHODS:
            out=ROOT/'results/online'/case/method;out.mkdir(parents=True,exist_ok=True)
            if (out/'ONLINE_SEALED.json').exists():continue
            if time.time()>cut:raise RuntimeError('No-new-heavy deadline reached before next job')
            command=[sys.executable,str(ROOT/'tools/guarded_online.py'),'--audit-log',str(out/'FILE_ACCESS.jsonl'),str(ROOT/'tools/run_r3.py'),'run','--measurement',str(inp/'observation.npz'),'--config',str(inp/'config.json'),'--output',str(out),'--method',method]
            started=time.time();row=dict(event='job_start',case=case,method=method,started_unix=started,command=command)
            with log.open('a') as f:f.write(json.dumps(row)+'\n')
            print('START',case,method,flush=True)
            with (out/'OUTPUT.log').open('w') as stdout:result=subprocess.run(command,cwd=ROOT,stdout=stdout,stderr=subprocess.STDOUT,env={**os.environ,'R3_SOURCE_ROOT':str(ROOT/'frozen/src'),'PYTHONDONTWRITEBYTECODE':'1'})
            row=dict(event='job_finish',case=case,method=method,returncode=result.returncode,wall_seconds=time.time()-started,sealed=(out/'ONLINE_SEALED.json').exists())
            with log.open('a') as f:f.write(json.dumps(row)+'\n')
            print('FINISH',case,method,row,flush=True)
            if result.returncode:print('FAILED retained',out,flush=True)
    (ROOT/'results'/f'CORE_QUEUE_r{a.noise}_FINISHED.json').write_text(json.dumps(dict(time_unix=time.time(),noise=a.noise,methods=METHODS,cases=['blind_normal','blind_cavity','blind_strong']),indent=2))
if __name__=='__main__':main()

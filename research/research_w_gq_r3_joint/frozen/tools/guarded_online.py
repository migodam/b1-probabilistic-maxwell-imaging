"""Run an online entry point with a durable file-access audit and truth-read guard.

The Python audit hook records open events, and rejects private generator/scoring
reads. It covers Python file opens in this process, not arbitrary native code or
shell children. The command must be a Python entry-point script.
"""
from pathlib import Path
import argparse,os,sys,json,time,runpy
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--audit-log',required=True,type=Path);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    if not a.command:raise ValueError('entry-point command required')
    a.audit_log.parent.mkdir(parents=True,exist_ok=True);fd=os.open(a.audit_log,os.O_WRONLY|os.O_CREAT|os.O_EXCL,0o600)
    forbidden=[str(ROOT/'inputs/private'),str(ROOT/'inputs/invalid_sampler_attempt'),str(ROOT/'tools/generate.py'),str(ROOT/'tools/score.py')]
    own=Path(a.command[0]).resolve();t=time.time()
    def log(row):os.write(fd,(json.dumps(dict(time_unix=time.time(),pid=os.getpid(),**row))+'\n').encode())
    def hook(event,args):
        if event!='open':return
        value=args[0]
        if not isinstance(value,(str,bytes,os.PathLike)):return
        path=os.path.abspath(os.fsdecode(value));mode=args[1];flags=args[2]
        write=(isinstance(mode,str) and any(x in mode for x in 'wax+')) or bool(flags&(os.O_WRONLY|os.O_RDWR))
        denied=not write and (any(path==x or path.startswith(x+os.sep) for x in forbidden) or path.endswith('/SCORES.json') or path.endswith('/observation.json'))
        log(dict(event='file_open',path=path,mode=mode,flags=flags,write=write,denied=denied))
        if denied:raise PermissionError('online read guard rejected private generator/scoring path')
    log(dict(event='guard_started',entrypoint=str(own),argv=a.command,forbidden_roots=forbidden,scope='Python open audit only; no OS sandbox or external process visibility'))
    sys.addaudithook(hook);sys.argv=a.command
    try:runpy.run_path(str(own),run_name='__main__')
    finally:log(dict(event='guard_finished',elapsed_seconds=time.time()-t)) # OS closes audit fd at process exit
if __name__=='__main__':main()

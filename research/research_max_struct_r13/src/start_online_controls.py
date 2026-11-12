"""Run registered extra controls after the finite base sequence releases cores."""
import argparse,json,time,sys,subprocess
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--legacy-repo',type=Path,required=True);p.add_argument('--base-root',type=Path,required=True);a=p.parse_args();r=a.root
    dependency=a.base_root/'FINITE_SEQUENCE_RECEIPT.json'
    print('Waiting for the immutable base confirmation sequence; no fresh scoring is performed.',flush=True)
    while not dependency.exists():time.sleep(1)
    d=json.loads(dependency.read_text())
    if d['status']!='FINITE_SEQUENCE_COMPLETED':raise RuntimeError('Base sequence not complete')
    name='confirmation_online_controls16'
    argv=[sys.executable,str(r/'src/finite_queue.py'),'--root',str(r),'--legacy-repo',str(a.legacy_repo),'--manifest',str(r/'configs'/(name+'.json')),'--cores','0,2,4,6,8,10,12,14']
    code=subprocess.call(argv,cwd=r)
    if code:raise RuntimeError('Technical queue failure; preserve evidence for parent recovery')
    raise SystemExit(subprocess.call([sys.executable,str(r/'src/archive_campaign.py'),'--root',str(r),'--campaign',name],cwd=r))


if __name__=='__main__':main()

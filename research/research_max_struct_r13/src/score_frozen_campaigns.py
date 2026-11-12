"""Finite sequential post-seal scoring with whole-child resource receipts.

Requires the parent's joint authorization before any scorer is launched. This
script never changes predictions, authorizes seals, chooses algorithms or runs
an inverse/forward solver. Separate physics validation is not part of scoring.
"""
import argparse,json,os,subprocess,sys,time
from pathlib import Path
from seal_results import validate_seal


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--authorization',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();auth=json.loads(a.authorization.read_text())
    if auth.get('authorized_joint_offline_scoring') is not True:
        raise PermissionError('Parent joint offline scoring authorization is required')
    campaigns=auth['score_campaigns']
    for name in auth['required_campaigns']:
        seal=json.loads((root/'audit'/(name+'_seal.json')).read_text())
        validate_seal(seal,require_authorization=True)
    for item in auth['public_selection_files']:
        path=root/item['relative_path'];actual=path.stat()
        if actual.st_size!=item['bytes'] or actual.st_mtime_ns!=item['mtime_ns']:
            raise RuntimeError('Public selection index changed after joint freeze: '+str(path))
    logs=root/'audit/scoring_execution';logs.mkdir(parents=True,exist_ok=True)
    env=os.environ.copy()
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS'):
        env[key]='1'
    for campaign in campaigns:
        output=root/'analysis'/campaign
        if output.exists():raise FileExistsError('Scored campaign already exists; do not overwrite: '+campaign)
        argv=[sys.executable,str(root/'src/score_results.py'),'--root',str(root),
            '--seal',str(root/'audit'/(campaign+'_seal.json')),'--out',str(output),
            '--reference','dr','--structure-gallery']
        started=time.perf_counter()
        with (logs/(campaign+'.log')).open('xb') as log:
            child=subprocess.Popen(argv,env=env,cwd=root,stdout=log,stderr=subprocess.STDOUT)
            pid,status,usage=os.wait4(child.pid,0)
            child.returncode=os.waitstatus_to_exitcode(status)
        receipt=dict(campaign=campaign,returncode=child.returncode,
            CPU_seconds=usage.ru_utime+usage.ru_stime,wall_seconds=time.perf_counter()-started,
            CPU_scope='exclusive wait4 whole scorer child, including imports, serialization and exit',
            peak_RSS_bytes=int(usage.ru_maxrss) if sys.platform=='darwin' else int(usage.ru_maxrss*1024),
            new_physical_forward=0,new_optimization=0,host='Apple_M4',
            truth_access='OFFLINE_AFTER_JOINT_AUTHORIZATION')
        (logs/(campaign+'_RESOURCE.json')).write_text(json.dumps(receipt,indent=2)+'\n')
        if child.returncode:
            print(json.dumps(receipt),flush=True)
            raise RuntimeError('Scoring failed; partial output/log/cost retained, no automatic rerun')
        target=output/'SCORING_RECEIPT.json';scored=json.loads(target.read_text())
        scored.update(receipt);target.write_text(json.dumps(scored,indent=2)+'\n')
        print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()

"""Record a campaign's public predictions without reading truth or authorizing scoring.

Only the parent orchestrator may set authorized_scoring=true after inspecting the
seal. File sizes/mtime are provenance snapshots, not cryptographic certificates.
"""
import argparse,json
from pathlib import Path
from datetime import datetime,timezone


def read(path):
    return json.loads(Path(path).read_text())


def snapshot(path):
    p=Path(path).resolve();s=p.stat()
    return {'path':str(p),'bytes':s.st_size,'mtime_ns':s.st_mtime_ns}


def inventory(path):
    p=Path(path)
    return [snapshot(f) for f in sorted(p.rglob('*')) if f.is_file()] if p.exists() else []


def validate_seal(seal,require_authorization=True):
    if require_authorization and seal.get('authorized_scoring') is not True:
        raise PermissionError('Parent has not authorized scoring: authorized_scoring must be true in the reviewed seal.')
    terminal=['COMPLETED','COMPLETED_WITH_ERRORS','COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY']
    if seal.get('queue_status') not in terminal:
        raise RuntimeError('A final terminal queue receipt is required; missing/unknown/running queue status cannot authorize scoring.')
    if seal.get('queue_running') or seal.get('queue_pending'):
        raise RuntimeError('Queue was still running/pending at seal time; predictions are not final.')
    errors=[]
    for item in seal['files']:
        p=Path(item['path'])
        if not p.is_file():errors.append('MISSING '+str(p));continue
        s=p.stat()
        if s.st_size!=item['bytes'] or s.st_mtime_ns!=item['mtime_ns']:errors.append('CHANGED '+str(p))
    for directory in seal['directories']:
        actual={r['path'] for r in inventory(directory['path'])}
        if actual!=set(directory['files']):errors.append('DIRECTORY_CONTENTS_CHANGED '+directory['path'])
    if errors:raise RuntimeError('Prediction seal no longer matches:\n'+'\n'.join(errors))
    return True


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    root=a.root.resolve();manifest=a.manifest.resolve();cfg=read(manifest);campaign=cfg['campaign'];results=root/'RESULTS'/campaign;queue=root/'queues'/campaign
    files={};directories=[]
    for d in [results,queue]:
        inv=inventory(d);files.update({x['path']:x for x in inv});directories.append({'path':str(d),'files':[x['path'] for x in inv]})
    for f in [manifest]+[root/t['input_file'].replace('\\','/') for t in cfg['tasks']]:
        if f.is_file():files[str(f.resolve())]=snapshot(f)
    tasks=[]
    for t in cfg['tasks']:
        folder=results/t['task_id'];rp=folder/'RECEIPT.json';reason=[];receipt=None
        if rp.is_file():
            try:receipt=read(rp)
            except Exception as ex:reason.append('UNREADABLE_RECEIPT:'+repr(ex))
        else:reason.append('MISSING_RECEIPT')
        if not(root/t['input_file'].replace('\\','/')).is_file():reason.append('MISSING_PUBLIC_INPUT')
        for name in ['CANDIDATES.npz','HISTORY.json','PREFIXES.json']:
            if not(folder/name).is_file():reason.append('MISSING_'+name)
        if receipt and receipt.get('truth_access')!='NONE':reason.append('UNEXPECTED_ONLINE_TRUTH_DECLARATION')
        tasks.append({'task_id':t['task_id'],'case_id':t['case_id'],'method':t['method'],'receipt_status':receipt.get('status') if receipt else None,'primary_candidate':receipt.get('primary_candidate') if receipt else None,'issues':reason,'directory_exists':folder.exists()})
    q=read(queue/'RECEIPT.json') if(queue/'RECEIPT.json').is_file() else None
    seal={'schema_version':'r13_prediction_seal_v1','created_at_utc':datetime.now(timezone.utc).isoformat(),'root':str(root),'campaign':campaign,'manifest_path':str(manifest),'manifest_snapshot':cfg,'authorized_scoring':False,'authorization_note':'Only parent may change authorized_scoring to true after reviewing complete expected-task accounting. This script never authorizes.','expected_tasks':len(tasks),'task_accounting':tasks,'missing_or_invalid_tasks':sum(bool(t['issues']) for t in tasks),'failed_tasks':sum(t['receipt_status']=='FAILED' for t in tasks),'queue_status':q.get('status') if q else None,'queue_pending':q.get('pending') if q else None,'queue_running':q.get('running') if q else None,'files':list(files.values()),'directories':directories,'truth_access':'NONE','provenance_scope':'mtime_ns and byte counts only; no SHA checks; missing/failed tasks retained in expected denominator'}
    if a.output.resolve().is_relative_to(results) or a.output.resolve().is_relative_to(queue):
        raise ValueError('Place seal outside frozen result/queue directories, for example root/audit/.')
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open('x') as f:json.dump(seal,f,indent=2);f.write('\n')
    print(json.dumps({k:seal[k] for k in ['campaign','expected_tasks','missing_or_invalid_tasks','failed_tasks','queue_status','authorized_scoring']}))
if __name__=='__main__':main()

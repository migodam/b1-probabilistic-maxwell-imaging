"""Explicitly root-authorized finite confirmation staging/launch/status transport."""
from pathlib import Path
import argparse,json,subprocess,resource,time
ROOT=Path(__file__).resolve().parent
HOST='jayzh@100.121.97.39';KEY='/Users/migodam/.ssh/heatcryst_xinan_ed25519'
REMOTE='D:/AI/MAX_STRUCT_R12_20261111';LEGACY='D:/AI/MAX_STRUCT_R11_DEV_20261107'
def run(argv,label,stream=False):
    before=resource.getrusage(resource.RUSAGE_CHILDREN);cpu,t=time.process_time(),time.perf_counter()
    if stream:
        with (ROOT/(label+'.log')).open('x') as log:
            p=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT)
        output='Output retained in '+str(ROOT/(label+'.log'))
    else:
        p=subprocess.run(argv,capture_output=True,text=True,encoding='utf-8',errors='replace');output=p.stdout
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (ROOT/'CONFIRMATION_TRANSPORT_COSTS.jsonl').open('a') as f:f.write(json.dumps(dict(label=label,exit=p.returncode,
        child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        controller_CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-t))+'\n')
    if p.returncode:raise RuntimeError(output+('' if stream else p.stderr))
    return output
def ssh(command,label,stream=False):
    return run(['ssh','-i',KEY,'-o','BatchMode=yes','-o','ServerAliveInterval=30',HOST,command],label,stream)
def python(code,label):return ssh('D:/python/python.exe -c "'+code+'"',label)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=('deploy','launch','status','stop','fetch'))
    ap.add_argument('--compact',action='store_true')
    ap.add_argument('--queue',default='recovery02_secondary_registered_dev_then_confirm');a=ap.parse_args()
    if a.mode=='deploy':
        current=ROOT/'CURRENT_PUBLIC_PACKAGE.json'
        bundle=Path(json.loads(current.read_text())['current']).name if current.exists() else 'R12_EIGHT_METHOD_CONFIRMATION_PUBLIC_20261111.zip'
        code="from pathlib import Path;p=Path('"+REMOTE+"/confirmation/RUNTIME_DEV6_EIGHT_METHOD_TASKS.json');assert not p.exists(),'Existing runtime preserved: use explicit targeted source patch, not a repeat bundle deployment'"
        python(code,'refuse_repeat_runtime_deployment')
        run(['scp','-i',KEY,str(ROOT/bundle),HOST+':'+REMOTE+'/'+bundle],'upload_eightmethod_bundle')
        print(python("import zipfile;z=zipfile.ZipFile('"+REMOTE+'/'+bundle+"');z.extractall('"+REMOTE+"');print(len(z.namelist()))",'extract_public_eightmethod_bundle'))
        command='set R12_LEGACY_REPO='+LEGACY+'&& D:/python/python.exe '+REMOTE+'/confirmation/prepare_runtime.py --legacy-repo '+LEGACY+' --stage all'
        print(ssh(command,'prepare_eightmethod_runtime'))
    elif a.mode=='launch':
        # Foreground SSH retains this finite controller; no daemon or unrelated service.
        command='set R12_LEGACY_REPO='+LEGACY+'&& D:/python/python.exe '+REMOTE+'/confirmation/finite_fifo.py --execute --legacy-repo '+LEGACY+' --receipt-dir '+REMOTE+'/confirmation/queues/dev48_then_confirm192 --manifest '+REMOTE+'/confirmation/RUNTIME_DEV6_EIGHT_METHOD_TASKS.json '+REMOTE+'/confirmation/RUNTIME_CONFIRM24_EIGHT_METHOD_TASKS.json --cpus 0 2 4 6 8 10 12 14'
        print(ssh(command,'finite_fifo_remote_stdout',True))
    elif a.mode=='status':
        code="import json,pathlib;d=pathlib.Path('"+REMOTE+"/confirmation/queues');qs=[(p.parent.name,json.loads(p.read_text())) for p in d.glob('*/RECEIPT.json')];print(json.dumps(dict(queues=[dict(queue=n,status=r['status'],stage=r.get('current_stage'),done=sum(j['status']=='FINISHED' for j in r['jobs']),failed=sum(j['status']!='FINISHED' for j in r['jobs']),running=r['running'],pending=r.get('pending'),childCPU=r.get('completed_child_CPU_seconds')) for n,r in qs],canonical_valid_tasks=len({j['task_id'] for n,r in qs for j in r['jobs'] if j['status']=='FINISHED'}),failed_attempt_CPU=sum(j.get('child_process_CPU_seconds',0) for n,r in qs for j in r['jobs'] if j['status']!='FINISHED'))))"
        if a.compact:
            if not a.queue.replace('_','').isalnum():raise ValueError('Invalid own queue name')
            code="import json,pathlib;r=json.loads(pathlib.Path('"+REMOTE+"/confirmation/queues/"+a.queue+"/RECEIPT.json').read_text());print(json.dumps(dict(status=r['status'],stage=r.get('current_stage'),confirmation_done=sum(j['status']=='FINISHED' and j['stage']=='confirm24' for j in r['jobs']),failed=[dict(task=j['task_id'],error=j.get('error_log_tail','')[-1600:]) for j in r['jobs'] if j['status']!='FINISHED'],running=len(r['running']),pending=r.get('pending'),external_child_CPU=r.get('completed_child_CPU_seconds'))))"
        print(python(code,'status_'+str(time.time_ns())))
    elif a.mode=='stop':
        if not a.queue.replace('_','').isalnum():raise ValueError('Invalid own queue name')
        print(python("from pathlib import Path;Path('"+REMOTE+"/confirmation/queues/"+a.queue+"/STOP_DISPATCH').touch();print('OWN_NEW_DISPATCH_STOPPED_ACTIVE_CHILDREN_FINISH')",'stop_own_dispatch_'+a.queue))
    else:
        archive='R12_CONFIRMATION_PUBLIC_OUTPUTS.zip'
        code="import zipfile,pathlib;root=pathlib.Path('"+REMOTE+"/confirmation');z=zipfile.ZipFile('"+REMOTE+'/'+archive+"','x',zipfile.ZIP_DEFLATED);[z.write(p,str(p.relative_to(root))) for d in ['RESULTS','queues','engineering'] for p in (root/d).rglob('*') if p.is_file()];z.close();print('PACKED_PUBLIC_ONLY')"
        print(python(code,'pack_public_completed_outputs'))
        run(['scp','-i',KEY,HOST+':'+REMOTE+'/'+archive,str(ROOT/archive)],'fetch_public_completed_outputs')
if __name__=='__main__':main()

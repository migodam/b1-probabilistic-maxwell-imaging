"""Foreground SSH launch/status helper for this finite CPU cohort; no GT."""
from pathlib import Path
import argparse,json,subprocess,time,os
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=['launch','status'])
    ap.add_argument('--host',default='jayzh@100.121.97.39')
    ap.add_argument('--key',default='/Users/migodam/.ssh/heatcryst_xinan_ed25519')
    ap.add_argument('--remote-root',default='D:/AI/MAX_STRUCT_R10_2D_20261107')
    ap.add_argument('--tag',default='REMOTE',help='Separate local receipts for a separate frozen cohort')
    ap.add_argument('--manifest-name',default='twod_tasks.json')
    ap.add_argument('--config-name',default='twod_solver.json')
    ap.add_argument('--cpus',nargs='+',type=int,default=[0,2,4,6])
    ap.add_argument('--conditional-pipeline',action='store_true')
    ap.add_argument('--resume-startup-guard',action='store_true')
    args=ap.parse_args();out=ROOT/'delegated/twod'
    if not args.tag.replace('_','').isalnum():raise ValueError('Invalid receipt tag')
    remote=args.remote_root+'/research/research_max_struct_r10'
    monitor=remote+('/delegated/twod/conditional_confirmation/RECEIPT.json' if args.conditional_pipeline else '/delegated/twod/queue/RECEIPT.json')
    if args.mode=='status':
        import resource
        target=out/(args.tag+'_QUEUE_RECEIPT.json')
        command=['scp','-i',args.key,args.host+':'+monitor,str(target)]
        before=time.process_time();wall=time.perf_counter();prior=resource.getrusage(resource.RUSAGE_CHILDREN)
        p=subprocess.run(command,capture_output=True,text=True)
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        with (out/('STATUS_FETCH_COSTS.jsonl' if args.tag=='REMOTE' else args.tag+'_STATUS_FETCH_COSTS.jsonl')).open('a') as fh:
            fh.write(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),returncode=p.returncode,
                local_controller_CPU_seconds=time.process_time()-before,
                wall_seconds=time.perf_counter()-wall,
                child_transfer_CPU_seconds=after.ru_utime+after.ru_stime-prior.ru_utime-prior.ru_stime,
                scope='Metadata SCP; exact local child transfer CPU, no physics'))+'\n')
        if p.returncode:raise RuntimeError(p.stderr)
        receipt=json.loads(target.read_text())
        if args.conditional_pipeline:
            print(json.dumps(dict(status=receipt['status'],queues=[dict(name=q['name'],terminal=len(q['jobs']),status=q['status']) for q in receipt['queues']],triggers=receipt.get('triggers',[]))))
            return
        print(json.dumps(dict(status=receipt['status'],completed=len(receipt['jobs']),
            failed=sum(j['status']!='FINISHED' for j in receipt['jobs']),running=receipt['running'],
            CPU_seconds=receipt['completed_child_CPU_seconds'],wall_seconds=receipt['wall_seconds'])))
        return
    active=out/(args.tag+'_ACTIVE.json')
    log=out/('remote_queue_stdout.log' if args.tag=='REMOTE' else args.tag.lower()+'_queue_stdout.log')
    if active.exists() or log.exists():raise FileExistsError('Existing launch is preserved')
    command=['ssh','-i',args.key,'-o','BatchMode=yes','-o','ServerAliveInterval=30',args.host,
        'D:/python/python.exe '+remote+'/experiments/twod_queue.py --execute --receipt-dir '+
        remote+'/delegated/twod/queue --vendor-path '+args.remote_root+'/vendor --manifest '+remote+
        '/configs/'+args.manifest_name+' --config '+remote+'/configs/'+args.config_name+
        ' --cpus '+' '.join(map(str,args.cpus))]
    if args.conditional_pipeline:
        command[-1]='D:/python/python.exe '+remote+'/experiments/twod_conditional_pipeline.py --execute --vendor-path '+args.remote_root+'/vendor --cpus '+' '.join(map(str,args.cpus))
        if args.resume_startup_guard:command[-1]+=' --resume-startup-guard'
    start,cpu=time.perf_counter(),time.process_time()
    import resource
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    state=dict(status='LAUNCHING',started_utc=datetime.now(timezone.utc).isoformat(),command=command,
        remote_monitor=monitor,local_stdout=str(log),
        truth_access='NONE',GPU=False,local_controller_pid=os.getpid())
    with log.open('x') as fh:
        process=subprocess.Popen(command,stdout=fh,stderr=subprocess.STDOUT)
        state.update(status='REMOTE_FOREGROUND_ACTIVE',SSH_pid=process.pid)
        active.write_text(json.dumps(state,indent=2)+'\n')
        print(json.dumps(state),flush=True)
        code=process.wait()
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    state.update(status='REMOTE_PROCESS_EXITED',returncode=code,
        finished_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.perf_counter()-start,
        local_controller_CPU_seconds=time.process_time()-cpu,
        local_SSH_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime)
    active.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps(state),flush=True)
    return code


if __name__=='__main__':raise SystemExit(main())

"""Fetch terminal public predictions and metadata, never any truth container."""
from pathlib import Path
import argparse,json,subprocess,tarfile,time,resource
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cohort',choices=['baseline','denominator','right','conditional'],default='baseline')
    ap.add_argument('--host',default='jayzh@100.121.97.39')
    ap.add_argument('--key',default='/Users/migodam/.ssh/heatcryst_xinan_ed25519')
    ap.add_argument('--remote-root')
    args=ap.parse_args()
    remote=args.remote_root or {'baseline':'D:/AI/MAX_STRUCT_R10_2D_20261107','denominator':'D:/AI/MAX_STRUCT_R10_DENOM_20261107','right':'D:/AI/MAX_STRUCT_R10_RIGHT_DIAG_20261107','conditional':'D:/AI/MAX_STRUCT_R10_CONDITIONAL_20261107'}[args.cohort]
    result_subdir={'right':'right_projection_diagnostic','conditional':'conditional_confirmation'}.get(args.cohort,'twod')
    queue_subdir='delegated/twod/conditional_confirmation' if args.cohort=='conditional' else 'delegated/twod/queue'
    out=ROOT/'delegated/twod'/('fetch_'+args.cohort)
    out.mkdir(parents=True,exist_ok=False)
    cpu,wall=time.process_time(),time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    commands=[];receipt=dict(status='RUNNING',cohort=args.cohort,truth_access='NONE',new_physical_calls=0,commands=commands)
    remote_archive=remote+'/completed_'+args.cohort+'.tgz'
    arc_receipt=remote+'/completed_'+args.cohort+'_ARCHIVE_RECEIPT.json'
    code="""import os,time,json,tarfile
from pathlib import Path
import psutil
psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
cpu,wall=time.process_time(),time.perf_counter()
root=Path(__REMOTE__)/'research/research_max_struct_r10'
queue=root/'__QUEUE__'
q=json.loads((queue/'RECEIPT.json').read_text())
if q['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS','STOPPED_DISPATCH_CHILDREN_FINISHED','ALL_ACTUAL_CONDITIONAL_PREDICTIONS_SEALED'):raise RuntimeError('Queue not terminal')
archive=Path(__ARCHIVE__)
if archive.exists():raise FileExistsError(archive)
with tarfile.open(archive,'w:gz',compresslevel=1) as tar:
 tar.add(root/'RESULTS/__RESULT__',arcname='research/research_max_struct_r10/RESULTS/__RESULT__')
 tar.add(queue,arcname='research/research_max_struct_r10/delegated/twod/fetched_'+__COHORT__+'_queue')
Path(__RECEIPT__).write_text(json.dumps(dict(status='PUBLIC_PREDICTIONS_ARCHIVED',cohort=__COHORT__,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,bytes=archive.stat().st_size,truth_access='NONE',physical_calls=0),indent=2))
""".replace('__REMOTE__',repr(remote)).replace('__ARCHIVE__',repr(remote_archive)).replace('__COHORT__',repr(args.cohort)).replace('__RECEIPT__',repr(arc_receipt)).replace('__RESULT__',result_subdir).replace('__QUEUE__',queue_subdir)
    # Script is sent through stdin; avoid shell interpolation of Python source.
    try:
        command=['ssh','-i',args.key,args.host,'D:/python/python.exe -']
        p=subprocess.run(command,input=code,capture_output=True,text=True)
        commands.append(dict(command=command,returncode=p.returncode))
        (out/'archive_stdout.log').write_text(p.stdout+p.stderr)
        if p.returncode:raise RuntimeError(p.stderr)
        for source,target in [(remote_archive,out/'completed.tgz'),(arc_receipt,out/'ARCHIVE_RECEIPT.json')]:
            command=['scp','-i',args.key,args.host+':'+source,str(target)]
            p=subprocess.run(command,capture_output=True,text=True)
            commands.append(dict(command=command,returncode=p.returncode))
            if p.returncode:raise RuntimeError(p.stderr)
        archive=out/'completed.tgz'
        with tarfile.open(archive) as tar:
            for member in tar.getmembers():
                path=Path(member.name)
                target=(REPO/path).resolve()
                allowed=[ROOT/'RESULTS'/result_subdir,ROOT/('delegated/twod/fetched_'+args.cohort+'_queue')]
                if not any(target.is_relative_to(p) for p in allowed):raise ValueError('Unexpected archive path')
                if member.issym() or member.islnk():raise ValueError('Archive contains link')
                if member.isfile() and target.exists():raise FileExistsError('Existing result preserved: '+str(target))
            tar.extractall(REPO,filter='data')
        receipt['status']='FETCHED_PUBLIC_PREDICTIONS_NO_SCORING'
    except BaseException as exc:
        receipt.update(status='FAILED_PRESERVED',error=repr(exc));raise
    finally:
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt.update(local_controller_CPU_seconds=time.process_time()-cpu,
            local_SSH_SCP_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            wall_seconds=time.perf_counter()-wall,
            remote_archive_CPU='ARCHIVE_RECEIPT.json; additive separately from local controller/SSH CPU')
        (out/'FETCH_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))


if __name__=='__main__':main()

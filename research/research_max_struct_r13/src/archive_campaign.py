"""Archive terminal campaign evidence without deleting or rerunning science."""
import argparse,json,time,zipfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--campaign',required=True);a=p.parse_args()
    r=a.root;name=a.campaign;receipt=json.loads((r/'queues'/name/'RECEIPT.json').read_text())
    if receipt['pending'] or receipt['running']:raise RuntimeError('Queue not terminal')
    dest=r/'archives';dest.mkdir(exist_ok=True);target=dest/(name+'_full_results.zip')
    cpu,wall=time.process_time(),time.perf_counter();items=[]
    dirs=[r/'RESULTS'/name,r/'queues'/name]
    if name=='development_v1':dirs += [r/'queues/development_v1_recovery',r/'audit/aborted_controller_children']
    with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=3) as z:
        for d in dirs:
            for f in sorted(d.rglob('*')):
                if f.is_file():z.write(f,f.relative_to(r));items.append(dict(name=str(f.relative_to(r)),bytes=f.stat().st_size))
        p=r/'audit/CONTROLLER_FAILURE_AND_RECOVERY.json'
        if name=='development_v1' and p.exists():z.write(p,p.relative_to(r))
    record=dict(campaign=name,archive=str(target.relative_to(r)),bytes=target.stat().st_size,files=len(items),entries=items,
        CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_physics=0,source_deleted=False)
    (dest/(name+'_ARCHIVE_RECEIPT.json')).write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='entries'}),flush=True)


if __name__=='__main__':main()

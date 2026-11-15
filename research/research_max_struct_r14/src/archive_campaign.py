"""Archive a terminal finite campaign, in bounded files, without deleting it."""
import argparse,json,time,zipfile
from pathlib import Path


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--campaign',required=True)
    p.add_argument('--part-bytes',type=int,default=650000000);a=p.parse_args();r=a.root
    q=r/'queues'/a.campaign;receipt=json.loads((q/'RECEIPT.json').read_text())
    if receipt.get('running') or receipt.get('pending_tasks',0):raise RuntimeError('Campaign still has active tasks')
    dirs=[r/'RESULTS'/a.campaign,q];files=[f for d in dirs for f in sorted(d.rglob('*')) if f.is_file()]
    config=r/'configs'/(a.campaign+'.json')
    if config.exists():files.append(config)
    dest=r/'archives';dest.mkdir(exist_ok=True);parts=[];archive=None;raw_bytes=0;cpu=time.process_time();wall=time.perf_counter()
    try:
        for f in files:
            size=f.stat().st_size
            if archive is None or raw_bytes+size>a.part_bytes:
                if archive is not None:
                    archive.close();parts[-1]['archive_bytes']=(r/parts[-1]['path']).stat().st_size
                path=dest/f'{a.campaign}_part{len(parts)+1:02d}.zip';archive=zipfile.ZipFile(path,'x',zipfile.ZIP_DEFLATED,compresslevel=1)
                raw_bytes=0;parts.append(dict(path=str(path.relative_to(r)),members=[],archive_bytes=None))
            compression=zipfile.ZIP_STORED if f.suffix.lower() in ('.npz','.npy','.zip','.png','.pdf') else zipfile.ZIP_DEFLATED
            archive.write(f,str(f.relative_to(r)),compress_type=compression);raw_bytes+=size
            parts[-1]['members'].append(dict(path=str(f.relative_to(r)),bytes=size))
    finally:
        if archive is not None:
            archive.close();parts[-1]['archive_bytes']=(r/parts[-1]['path']).stat().st_size
    result=dict(status='ARCHIVED_NO_DELETION',campaign=a.campaign,parts=parts,files=len(files),
                CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                new_physics=0,hash_checks='NONE',originals_retained=True)
    (dest/(a.campaign+'_ARCHIVE_RECEIPT.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='parts'}))
    print(json.dumps([dict(path=p['path'],bytes=p['archive_bytes']) for p in parts]))


if __name__=='__main__':main()

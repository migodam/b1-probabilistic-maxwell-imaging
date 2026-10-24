"""Upload only fsynced closed evidence parts to an existing draft Release.

Uses the existing gh account, never reads credentials, never hashes or changes
assets, and never starts an experiment. Resumes by exact asset name and size.
The delivery receipt is separate from the already closed research ledger.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

REPO='migodam/b1-probabilistic-maxwell-imaging'


def write(path,obj):
    tmp=path.with_suffix(path.suffix+'.pending')
    tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(path)


def remote_assets(tag):
    r=subprocess.run(['gh','release','view',tag,'--repo',REPO,'--json','assets'],capture_output=True,text=True,check=True,timeout=90)
    return {v['name']:v['size'] for v in json.loads(r.stdout)['assets']}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--tag',default='research-p3lr4-20261025');p.add_argument('--max-wall',type=float,default=43200)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();records=[];uploaded={};status='ERROR';error=None
    try:
        while time.perf_counter()-start<a.max_wall:
            markers=sorted(a.bundle.glob('closed_parts/*.json'))
            for marker in markers:
                row=json.loads(marker.read_text());name=row['file']
                if name in uploaded:continue
                part=a.bundle/name
                if part.stat().st_size!=row['bytes']:raise ValueError('closed part size mismatch: '+name)
                existing=remote_assets(a.tag)
                if name in existing:
                    if existing[name]!=row['bytes']:raise ValueError('existing asset size mismatch; no clobber: '+name)
                    result='ALREADY_UPLOADED_NAME_SIZE_MATCH'
                else:
                    t=time.perf_counter();result='FAILED'
                    for attempt in range(3):
                        remaining=a.max_wall-(time.perf_counter()-start)
                        if remaining<=0:raise TimeoutError('delivery wall cap')
                        r=subprocess.run(['gh','release','upload',a.tag,str(part),'--repo',REPO],capture_output=True,text=True,timeout=remaining)
                        records.append(dict(file=name,attempt=attempt+1,returncode=r.returncode,stderr=r.stderr[-2000:]))
                        if r.returncode==0:result='UPLOADED';break
                        existing=remote_assets(a.tag)
                        if existing.get(name)==row['bytes']:result='UPLOADED_VERIFIED_AFTER_UNCERTAIN_RETURN';break
                        if name in existing:raise ValueError('uncertain asset has wrong size; no overwrite')
                        time.sleep(min(30,remaining))
                    if result=='FAILED':raise RuntimeError('upload failed after 3 attempts: '+name)
                    records.append(dict(file=name,upload_wall_seconds=time.perf_counter()-t))
                uploaded[name]=dict(bytes=row['bytes'],status=result)
                write(a.out/'PROGRESS.json',dict(status='UPLOADING',assets=uploaded,uploaded_bytes=sum(v['bytes'] for v in uploaded.values()),wall_seconds=time.perf_counter()-start))
            receipt=a.bundle/'RECEIPT.json'
            if receipt.exists():
                done=json.loads(receipt.read_text())
                if done['status']!='COMPLETE':raise RuntimeError('packager closed unsuccessfully; preserve partial assets')
                expected={v['file']:v['bytes'] for v in done['parts']}
                if set(expected)==set(uploaded):
                    existing=remote_assets(a.tag)
                    if any(existing.get(k)!=n for k,n in expected.items()):raise ValueError('final remote name/size audit failed')
                    status='ALL_PARTS_UPLOADED_NAME_SIZE_VERIFIED';break
            time.sleep(15)
        else:raise TimeoutError('delivery wall cap')
    except BaseException as exc:error=repr(exc);raise
    finally:
        write(a.out/'UPLOAD_RECEIPT.json',dict(status=status,error=error,tag=a.tag,assets=uploaded,attempts=records,
            wall_seconds=time.perf_counter()-start,RHS=0,training=0,hashes_used=False,no_clobber=True,
            scope='delivery only; does not alter scientific qualification'))


if __name__=='__main__':main()

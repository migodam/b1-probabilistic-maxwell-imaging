"""Windows/Python-stdlib streaming upload to an EXISTING GitHub draft release.

Token: one line from stdin only, in memory. No create/delete/publish endpoints,
asset replacement, retries, hashes, scientific computation, or local GB copy.
"""
from __future__ import annotations
import argparse,http.client,json,mimetypes,re,sys,time
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode

REPO='migodam/b1-probabilistic-maxwell-imaging'
VERSION='2026-03-10'


def append(path,row):
    with path.open('a',encoding='utf8') as stream:stream.write(json.dumps(row,ensure_ascii=False)+'\n')
    print(json.dumps(row,ensure_ascii=False),flush=True)


def headers(token):
    return {'Authorization':'Bearer '+token,'Accept':'application/vnd.github+json',
        'X-GitHub-Api-Version':VERSION,'User-Agent':'MAX-STRUCT-R14-stream-release'}


def request_json(path,token,timeout):
    connection=http.client.HTTPSConnection('api.github.com',timeout=timeout)
    try:
        connection.request('GET',path,headers=headers(token));response=connection.getresponse()
        status=response.status;payload=response.read(8*1024*1024)
        if status!=200:raise RuntimeError('GitHubmetadataHTTP'+str(status))
        return json.loads(payload),status
    finally:connection.close()


def asset_list(data,archive_dir):
    entries=data.get('parts',data.get('assets',[]));result=[];names=set()
    if not entries:raise ValueError('No actual archive parts/assets in manifest')
    for item in entries:
        name=item.get('archive',item.get('name'))
        if not name or Path(name).name!=name or any(c in name for c in '/\\\r\n') or name in names:
            raise ValueError('Require unique plain asset filenames')
        names.add(name);path=(archive_dir/name).resolve()
        if not path.is_relative_to(archive_dir.resolve()) or not path.is_file():raise FileNotFoundError('Declaredarchiveabsent')
        expected=int(item['bytes']);actual=path.stat().st_size
        if expected!=actual:raise ValueError('Assetsize differs from executed packaging receipt: '+name)
        result.append(dict(name=name,path=path,bytes=actual))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo',default=REPO);p.add_argument('--release-id',required=True,type=int)
    p.add_argument('--archive-dir',type=Path,required=True);p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True);p.add_argument('--timeout',type=float,default=900.)
    a=p.parse_args()
    if a.repo!=REPO:raise ValueError('This helper is scoped to the existing project repository')
    if a.release_id<=0:raise ValueError('Existing positive release ID required')
    # Validate all local metadata before touching credentials or the network.
    assets=asset_list(json.loads(a.manifest.read_text()),a.archive_dir)
    if a.receipt.exists():raise FileExistsError('Use a new uploadattempt receipt filename')
    a.receipt.parent.mkdir(parents=True,exist_ok=True)
    raw=sys.stdin.buffer.readline(8192).strip()
    if not raw or not re.fullmatch(rb'[A-Za-z0-9_]+',raw):raise ValueError('Invalid or absent stdin credential; value suppressed')
    token=raw.decode('ascii');del raw
    start=time.perf_counter();cpu=time.process_time();failed=0;rows=[]
    try:
        release,status=request_json(f'/repos/{a.repo}/releases/{a.release_id}',token,a.timeout)
        if release.get('draft') is not True:raise ValueError('Target release is not a draft; refusing any upload')
        existing={};page=1
        while True:
            batch,_=request_json(f'/repos/{a.repo}/releases/{a.release_id}/assets?per_page=100&page={page}',token,a.timeout)
            for item in batch:existing[item['name']]=item
            if len(batch)<100:break
            page+=1
        append(a.receipt,dict(event='DRAFT_VERIFIED',release_id=a.release_id,repo=a.repo,HTTP_status=status,
            draft=True,asset_count=len(assets),existing_assets=len(existing),created_utc=datetime.now(timezone.utc).isoformat(),
            token_transport='stdin memoryonly;neverarguments/environment/files/logs',retry_policy='NONE',publication='NOT_PERFORMED'))
        for asset in assets:
            name=asset['name'];t=time.perf_counter();c=time.process_time();sent=0;connection=None;http_status=None
            row=dict(event='ASSET',name=name,expected_bytes=asset['bytes'],actual_file_bytes=asset['bytes'],
                uploaded_bytes_sent=0,HTTP_status=None,url=None,repo=a.repo,release_id=a.release_id)
            try:
                if name in existing:
                    old=existing[name];same=int(old.get('size',-1))==asset['bytes'] and old.get('state')=='uploaded'
                    row.update(status='SKIPPED_EXISTING_SAME_SIZE_NOT_CONTENT_VERIFIED' if same else 'FAILED_EXISTING_NAME_COLLISION',
                        existing_asset_id=old.get('id'),existing_bytes=old.get('size'),url=old.get('browser_download_url'))
                    failed+=not same
                else:
                    connection=http.client.HTTPSConnection('uploads.github.com',timeout=a.timeout)
                    path=f'/repos/{a.repo}/releases/{a.release_id}/assets?'+urlencode({'name':name})
                    connection.putrequest('POST',path)
                    for key,value in headers(token).items():connection.putheader(key,value)
                    connection.putheader('Content-Type',mimetypes.guess_type(name)[0] or 'application/octet-stream')
                    connection.putheader('Content-Length',str(asset['bytes']));connection.endheaders()
                    with asset['path'].open('rb') as stream:
                        while block:=stream.read(1024*1024):connection.send(block);sent+=len(block)
                    response=connection.getresponse();http_status=response.status
                    payload=response.read(8*1024*1024)
                    if http_status==201:
                        published=json.loads(payload)
                        good=published.get('name')==name and published.get('size')==asset['bytes'] and published.get('state')=='uploaded' and sent==asset['bytes']
                        row.update(status='UPLOADED_TO_DRAFT' if good else 'FAILED_RESPONSE_METADATA_MISMATCH',
                            server_bytes=published.get('size'),asset_id=published.get('id'),url=published.get('browser_download_url'))
                        failed+=not good
                        if good:existing[name]=published
                    else:
                        row.update(status='FAILED_HTTP_NO_RETRY');failed+=1
            except Exception as exc:
                # Do not serialize raw HTTP bodies, headers or exception values.
                row.update(status='FAILED_EXCEPTION_NO_RETRY',error_type=type(exc).__name__);failed+=1
            finally:
                if connection is not None:connection.close()
                row.update(uploaded_bytes_sent=sent,sent_byte_scope='Completed HTTPS send calls only; failed partial-call bytes unknown',HTTP_status=http_status,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t)
                rows.append(row);append(a.receipt,row)
        append(a.receipt,dict(event='BATCH_END',status='COMPLETE_DRAFT_UPLOADS' if not failed else 'COMPLETE_WITH_ERRORS',
            assets=len(rows),failed=failed,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start,
            process_CPU_seconds_including_imports=time.process_time(),published=False,deleted_assets=0,overwritten_assets=0,retries=0,
            original_archives_retained=True,hash_operations='NONE'))
        return 0 if not failed else 1
    except Exception as exc:
        append(a.receipt,dict(event='BATCH_FAILED',error_type=type(exc).__name__,
            status='NO_RELEASE_CREATION_OR_OVERWRITE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-start))
        return 1
    finally:
        token=None # Lifetime ends with this finite process; no disk cache.


if __name__=='__main__':raise SystemExit(main())

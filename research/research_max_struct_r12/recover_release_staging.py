"""Resume authorized R12 evidence packaging after physical SSD reconnection.

All newly created files are on the verified external volume. Internal copies,
original ZIPs, interrupted wrappers and scientific sources are retained.
"""
from pathlib import Path
import argparse,filecmp,json,os,resource,shutil,subprocess,time,zipfile

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'release_assets';AUDIT=ROOT/'audit/release'
REC=AUDIT/'recovery_20261112'
INTERNAL=Path('/tmp/r12_release_staging')

def guard():
    assert ROOT.is_dir()and(ROOT/'package_release.py').is_file()
    mount=subprocess.check_output(['mount'],text=True)
    assert " on /Volumes/migodam's-external-brain (apfs"in mount,'Expected APFS mount missing; never write placeholder mount path'
def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True);temp=path.with_name(path.name+'.tmp');temp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n');os.replace(temp,path)
def log(row):
    REC.mkdir(parents=True,exist_ok=True)
    with(REC/'OPERATIONS.jsonl').open('a')as f:f.write(json.dumps(row)+'\n')
def run(argv,label):
    before=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=time.process_time();wall=time.perf_counter();p=subprocess.run(argv,capture_output=True,text=True);after=resource.getrusage(resource.RUSAGE_CHILDREN)
    row=dict(label=label,returncode=p.returncode,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    log(row)
    if p.returncode:raise RuntimeError(p.stderr[-5000:]+p.stdout[-5000:])
    return p.stdout,row
def copy_equal(src,dest,label):
    guard();cpu=time.process_time();wall=time.perf_counter()
    if dest.exists():
        if not filecmp.cmp(src,dest,shallow=False):raise ValueError('Existing destination differs; preserve both and inspect: '+str(dest))
        action='REUSED_BYTE_EQUAL'
    else:
        dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
        if not filecmp.cmp(src,dest,shallow=False):raise ValueError('Copy byte mismatch')
        action='COPIED_AND_BYTE_COMPARED'
    log(dict(label=label,action=action,source=str(src),destination=str(dest),bytes=src.stat().st_size,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,hashes='NOT_COMPUTED',source_deleted=False))
def migrate():
    guard();REC.mkdir(parents=True,exist_ok=True)
    for src in sorted(INTERNAL.iterdir()):
        if not src.is_file():continue
        if src.name in {'max-struct-r12-remote-01.zip.zst','max-struct-r12-remote-02.zip.zst','max-struct-r12-remote-01.zip','max-struct-r12-remote-03.zip.zst.partial'}:
            copy_equal(src,OUT/src.name,'reuse_internal_'+src.name)
        elif src.suffix in {'.json','.jsonl'}:
            dest=REC/src.name
            if not dest.exists():copy_equal(src,dest,'retain_internal_receipt_'+src.name)
        elif src.suffix in {'.py','.md'}:
            dest=REC/'operational_handoff'/src.name
            if not dest.exists():copy_equal(src,dest,'retain_operational_handoff_'+src.name)
    # Preserve prior source/clock receipts under distinct names, do not replace.
    for name,target in [('PLAN.json','REMOTE_PLAN_RECOVERED.json'),('PACKAGING_RECEIPT.json','REMOTE_PACKAGING_RECEIPT_RECOVERED.json')]:
        if(INTERNAL/name).exists()and not(AUDIT/target).exists():copy_equal(INTERNAL/name,AUDIT/target,'recover_'+target)
    dump(REC/'MIGRATION_COMPLETE.json',dict(status='COMPLETE',all_new_files_external=True,internal_originals_preserved=True,hashes='NOT_COMPUTED',new_scientific_forward_calls=0))
def validate_zip(path,index):
    data=json.loads(index.read_text())
    with zipfile.ZipFile(path)as z:actual=[(i.filename,i.file_size)for i in z.infolist()]
    expected=[(e['archive_path'],e['bytes'])for e in data['entries']]
    if actual!=expected:raise ValueError('ZIP directory/index mismatch')
    return len(actual)
def ready_manifest():
    assets=[]
    # Migrated remote01 was decoded and directory-verified before interruption;
    # copy_equal additionally compared source/destination bytes after remount.
    prior=REC/'VERIFICATION_RECEIPT.json'
    remote=[]
    if prior.exists():remote+=json.loads(prior.read_text())['assets']
    extra=REC/'REMOTE_VERIFIED_ASSETS.json'
    if extra.exists():remote+=json.loads(extra.read_text())['assets']
    if remote:
        for x in {a['name']:a for a in remote}.values():
            name=x['name'];p=OUT/name;index=AUDIT/(name[:-8]+'.index.json');zpath=OUT/name[:-4]
            if p.exists()and zpath.exists():
                count=validate_zip(zpath,index);assets.append(dict(name=name,path='release_assets/'+name,upload_local_path=str(p),bytes=p.stat().st_size,format='zstandard-wrapped ZIP64',content_index=str(index.relative_to(ROOT)),file_count=count,ready=True))
    local=REC/'LOCAL_ZSTD_RECOVERED_RECEIPT.json'
    if local.exists():
        for x in json.loads(local.read_text())['archives']:
            if x.get('selected')and not x.get('decoded_byte_equality'):continue
            name=x['name']if x['selected']else x['zip_name'];p=OUT/name;index=AUDIT/(x['zip_name'][:-4]+'.index.json');count=json.loads(index.read_text())['file_count']
            if not x['selected']:validate_zip(p,index)
            assets.append(dict(name=name,path='release_assets/'+name,upload_local_path=str(p),bytes=p.stat().st_size,format='zstandard-wrapped ZIP64'if x['selected']else'ZIP64',content_index=str(index.relative_to(ROOT)),file_count=count,ready=True))
    assets.sort(key=lambda x:x['name']);expected=len(list(AUDIT.glob('max-struct-r12-*.index.json')));complete=len(assets)==expected
    data=dict(status='READY_ALL_VERIFIED'if complete else'PARTIAL_READY_VERIFIED',assets=assets,pending=None if complete else'Other logical assets still downloading/validating; only ready=true entries are accepted',publication_transport_manifest='PUBLICATION_PARTS_MANIFEST.json',publication_scope='These are logical evidence wrappers; root splits/uploads128MiB transport parts without altering them',source_files_deleted=0,published=False,hashes='NOT_COMPUTED')
    dump(ROOT/'UPLOAD_MANIFEST_READY.json',data)
    lines=['# R12已验收逻辑证据包（逐步更新）','',f"状态{data['status']}；下表每包已验收。未列包仍pending，不能当作整组完成。平台上传采用root生成的128MiB运输parts，另见PUBLICATION_PARTS_MANIFEST.json，逻辑包不因split而改变。",'', '|逻辑包|字节|文件条目|索引|','|---|---:|---:|---|']
    for x in assets:lines.append(f"|`{x['name']}`|{x['bytes']}|{x['file_count']}|[{Path(x['content_index']).name}]({x['content_index']})|")
    lines+=['','只上传UPLOAD_MANIFEST_READY.json列出的正式文件，勿上传ZIP副本、partial或操作临时文件。完整源文件与历史中断产物均保留，尚未发布。Git上的最新科学文稿/图/代码优先；大包中的历史文稿仅是执行时快照，不用为小文稿修订重打原始证据包。']
    (ROOT/'ARTIFACTS.md').write_text('\n'.join(lines)+'\n')
    dump(ROOT/'ARTIFACTS.json',dict(status=data['status'],assets=assets,remaining_status=data['pending'],full_indices='audit/release/max-struct-r12-*.index.json',source_files_deleted=0,hashes='NOT_COMPUTED'))
    return data
def verify_available_remote():
    guard();meta=json.loads((AUDIT/'REMOTE_ZSTD_RECEIPT.json').read_text());rows=[]
    for a in meta['archives']:
        source=OUT/a['name'];zpath=OUT/a['zip_name'];index=AUDIT/(a['zip_name'][:-4]+'.index.json')
        if not source.exists():continue
        if source.stat().st_size!=a['bytes']:raise ValueError('Formal compressed source size mismatch')
        if not zpath.exists():run(['zstd','-q','-d',str(source),'-o',str(zpath)],'decode_available_remote_'+source.name)
        count=validate_zip(zpath,index);rows.append(dict(name=a['name'],compressed_bytes=source.stat().st_size,decoded_zip_bytes=zpath.stat().st_size,entries=count,central_directory_names_sizes_exact=True,hashes='NOT_COMPUTED'))
        dump(REC/'REMOTE_VERIFIED_ASSETS.json',dict(status='PARTIAL_READY'if len(rows)<len(meta['archives'])else'COMPLETE',assets=rows));ready_manifest();print(json.dumps(rows[-1]),flush=True)
def local_wrappers():
    guard();old=AUDIT/'LOCAL_ZSTD_RECEIPT.json';prior={x['zip_name']:x for x in json.loads(old.read_text())['archives']}if old.exists()else{};rows=[]
    for index in sorted(AUDIT.glob('max-struct-r12-local*.index.json')):
        guard();a=json.loads(index.read_text());src=OUT/a['name'];target=OUT/(a['name']+'.zst');validate_zip(src,index)
        if target.exists():
            # Existing wrapper may have been interrupted at USB termination.
            try:_,check=run(['zstd','-q','-t',str(target)],'existing_wrapper_decode_check_'+target.name)
            except RuntimeError:
                preserved=target.with_name(target.name+'.interrupted_20261111')
                if preserved.exists():raise ValueError('Prior interruption copy exists; inspect, never overwrite')
                target.rename(preserved);log(dict(label='preserve_interrupted_wrapper',path=str(preserved),bytes=preserved.stat().st_size,scientific_source_changed=False))
        if not target.exists():
            _,cost=run(['zstd','-q','-6','--long=27','--no-check','-T1',str(src),'-o',str(target)],'sequential_singlethread_local_zstd_'+src.name)
            row=dict(zip_name=src.name,name=target.name,bytes=target.stat().st_size,zip_bytes=src.stat().st_size,CPU_seconds=cost['child_CPU_seconds'],wall_seconds=cost['wall_seconds'],threads=1,window_log=27,selected=target.stat().st_size<src.stat().st_size*.95,recovered_after_USB_disconnect=True)
        else:
            row=prior.get(src.name,dict(zip_name=src.name,name=target.name,bytes=target.stat().st_size,zip_bytes=src.stat().st_size,CPU_seconds=None,wall_seconds=None,selected=target.stat().st_size<src.stat().st_size*.95))
        if row['bytes']!=target.stat().st_size:raise ValueError('Existing recorded wrapper size mismatch')
        if row['selected']:
            # Full decoded byte comparison to the original ZIP; no hashes.
            cpu=time.process_time();wall=time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN)
            p=subprocess.Popen(['zstd','-q','-d','-c',str(target)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            equal=True
            with src.open('rb')as f:
                while True:
                    buf=p.stdout.read(1024**2)
                    if buf!=f.read(len(buf)):equal=False;break
                    if not buf:
                        if f.read(1):equal=False
                        break
            if not equal:p.terminate()
            error=p.stderr.read();rc=p.wait();after=resource.getrusage(resource.RUSAGE_CHILDREN)
            if rc or not equal:raise ValueError('Decoded wrapper equality failed: '+error.decode('utf8','replace'))
            row.update(decoded_byte_equality=True,verification_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime+time.process_time()-cpu,verification_wall_seconds=time.perf_counter()-wall,hashes='NOT_COMPUTED')
        rows.append(row);dump(REC/'LOCAL_ZSTD_RECOVERED_RECEIPT.json',dict(status='RUNNING',archives=rows,source_files_deleted=0));ready_manifest()
        print(json.dumps(dict(name=target.name,bytes=row['bytes'],selected=row['selected'],decoded_byte_equality=row.get('decoded_byte_equality'))),flush=True)
    if old.exists()and not(REC/'LOCAL_ZSTD_BEFORE_RECOVERY.json').exists():shutil.copy2(old,REC/'LOCAL_ZSTD_BEFORE_RECOVERY.json')
    receipt=dict(status='COMPLETE',archives=rows,source_files_deleted=0,scientific_forward_calls=0,scope='Sequential, at most1compression thread; complete decoded bytes compared to original ZIPs; interruption copies preserved')
    dump(old,receipt);dump(REC/'LOCAL_ZSTD_RECOVERED_RECEIPT.json',receipt)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['migrate','local-wrappers','ready','verify-remote']);a=p.parse_args()
    if a.mode=='migrate':migrate()
    elif a.mode=='ready':print(json.dumps(ready_manifest()))
    elif a.mode=='verify-remote':verify_available_remote()
    else:local_wrappers()

"""Prepare lossless R12 release evidence; never publishes or deletes sources.

No content hashes. Reuses R11 paths/byte-size indices, ZIP64 and explicit costs.
Remote mode only packages a frozen file plan, reads finished scientific outputs,
and writes archives under a dedicated export directory. SSH credentials are
arguments/environment configuration, never archived key contents.
"""
from pathlib import Path
import argparse,json,os,subprocess,resource,sys,time,zipfile,tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT=Path(__file__).resolve().parent
PREFIX=Path('research/research_max_struct_r12')
MAX_ASSET=1_800_000_000
GROUP_LIMIT=1_650_000_000
OUT=ROOT/'release_assets'
AUDIT=ROOT/'audit/release'
DENIED_PARTS={'.git','.ssh','.codex','__pycache__','node_modules','.venv','.venv_nn','primary_texts'}
DENIED_NAMES={'.env','id_rsa','id_ed25519','authorized_keys','known_hosts','chat_history.json','conversation.json','conversations.json','credentials.json'}
TRANSPORT_DUPLICATES={'R12_SELECTED_PUBLIC_OUTPUTS_20261111.zip','R12_CONFIRMATION_PUBLIC_OUTPUTS.zip'}

def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def denied(p):
    return bool(set(p.parts)&DENIED_PARTS)or p.name.lower()in DENIED_NAMES or p.suffix in {'.pyc','.pem','.key','.exp'}
def groups(entries,stem):
    result=[];current=[];size=0
    for e in sorted(entries,key=lambda x:x['archive_path']):
        if e['bytes']>GROUP_LIMIT:raise ValueError('Single file exceeds group cap; explicit lossless chunk plan required: '+e['archive_path'])
        if current and size+e['bytes']>GROUP_LIMIT:
            result.append(current);current=[];size=0
        current.append(e);size+=e['bytes']
    if current:result.append(current)
    return [dict(name=f'max-struct-r12-{stem}-{i:02d}.zip',entries=g,source_bytes=sum(x['bytes']for x in g))for i,g in enumerate(result,1)]
def run(argv,label):
    before=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=time.process_time();wall=time.perf_counter();p=subprocess.run(argv,capture_output=True,text=True,encoding='utf8',errors='replace');after=resource.getrusage(resource.RUSAGE_CHILDREN)
    AUDIT.mkdir(parents=True,exist_ok=True)
    with (AUDIT/'TRANSPORT_COSTS.jsonl').open('a')as f:f.write(json.dumps(dict(label=label,returncode=p.returncode,child_CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))+'\n')
    if p.returncode:raise RuntimeError(label+': '+p.stderr[-5000:]+p.stdout[-5000:])
    return p.stdout

REMOTE_WORKER=r'''from pathlib import Path
import json,time,zipfile,sys
folder=Path(sys.argv[1]);plan=json.loads((folder/'PLAN.json').read_text());rows=[];cpu=time.process_time();wall=time.perf_counter()
for asset in plan['archives']:
 out=folder/asset['name'];index=out.with_suffix('.index.json')
 if out.exists()and index.exists():
  old=json.loads(index.read_text())
  if old['archive_bytes']!=out.stat().st_size or old['entries']!=asset['entries']:raise RuntimeError('Existing archive does not match frozen plan')
  rows.append(old);print(json.dumps(dict(asset=out.name,resumed=True,bytes=out.stat().st_size)),flush=True);continue
 if out.exists():raise RuntimeError('Partial archive preserved; inspect and choose explicit new name, never overwrite')
 c,t=time.process_time(),time.perf_counter()
 with zipfile.ZipFile(out,'x',compression=zipfile.ZIP_STORED,allowZip64=True)as z:
  for e in asset['entries']:
   p=Path(e['source_root'])/e['source_relative']
   if not p.is_file()or p.stat().st_size!=e['bytes']:raise RuntimeError('Source size changed: '+str(p))
   z.write(p,e['archive_path'])
 if out.stat().st_size>=plan['max_asset_bytes']:raise RuntimeError('Asset cap exceeded')
 with zipfile.ZipFile(out)as z:
  actual=[(x.filename,x.file_size)for x in z.infolist()]
  expected=[(e['archive_path'],e['bytes'])for e in asset['entries']]
  if actual!=expected:raise RuntimeError('Central directory index mismatch')
 row=dict(name=out.name,archive_bytes=out.stat().st_size,file_count=len(asset['entries']),source_bytes=asset['source_bytes'],entries=asset['entries'],CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,hashes='NOT_COMPUTED',verification='ZIP central-directory names/counts/byte sizes')
 index.write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(json.dumps(dict(asset=out.name,bytes=row['archive_bytes'],files=row['file_count'])),flush=True)
receipt=dict(status='COMPLETE',archives=[{k:v for k,v in x.items()if k!='entries'}for x in rows],archive_phase_CPU_seconds_sum=sum(x['CPU_seconds']for x in rows),CPU_seconds=time.process_time()-cpu,process_CPU_seconds_near_exit=time.process_time(),wall_seconds=time.perf_counter()-wall,scientific_forward_calls=0,source_files_deleted=0,hashes='NOT_COMPUTED')
if(folder/'PACKAGING_RECEIPT.json').exists():
 with(folder/'PACKAGING_RESUMES.jsonl').open('a')as f:f.write(json.dumps(receipt)+'\n')
else:(folder/'PACKAGING_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
'''

ZSTD_WORKER=r'''from pathlib import Path
import json,time,sys,zstandard
folder=Path(sys.argv[1]);plan=json.loads((folder/'PLAN.json').read_text());receipt=folder/'ZSTD_RECEIPT.json';old=json.loads(receipt.read_text())if receipt.exists()else dict(archives=[]);prior={x['zip_name']:x for x in old['archives']};rows=[]
params=zstandard.ZstdCompressionParameters.from_level(6,window_log=27,enable_ldm=True,threads=2,write_checksum=0)
for asset in plan['archives']:
 src=folder/asset['name'];out=folder/(asset['name']+'.zst')
 if out.exists()and asset['name']in prior:
  row=prior[asset['name']]
  if out.stat().st_size!=row['bytes']:raise RuntimeError('Existing zstd size mismatch')
  rows.append(row);continue
 if out.exists():raise RuntimeError('Partial zstd preserved; explicit inspection required')
 c,t=time.process_time(),time.perf_counter()
 with src.open('rb')as f,out.open('xb')as target:zstandard.ZstdCompressor(compression_params=params).copy_stream(f,target)
 row=dict(zip_name=src.name,name=out.name,bytes=out.stat().st_size,zip_bytes=src.stat().st_size,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,version=zstandard.__version__,window_log=27,enable_ldm=True,threads=2,frame_checksum=False)
 rows.append(row);receipt.write_text(json.dumps(dict(status='COMPLETE'if len(rows)==len(plan['archives'])else'RUNNING',archives=rows,scientific_forward_calls=0,source_files_deleted=0),indent=2)+'\n');print(json.dumps(row),flush=True)
receipt.write_text(json.dumps(dict(status='COMPLETE',archives=rows,scientific_forward_calls=0,source_files_deleted=0),indent=2)+'\n')
'''

def prepare_remote(inventory,legacy_inventory):
    x=json.loads(Path(inventory).read_text());entries=[];excluded=[]
    for e in x['files']:
        p=Path(e['path'])
        if denied(p)or p.name in TRANSPORT_DUPLICATES:
            excluded.append(dict(path=e['path'],reason='cache/private-key/raw-measurement/duplicate-derived-transport'));continue
        if p.parts[0]=='confirmation':
            target=PREFIX/'delegated/confirmation/remote_complete'/Path(*p.parts[1:])
        elif str(p).startswith('research/research_max_struct_r12/RESULTS/basin_v1/'):
            target=PREFIX/'RESULTS/basin_remote'/Path(*p.parts[4:])
        elif str(p).startswith('research/research_max_struct_r12/'):
            target=PREFIX/'remote_runtime'/p
        elif p.suffix=='.zip':
            # Original staged input is retained once remotely, local matching
            # deployment copy is listed as a logical transport duplicate.
            target=PREFIX/'remote_runtime/input_bundles'/p.name
        else:
            excluded.append(dict(path=e['path'],reason='outside scientific runtime allowlist'));continue
        entries.append(dict(source_root=x['root'],source_relative=e['path'],archive_path=target.as_posix(),bytes=e['bytes'],origin='finished_remote_R12'))
    if legacy_inventory:
        y=json.loads(Path(legacy_inventory).read_text())
        for e in y['files']:
            p=Path(e['path'])
            if denied(p):continue
            entries.append(dict(source_root=y['root'],source_relative=e['path'],archive_path=(PREFIX/'remote_runtime/shared_legacy'/p).as_posix(),bytes=e['bytes'],origin='actual_shared_legacy_runtime_source_once'))
    plan=dict(format=1,max_asset_bytes=MAX_ASSET,group_source_limit_bytes=GROUP_LIMIT,compression='ZIP_STORED (NPZ already compressed)',archives=groups(entries,'remote'),excluded=excluded,hashes='NOT_COMPUTED',scientific_writers='Root explicitly confirmed finished remote scientific processes',private_chats_keys_repo_history=False)
    dump(AUDIT/'REMOTE_PLAN.json',plan);return plan

def remote_fetch(args):
    plan=prepare_remote(args.inventory,args.legacy_inventory);OUT.mkdir(exist_ok=True);AUDIT.mkdir(parents=True,exist_ok=True)
    worker=AUDIT/'remote_pack_worker.py';worker.write_text(REMOTE_WORKER)
    ssh=['ssh','-i',args.key,'-o','BatchMode=yes','-o','ServerAliveInterval=30',args.host]
    scp=['scp','-i',args.key,'-o','BatchMode=yes']
    folder=args.remote_export.rstrip('/')
    run(ssh+[args.remote_python+' -c "from pathlib import Path;Path(\''+folder+'\').mkdir(parents=True,exist_ok=True)"'],'remote_export_directory')
    run(scp+[str(worker),args.host+':'+folder+'/remote_pack_worker.py'],'upload_pack_worker')
    run(scp+[str(AUDIT/'REMOTE_PLAN.json'),args.host+':'+folder+'/PLAN.json'],'upload_frozen_remote_plan')
    output=run(ssh+[args.remote_python+' '+folder+'/remote_pack_worker.py '+folder],'remote_lossless_pack_or_resume')
    with (AUDIT/'REMOTE_PACK_STDOUT.log').open('a')as f:f.write(output)
    zstd_meta={}
    if args.zstd:
        worker=AUDIT/'remote_zstd_worker.py';worker.write_text(ZSTD_WORKER)
        run(scp+[str(worker),args.host+':'+folder+'/remote_zstd_worker.py'],'upload_existing_zstd_worker')
        text=run(ssh+[args.remote_python+' '+folder+'/remote_zstd_worker.py '+folder],'remote_long_window_lossless_zstd')
        with(AUDIT/'REMOTE_ZSTD_STDOUT.log').open('a')as f:f.write(text)
        run(scp+[args.host+':'+folder+'/ZSTD_RECEIPT.json',str(AUDIT/'REMOTE_ZSTD_RECEIPT.json')],'fetch_zstd_receipt')
        zstd_meta={x['zip_name']:x for x in json.loads((AUDIT/'REMOTE_ZSTD_RECEIPT.json').read_text())['archives']}
    jobs=[]
    for a in plan['archives']:
        indexname=Path(a['name']).with_suffix('.index.json').name;index=AUDIT/indexname
        if not index.exists():run(scp+[args.host+':'+folder+'/'+indexname,str(index)],'fetch_'+indexname)
        metadata=json.loads(index.read_text());jobs.append((a,metadata,zstd_meta.get(a['name'])))
    with ThreadPoolExecutor(max_workers=args.download_workers)as pool:
        futures=[pool.submit(stream_download,args,folder+'/'+(z['name']if z else a['name']),OUT/(z['name']if z else a['name']),z['bytes']if z else m['archive_bytes'])for a,m,z in jobs]
        for future in futures:future.result()
    for a,m,z in jobs:
        if z and not(OUT/a['name']).exists():run(['zstd','-q','-d',str(OUT/z['name']),'-o',str(OUT/a['name'])],'local_lossless_decode_'+a['name'])
        verify_archive(OUT/a['name'],a['entries'])
    run(scp+[args.host+':'+folder+'/PACKAGING_RECEIPT.json',str(AUDIT/'REMOTE_PACKAGING_RECEIPT.json')],'fetch_remote_packaging_receipt')

def stream_download(args,remote_path,local,expected):
    """Resumable SSH binary stream avoids Windows SFTP request latency.

    Existing partial bytes are retained. SSH transport provides channel integrity;
    no content hashes are computed. Final size and ZIP directory are checked.
    """
    if local.exists()and local.stat().st_size==expected:return
    partial=local.with_suffix(local.suffix+'.partial')
    if local.exists():
        if partial.exists():raise RuntimeError('Two partial downloads preserved; inspect explicitly')
        local.rename(partial)
    offset=partial.stat().st_size if partial.exists()else 0
    if offset>expected:raise ValueError('Partial exceeds remote archive size')
    before=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=time.process_time();wall=time.perf_counter()
    code="import sys,shutil;f=open('"+remote_path+"','rb');f.seek("+str(offset)+");shutil.copyfileobj(f,sys.stdout.buffer,1048576)"
    argv=['ssh','-T','-i',args.key,'-o','BatchMode=yes','-o','ServerAliveInterval=30','-o','Compression=no',args.host,args.remote_python+' -c "'+code+'"']
    with partial.open('ab')as target,tempfile.TemporaryFile()as err:
        p=subprocess.Popen(argv,stdout=target,stderr=err)
        if hasattr(os,'wait4'):
            _,status,usage=os.wait4(p.pid,0);p.returncode=os.waitstatus_to_exitcode(status);child_cpu=usage.ru_utime+usage.ru_stime
        else:p.wait();child_cpu=None
        err.seek(0);error=err.read()
    row=dict(label='binary_stream_'+local.name,returncode=p.returncode,start_offset=offset,received_bytes=partial.stat().st_size-offset,expected_total_bytes=expected,child_CPU_seconds=child_cpu,controller_CPU_seconds=None,controller_CPU_scope='Concurrent controller clock overlaps other downloads; not independently additive',wall_seconds=time.perf_counter()-wall)
    with (AUDIT/'TRANSPORT_COSTS.jsonl').open('a')as f:f.write(json.dumps(row)+'\n')
    if p.returncode:raise RuntimeError(error.decode('utf8','replace'))
    if partial.stat().st_size!=expected:raise ValueError('Stream byte count mismatch; partial retained')
    partial.rename(local);print(json.dumps(dict(downloaded=local.name,bytes=expected,resumed_bytes=offset,wall_seconds=row['wall_seconds'])),flush=True)

def verify_archive(path,entries):
    if path.stat().st_size>=MAX_ASSET:raise ValueError('Asset >=1.8GB')
    with zipfile.ZipFile(path)as z:
        actual=[(i.filename,i.file_size)for i in z.infolist()]
    expected=[(e['archive_path'],e['bytes'])for e in entries]
    if actual!=expected:raise ValueError('Archive central-directory mismatch: '+str(path))

def local_pack():
    OUT.mkdir(exist_ok=True);AUDIT.mkdir(parents=True,exist_ok=True);entries=[];retained=[];excluded=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(ROOT)
        if rel.parts[0]=='release_assets'or rel.parts[:2]==('audit','release')or denied(rel):continue
        if rel.parts[:2]==('audit','remote_snapshot'):
            excluded.append(dict(path=rel.as_posix(),reason='transport mirror: canonical basin_remote/remote assets retained'));continue
        if p.name in TRANSPORT_DUPLICATES:
            excluded.append(dict(path=rel.as_posix(),reason='derived selected-output transport duplicate, expanded canonical outputs retained'));continue
        n=p.stat().st_size;is_bundle=p.suffix in {'.zip','.tgz','.tar','.gz','.7z'}
        release=(p.suffix in {'.npz','.npy'}and(n>8*1024**2 or p.name=='states.npz'))or is_bundle or n>50*1024**2
        row=dict(source_root=str(ROOT),source_relative=rel.as_posix(),archive_path=(PREFIX/rel).as_posix(),bytes=n,origin='local_scientific_evidence')
        if release:entries.append(row)
        else:retained.append(dict(path=(PREFIX/rel).as_posix(),bytes=n,placement='git_code_metadata_figures_small_predictions'))
    plan=dict(format=1,max_asset_bytes=MAX_ASSET,archives=groups(entries,'local'),git_retained=retained,excluded=excluded,hashes='NOT_COMPUTED')
    dump(AUDIT/'LOCAL_PLAN.json',plan);cpu=time.process_time();wall=time.perf_counter();rows=[]
    for a in plan['archives']:
        target=OUT/a['name'];index=AUDIT/target.with_suffix('.index.json').name
        if target.exists():verify_archive(target,a['entries']);rows.append(json.loads(index.read_text()));continue
        c,t=time.process_time(),time.perf_counter()
        with zipfile.ZipFile(target,'x',compression=zipfile.ZIP_STORED,allowZip64=True)as z:
            for e in a['entries']:
                p=Path(e['source_root'])/e['source_relative']
                if p.stat().st_size!=e['bytes']:raise ValueError('Source size changed')
                z.write(p,e['archive_path'])
        verify_archive(target,a['entries'])
        row=dict(name=target.name,archive_bytes=target.stat().st_size,source_bytes=a['source_bytes'],file_count=len(a['entries']),entries=a['entries'],CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,hashes='NOT_COMPUTED',verification='ZIP central-directory names/counts/byte sizes')
        dump(index,row);rows.append(row);print(json.dumps(dict(asset=target.name,bytes=target.stat().st_size,files=len(a['entries']))),flush=True)
    receipt=dict(status='COMPLETE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,archives=[{k:v for k,v in x.items()if k!='entries'}for x in rows],scientific_forward_calls=0,source_files_deleted=0)
    if (AUDIT/'LOCAL_PACKAGING_RECEIPT.json').exists():
        with (AUDIT/'LOCAL_PACKAGING_RESUMES.jsonl').open('a')as f:f.write(json.dumps(receipt)+'\n')
    else:dump(AUDIT/'LOCAL_PACKAGING_RECEIPT.json',receipt)
    marker='# BEGIN R12 GENERATED RELEASE ASSETS';end='# END R12 GENERATED RELEASE ASSETS';p=ROOT/'.gitignore';old=p.read_text()if p.exists()else''
    if marker in old:old=old.split(marker)[0]+old.split(end,1)[1]
    block=[marker,'/release_assets/','**/__pycache__/','*.pyc']+['/'+e['source_relative']for e in entries]+['/'+e['path']+'/'for e in []]+['/audit/remote_snapshot/',end]
    p.write_text(old.rstrip()+'\n\n'+'\n'.join(block)+'\n')

def compact_predictions():
    import numpy as np
    cpu=time.process_time();wall=time.perf_counter();plans=[json.loads(p.read_text())for p in AUDIT.glob('LOCAL*PLAN.json')];rows=[]
    for a in [a for plan in plans for a in plan['archives']]:
        for e in a['entries']:
            p=Path(e['source_root'])/e['source_relative']
            if p.name.lower()not in {'prediction.npz','derived.npz'}:continue
            out=p.with_name(p.stem+'_compact.npz')
            with np.load(p,allow_pickle=False)as z:
                keys=[k for k in ['m','prediction','raw','initial','Y_prediction','gradient','centers','train_mask','withheld_mask','held_sources']if k in z.files]
                if not out.exists():np.savez_compressed(out,**{k:z[k]for k in keys})
            rows.append(dict(original=e['archive_path'],compact=(PREFIX/out.relative_to(ROOT)).as_posix(),keys=keys,bytes=out.stat().st_size,selection='exact saved endpoint arrays; candidates remain in archived original; no numerical solve'))
    dump(AUDIT/'COMPACT_PREDICTIONS_RECEIPT.json',dict(status='COMPLETE',rows=rows,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,scientific_forward_calls=0))

def supplement():
    """Small safe Git reduction: preserve all original receipts and document links."""
    import re
    existing={e['source_relative']for p in AUDIT.glob('max-struct-r12-local*.index.json')for e in json.loads(p.read_text())['entries']};linked=set()
    for md in ROOT.rglob('*.md'):
        if 'release_assets'in md.parts:continue
        for value in re.findall(r'\]\(([^)]+)\)',md.read_text(errors='replace')):
            if '://'not in value:linked.add((md.parent/value.split('#')[0]).resolve())
    entries=[];kept=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(ROOT)
        if rel.parts[0]=='release_assets'or rel.parts[:2]in {('audit','release'),('audit','remote_snapshot')}or denied(rel)or rel.as_posix()in existing:continue
        n=p.stat().st_size;name=p.name.lower();candidate=p.suffix=='.npz'and name in {'paid_materials.npz','candidates.npz','derived.npz','prediction.npz'}and n>1_000_000
        history=name in {'history.json','corrector_history.json','auxiliary_history.json'}and n>1_000_000
        transport=p.suffix in {'.zip','.tgz','.7z'}
        if not(candidate or history or transport):continue
        if p.resolve()in linked:
            kept.append(dict(path=rel.as_posix(),bytes=n,reason='Direct existing Markdown link preserved; no original overwritten'));continue
        entries.append(dict(source_root=str(ROOT),source_relative=rel.as_posix(),archive_path=(PREFIX/rel).as_posix(),bytes=n,origin='local_scientific_evidence'))
    plan=dict(format=1,archives=groups(entries,'local-supplement'),kept_direct_links=kept,original_receipts_kept_in_git=True);dump(AUDIT/'LOCAL_SUPPLEMENT_PLAN.json',plan);cpu=time.process_time();wall=time.perf_counter()
    for a in plan['archives']:
        out=OUT/a['name'];index=AUDIT/out.with_suffix('.index.json').name
        if out.exists():verify_archive(out,a['entries']);continue
        c,t=time.process_time(),time.perf_counter()
        with zipfile.ZipFile(out,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True)as z:
            for e in a['entries']:z.write(ROOT/e['source_relative'],e['archive_path'])
        verify_archive(out,a['entries']);dump(index,dict(name=out.name,archive_bytes=out.stat().st_size,source_bytes=a['source_bytes'],file_count=len(a['entries']),entries=a['entries'],CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,hashes='NOT_COMPUTED'))
        for e in a['entries']:
            p=ROOT/e['source_relative'];locator=os.path.relpath(ROOT/'ARTIFACTS.md',p.parent).replace(os.sep,'/');idx=os.path.relpath(index,p.parent).replace(os.sep,'/')
            note=f"# 原始轨迹的归档入口\n\n原件保持本地原位，未覆盖。原路径：`{e['archive_path']}`；字节：{e['bytes']}。\n\nRelease资产：`{out.name}`，同名仓库相对路径可无损恢复。见[资产总索引]({locator})及[逐文件索引]({idx})。大prediction/derived另有邻近`*_compact.npz`小终点；本入口不代替原始receipt。\n"
            p.with_name(p.name+'.archive-index.md').write_text(note)
    with(ROOT/'.gitignore').open('a')as f:f.write('\n# Supplemental paid candidates/history/transports; original bytes remain in release assets.\n'+'\n'.join('/'+e['source_relative']for e in entries)+'\n')
    dump(AUDIT/'SUPPLEMENT_PACKAGING_RECEIPT.json',dict(status='COMPLETE',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,source_bytes=sum(e['bytes']for e in entries),files=len(entries),original_receipts_preserved=True,kept_direct_links=kept,scientific_forward_calls=0))
    print(json.dumps(dict(supplement_files=len(entries),source_bytes=sum(e['bytes']for e in entries),kept_direct_links=len(kept))))

def compress_local():
    """Existing zstd CLI, lossless long-window wrapper; original ZIPs retained."""
    rows=[]
    for index in sorted(AUDIT.glob('max-struct-r12-local*.index.json')):
        a=json.loads(index.read_text());src=OUT/a['name'];out=OUT/(a['name']+'.zst')
        if not out.exists():
            before=resource.getrusage(resource.RUSAGE_CHILDREN);cpu=time.process_time();wall=time.perf_counter()
            run(['zstd','-q','-6','--long=27','--no-check','-T2',str(src),'-o',str(out)],'local_long_window_zstd_'+src.name)
            after=resource.getrusage(resource.RUSAGE_CHILDREN)
            row=dict(zip_name=src.name,name=out.name,bytes=out.stat().st_size,zip_bytes=src.stat().st_size,CPU_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,controller_CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,selected=out.stat().st_size<src.stat().st_size*.95,window_log=27,threads=2,frame_checksum=False)
        else:
            old=AUDIT/'LOCAL_ZSTD_RECEIPT.json';prior={x['zip_name']:x for x in json.loads(old.read_text())['archives']}if old.exists()else{}
            row=prior.get(src.name,dict(zip_name=src.name,name=out.name,bytes=out.stat().st_size,zip_bytes=src.stat().st_size,CPU_seconds=None,wall_seconds=None,selected=out.stat().st_size<src.stat().st_size*.95))
        rows.append(row);dump(AUDIT/'LOCAL_ZSTD_RECEIPT.json',dict(status='RUNNING',archives=rows,source_files_deleted=0,scientific_forward_calls=0));print(json.dumps(row),flush=True)
    dump(AUDIT/'LOCAL_ZSTD_RECEIPT.json',dict(status='COMPLETE',archives=rows,source_files_deleted=0,scientific_forward_calls=0))

def manifest():
    assets=[]
    zstdpath=AUDIT/'REMOTE_ZSTD_RECEIPT.json';zstd_meta={x['zip_name']:x for x in json.loads(zstdpath.read_text())['archives']}if zstdpath.exists()else{}
    local_zstd=AUDIT/'LOCAL_ZSTD_RECEIPT.json'
    if local_zstd.exists():zstd_meta.update({x['zip_name']:x for x in json.loads(local_zstd.read_text())['archives']if x.get('selected')})
    for p in sorted(AUDIT.glob('max-struct-r12-*.index.json')):
        a=json.loads(p.read_text());local=OUT/a['name']
        if not local.exists():raise FileNotFoundError(local)
        verify_archive(local,a['entries']);z=zstd_meta.get(a['name']);asset=OUT/z['name']if z else local
        if not asset.exists()or(z and asset.stat().st_size!=z['bytes']):raise ValueError('Compressed asset missing/size mismatch')
        assets.append(dict(name=asset.name,relative_asset_path='release_assets/'+asset.name,bytes=asset.stat().st_size,file_count=a['file_count'],source_bytes=a['source_bytes'],index='audit/release/'+p.name,entries=a['entries'],format='zstandard-wrapped ZIP64'if z else'ZIP64',embedded_zip_name=a['name']if z else None,uncompressed_zip_bytes=local.stat().st_size))
    localplan=json.loads((AUDIT/'LOCAL_PLAN.json').read_text())if(AUDIT/'LOCAL_PLAN.json').exists()else{}
    remoteplan=json.loads((AUDIT/'REMOTE_PLAN.json').read_text())if(AUDIT/'REMOTE_PLAN.json').exists()else{}
    archived_local={e['source_relative']for a in assets for e in a['entries']if e['origin']=='local_scientific_evidence'};retained=[];unplanned_large=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(ROOT)
        if rel.parts[0]=='release_assets'or rel.parts[:2]==('audit','remote_snapshot')or denied(rel)or rel.as_posix()in archived_local or p.name in TRANSPORT_DUPLICATES:continue
        n=p.stat().st_size
        if p.suffix in {'.npz','.npy'}and n>8*1024**2:unplanned_large.append(dict(path=rel.as_posix(),bytes=n));continue
        retained.append(dict(path=(PREFIX/rel).as_posix(),bytes=n,placement='git_code_metadata_figures_small_predictions'))
    m=dict(format=1,status='PREPARED_NOT_PUBLISHED',asset_limit_bytes=MAX_ASSET,assets=assets,asset_bytes=sum(a['bytes']for a in assets),archived_file_entries=sum(a['file_count']for a in assets),git_retained=retained,unplanned_large_local=unplanned_large,excluded_local=localplan.get('excluded',[]),excluded_remote=remoteplan.get('excluded',[]),hashes='NOT_COMPUTED_USER_INSTRUCTION',source_files_deleted=0,private_chat_key_repo_history_archived=False,remote_runtime_shared_sources='Canonical source snapshot once under remote_runtime/shared_legacy; no git history or installed libraries',restore='ZIP entries use repo-relative paths; extract into an empty review checkout/directory. remote_complete retains full remote confirmation outputs without overwriting compact Git selections. Source roots in indices record actual provenance.')
    dump(ROOT/'ARTIFACTS.json',m)
    upload=dict(status='PREPARED_NOT_UPLOADED',assets=[dict(name=a['name'],path=a['relative_asset_path'],bytes=a['bytes'],format=a['format'],content_index=a['index'])for a in assets],metadata=[],exclusions='Do not upload ZIP originals behind selectedzip.zst, .partial, incomplete_transfers/, workers, plans, temporary files or source duplicates',original_zip_files_preserved_locally=True,published=False)
    upload['metadata']=[dict(path='ARTIFACTS.json',purpose='full sourcepath/byte/content mapping'),dict(path='ARTIFACTS.md',purpose='reader restoration index')]+[dict(path=a['index'],purpose='perasset content index')for a in assets]
    dump(ROOT/'UPLOAD_MANIFEST.json',upload)
    lines=['# R12 原始证据归档（准备完成，尚未发布）','', '所有源文件保持原位。代码、JSON/CSV/MD、图片与小预测留Git；巨大数组及运行bundle保存在分包资产。单包严格小于1.8GB。未计算内容hash；核对ZIP中央目录的路径、条数、字节和实际包体大小。','', '|资产相对路径|字节|文件条目|内容索引|','|---|---:|---:|---|']
    for a in assets:lines.append(f"|`{a['relative_asset_path']}`|{a['bytes']}|{a['file_count']}|[{Path(a['index']).name}]({a['index']})|")
    lines+=['',f"包体总计{m['asset_bytes']}字节，{m['archived_file_entries']}个条目。完整映射在[ARTIFACTS.json](ARTIFACTS.json)。",'', '远端basin原始states恢复到RESULTS/basin_remote；confirmation完整RESULTS/queues/engineering在delegated/confirmation/remote_complete；实际部署的R12源与共享legacy源在remote_runtime中。单个源文件只在对应canonical档案计一次；展开后的运行bundle本身作为历史输入保留，其内部共享源不是又一次科学计算。没有私人聊天、密钥文件、原始Fresnel测量或仓库全history。','', '恢复：.zip.zst先用`zstd -d <asset>.zip.zst -o <asset>.zip`无损解码；然后在空目录运行`python -m zipfile -e <asset>.zip <empty-destination>`。普通.zip直接第二步。归档路径均相对仓库，不要解到正在写入的科学工作树。','', '费用和运输见audit/release/{LOCAL_PACKAGING_RECEIPT,REMOTE_PACKAGING_RECEIPT,REMOTE_ZSTD_RECEIPT}.json及TRANSPORT_COSTS.jsonl。远端压缩只读取已结束任务，写入专用export目录，不触碰其他会话。没有发布、push或删除源。']
    (ROOT/'ARTIFACTS.md').write_text('\n'.join(lines)+'\n');print(json.dumps(dict(assets=len(assets),bytes=m['asset_bytes'],entries=m['archived_file_entries'],status=m['status'])))

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='mode',required=True);sub.add_parser('local');sub.add_parser('manifest');sub.add_parser('compact');sub.add_parser('supplement');sub.add_parser('compress-local');q=sub.add_parser('remote');q.add_argument('--inventory',required=True);q.add_argument('--legacy-inventory');q.add_argument('--host',default=os.getenv('R12_RELEASE_SSH_HOST'));q.add_argument('--key',default=os.getenv('R12_RELEASE_SSH_KEY'));q.add_argument('--remote-export',required=True);q.add_argument('--remote-python',default='D:/python/python.exe');q.add_argument('--download-workers',type=int,default=3);q.add_argument('--zstd',action='store_true');a=p.parse_args()
    if a.mode=='local':local_pack()
    elif a.mode=='manifest':manifest()
    elif a.mode=='compact':compact_predictions()
    elif a.mode=='supplement':supplement()
    elif a.mode=='compress-local':compress_local()
    else:
        if not a.host or not a.key:p.error('Remote mode needs explicit SSH host/key path; no key contents are read into reports')
        remote_fetch(a)
if __name__=='__main__':main()

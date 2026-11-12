"""Plan/package R13-owned release artifacts; never upload, publish or run science.

Default mode=plan creates metadata only. --mode prepare may split existing terminal
archives/runtime bundles. Private post-seal inputs require --include-postseal and
explicit parent authorization covering every declared independent campaign seal.
--mode downloads writes active DOWNLOADS.json only with --published and actual
uploaded URLs/bytes from a verification receipt, followed by public size readback.
"""
import argparse,csv,json,os,re,stat,time,urllib.request,zipfile
from pathlib import Path,PurePosixPath
from datetime import datetime,timezone

PART_BYTES=1_000_000_000
TERMINAL={'COMPLETED','COMPLETED_WITH_ERRORS','COMPLETED_AFTER_EXPLICIT_CONTROLLER_RECOVERY'}
BLOCKED_PARTS={'.git','.venv','venv','__pycache__','.cache','node_modules','scorer_fixture','private_operational'}
BLOCKED_FILES={'PRO_CHAT_READBACK.json','CHAT_READBACK.json','.DS_Store','.env','credentials.json','token.json','id_rsa','id_ed25519'}


def load(path):return json.loads(Path(path).read_text())


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def exclusion(name):
    path=PurePosixPath(str(name).replace('\\','/'));parts=path.parts
    if path.is_absolute() or '..' in parts:return 'absolute/traversal path'
    if any(parts[i:i+2] in [('audit','transfers'),('audit','remote_metadata')] for i in range(max(0,len(parts)-1))) or path.name=='REMOTE_TRANSFER_PROGRESS.json':return 'Private operational SSH/controller/transfer metadata excluded from all publication assets'
    if any(p in BLOCKED_PARTS or 'scorer_fixture' in p.lower() for p in parts):return 'git/environment/cache/fake fixture'
    if path.name in BLOCKED_FILES or any(x in path.name.lower() for x in ['chat_readback','private_chat','conversation_readback','thread_readback']) or path.name.lower().startswith(('.env.','credentials.','secret.')):return 'private chat payload or sensitive filename'
    if path.suffix.lower() in ['.pem','.key','.p12','.pfx','.exp']:return 'sensitive key/original third-party raw measurement'
    if path.suffix.lower()=='.pdf' and not('FIGURES' in parts or 'figures' in parts or path.name=='MAX_STRUCT_Research_20261112.pdf'):return 'third-party/full-paper PDF excluded unless explicit own figure/Pro research report'
    return None


def inspect_zip(path):
    bad=[]
    with zipfile.ZipFile(path) as z:
        for item in z.infolist():
            why=exclusion(item.filename)
            if stat.S_ISLNK(item.external_attr>>16):why='symlink member'
            if why:bad.append(dict(name=item.filename,reason=why))
    return bad


def terminal_assets(root):
    candidates=[];missing=[];recovery=set()
    for receipt in (root/'queues').glob('*/RECEIPT.json'):
        q=load(receipt)
        if q.get('recovery_source'):recovery.add((receipt.parent/str(q['recovery_source']).replace('\\','/')).resolve())
    for receipt in (root/'queues').glob('*/RECEIPT.json'):
        if receipt.resolve() in recovery:continue
        q=load(receipt);campaign=receipt.parent.name
        if q.get('status') not in TERMINAL or q.get('pending') or q.get('running'):continue
        archive_name=q.get('archive_filename',campaign+'_full_results.zip')
        if Path(archive_name).name!=archive_name or not archive_name.endswith('.zip'):raise ValueError('Invalid explicit campaign archive filename')
        path=root/'archives'/archive_name
        if path.is_file():candidates.append(dict(path=path,kind='terminal_campaign',campaign=campaign,terminal_receipt=str(receipt)))
        else:missing.append(dict(campaign=campaign,reason='Terminal queue exists but complete campaign archive absent',expected_path=str(path)))
        # Composed physics retains its immutable failed-initial archive as a
        # separate asset; never replace it with the recovered final ZIP.
        if campaign=='physics_v1' and q.get('recovery_source') and archive_name!='physics_v1_full_results.zip':
            initial=root/'archives/physics_v1_full_results.zip'
            attempt=receipt.parent/'attempts/initial_RECEIPT.json'
            if initial.is_file() and attempt.is_file():
                original=load(attempt)
                if original.get('status') not in TERMINAL:raise ValueError('Initial physics attempt receipt is not terminal')
                candidates.append(dict(path=initial,kind='immutable_failed_initial_campaign',campaign='physics_v1_initial_attempt',terminal_receipt=str(attempt)))
            else:missing.append(dict(campaign='physics_v1_initial_attempt',reason='Composed physics requires original immutable archive AND initial attempt receipt for accessible failed-attempt logs',expected_path=str(initial),expected_receipt=str(attempt)))
    candidates += [dict(path=f,kind='frozen_runtime',campaign=None,terminal_receipt=None) for f in sorted(root.glob('runtime*.zip'))]
    return candidates,missing


def postseal_authorization(root,path):
    if path is None:raise PermissionError('--include-postseal requires explicit parent packaging authorization')
    auth=load(path)
    if auth.get('authorized_release_packaging') is not True:raise PermissionError('Parent has not authorized release packaging of post-seal data')
    declared=set()
    for f in (root/'configs').glob('*.json'):
        d=load(f);name=d.get('campaign','')
        if name.startswith('confirmation_') or name=='physics_v1':declared.add(name)
    required=set(auth.get('required_campaigns',[]))
    if not required or not declared<=required:raise PermissionError('Authorization must cover every declared independent/physics campaign: '+', '.join(sorted(declared-required)))
    from seal_results import validate_seal
    verified=[]
    for name in sorted(required):
        p=Path(auth.get('seals',{}).get(name,root/'audit'/(name+'_seal.json')))
        if not p.is_absolute():p=root/p
        seal=load(p)
        if seal.get('campaign')!=name:raise ValueError('Seal campaign mismatch '+name)
        validate_seal(seal,require_authorization=True);verified.append(dict(campaign=name,seal=str(p),expected_tasks=seal.get('expected_tasks')))
    return verified


def gather(root,directories,extensions=None):
    files=[];excluded=[]
    for rel in directories:
        start=root/rel
        if not start.exists():continue
        paths=[start] if start.is_file() else sorted(start.rglob('*'))
        for f in paths:
            if not f.is_file():continue
            name=f.relative_to(root).as_posix();why='symlink' if f.is_symlink() else 'nested transport/archive excluded from source snapshots' if f.suffix.lower() in ['.zip','.zst'] or re.search(r'\.part\d+$',f.name) else exclusion(name)
            if why:excluded.append(dict(path=name,reason=why));continue
            if extensions and f.suffix.lower() not in extensions:continue
            files.append(f)
    return sorted(set(files)),excluded


def create_zip(root,out,name,files,extra_json=None):
    path=out/name;entries=[];cpu,wall=time.process_time(),time.perf_counter()
    with zipfile.ZipFile(path,'x',zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
        for f in files:
            before=f.stat();name=f.relative_to(root).as_posix();z.write(f,name);after=f.stat()
            if before.st_size!=after.st_size or before.st_mtime_ns!=after.st_mtime_ns:raise RuntimeError('Source changed during archive: '+str(f))
            entries.append(dict(name=name,bytes=before.st_size,mtime_ns=before.st_mtime_ns))
        for member,value in (extra_json or{}).items():
            encoded=(json.dumps(value,indent=2)+'\n').encode();z.writestr(member,encoded);entries.append(dict(name=member,bytes=len(encoded),generated_metadata=True))
    index=out/(path.name+'.index.json');dump(index,dict(archive=path.name,bytes=path.stat().st_size,files=len(entries),entries=entries,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,new_physics=0,source_deleted=False))
    return path,index


def transport(path,out,max_bytes):
    before=path.stat();parts=[]
    if before.st_size<=max_bytes:return [dict(name=path.name,bytes=before.st_size,order=1,local_path=str(path.resolve()))]
    directory=out/'parts';directory.mkdir(exist_ok=True)
    with path.open('rb') as source:
        index=1
        while source.tell()<before.st_size:
            target=directory/(path.name+f'.part{index:03d}');written=0
            with target.open('xb') as dest:
                while written<max_bytes:
                    block=source.read(min(8*1024*1024,max_bytes-written))
                    if not block:break
                    dest.write(block);written+=len(block)
            parts.append(dict(name=target.name,bytes=written,order=index,local_path=str(target.resolve())));index+=1
    after=path.stat()
    if sum(p['bytes'] for p in parts)!=before.st_size or after.st_size!=before.st_size or after.st_mtime_ns!=before.st_mtime_ns:raise RuntimeError('Incomplete split or changed source '+str(path))
    return parts


def verified_downloads(root,args):
    if not args.published:raise PermissionError('Active DOWNLOADS.json requires explicit --published after actual upload')
    if not args.transport_manifest or not args.upload_verification:raise ValueError('Transport manifest and actual upload verification receipt required')
    manifest=load(args.transport_manifest);verification=load(args.upload_verification);uploaded={r['name']:r for r in verification.get('parts',[])};logical=[]
    for asset in manifest['logical_assets']:
        public=dict(name=asset['name'],bytes=asset['bytes'],parts=[])
        for part in asset['parts']:
            record=uploaded.get(part['name'],{});url=record.get('url','')
            if record.get('verified') is not True or record.get('bytes')!=part['bytes']:raise ValueError('Missing verified actual upload '+part['name'])
            if not url.startswith('https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/') or '?' in url:raise ValueError('Use stable public existing-project Release URLs, not guesses or temporary credential URLs')
            request=urllib.request.Request(url,headers={'Range':'bytes=0-0','User-Agent':'R13-artifact-size-readback'})
            with urllib.request.urlopen(request) as response:
                content_range=response.headers.get('Content-Range','');match=re.search(r'/(\d+)$',content_range);size=int(match.group(1)) if match else int(response.headers.get('Content-Length','-1'));response.read(1)
            if size!=part['bytes']:raise RuntimeError('Public readback size mismatch '+part['name'])
            public['parts'].append(dict(name=part['name'],bytes=part['bytes'],url=url,order=part['order'],url_readback_verified=True))
        logical.append(public)
    target=(args.downloads or root/'DOWNLOADS.json').resolve()
    if not target.is_relative_to(root):raise ValueError('Active R13 download index must stay under the R13 study root')
    dump(target,dict(status='PUBLISHED_PUBLIC_SIZE_READBACK_VERIFIED',published=True,created_at_utc=datetime.now(timezone.utc).isoformat(),logical_assets=logical,part_size_bytes=manifest['part_size_bytes'],order='parts listed in binary concatenation order',hashes='not computed; exact byte sizes and public size readback only',upload_verification_source=str(args.upload_verification),uploaded_by_this_script=False))
    print(json.dumps(dict(status='ACTIVE_DOWNLOAD_INDEX_CREATED',path=str(target),logical_assets=len(logical),uploaded_by_this_script=False)))


def main():
    process_cpu0,process_wall0=time.process_time(),time.perf_counter()
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path);p.add_argument('--mode',choices=['plan','prepare','downloads'],default='plan');p.add_argument('--include-postseal',action='store_true');p.add_argument('--authorization',type=Path);p.add_argument('--part-bytes',type=int,default=PART_BYTES);p.add_argument('--transport-manifest',type=Path);p.add_argument('--upload-verification',type=Path);p.add_argument('--published',action='store_true');p.add_argument('--downloads',type=Path);a=p.parse_args();root=a.root.resolve()
    if a.mode=='downloads':return verified_downloads(root,a)
    if not a.out:raise ValueError('--out NEW_DIRECTORY required')
    if not 1<=a.part_bytes<=PART_BYTES:raise ValueError('Transport part limit must be positive and at most1,000,000,000bytes')
    out=a.out.resolve()
    if out.is_relative_to(root.parent/'research_max_struct_r12'):raise ValueError('R13 tooling never writes into R12 directories')
    out.mkdir(parents=True,exist_ok=False);assets,missing=terminal_assets(root)
    if len({x['path'].name for x in assets})!=len(assets):raise ValueError('Logical archive filenames must be unique')
    blocked=[]
    for asset in assets:
        bad=inspect_zip(asset['path'])
        if bad:blocked.append(dict(asset=asset['path'].name,excluded_entries=bad))
    plan=dict(status='PLAN_ONLY_NO_ARCHIVES_OR_PARTS_WRITTEN',existing_assets=[dict(name=x['path'].name,bytes=x['path'].stat().st_size,kind=x['kind'],campaign=x['campaign'],local_path=str(x['path'])) for x in assets],missing_terminal_archives=missing,policy_blocked_existing_archives=blocked,postseal_packages=['max-struct-r13-postseal-inputs-sources.zip','max-struct-r13-postseal-analysis-tables.zip'],part_size_bytes=a.part_bytes,public_urls_created=False,private_inputs_archived=False,new_physics=0)
    if a.mode=='plan':dump(out/'PACKAGE_PLAN.json',plan);print(json.dumps(dict(status=plan['status'],existing_assets=len(assets),missing_terminal_archives=len(missing),blocked_archives=len(blocked))));return
    if missing or blocked:raise RuntimeError('Prepare refused: missing terminal archives or excluded/private members; inspect plan-mode output first')
    verified=[];exclusions=[];indices={}
    source=dict(repository='https://github.com/migodam/b1-probabilistic-maxwell-imaging',R13_directory='research/research_max_struct_r13',base_checkout='ec3a37da6f3ee90beb61e753067d4f639f4118f5',source_freeze_metadata=[str(p.relative_to(root)) for p in(root/'audit').glob('*FREEZE*.json')],runtime_bundles=[p.name for p in root.glob('runtime*.zip')],NN_weights_reference='research/research_max_struct_r8/RESULTS/training/physics/weights.pt; original R8 dependency, not a new R13 training',legacy_reference='D:/AI/MAX_STRUCT_R11_DEV_20261107 was the execution dependency; restore equivalent repository sources rather than reproducing that absolute path')
    if a.include_postseal:
        verified=postseal_authorization(root,a.authorization)
        inputs,excluded=gather(root,['data_public','sealed_truth','sources/pro/MAX_STRUCT_20261112','sources/user','src','runtime','experiments','configs','audit','delegated/full_current','delegated/physics','delegated/evidence','THEORY_DERIVATIONS.md','TASK_RESUME.md','REPRODUCE.md','ARTIFACTS.md','RESEARCH_REPORT_ZH.md','PAPER_DRAFT.md']);exclusions+=excluded
        selection_indices,excluded=gather(root,['predictions/min_change_paid_pool'],extensions={'.json'});exclusions+=excluded
        for f in selection_indices:
            if f.stat().st_size<=2_000_000:inputs.append(f)
            else:exclusions.append(dict(path=f.relative_to(root).as_posix(),reason='Selection JSON exceeds compact public-index limit; selected NPZ never included by this gather'))
        analysis,excluded=gather(root,['analysis'],extensions={'.json','.csv','.jsonl','.tsv','.npz'});exclusions+=excluded
        for name,files in [('max-struct-r13-postseal-inputs-sources.zip',inputs),('max-struct-r13-postseal-analysis-tables.zip',analysis)]:
            path,index=create_zip(root,out,name,files,extra_json={'RELEASE_REPRODUCIBILITY_SOURCE_INDEX.json':source} if 'inputs-sources' in name else None);assets.append(dict(path=path,kind='postseal_own_study',campaign=None,terminal_receipt=None));indices[path.name]=str(index)
    cpu,wall=time.process_time(),time.perf_counter();logical=[]
    for asset in assets:
        path=asset['path'];logical.append(dict(name=path.name,bytes=path.stat().st_size,kind=asset['kind'],original_zip_retained_at=str(path.resolve()),content_index=indices.get(path.name),parts=transport(path,out,a.part_bytes)))
    dump(out/'REPRODUCIBILITY_SOURCE_INDEX.json',source);dump(out/'TRANSPORT_MANIFEST.json',dict(status='PREPARED_LOCAL_NOT_PUBLISHED',published=False,logical_assets=logical,part_size_bytes=a.part_bytes,order='parts listed in binary concatenation order',hashes='not computed; exact byte counts only',verified_required_seals=verified,excluded_entries=exclusions,source_index='REPRODUCIBILITY_SOURCE_INDEX.json',private_inputs_archived=a.include_postseal,split_CPU_seconds=time.process_time()-cpu,split_wall_seconds=time.perf_counter()-wall,new_physics=0,source_deleted=False,public_urls_created=False))
    dump(out/'PACKAGING_RECEIPT.json',dict(status='LOCAL_PACKAGING_COMPLETE_NOT_PUBLISHED',CPU_seconds=time.process_time()-process_cpu0,wall_seconds=time.perf_counter()-process_wall0,logical_assets=len(logical),transport_parts=sum(len(x['parts']) for x in logical),timer_scope='Whole prepare process through manifest write; per-archive index and split timers nested, not additive; final receipt/exit outside interval',new_physics=0,source_deleted=False))
    print(json.dumps(dict(status='PREPARED_LOCAL_NOT_PUBLISHED',logical_assets=len(logical),transport_parts=sum(len(x['parts']) for x in logical),private_inputs_archived=a.include_postseal,new_physics=0)))
if __name__=='__main__':main()

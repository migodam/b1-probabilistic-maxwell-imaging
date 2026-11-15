"""Plan/store full evidence in bounded release ZIPs; no array/scoring operations.

Defaults to a small manifest-only plan. Execution is explicit and never deletes
originals. Text uses lossless DEFLATED level1; already-compressed/binary files stay STORED.
"""
from __future__ import annotations
import argparse,csv,json,os,re,shutil,sys,time,zipfile
from datetime import datetime,timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TARGET_BYTES=600_000_000
HARD_BYTES=650_000_000
EXCLUDE_PARTS={'private_operational','sealed_truth','fulltext','__pycache__','.git','.venv','node_modules'}
EXCLUDE_NAMES={'CHAT_READBACK.json','TASK_RESUME.md','.DS_Store'}
SECRET_RE=re.compile(rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{30,}|(?<![A-Za-z0-9])sk-[A-Za-z0-9_-]{30,}|Bearer\s+[A-Za-z0-9._~-]{24,}|(?:api[_-]?key|access[_-]?token|refresh[_-]?token|client[_-]?secret|password)\s*["\x27]?\s*[:=]\s*["\x27][A-Za-z0-9_./+=-]{16,})',re.I)
TEXT_SUFFIX={'.py','.json','.jsonl','.csv','.txt','.md','.html','.log','.yml','.yaml','.toml','.ini','.cfg','.ps1','.sh'}
PUBLIC_EXT={'.npz','.npy','.json','.csv','.md','.txt'}


def now():return datetime.now(timezone.utc).isoformat()

def safe_name(value):
    if not value or Path(value).name!=value or value in ('.','..') or any(c in value for c in '/\\'):
        raise ValueError('Expected one safe campaign/package name')
    return value


def excluded(path,root):
    rel=path.relative_to(root)
    return any(x.lower() in EXCLUDE_PARTS for x in rel.parts) or path.name in EXCLUDE_NAMES or any(x in path.name.lower() for x in ('credential','secret','token_cache','cached_chat')) or path.is_symlink()


def files_under(directory,root):
    if not directory.is_dir() or directory.is_symlink():return
    for base,dirs,names in os.walk(directory,followlinks=False):
        here=Path(base);dirs[:]=[d for d in dirs if not excluded(here/d,root)]
        for name in sorted(names):
            path=here/name
            if path.is_file() and not excluded(path,root) and path.resolve().is_relative_to(root):yield path


def scan_text(path):
    # These are generated numerical evidence trees, not credential-bearing
    # source/configuration. Private directories are excluded by inventory.
    if any(part in ('RESULTS','queues','analysis') for part in path.parts):return
    if path.suffix.lower() in ('.json','.jsonl','.csv'):return
    if path.suffix.lower() in ('.json','.jsonl','.csv') and 'audit' in path.parts and any(part in ('shared_readout4_v1','vector4_v1','cost_ledger_partial_v3') for part in path.parts):return
    if path.suffix.lower() not in TEXT_SUFFIX:return
    with path.open('rb') as stream:
        tail=b''
        while block:=stream.read(1024*1024):
            combined=tail+block
            if SECRET_RE.search(combined):raise ValueError('Possiblecredential/token literal; omitnothing silently, review file: '+str(path))
            tail=combined[-1024:]


def inventory(paths,root):
    rows=[]
    for path in sorted(set(paths)):
        if excluded(path,root):raise ValueError('Excluded file cannot be admitted: '+str(path))
        scan_text(path);stat=path.stat()
        rows.append(dict(path=path.relative_to(root).as_posix(),bytes=stat.st_size,mtime_ns=stat.st_mtime_ns,
            archive_compression='ZIP_DEFLATED' if path.suffix.lower() in TEXT_SUFFIX else 'ZIP_STORED',compression_level=1 if path.suffix.lower() in TEXT_SUFFIX else None,numerical_transform='NONE'))
    return rows


def volumes(rows,prefix):
    parts=[];current=[];total=0
    for row in rows:
        if row['bytes']>TARGET_BYTES and row['archive_compression']=='ZIP_STORED':
            raise ValueError('Single STORED file exceeds600MB; no fragmentation implemented: '+row['path'])
        # ZIP headers counted conservatively in addition to the 50MB margin.
        charge=row['bytes']+512+2*len(row['path'].encode('utf8'))
        if current and total+charge>TARGET_BYTES:
            parts.append(dict(archive=f'{prefix}_part{len(parts)+1:02d}.zip',members=current,raw_bytes_plus_header_allowance=total+4096,actual_compressed_bytes='UNKNOWN_UNTIL_EXECUTE'));current=[];total=0
        current.append(row);total+=charge
    if current:parts.append(dict(archive=f'{prefix}_part{len(parts)+1:02d}.zip',members=current,raw_bytes_plus_header_allowance=total+4096,actual_compressed_bytes='UNKNOWN_UNTIL_EXECUTE'))
    return parts


def campaign_package(root,campaign):
    campaign=safe_name(campaign);queue=root/'queues'/campaign;result=root/'RESULTS'/campaign
    terminal=json.loads((queue/'RECEIPT.json').read_text())
    if terminal.get('status') not in ('COMPLETED','COMPLETED_WITH_ERRORS') or terminal.get('running') or terminal.get('pending_tasks',0) or terminal.get('pending_cases',0):
        raise RuntimeError('Require a terminal campaign; running data cannot be frozen: '+campaign)
    frozen=json.loads((queue/'FROZEN_MANIFEST.json').read_text());tasks=frozen['tasks'];missing=[]
    for task in tasks:
        folder=result/task['task_id']
        if not folder.is_dir():missing.append(task['task_id'])
    if missing:raise FileNotFoundError('Taskoutputs missing, including failed tasks: '+','.join(missing))
    paths=[*files_under(result,root),*files_under(queue,root)]
    config=root/'configs'/(campaign+'.json')
    if config.is_file():paths.append(config)
    rows=inventory(paths,root)
    return dict(name='raw_'+campaign,kind='FULL_RAW_CAMPAIGN',campaign=campaign,
        queue_status=terminal['status'],tasks=len(tasks),files=rows,total_original_bytes=sum(x['bytes'] for x in rows),
        max_file_bytes=max((x['bytes'] for x in rows),default=0),
        completeness='All admitted RESULTS and queues files, every endpoint/full-field NPZ, checkpoint, candidatepool, JSON/JSONL and log; no material-only reduction',
        parts=volumes(rows,'r14_raw_'+campaign))


def local_group_package(root,campaign):
    """Only these two known terminal local groups; preserve every attempt."""
    if campaign not in ('shared_readout4_v1','vector4_v1'):
        raise ValueError('Only the two explicitly authorized local groups are supported')
    attempt='attempt_02' if campaign=='shared_readout4_v1' else 'attempt_01'
    queue=root/'queues'/campaign;success_path=queue/attempt/'RECEIPT.json'
    success=json.loads(success_path.read_text())
    expected=16 if campaign=='shared_readout4_v1' else 24
    allowed=('ALL_EXPECTED_TASKS_TERMINAL','ALL_EXPECTED_CHILDREN_TERMINAL')
    if success.get('status') not in allowed or success.get('expected')!=expected or success.get('terminal')!=expected or success.get('successful')!=expected or success.get('failed_or_missing'):
        raise RuntimeError('Explicit successful attempt is not complete: '+str(success_path))
    jobs=success.get('jobs',[])
    if len(jobs)!=expected or any(j.get('returncode')!=0 or j.get('child_status')!='COMPLETE' for j in jobs):
        raise RuntimeError('Successful attempt must explicitly account for every complete child')
    config=root/'configs'/(campaign+'.json');cfg=json.loads(config.read_text())
    if len(cfg['tasks'])!=expected or {t['task_id'] for t in cfg['tasks']}!={j['task_id'] for j in jobs}:
        raise ValueError('Successful attempt differs from actual frozen config task IDs')
    result=root/'RESULTS'/campaign
    for task in cfg['tasks']:
        folder=result/task['task_id']
        receipt=folder/'RECEIPT.json' if campaign=='shared_readout4_v1' else folder/'attempt_01/RECEIPT.json'
        if not receipt.is_file() or json.loads(receipt.read_text()).get('status')!='COMPLETE':
            raise FileNotFoundError('Missing successful original result receipt: '+str(receipt))
    paths=[*files_under(result,root),*files_under(queue,root),*files_under(root/'predictions'/campaign,root),config]
    audit=root/'audit'/campaign
    if campaign=='shared_readout4_v1':
        paths.extend(files_under(audit,root))
    else:
        # No GT-derived support arrays/tables or scored figures are opened.
        paths.extend(p for p in files_under(audit,root)
            if p.parent==audit or 'RECEIPT' in p.name)
    rows=inventory(paths,root)
    old_root=queue/'RECEIPT.json'
    return dict(name='local_fullstates_'+campaign,kind='LOCAL_KNOWN_TERMINAL_FULLSTATES',campaign=campaign,
        successful_attempt=attempt,successful_receipt=success_path.relative_to(root).as_posix(),
        successful_status=success['status'],tasks=expected,all_queue_attempts_included=True,
        old_root_receipt_included=old_root.is_file(),
        root_receipt_not_substituted_for_success=campaign=='shared_readout4_v1',
        files=rows,total_original_bytes=sum(r['bytes'] for r in rows),max_file_bytes=max(r['bytes'] for r in rows),
        completeness='All admitted original RESULTS andqueues/predictions/config bytes,all candidate/fullfield/checkpoint arrays;all attempts retained,successful attempt explicitly identified',
        audit_scope='Shared original source/runtime/tag/cost receipts;vector top-level audit plus nested receipt files only. GT-derived support arrays/tables and scoredfigures outside rawstate scope excluded without opening.',
        parts=volumes(rows,'r14_local_fullstates_'+campaign))


def local_summary_package(root):
    directory=root/'audit/cost_ledger_partial_v3'
    names=('REPORT.md','FULL_STATE_CATALOG.md','COST_TOTALS.json','COST_TOTALS.csv','COVERAGE.json','COVERAGE.csv','RECEIPT.json','READBACK_VERIFICATION.json')
    paths=[directory/n for n in names]
    if any(not p.is_file() for p in paths):raise FileNotFoundError('Accepted partial_v3 summary source missing')
    rows=inventory(paths,root)
    return dict(name='cost_ledger_partial_v3_summary',kind='EXISTING_COST_SUMMARY_PROVENANCE',files=rows,
        total_original_bytes=sum(r['bytes'] for r in rows),max_file_bytes=max(r['bytes'] for r in rows),
        scope='Current accepted summary provenance only;does not rearchive main192 data anddoesnot include this later archive cost',
        parts=volumes(rows,'r14_local_fullstates_cost_summary'))


def write_local_plan(root,destination):
    cpu,wall=time.process_time(),time.perf_counter();destination.mkdir(parents=True,exist_ok=False)
    packages=[local_group_package(root,c) for c in ('shared_readout4_v1','vector4_v1')]+[local_summary_package(root)]
    plan=dict(schema='R14_FULL_EVIDENCE_PLAN_V1',status='PLANNED_NO_ARCHIVES_WRITTEN',created_utc=now(),source_root=str(root),
        shard_target_bytes=TARGET_BYTES,actual_archive_hard_limit_bytes=HARD_BYTES,packages=packages,
        files=sum(len(p['files']) for p in packages),original_bytes=sum(p['total_original_bytes'] for p in packages),
        exclusions=dict(sealed_truth_access='NONE',main192_rearchived=False,
            known_generated_numeric_scan='RESULTS/queues andknown numericJSON/CSV auditgroups are not credential-scanned;small source/config text still uses existing scanner',
            vector_GT_diagnostics='Nested support arrays,tables andscoredfigures not requested/not traversed for content'),
        transformations='NONE;lossless textDEFLATE1,binarySTORED,original field/candidate/state bytes retained',
        hash_operations='NONE',archives_written=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
    for package in packages:write_new(destination/(package['name']+'_MANIFEST.json'),package)
    write_new(destination/'PACKAGE_PLAN.json',plan)
    with (destination/'FILE_INVENTORY.csv').open('x',newline='') as f:
        writer=csv.writer(f);writer.writerow(['package','path','bytes','mtime_ns','archive_compression'])
        for package in packages:writer.writerows([package['name'],r['path'],r['bytes'],r['mtime_ns'],r['archive_compression']] for r in package['files'])
    return plan


def source_package(root):
    paths=[]
    for folder in ('src','runtime','configs','data_public','sources'):
        for path in files_under(root/folder,root):
            if folder=='data_public' and path.suffix.lower() not in PUBLIC_EXT:continue
            if folder in ('src','runtime') and path.suffix.lower() not in ('.py','.md','.txt','.json','.ps1','.sh','.toml','.yml','.yaml'):continue
            paths.append(path)
    for path in root.iterdir():
        if path.is_file() and path.suffix.lower() in ('.md','.json','.html','.css') and not excluded(path,root):paths.append(path)
    # The self-contained source bundle deliberately omits cached conversations,
    # downloaded literature fulltexts and every private/sealed directory.
    rows=inventory(paths,root)
    return dict(name='sources_configs_public',kind='SOURCE_CONFIG_PUBLIC_NPZ',files=rows,
        total_original_bytes=sum(x['bytes'] for x in rows),max_file_bytes=max((x['bytes'] for x in rows),default=0),
        parts=volumes(rows,'r14_sources_configs_public'))


def diagnostics_package(root):
    paths=[]
    for base in ('audit','analysis','FIGURES','predictions','release'):
        for path in files_under(root/base,root):
            if path.is_relative_to(root/'audit/packaging_plan'):continue
            paths.append(path)
    rows=inventory(paths,root)
    return dict(name='analysis_diagnostics',kind='FULL_ANALYSIS_AUDIT_INCLUDING_NUMERIC_STATES',files=rows,
        total_original_bytes=sum(x['bytes'] for x in rows),max_file_bytes=max((x['bytes'] for x in rows),default=0),
        parts=volumes(rows,'r14_analysis_diagnostics'))


def write_new(path,data):
    with path.open('x',encoding='utf8') as stream:json.dump(data,stream,indent=2,ensure_ascii=False);stream.write('\n')


def write_plan(root,destination,campaigns,include_sources,include_diagnostics):
    cpu=time.process_time();wall=time.perf_counter();packages=[]
    # Only tiny JSON manifests/restore instructions are written in plan mode.
    destination.mkdir(parents=True,exist_ok=False)
    try:
        if include_sources:packages.append(source_package(root))
        for campaign in campaigns:packages.append(campaign_package(root,campaign))
        if include_diagnostics:packages.append(diagnostics_package(root))
        plan=dict(schema='R14_FULL_EVIDENCE_PLAN_V1',status='PLANNED_NO_ARCHIVES_WRITTEN',created_utc=now(),source_root=str(root),
            shard_target_bytes=TARGET_BYTES,actual_archive_hard_limit_bytes=HARD_BYTES,
            packages=packages,files=sum(len(x['files']) for x in packages),
            original_bytes=sum(x['total_original_bytes'] for x in packages),
            exclusions=dict(directory_components=sorted(EXCLUDE_PARTS),filenames=sorted(EXCLUDE_NAMES),
                credentials='Textstream conservative literal scan; positives block rather than remove evidence silently',
                third_party_fulltexts='fulltext directories never traversed',sealed_truth_access='NONE'),
            transformations='NONE; original bytes restored; text lossless DEFLATED level1, NPZ/NPY/PNG/PDF STORED; no NPZ load, cast, recompute or material-only substitution',
            hash_operations='NONE',archives_written=0,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
        for package in packages:write_new(destination/(package['name']+'_MANIFEST.json'),package)
        write_new(destination/'PACKAGE_PLAN.json',plan)
        with (destination/'FILE_INVENTORY.csv').open('x',newline='') as stream:
            writer=csv.writer(stream);writer.writerow(['package','path','bytes','mtime_ns','archive_compression'])
            for package in packages:
                writer.writerows([package['name'],r['path'],r['bytes'],r['mtime_ns'],r['archive_compression']] for r in package['files'])
        (destination/'RESTORE.md').write_text('''# Full raw evidence restore\n\nThese manifests describe ORIGINAL full-field files, all candidates/checkpoints and queue receipts, not export_scoring material-only derivatives. Plans contain no release URLs until actual upload/readback.\n\nText JSON/CSV/MD/source/log uses lossless ZIP_DEFLATED level1; NPZ/NPY/PNG/PDF and other binary assets use ZIP_STORED. All original file bytes restore without numerical transformation or lossy reduction. Manifest records each file compression and restore checks it. Volumes target600MB by raw input size; a single larger text file occupies its own compressed volume, whose actual ZIP must still stay≤650MB or execution fails. No SHA operations.\n\nAfter publication download every listed part plus PACKAGE_PLAN.json to one external-disk folder. Run:\n\n```text\npython package_evidence.py restore --plan /external/PACKAGE_PLAN.json --archive-dir /external/downloads --destination /external/restored_r14\n```\n\nRestore refuses an existing destination, checks exact member path sets and byte sizes, and opens each ZIP member as a stream. Python zipfile performs ordinary ZIP CRC checks; no scientific rerun, SHA, optimizer, scorer or pickle is used. Repeated config files across source/raw packages must compare byte-identical without hashes. Paths reconstruct RESULTS/<campaign>/..., queues/<campaign>/..., configs/..., data_public/... and original source layout. Missing parts are fatal; originals are never deleted.\n\nTo produce actual archives later:\n\n```text\npython package_evidence.py execute --root /external/r14 --plan /external/PACKAGE_PLAN.json --output /external/release_shards\n```\n\nBefore execution terminal campaign status and source sizes/mtime must still match the plan; stale plans must be regenerated. Successful creation reports archive actual sizes and per-file membership. Upload, anonymous access verification, data closure and scientific acceptance are parent tasks. Never substitute *_scoring_only.zip for these full raw archives.\n''')
        return plan
    except Exception as exc:
        write_new(destination/'PLAN_ERROR.json',dict(status='FAILED_PLAN_NO_ARCHIVES',error=repr(exc),created_utc=now()))
        raise


def external_output(path):
    if sys.platform=='darwin' and not path.resolve().is_relative_to('/Volumes'):
        raise ValueError('Large archive/restore writes must use external /Volumes disk')


def execute(root,plan_path,output):
    external_output(output);plan=json.loads(plan_path.read_text());output.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter();cpu=time.process_time();done=[]
    try:
        for package in plan['packages']:
            if package['kind']=='LOCAL_KNOWN_TERMINAL_FULLSTATES':
                current=local_group_package(root,package['campaign'])
                if current['files']!=package['files'] or current['successful_attempt']!=package['successful_attempt']:
                    raise ValueError('Local terminal group changed;regenerate explicit plan')
            if package['kind']=='FULL_RAW_CAMPAIGN':
                campaign=package['campaign'];queue=root/'queues'/campaign
                terminal=json.loads((queue/'RECEIPT.json').read_text())
                if terminal.get('status') not in ('COMPLETED','COMPLETED_WITH_ERRORS') or terminal.get('running') or terminal.get('pending_tasks',0):
                    raise RuntimeError('Campaign is no longer terminal: '+campaign)
                current={p.relative_to(root).as_posix() for d in (root/'RESULTS'/campaign,queue) for p in files_under(d,root)}
                config=root/'configs'/(campaign+'.json')
                if config.is_file():current.add(config.relative_to(root).as_posix())
                if current!={x['path'] for x in package['files']}:raise ValueError('Campaign file inventory changed; regenerate plan')
            for part in package['parts']:
                archive=output/part['archive']
                with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_STORED,allowZip64=True) as z:
                    for entry in part['members']:
                        path=root/entry['path']
                        if excluded(path,root) or not path.resolve().is_relative_to(root):raise ValueError('Excluded/unsafe original path')
                        st=path.stat()
                        if st.st_size!=entry['bytes'] or st.st_mtime_ns!=entry['mtime_ns']:raise ValueError('Stale source plan: '+entry['path'])
                        scan_text(path);z.write(path,entry['path'],compress_type=getattr(zipfile,entry['archive_compression']),compresslevel=entry.get('compression_level'))
                        after=path.stat()
                        if (after.st_size,after.st_mtime_ns)!=(entry['bytes'],entry['mtime_ns']):raise ValueError('Original changed while archiving: '+entry['path'])
                actual=archive.stat().st_size
                if actual>HARD_BYTES:raise ValueError('Actual ZIP exceeds650MB: '+archive.name)
                with zipfile.ZipFile(archive) as z:
                    infos=z.infolist()
                    if [x.filename for x in infos]!=[x['path'] for x in part['members']]:raise ValueError('Archive member mismatch')
                    if any(x.file_size!=e['bytes'] or x.compress_type!=getattr(zipfile,e['archive_compression']) for x,e in zip(infos,part['members'])):raise ValueError('Archive member size/compression mismatch')
                done.append(dict(archive=archive.name,bytes=actual,members=len(part['members']),status='WRITTEN_SIZE_AND_DIRECTORY_CHECKED'))
        shutil.copy2(plan_path,output/'PACKAGE_PLAN.json')
        receipt=dict(status='ARCHIVED_NO_DELETION_NOT_UPLOADED',parts=done,originals_retained=True,hash_operations='NONE',
            CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-started)
        write_new(output/'PACKAGE_RECEIPT.json',receipt);return receipt
    except Exception as exc:
        write_new(output/'PACKAGE_ERROR.json',dict(status='FAILED_RETAINED_ALL_ATTEMPTS',error=repr(exc),completed_parts=done))
        raise


def restore(plan_path,archive_dir,destination):
    external_output(destination);plan=json.loads(plan_path.read_text());destination.mkdir(parents=True,exist_ok=False)
    seen=set();count=0
    for package in plan['packages']:
        for part in package['parts']:
            with zipfile.ZipFile(archive_dir/part['archive']) as z:
                expected={x['path']:x for x in part['members']}
                if set(z.namelist())!=set(expected) or len(z.namelist())!=len(expected):raise ValueError('Missing/extra ZIP members')
                for info in z.infolist():
                    target=destination/info.filename
                    if Path(info.filename).is_absolute() or '\\' in info.filename or not target.resolve().is_relative_to(destination.resolve()):raise ValueError('Unsafe ZIP memberpath')
                    if info.file_size!=expected[info.filename]['bytes'] or info.compress_type!=getattr(zipfile,expected[info.filename]['archive_compression']):raise ValueError('Member bytes/compression differ from manifest')
                    if info.filename in seen:
                        # Campaign manifests deliberately carry their config
                        # again so raw releases also restore independently.
                        with z.open(info) as incoming,target.open('rb') as existing:
                            while chunk:=incoming.read(1024*1024):
                                if existing.read(len(chunk))!=chunk:raise ValueError('Conflicting repeated original file: '+info.filename)
                            if existing.read(1):raise ValueError('Repeated original length differs')
                        continue
                    seen.add(info.filename)
                    target.parent.mkdir(parents=True,exist_ok=True)
                    with z.open(info) as source,target.open('xb') as out:shutil.copyfileobj(source,out,length=1024*1024)
                    count+=1
    return dict(status='RESTORED_ORIGINAL_FILE_BYTES',files=count,hash_operations='NONE',new_physics=0)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    q=sub.add_parser('plan');q.add_argument('--root',type=Path,default=ROOT);q.add_argument('--campaign',action='append',default=[])
    q.add_argument('--output',type=Path);q.add_argument('--without-sources',action='store_true');q.add_argument('--include-diagnostics',action='store_true')
    q=sub.add_parser('plan-local-fullstates');q.add_argument('--root',type=Path,default=ROOT);q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('execute');q.add_argument('--root',type=Path,default=ROOT);q.add_argument('--plan',type=Path,required=True);q.add_argument('--output',type=Path,required=True)
    q=sub.add_parser('restore');q.add_argument('--plan',type=Path,required=True);q.add_argument('--archive-dir',type=Path,required=True);q.add_argument('--destination',type=Path,required=True)
    a=p.parse_args()
    if a.command=='plan':
        root=a.root.resolve();out=a.output or root/'audit/packaging_plan'/datetime.now().strftime('%Y%m%d_%H%M%S')
        result=write_plan(root,out.resolve(),a.campaign,not a.without_sources,a.include_diagnostics)
        print(json.dumps({k:result[k] for k in ('status','files','original_bytes','archives_written')}))
    elif a.command=='plan-local-fullstates':
        result=write_local_plan(a.root.resolve(),a.output.resolve());print(json.dumps({k:result[k] for k in ('status','files','original_bytes','archives_written')}))
    elif a.command=='execute':print(json.dumps(execute(a.root.resolve(),a.plan.resolve(),a.output.resolve())))
    else:print(json.dumps(restore(a.plan.resolve(),a.archive_dir.resolve(),a.destination.resolve())))


if __name__=='__main__':main()

"""Finite metadata planning/ZIP packaging for the existing B1 GitHub Release.

No NumPy/torch import, array decoding, scoring, network, hashes, or installation.
Planning is the default. Only an explicit pack action writes ZIP payloads.
"""
from pathlib import Path,PurePosixPath
from collections import Counter,deque
from datetime import datetime,timezone
import argparse,ast,json,os,sys,time,zipfile

STUDY=Path(__file__).resolve().parents[1]
REPO=STUDY.parents[1]
MAX_BYTES=1_500_000_000  # Decimal1.5GB; safely below the per-asset request.
MAX_MEMBERS=60000       # Avoid ZIP64 central-directory count overhead.
BLOCK=8*1024*1024
PREFIX='gn-nn-r4'
OLD_PY_ROOTS=('research/research_gn_nn_completion_r3','research/research_gn_search_space_r2','research/research_z_alg_r4_diversity')
FIXED_SPECTRA='research/research_gn_green_error_20261106/results/fixed_spectra.npz'
OLD_E0_METADATA=('research/research_gn_nn_completion_r3/audit/CONFIRMATION_PLAN_FROZEN.json','research/research_gn_nn_completion_r3/evaluation/CONFIRMATION_SEALED.json','research/research_gn_nn_completion_r3/evaluation/CONFIRMATION_SEALED_MANIFEST.json')
EXCLUDE_PARTS={'.git','__pycache__','.venv','node_modules'}


def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def utc():return datetime.now(timezone.utc).isoformat()
def within(path,root):
    try:path.resolve().relative_to(root.resolve());return True
    except ValueError:return False

def walk_files(root,suffixes=None):
    """Walk one explicit study subtree, never following directory/file symlinks."""
    root=Path(root)
    if not root.exists():return []
    if root.is_symlink():return []
    if root.is_file():return [root] if suffixes is None or root.suffix.lower()in suffixes else []
    result=[]
    for current,dirs,files in os.walk(root,followlinks=False):
        dirs[:]=sorted(d for d in dirs if d not in EXCLUDE_PARTS and not d.startswith('.') and not (Path(current)/d).is_symlink())
        for name in sorted(files):
            p=Path(current)/name
            if name.startswith('.') or p.is_symlink() or not p.is_file() or name.endswith(('.tmp','.partial')):continue
            if suffixes is None or p.suffix.lower()in suffixes:result.append(p)
    return result

def item(path,scope):
    path=Path(path)
    if path.is_symlink() or not within(path,REPO):raise ValueError('Payload must be a real file within the declared repository layout')
    st=path.stat();return dict(path=path.relative_to(REPO).as_posix(),size_bytes=st.st_size,mtime_ns=st.st_mtime_ns,provenance_scope=scope)

def old_python_dependencies():
    """AST-only import closure, restricted to three explicitly approved old studies."""
    source=[]
    for folder in ('src','experiments','vendor'):source+=walk_files(STUDY/folder,{'.py'})
    queue=deque(source);seen=set();found={};external=set();unresolved=[]
    while queue:
        path=queue.popleft()
        if path in seen:continue
        seen.add(path)
        try:tree=ast.parse(path.read_text(encoding='utf-8'),filename=str(path))
        except (OSError,SyntaxError)as exc:unresolved.append(dict(path=str(path.relative_to(REPO)),reason=str(exc)));continue
        for node in ast.walk(tree):
            modules=[]
            if isinstance(node,ast.Import):modules=[alias.name for alias in node.names]
            elif isinstance(node,ast.ImportFrom):
                if node.level:
                    package=path.parent
                    for _ in range(node.level-1):package=package.parent
                    prefix='.'.join(package.relative_to(REPO).parts)
                    modules=[prefix+'.'+node.module]if node.module else[prefix+'.'+a.name for a in node.names if a.name!='*']
                elif node.module:modules=[node.module]
            for module in modules:
                if not module.startswith('research.'):
                    external.add(module.split('.')[0]);continue
                relative=module.replace('.','/')
                if relative.startswith(STUDY.relative_to(REPO).as_posix()+'/'):continue
                if not any(relative==r or relative.startswith(r+'/')for r in OLD_PY_ROOTS):
                    unresolved.append(dict(imported_by=path.relative_to(REPO).as_posix(),module=module,status='OUTSIDE_APPROVED_OLD_SOURCE_SCOPE'));continue
                candidate=REPO/(relative+'.py')
                if not candidate.is_file():candidate=REPO/relative/'__init__.py'
                if not candidate.is_file():unresolved.append(dict(module=module,status='MISSING_IMPORTED_PYTHON'));continue
                if candidate.is_symlink():raise ValueError('Dependency symlink rejected')
                if candidate not in found:
                    found[candidate]=item(candidate,'frozen imported old-study Python');queue.append(candidate)
                parent=candidate.parent
                while parent!=REPO and within(parent,REPO):
                    init=parent/'__init__.py'
                    if init.is_file() and not init.is_symlink() and any(within(init,REPO/r)for r in OLD_PY_ROOTS)and init not in found:
                        found[init]=item(init,'old-study package initializer');queue.append(init)
                    parent=parent.parent
    # R2 Search loads this physics module through importlib; the explicit runtime
    # dependency is retained even if a future AST closure no longer sees a from-import.
    physical=REPO/'research/research_z_alg_r4_diversity/src/physics.py'
    if physical.is_file():found.setdefault(physical,item(physical,'explicit R2 importlib physics dependency'))
    return sorted(found.values(),key=lambda v:v['path'].encode()),sorted(external),unresolved

def gates():
    method=STUDY/'audit/METHOD_FREEZE.json';sealed=STUDY/'evaluation/TEST_SEALED.json';issues=[];seal=None;freeze=None
    for p in (method,sealed):
        if not p.is_file():issues.append(dict(path=p.relative_to(STUDY).as_posix(),status='NOT_RUN_MISSING'))
    if method.is_file():freeze=read_json(method)
    if sealed.is_file():
        seal=read_json(sealed)
        if seal.get('status')!='ALL_PLANNED_PREDICTIONS_SEALED' or set(seal.get('splits',[]))!={'test'}:issues.append(dict(path='evaluation/TEST_SEALED.json',status='HOLD_INVALID_TEST_SEAL'))
    return dict(pack_ready=not issues,issues=issues,method_freeze_present=method.is_file(),test_sealed_present=sealed.is_file()),freeze,seal

def metadata_coverage(seal):
    public=STUDY/'data/manifests/online.json';final=STUDY/'audit/FINAL_RUNS.json'
    result=dict(seal_disposition_counts=dict(Counter(r.get('disposition','UNKNOWN')for r in (seal or{}).get('receipts',[]))),sealed_NOT_RUN=[dict(run=r['run'],obs_id=r['obs_id'],status='NOT_RUN')for r in (seal or{}).get('receipts',[])if r.get('disposition')=='NOT_RUN'])
    if public.is_file():
        rows=read_json(public)['records'];valtest=[r for r in rows if r['split']in ('val','test')]
        specs=[('strict',rows,''),('discrepancy_final_only',valtest,'_discrepancy'),('weakreg1e-5',valtest,'_weakreg1e-5'),('discrepancy_all_scales',valtest,'_discrepancy_allscales')]
        result['anchor_coverage']={}
        for name,records,suffix in specs:
            missing=[r['obs_id']for r in records if not all((STUDY/'anchors'/(r['obs_id']+suffix+ext)).is_file()for ext in ('.npz','.json'))]
            result['anchor_coverage'][name]=dict(planned=len(records),complete_pairs=len(records)-len(missing),missing_NOT_RUN=missing)
        result['public_observations']=dict(planned=len(rows),missing_NOT_RUN=[r['obs_id']for r in rows if not (STUDY/'data'/r['public_path']).is_file()])
    if final.is_file():
        frozen=read_json(final);result['declared_runs']=frozen['runs'];result['declared_checkpoint_files']=[]
        for plan in frozen['plans']:
            if plan['kind']!='NN':continue
            checkpoint=plan.get('checkpoint')
            if checkpoint:
                p=Path(checkpoint);p=p if p.is_absolute()else STUDY/p
                result['declared_checkpoint_files'].append(dict(run=plan['run'],path=str(checkpoint),status='PRESENT'if within(p,STUDY)and p.is_file()else'NOT_RUN_MISSING_OR_OUTSIDE_STUDY'))
            else:result['declared_checkpoint_files'].append(dict(run=plan['run'],status='NOT_DECLARED_IN_FINAL_RUNS'))
    else:result['final_runs_status']='NOT_RUN_MISSING_FINAL_RUNS'
    return result

def member_cost(member):return member['length_bytes']+76+2*len(member['member'].encode('utf-8'))
def part_members(files,limit):
    members=[]
    for f in files:
        direct=dict(source=f['path'],member=f['path'],offset=0,length_bytes=f['size_bytes'])
        if member_cost(direct)+22<=limit:members.append(direct);continue
        offset=0;chunk=1
        while offset<f['size_bytes']:
            name='__release_chunks__/'+f['path']+f'.chunk{chunk:06d}'
            length=min(f['size_bytes']-offset,limit-22-76-2*len(name.encode('utf-8')))
            if length<=0:raise ValueError('Asset byte limit too small for a ZIP member name')
            members.append(dict(source=f['path'],member=name,offset=offset,length_bytes=length));offset+=length;chunk+=1
    parts=[];current=[];size=22
    for m in members:
        cost=member_cost(m)
        if current and (size+cost>limit or len(current)>=MAX_MEMBERS):parts.append(current);current=[];size=22
        current.append(m);size+=cost
    if current:parts.append(current)
    return parts

def build_plan(args):
    groups=[];used=set()
    def group(name,roots,suffixes=None,extra=(),optional=False):
        paths=[];components=[]
        for root in roots:
            root=STUDY/root;present=root.exists();components.append(dict(path=root.relative_to(STUDY).as_posix(),status='PRESENT'if present else'NOT_RUN_MISSING'))
            paths+=walk_files(root,suffixes)
        paths+=list(extra);files=[]
        for p in sorted(set(paths),key=lambda x:x.relative_to(REPO).as_posix().encode()):
            if p in used:continue
            used.add(p);files.append(item(p,'R4 study artifact'if within(p,STUDY)else'explicit frozen dependency'))
        groups.append(dict(name=name,status='PLANNED_NOT_PACKED'if files else('NOT_RUN_MISSING_OPTIONAL'if optional else'NOT_RUN_MISSING'),components=components,file_count=len(files),source_bytes=sum(f['size_bytes']for f in files),files=files))
    group('data-public',['data/public','data/operators.npz','data/reference_jacobian.npz','data/reference_jacobian.json','data/manifests/online.json','data/GENERATION_RECEIPT.json'],{'.npz','.json'})
    group('data-evaluator',['data/truth','data/manifests/evaluation.json'],{'.npz','.json'})
    group('anchors',['anchors'],{'.npz','.json'})
    group('models',['models'])  # Includes ignored best.pt files, configs and histories.
    group('predictions',['predictions'],{'.npz','.json'})
    group('evaluation',['evaluation'],{'.npz','.json'})
    spectra=REPO/FIXED_SPECTRA
    group('fixed-projectors',['cache'],{'.npz','.json'},extra=[spectra]if spectra.is_file()else[])
    group('analysis-tables',['results/e0','results/test_energy','results'],{'.csv','.tsv','.json','.npz'})
    group('canonical-receipts',['audit','delegated'],{'.json','.jsonl','.csv','.log','.err','.out','.txt'})
    # The plan is output metadata, not recursively its own input payload.
    for g in groups:
        g['files']=[f for f in g['files']if f['path']!=(STUDY/'audit/RELEASE_PACKAGE_PLAN.json').relative_to(REPO).as_posix()]
        g['file_count']=len(g['files']);g['source_bytes']=sum(f['size_bytes']for f in g['files'])
    deps,external,unresolved=old_python_dependencies()
    old_meta=[REPO/p for p in OLD_E0_METADATA if(REPO/p).is_file()]
    dep_files=deps+[item(p,'explicit E0 historical seal/plan metadata; no old raw arrays')for p in old_meta]
    groups.append(dict(name='frozen-python-dependencies',status='PLANNED_NOT_PACKED'if dep_files else'NOT_RUN_MISSING',file_count=len(dep_files),source_bytes=sum(f['size_bytes']for f in dep_files),files=sorted(dep_files,key=lambda f:f['path'].encode()),components=[dict(path=r,scope='AST import closure only; no directory-wide old-study payload')for r in OLD_PY_ROOTS]))
    atlas=Path(args.atlas_dir);atlas=atlas if atlas.is_absolute()else STUDY/atlas
    if not within(atlas,STUDY):raise ValueError('Optional atlas must stay within this study')
    if args.include_atlas:group('atlas',[atlas.relative_to(STUDY).as_posix()],optional=True)
    else:groups.append(dict(name='atlas',status='NOT_INCLUDED_OPTIONAL',file_count=0,source_bytes=0,files=[],components=[dict(path=atlas.relative_to(STUDY).as_posix(),status='PRESENT'if atlas.exists()else'NOT_RUN_MISSING')]))
    assets=[]
    for g in groups:
        for index,members in enumerate(part_members(g['files'],args.max_bytes),1):
            assets.append(dict(name=f'{PREFIX}-{g["name"]}.part{index:03d}.zip',function=g['name'],status='PLANNED_NOT_CREATED',planned_size_bytes=22+sum(member_cost(m)for m in members),member_count=len(members),logical_file_count=len({m['source']for m in members}),members=members))
    gate,freeze,seal=gates();protocol=read_json(STUDY/'audit/PROTOCOL.json')if(STUDY/'audit/PROTOCOL.json').is_file()else{}
    regular=[]
    for folder in ('src','experiments','vendor','theory','figures','reports'):
        regular += [p.relative_to(REPO).as_posix()for p in walk_files(STUDY/folder)if not(args.include_atlas and within(p,atlas))]
    regular += [p.relative_to(REPO).as_posix()for p in STUDY.glob('*.md')if not p.is_symlink()]
    regular += [p.relative_to(REPO).as_posix()for folder in ('audit','delegated')for p in walk_files(STUDY/folder,{'.md'})]
    return dict(schema='gn-nn-r4-release-package-v1',status='PLAN_ONLY_NO_ASSETS_CREATED',created_utc=utc(),target='GitHub Release assets on the existing B1 repository; no new repository and no upload by this script',repository_url=args.repository_url or'NOT_DECLARED_PARENT_PROVIDES_EXISTING_B1_URL',publication_ref=args.publication_ref or'NOT_DECLARED',provenance=dict(study=STUDY.relative_to(REPO).as_posix(),inherited_source_commit_recorded=protocol.get('source_commit'),current_R4_publication_commit='parent supplies publication_ref; no invented commit/URL',inventory='file size and mtime only; no SHA256 or other digest',array_access='metadata only during planning; opaque byte streaming during explicit packing, never array decoding'),max_asset_bytes=args.max_bytes,zip_storage='ZIP_STORED; fixed1980 metadata and sorted paths; deterministic offsets/chunks; ZIP intrinsic CRC only',gates=gate,groups=groups,assets=assets,file_count=sum(g['file_count']for g in groups),source_bytes=sum(g['source_bytes']for g in groups),coverage=metadata_coverage(seal),external_import_names=external,dependency_issues=unresolved,external_historical_E0_data=dict(status='EXPLICIT_EXTERNAL_PUBLIC_STUDY_DEPENDENCY',release_url=args.historical_e0_release_url or'NOT_PROVIDED_PARENT_MUST_SUPPLY_VERIFIED_DIRECT_RELEASE_LINK',link_verification='parent independent access check; this packaging script performs no network request',reason='E0 CSVs and canonical receipts are included; rerunning original historical E0 additionally needs the explicitly sealed R3 raw dataset/evaluation from its existing public study; this script does not sweep that old tree'),regular_git=dict(policy='R4 source/reports/figures remain normal git files; optional full atlas can be an asset',paths=sorted(set(regular))),packing_executed=False,upload_executed=False)

def validate_seal_metadata():
    gate,freeze,seal=gates()
    if not gate['pack_ready']:raise RuntimeError('METHOD_FREEZE and valid evaluation/TEST_SEALED are required for explicit packing')
    old_root=Path(seal['study']);issues=[]
    for row in seal.get('receipts',[]):
        if row.get('disposition')!='SEALED':continue
        for key in ('json','npz'):
            stamp=row.get(key)
            if not stamp:issues.append(dict(run=row['run'],obs_id=row['obs_id'],file=key,reason='MISSING_SEALED_STAMP'));continue
            original=Path(stamp['path'])
            try:relative=original.relative_to(old_root)
            except ValueError:issues.append(dict(reason='SEALED_PATH_OUTSIDE_STUDY',path=str(original)));continue
            current=STUDY/relative
            if not current.is_file():issues.append(dict(path=relative.as_posix(),reason='MISSING_SEALED_FILE'));continue
            st=current.stat()
            if (st.st_size,st.st_mtime_ns)!=(stamp['size'],stamp['mtime_ns']):issues.append(dict(path=relative.as_posix(),reason='SEALED_SIZE_OR_MTIME_CHANGED'))
    if issues:raise RuntimeError('Sealed metadata mismatch; preserve existing evidence: '+json.dumps(issues[:10]))
    return gate

def copy_exact(source,destination,length):
    left=length
    while left:
        block=source.read(min(BLOCK,left))
        if not block:raise IOError('Source ended before its planned opaque-byte length')
        destination.write(block);left-=len(block)

def pack(args):
    started=time.perf_counter();validate_seal_metadata();plan=build_plan(args)
    if plan['dependency_issues']:raise RuntimeError('Resolve declared Python dependency issues before packing')
    output=Path(args.output).resolve()
    if output.exists() and any(output.iterdir()):raise FileExistsError('Choose a fresh empty output directory; existing assets/partial attempts are preserved')
    output.mkdir(parents=True,exist_ok=True);index_path=output/f'{PREFIX}-release-index.json';receipt_path=output/'PACKAGING_RECEIPT.json'
    by_path={f['path']:f for g in plan['groups']for f in g['files']};done=[];error=None
    write_json(STUDY/'audit/RELEASE_PACKAGE_PLAN.json',plan)
    try:
        for asset in plan['assets']:
            path=output/asset['name'];partial=path.with_suffix('.zip.partial');tick=time.perf_counter()
            with zipfile.ZipFile(partial,'w',compression=zipfile.ZIP_STORED,allowZip64=False)as archive:
                for member in asset['members']:
                    source=REPO/member['source'];expected=by_path[member['source']];before=source.stat()
                    if (before.st_size,before.st_mtime_ns)!=(expected['size_bytes'],expected['mtime_ns']):raise RuntimeError('Source changed after inventory: '+member['source'])
                    info=zipfile.ZipInfo(member['member'],date_time=(1980,1,1,0,0,0));info.compress_type=zipfile.ZIP_STORED;info.create_system=3;info.external_attr=0o100644<<16;info.file_size=member['length_bytes']
                    with source.open('rb')as reader,archive.open(info,'w',force_zip64=False)as writer:reader.seek(member['offset']);copy_exact(reader,writer,member['length_bytes'])
                    after=source.stat()
                    if (after.st_size,after.st_mtime_ns)!=(before.st_size,before.st_mtime_ns):raise RuntimeError('Source mutated during packing: '+member['source'])
            actual=partial.stat().st_size
            if actual>args.max_bytes or actual!=asset['planned_size_bytes']:raise RuntimeError('ZIP size differs from deterministic header/payload plan; partial retained')
            partial.replace(path);asset.update(status='CREATED',size_bytes=actual,pack_wall_seconds=time.perf_counter()-tick);done.append(asset['name'])
            write_json(receipt_path,dict(status='PACKING',created_assets=done,wall_seconds=time.perf_counter()-started,upload_executed=False))
        plan.update(status='PACKAGED_NOT_UPLOADED',packing_executed=True,packed_utc=utc(),pack_wall_seconds=time.perf_counter()-started)
        write_json(index_path,plan)
    except Exception as exc:error=repr(exc);raise
    finally:write_json(receipt_path,dict(status='COMPLETE_NOT_UPLOADED'if error is None else'FAILED_PARTIAL_EVIDENCE_PRESERVED',error=error,created_assets=done,wall_seconds=time.perf_counter()-started,index=str(index_path)if index_path.exists()else None,digests_used=False,numeric_or_scoring_work=False,upload_executed=False))
    print(json.dumps(dict(status=plan['status'],assets=len(done),index=str(index_path),wall_seconds=time.perf_counter()-started)))

def restore(args):
    index_path=Path(args.index).resolve();index=read_json(index_path);root=Path(args.target).resolve();root.mkdir(parents=True,exist_ok=True)
    mapping={}
    for asset in index['assets']:
        if asset.get('status')!='CREATED':raise ValueError('Index includes an asset that was not created')
        asset_name=PurePosixPath(asset['name'])
        if asset_name.is_absolute() or len(asset_name.parts)!=1:raise ValueError('Unsafe asset filename')
        payload=index_path.parent/asset['name']
        if not payload.is_file()or payload.stat().st_size!=asset['size_bytes']:raise ValueError('Missing asset or incorrect published file size: '+asset['name'])
        for member in asset['members']:mapping.setdefault(member['source'],[]).append((asset['name'],member))
    inventory={f['path']:f for g in index['groups'] for f in g['files']}
    restored=0;reused=0;timestamp_limits=[]
    for logical,chunks in sorted(mapping.items()):
        rel=PurePosixPath(logical)
        if rel.is_absolute()or'..'in rel.parts:raise ValueError('Unsafe restoration path')
        destination=root/rel
        if not within(destination,root):raise ValueError('Restore path escapes target')
        if logical not in inventory:raise ValueError('Missing logical-file inventory')
        exists=destination.exists()
        if exists and (not args.reuse_identical or destination.is_symlink() or not destination.is_file()):raise FileExistsError('Restore preserves existing files; use empty staging or explicit --reuse-identical: '+str(destination))
        destination.parent.mkdir(parents=True,exist_ok=True);partial=destination.with_suffix(destination.suffix+'.release-partial')
        if exists and destination.stat().st_size!=inventory[logical]['size_bytes']:raise FileExistsError('Existing file differs in size; preserved: '+str(destination))
        with (destination.open('rb')if exists else partial.open('xb'))as out:
            offset=0
            for name,member in sorted(chunks,key=lambda x:x[1]['offset']):
                if member['offset']!=offset:raise ValueError('Chunk offset gap/overlap')
                with zipfile.ZipFile(index_path.parent/name)as archive,archive.open(member['member'])as reader:
                    if archive.getinfo(member['member']).file_size!=member['length_bytes']:raise ValueError('ZIP member length differs from index')
                    if exists:
                        left=member['length_bytes']
                        while left:
                            block=reader.read(min(BLOCK,left))
                            if not block or out.read(len(block))!=block:raise FileExistsError('Existing bytes differ; preserved without overwrite: '+str(destination))
                            left-=len(block)
                    else:copy_exact(reader,out,member['length_bytes'])
                offset+=member['length_bytes']
            if offset!=inventory[logical]['size_bytes']:raise ValueError('Restored length differs from canonical inventory')
        if not exists:partial.replace(destination)
        else:reused+=1
        recorded=inventory[logical]['mtime_ns']
        try:os.utime(destination,ns=(recorded,recorded))
        except OSError as exc:timestamp_limits.append(dict(path=logical,recorded_mtime_ns=recorded,error=repr(exc)))
        else:
            actual=destination.stat().st_mtime_ns
            if actual!=recorded:timestamp_limits.append(dict(path=logical,recorded_mtime_ns=recorded,actual_mtime_ns=actual))
        restored+=1
    receipt=dict(status='RESTORED_OPAQUE_BYTES_NO_ARRAY_DECODING',files=restored,reused_identical_files=reused,target=str(root),timestamp_limits=timestamp_limits,original_seals_rewritten=False,reproduction_scoring='create a separate new local seal from canonical pairs or perform a new public replay; retain original seals as historical provenance')
    write_json(root/f'R4_RESTORE_RECEIPT_{time.time_ns()}.json',receipt)
    print(json.dumps(receipt))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',nargs='?',choices=['plan','pack','restore'],default='plan');parser.add_argument('--output',type=Path);parser.add_argument('--max-bytes',type=int,default=MAX_BYTES);parser.add_argument('--include-atlas',action='store_true');parser.add_argument('--atlas-dir',default='figures/atlas');parser.add_argument('--repository-url');parser.add_argument('--publication-ref');parser.add_argument('--historical-e0-release-url');parser.add_argument('--reuse-identical',action='store_true');parser.add_argument('--index',type=Path);parser.add_argument('--target',type=Path);args=parser.parse_args()
    if args.historical_e0_release_url and not(args.historical_e0_release_url.startswith('https://') and '/releases/'in args.historical_e0_release_url):parser.error('historical E0 requires an explicit direct HTTPS Release page/asset URL')
    if not 1024*1024<=args.max_bytes<=MAX_BYTES:parser.error('--max-bytes must be between1MiB and1500000000bytes')
    if args.action=='restore':
        if not args.index or not args.target:parser.error('restore requires --index and --target')
        restore(args)
    elif args.action=='pack':
        if not args.output:parser.error('pack requires --output fresh empty directory')
        pack(args)
    else:
        plan=build_plan(args);write_json(STUDY/'audit/RELEASE_PACKAGE_PLAN.json',plan)
        print(json.dumps(dict(status=plan['status'],pack_ready=plan['gates']['pack_ready'],logical_files=plan['file_count'],source_bytes=plan['source_bytes'],planned_assets=len(plan['assets']),groups=[dict(name=g['name'],files=g['file_count'],bytes=g['source_bytes'],status=g['status'])for g in plan['groups']],dependency_issues=plan['dependency_issues'],zip_assets_created=0)))
if __name__=='__main__':main()

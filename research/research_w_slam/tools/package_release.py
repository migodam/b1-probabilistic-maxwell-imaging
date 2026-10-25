#!/usr/bin/env python3
"""Package actual saved evidence; no solver, network, hashes or source mutation."""
import argparse
import json
import os
from pathlib import Path,PurePosixPath
import shutil
import tempfile
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
BOUNDED_TREES=('inputs','legacy_inputs','results','figures','source_package','research/delegated')
EXCLUDED={'.git','.aws','.codex','.agents','.venv','.deps','__pycache__','node_modules','release_assets'}
EXTENSIONS={'.npz','.npy','.png'}


def write_json(path,value):
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n');temporary.replace(path)


def file_group(relative):
    p=relative.lower()
    if p.startswith(('inputs/','legacy_inputs/')):return 'inputs'
    if p.startswith(('results/e0/','results/e1/')) or 'actual_e0_e1' in p:return 'E0E1'
    if p.startswith('results/a0/') or '/a0/' in p:return 'A0'
    if p.startswith('results/e2e3/'):return 'controlled_E2E3'
    if p.startswith(('results/e4/','results/e5/')):return 'E4E5'
    return 'auxiliary_saved_evidence'


def enumerate_bounded(root):
    files=[]
    for tree in BOUNDED_TREES:
        start=root/tree
        if not start.exists():continue
        for directory,dirs,names in os.walk(start,followlinks=False):
            dirs[:]=sorted(d for d in dirs if d not in EXCLUDED and not (Path(directory)/d).is_symlink())
            for name in sorted(names):
                path=Path(directory)/name
                if not path.is_symlink() and path.suffix.lower() in EXTENSIONS:
                    files.append(dict(path=path.relative_to(root).as_posix()))
    return files


def read_optional(path):
    return json.loads(path.read_text()) if path.exists() else None


def terminal_exception(section,status_sources):
    """Only explicit matching terminal failure/not-run declarations qualify."""
    accepted=('NOT_RUN','BLOCKED','FAILED','FAIL','STOPPED','NOT_EXECUTED')
    records=[]
    def visit(value,pointer=''):
        if isinstance(value,dict):
            identity=' '.join((pointer,str(value.get('stage','')),str(value.get('section','')))).lower()
            matches=any(tag in identity for tag in ('e2e3',section,'finite_reference')) or any('/'+tag+'/' in identity or identity.rstrip().endswith('/'+tag) for tag in ('e2','e3'))
            statuses=[str(value.get(k,'')).upper() for k in ('status','decision','execution_status')]
            if matches and (any(any(tag in s for tag in accepted) for s in statuses) or value.get('executed') is False and any(value.get(k) for k in ('reason','blocked_by','status'))):
                records.append(dict(pointer=pointer,declaration=value))
            for key,child in value.items():visit(child,pointer+'/'+str(key))
        elif isinstance(value,list):
            for i,child in enumerate(value):visit(child,pointer+'/'+str(i))
        elif isinstance(value,str) and any(tag in pointer.lower() for tag in ('e2e3','bank','inference','finite_reference','/e2','/e3')) and any(tag in value.upper() for tag in accepted):
            records.append(dict(pointer=pointer,declaration=value))
    result=[]
    for source,data in status_sources:
        records=[];visit(data)
        result.extend(dict(source=source,**row) for row in records)
    return result


def final_prerequisites(root,args,status_sources):
    bank=read_optional(root/'results/e2e3/finite_bank/COMPLETE.json')
    checkpoint=read_optional(root/'results/e2e3/finite_bank/CHECKPOINT.json')
    costs=read_optional(root/'results/e2e3/finite_bank/online/JOB_COST_FINAL.json')
    inference=read_optional(root/'results/e2e3/inference/gate.json')
    bank_ok=bool(bank and checkpoint and costs and bank.get('status')=='COMPLETE' and checkpoint.get('status')=='COMPLETE' and checkpoint.get('completed')==2950 and costs.get('complete') is True)
    inference_ok=bool(inference and inference.get('prior_id')=='FINITE_PRIOR_ENUM4_V1' and inference.get('decision') in ('PASS','FAIL') and 'reliable_actual_reference' in inference)
    reports={}
    for section,ok in [('bank',bank_ok),('inference',inference_ok)]:
        declarations=terminal_exception(section,status_sources) if args.allow_not_run else []
        if not ok and not declarations:
            raise RuntimeError(f'final mode requires terminal {section} evidence, or --allow-not-run with explicit matching terminal declarations; no missing arrays fabricated')
        reports[section]=dict(complete_saved_evidence=ok,explicit_terminal_exception=declarations if not ok else [],
            scope='Terminal saved evidence or declared failed/blocked/not-run publication, never inferred scientific PASS')
    ledger=read_optional(root/'results/budget_ledger.json')
    if ledger:
        events=ledger.get('events',[])
        for i,event in enumerate(events):
            if event.get('kind')!='reservation' or event.get('source')!='run_job':continue
            commands=[str(x) for x in event.get('command',[])]
            if not any(Path(x).name in ('finite_bank.py','e2_revisit.py') for x in commands):continue
            done=any(e.get('kind')=='job_completion' and e.get('label')==event.get('label') and e.get('stage')==event.get('stage') for e in events[i+1:])
            if not done:raise RuntimeError('bank/inference wrapper still lacks terminal receipt: '+event.get('label',''))
    return reports


def run(args):
    root=Path(args.root).resolve();output=root/'release_assets'
    status_paths=[root/'GATE_SUMMARY.json']+[Path(p).resolve() for p in args.status_manifest]
    statuses=[(str(p),read_optional(p)) for p in status_paths if p.exists()]
    if args.mode=='final':prerequisites=final_prerequisites(root,args,statuses)
    else:prerequisites=dict(scope='Running/partial snapshot; no terminal or completeness claim')
    layout_path=Path(args.layout).resolve() if args.layout else root/'RELEASE_LAYOUT.json'
    layout=read_optional(layout_path)
    if layout is not None:
        entries=layout.get('files')
        if not isinstance(entries,list):raise ValueError('RELEASE_LAYOUT must contain files[] with stage-relative path')
        listed={row['path'] if isinstance(row,dict) else str(row) for row in entries}
        extra_arrays=[row for row in enumerate_bounded(root) if Path(row['path']).suffix.lower() in ('.npz','.npy') and row['path'] not in listed]
        entries=entries+extra_arrays
        selection_source=str(layout_path)+' plus all actual bounded NPZ/NPY evidence (retains unlisted late/failed/full-field arrays)'
    else:entries=enumerate_bounded(root);selection_source='bounded scientific tree enumeration';extra_arrays=[]
    selected={};missing=[];changed=[]
    for raw in entries:
        relative=raw['path'] if isinstance(raw,dict) else str(raw)
        posix=PurePosixPath(relative)
        if posix.is_absolute() or '..' in posix.parts or '\\' in relative:raise ValueError('unsafe archive path: '+relative)
        if any(part in EXCLUDED for part in posix.parts):raise ValueError('excluded/non-scientific archive path: '+relative)
        path=root/relative
        if path.suffix.lower() not in EXTENSIONS:continue
        if path.is_symlink() or not path.resolve().is_relative_to(root):raise ValueError('archive member escapes configured root: '+relative)
        if not path.is_file():missing.append(relative);continue
        actual=path.stat().st_size
        expected=raw.get('bytes',raw.get('size_bytes')) if isinstance(raw,dict) else None
        if expected is not None and int(expected)!=actual:changed.append(dict(path=relative,layout_bytes=expected,actual_bytes=actual))
        selected[posix.as_posix()]=dict(source_path=posix.as_posix(),file_bytes=actual,group=file_group(relative))
    if args.mode=='final' and missing and not args.allow_not_run:raise RuntimeError('final layout has missing saved files; no replacement solves or invented completeness')
    if not selected:raise RuntimeError('no actual saved evidence selected')
    output.mkdir(parents=True,exist_ok=True)
    manifest=dict(schema='WSLAM_RELEASE_ASSETS_V1',mode=args.mode,
        status='TERMINAL_EVIDENCE_PACKAGE' if args.mode=='final' else 'PARTIAL_SNAPSHOT_PACKAGE',
        selection_source=selection_source,unlisted_actual_arrays_added=[row['path'] for row in extra_arrays],unpack_relative_to='research/research_w_slam/',
        archive_member_rule='Stage-relative path; no additional enclosing directory',
        scientific_completeness_claim=False,scientific_PASS_required=False,prerequisites=prerequisites,
        missing_layout_paths=sorted(missing),layout_size_changes=changed,archives=[],
        preserves_failed_raw_arrays=True,full_field_snapshots_selected_by_same_rules=True,
        raw_arrays_parsed=False,hashes_generated=False,source_deleted=False,
        git_exclusion='/release_assets/ in root .gitignore',
        terminal_not_run_publication_allowed=bool(args.allow_not_run))
    groups={}
    for relative,row in sorted(selected.items()):groups.setdefault(row['group'],[]).append(row)
    started=time.perf_counter();cpu=time.process_time()
    for group,members in sorted(groups.items()):
        archive=output/f'research-w-slam-{group}.zip'
        fd,tmp=tempfile.mkstemp(prefix=archive.name+'.',dir=output);os.close(fd)
        try:
            with zipfile.ZipFile(tmp,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as zipped:
                for row in members:
                    source=root/row['source_path']
                    if source.stat().st_size!=row['file_bytes']:raise RuntimeError('source changed during packaging: '+row['source_path'])
                    info=zipfile.ZipInfo(row['source_path'],date_time=(1980,1,1,0,0,0))
                    info.compress_type=zipfile.ZIP_DEFLATED;info.external_attr=0o100644<<16
                    with source.open('rb') as incoming,zipped.open(info,'w',force_zip64=True) as outgoing:
                        shutil.copyfileobj(incoming,outgoing,length=1024*1024)
                    if source.stat().st_size!=row['file_bytes']:raise RuntimeError('source changed during packaging: '+row['source_path'])
                saved=zipped.infolist()
                records=[dict(**row,archive_member_bytes=info.file_size,compressed_member_bytes=info.compress_size) for row,info in zip(members,saved)]
            os.replace(tmp,archive)
        finally:
            if os.path.exists(tmp):os.unlink(tmp)
        manifest['archives'].append(dict(file=archive.name,archive_bytes=archive.stat().st_size,members=records,
            uncompressed_source_bytes=sum(row['file_bytes'] for row in members),member_count=len(members)))
        write_json(output/'ASSET_MANIFEST.json',manifest)
    manifest.update(artifact_creation_complete=True,wall_seconds=time.perf_counter()-started,
        diagnostic_process_cpu_seconds=time.process_time()-cpu,
        CPU_accounting='Enclosing wrapper only; no manual or duplicate ledger charge')
    write_json(output/'ASSET_MANIFEST.json',manifest)
    print(json.dumps(dict(status=manifest['status'],manifest=str(output/'ASSET_MANIFEST.json'),archives=len(groups),members=len(selected))))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',default=str(ROOT))
    parser.add_argument('--mode',choices=('partial','final'),default='partial');parser.add_argument('--layout')
    parser.add_argument('--allow-not-run',action='store_true',help='Publish explicit terminal failed/blocked stages; requires matching status declarations')
    parser.add_argument('--status-manifest',action='append',default=[],help='Explicit terminal scientific status document, e.g. GATE_SUMMARY.json')
    run(parser.parse_args(argv))


if __name__=='__main__':main()

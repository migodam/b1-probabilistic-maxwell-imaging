"""Prepare/build the R3 offline audit ZIP; never publishes or modifies sources.

Default --plan does no archive I/O and does not parse observations/truths.
--build is reserved for final documents near18:20UTC on root's instruction.
"""
from pathlib import Path
import argparse,datetime,json,os,time,zipfile
ROOT=Path(__file__).resolve().parents[1]
SKIP_PARTS={'__pycache__','.pytest_cache','.mypy_cache','.git','.mpl','.cache'}
def exclusion(path,archive):
    rel=path.relative_to(ROOT);low=[p.lower() for p in rel.parts];name=path.name.lower()
    if 'PATH_ONLY_REPAIR_PLACEHOLDERS' in rel.parts:return 'metadata-only path preparation failure retained locally; manifest only'
    if path.is_symlink():return 'symlink not followed; source remains local'
    if path.resolve()==archive.resolve():return 'output ZIP excluded from its own archive'
    if any(p in SKIP_PARTS for p in low):return 'cache or repository metadata'
    if name in {'package_plan.json','package_manifest.json','package_build_receipt.json'}:return 'derived packaging metadata included separately, not recursively'
    if path.suffix.lower() in {'.pyc','.pyo','.lock','.tmp','.part'} or name.endswith('~') or name=='.ds_store':return 'cache, lock or temporary file'
    if path.suffix.lower()=='.zip' and 'delivery' in low:return 'prior delivery archive retained locally'
    if name.startswith('chen_') and path.suffix.lower() in {'.pdf','.txt','.png'}:return 'private CHEN book source/excerpt'
    if path.suffix.lower() in {'.pdf','.txt','.png'} and any(('book' in p or p.startswith('chen')) for p in low):return 'third-party book/excerpt material'
    return None
def license_sources():
    r2=ROOT.parent/'research_w_gq_r2';roots=[ROOT,ROOT.parents[1],r2,r2/'reuse/fixed_vendor',r2/'reuse/previous_gq']
    found=[]
    for base in roots:
        if not base.exists():continue
        for p in base.iterdir():
            if p.is_file() and (p.name.upper().startswith('LICENSE') or p.name.upper().startswith('COPYING')):found.append(p)
    return sorted(set(found))
def main():
    ap=argparse.ArgumentParser();mode=ap.add_mutually_exclusive_group();mode.add_argument('--plan',action='store_true');mode.add_argument('--build',action='store_true');ap.add_argument('--output',type=Path,default=ROOT/'delivery/RESEARCH_W_GQ_R3_OFFLINE_AUDIT.zip');a=ap.parse_args();archive=a.output.resolve();started=time.time();included=[];excluded=[]
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file() and not p.is_symlink():continue
        reason=exclusion(p,archive)
        if reason:excluded.append(dict(path=str(p.relative_to(ROOT)),reason=reason));continue
        included.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,role='offline truth: opaque bytes only' if 'private' in p.relative_to(ROOT).parts else 'R3 artifact/source'))
    licenses=license_sources();failures=sorted({str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if any(s in p.name.upper() for s in ['INVALID','REPAIR','FAILURE','ABLATION'])})
    manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='R3 offline audit package; not a live online input bundle',publication='NOT_PUBLISHED; ZIP only',files=included,excluded=excluded,total_uncompressed_source_bytes=sum(r['bytes'] for r in included),offline_truth_handling='Included only as opaque file bytes after online methods froze; packaging never loads or parses truth arrays',failure_retention='All original/repair/failed attempts remain locally in-place; included when within package scope',failure_markers=failures,license_files_found=[str(p) for p in licenses],license_source_records=[dict(path=str(p),kind='licensing status/provenance notice' if any(x in p.name.upper() for x in ['STATUS','PROVENANCE','NOTICE']) else 'existing LICENSE/COPYING file') for p in licenses],license_policy='Copy existing scoped LICENSE/COPYING files with source provenance only; no new MIT or inferred license assigned',third_party_history='No old R2 directory copied; portable contains the12actually imported vendor source files',scientific_source_unchanged=True,no_extra_physics=True)
    (ROOT/'PACKAGE_PLAN.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if not a.build:print(json.dumps(dict(mode='PLAN_ONLY',files=len(included),bytes=manifest['total_uncompressed_source_bytes'],excluded=len(excluded),licensing_source_files=len(licenses),archive_not_created=str(archive))));return
    archive.parent.mkdir(parents=True,exist_ok=True)
    if archive.exists():raise FileExistsError('Refuse replacing an existing delivery ZIP')
    temporary_archive=archive.with_name(archive.name+'.'+str(int(time.time()))+'.part')
    with zipfile.ZipFile(temporary_archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as z:
        for row in included:z.write(ROOT/row['path'],'research_w_gq_r3_joint/'+row['path'])
        z.writestr('research_w_gq_r3_joint/PACKAGE_MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
        for i,p in enumerate(licenses):
            if p.is_relative_to(ROOT):continue
            z.write(p,f'research_w_gq_r3_joint/INHERITED_LICENSE_SOURCES/{i:02d}_{p.name}')
    if archive.exists():raise FileExistsError('Final archive appeared during build; completed temporary archive retained')
    temporary_archive.replace(archive)
    with zipfile.ZipFile(archive) as z:entry_count=len(z.infolist());stored_bytes=sum(x.file_size for x in z.infolist())
    receipt=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),archive=str(archive),zip_bytes=archive.stat().st_size,entries=entry_count,stored_uncompressed_bytes=stored_bytes,wall_seconds=time.time()-started,source_files_unchanged=True,original_failures_retained=True,published=False,truth_arrays_parsed=False,validation='ZIP directory metadata readback only; no SHA checks or extra physics')
    (ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n');(ROOT/'PACKAGE_BUILD_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
if __name__=='__main__':main()

#!/usr/bin/env python3
"""Inventory or package scoped R7 evidence without reading array contents.

Dry inventory is the default. --write requires the parent's explicit seal file;
this script neither scores, publishes, runs physics, hashes, nor uses Git's ignore
rules. ZIP members are streamed, so gitignored raw arrays/traces are retained.
"""
import argparse
import datetime
import fnmatch
import json
import os
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_DEST=Path("/Volumes/migodam's-external-brain/Research/Inv_SLAM/research/delegated/max_struct_r7_artifacts")


def exclusion(rel, sealed):
    parts=rel.parts; name=rel.name.lower(); s=rel.as_posix()
    if any(x in parts for x in ['.git','__pycache__','.venv','node_modules']):return 'runtime/cache or Git internals'
    if name.endswith(('.pyc','.lock')):return 'runtime/cache'
    if name.endswith(('.zip','.tar','.tgz','.gz','.7z','.rar')):return 'input or copied archive; original third-party handoff archives are not redistributable'
    if s.startswith('data_public/fresnel/') and name.endswith('.exp'):return 'raw Fresnel measurements; redistribution rights not established'
    if fnmatch.fnmatch(s,'RESULTS/fresnel_*/endpoint.npz'):return 'contains measured_Y; use fresnel_public_derived.npz'
    if name.endswith('.pdf') and not s.startswith('FIGURES/'):return 'PDF outside original generated FIGURES; conservatively exclude third-party copied documents'
    if s.startswith('delegated/public_literature/') and name.endswith(('.pdf','.html','.htm','.txt')):return 'third-party publication/full text'
    if name.endswith(('.pdf','.html','.htm','.txt')) and any(x in s.lower() for x in ['fulltext','full_text','paper_source','download','handoff','pro_input']):return 'copied third-party full text or handoff'
    if not sealed and (any('truth' in x.lower() or x.lower() in ['private','sealed','gt'] for x in parts) or 'mapping' in name):return 'synthetic truth/private mapping withheld until parent seal'
    return None


def inventory(root,sealed):
    included=[];excluded=[]
    for base,dirs,files in os.walk(root,followlinks=False):
        # File-level exclusions are explicit in the manifest; avoid traversing caches.
        for directory in list(dirs):
            p=Path(base)/directory
            if p.is_symlink():
                excluded.append(dict(path=p.relative_to(root).as_posix(),reason='symlink excluded; no target traversal'));dirs.remove(directory)
            elif directory in ['.git','__pycache__','.venv','node_modules']:
                excluded.append(dict(path=p.relative_to(root).as_posix()+'/',reason='runtime/cache or Git internals'));dirs.remove(directory)
        for name in sorted(files):
            p=Path(base)/name;rel=p.relative_to(root)
            reason='symlink excluded; no target traversal' if p.is_symlink() else exclusion(rel,sealed)
            if reason:excluded.append(dict(path=rel.as_posix(),reason=reason));continue
            included.append(dict(path=rel.as_posix(),bytes=p.stat().st_size))
    return sorted(included,key=lambda r:r['path']),sorted(excluded,key=lambda r:r['path'])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,default=ROOT);ap.add_argument('--output-dir',type=Path,default=DEFAULT_DEST);ap.add_argument('--manifest-output',type=Path,required=True);ap.add_argument('--write',action='store_true');ap.add_argument('--publish-seal',type=Path);ap.add_argument('--max-part-bytes',type=int,default=1_900_000_000);args=ap.parse_args()
    root=args.root.resolve();seal=None
    if args.publish_seal:
        seal=json.loads(args.publish_seal.read_text())
        if seal.get('status')!='ALL_PRIMARY_SEALED' or seal.get('release_authorized') is not True:ap.error('seal must state status=ALL_PRIMARY_SEALED and release_authorized=true, written by parent')
    if args.write and seal is None:ap.error('--write requires explicit parent --publish-seal')
    if not 1_000_000 <= args.max_part_bytes < 2_000_000_000:ap.error('part limit must be >=1 MB and strictly below 2 GB')
    if args.output_dir.resolve().is_relative_to(root):ap.error('output must be outside R7 source tree')
    included,excluded=inventory(root,seal is not None)
    manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='PACKAGED' if args.write else 'INVENTORY_ONLY',root_scope='research/research_max_struct_r7',synthetic_truth_authorized=seal is not None,publish_seal=str(args.publish_seal) if seal else None,raw_arrays_opened=False,hashes_computed=False,included=included,excluded=excluded,release_files=[],max_part_bytes=args.max_part_bytes,notes=['All R7 files, including Git-ignored raw traces/states, are inventoried independently of Git.','No raw Fresnel endpoint or third-party source archive is included.','Self-generated figure PDFs remain included.','Parent must publish small Git files and upload these Release assets; this script does neither.'])
    if args.write:
        args.output_dir.mkdir(parents=True,exist_ok=True)
        # Group by uncompressed size, then deflate to reduce publication transfer. Budget
        # exact local/central-directory overhead per UTF-8 name plus a safety margin.
        groups=[];batch=[];size=4096
        for row in included:
            cost=row['bytes']+2*len(('research/research_max_struct_r7/'+row['path']).encode('utf-8'))+512
            if cost+4096>args.max_part_bytes:raise ValueError('single file exceeds part limit: '+row['path'])
            if size+cost>args.max_part_bytes:groups.append(batch);batch=[];size=4096
            batch.append(row);size+=cost
        if batch:groups.append(batch)
        for i,batch in enumerate(groups,1):
            dest=args.output_dir/('max-struct-r7-evidence-part%03d.zip'%i)
            if dest.exists():raise FileExistsError('refusing to overwrite '+str(dest))
            with zipfile.ZipFile(dest,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
                for row in batch:
                    p=root/row['path']
                    if p.is_symlink() or p.stat().st_size!=row['bytes']:raise RuntimeError('source changed after inventory: '+row['path'])
                    archive.write(p,'research/research_max_struct_r7/'+row['path'])
            if dest.stat().st_size>=args.max_part_bytes:raise RuntimeError('part exceeds declared size: '+str(dest))
            manifest['release_files'].append(dict(filename=dest.name,bytes=dest.stat().st_size,members=len(batch)))
    args.manifest_output.parent.mkdir(parents=True,exist_ok=True);args.manifest_output.write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(status=manifest['status'],included=len(included),excluded=len(excluded),release_files=manifest['release_files'])))
if __name__=='__main__':main()

"""Plan now; explicit metadata inventory/serial tar.gz only after formal closure.

No SHA/content hashing, source mutation, model import, solver, or truth analysis.
This script has deliberately NOT been executed in inventory/pack mode during D.
"""
from __future__ import annotations
import argparse,io,json,os,stat,tarfile,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
TRACK=HERE.parent
OUT=TRACK.parent


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')


def require_closed():
    if not (OUT/'CLOSURE.json').exists():raise RuntimeError('HOLD: formal CLOSURE.json absent; no raw inventory/archive during D')
    ledger=json.loads((OUT/'RESOURCE_LEDGER_SNAPSHOT.json').read_text())
    if ledger.get('active_jobs') or ledger.get('active_actions'):raise RuntimeError('HOLD: central job/action still active')


def roots(stage):
    selected=[OUT/f'track_{stage.lower()}',OUT/'configs',OUT/'shared',OUT/'src']
    for name in ('RESOURCE_LEDGER.jsonl','RESOURCE_LEDGER_SNAPSHOT.json','CLOSURE.json'):
        if (OUT/name).exists():selected.append(OUT/name)
    for start in (OUT/'raw/jobs').glob('*/START.json'):
        if json.loads(start.read_text()).get('stage')==stage:selected.append(start.parent)
    groups=['CDEV'] if stage=='C' else ([p.stem.removeprefix('parents_') for p in (OUT/'configs').glob('parents_D*.json')] if stage=='D' else [])
    for group in groups:
        manifest=OUT/f'configs/parents_{group}.json'
        for row in json.loads(manifest.read_text())['parents']:
            for directory in ('online','offline_truth'):
                path=OUT/'raw'/directory/f"{row['parent']}.npz"
                if path.exists():selected.append(path)
        receipt=OUT/f'raw/GENERATION_{group}.json'
        if receipt.exists():selected.append(receipt)
    return selected


def metadata(path):
    value=path.lstat();mode=value.st_mode
    kind='symlink' if stat.S_ISLNK(mode) else ('directory' if stat.S_ISDIR(mode) else 'file')
    if kind=='file' and not stat.S_ISREG(mode):raise ValueError(f'unsupported source type: {path}')
    row=dict(path=str(path.relative_to(OUT)),kind=kind,bytes=value.st_size,mtime_ns=value.st_mtime_ns,mode=stat.S_IMODE(mode))
    if kind=='symlink':row['link_target']=os.readlink(path)
    return row


def inventory(stage,output):
    require_closed();paths={}
    for root in roots(stage):
        if not root.exists() and not root.is_symlink():continue
        if root.is_dir() and not root.is_symlink():
            for directory,subdirs,files in os.walk(root,followlinks=False):
                directory=Path(directory)
                subdirs[:]=[name for name in subdirs if (directory/name).resolve()!=output.resolve()]
                for path in [directory,*[directory/name for name in files],*[directory/name for name in subdirs if (directory/name).is_symlink()]]:
                    paths[str(path.relative_to(OUT))]=metadata(path)
        else:paths[str(root.relative_to(OUT))]=metadata(root)
    rows=[paths[key] for key in sorted(paths)]
    result=dict(stage=stage,status='METADATA_INVENTORY_ONLY',created_unix=time.time(),entries=rows,
        regular_files=sum(row['kind']=='file' for row in rows),regular_file_bytes=sum(row['bytes'] for row in rows if row['kind']=='file'),
        SHA_checks=False,content_hashes=False,external_historical_assets='referenced by provenance; not automatically copied; no self-contained cold-rebuild claim')
    path=output/f'STAGE_{stage}_INVENTORY.json'
    if path.exists():raise FileExistsError('preserve previous inventory; choose a new output directory')
    dump(path,result);return result


def pack(stage,output):
    require_closed();manifest_path=output/f'STAGE_{stage}_INVENTORY.json'
    manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else inventory(stage,output)
    if manifest['stage']!=stage:raise ValueError('inventory stage mismatch')
    destination=output/f'STAGE_{stage}.tar.gz';partial=output/f'STAGE_{stage}.tar.gz.partial'
    if destination.exists() or partial.exists():raise FileExistsError('existing archive/partial preserved; no automatic overwrite')
    wall,cpu=time.perf_counter(),time.process_time();status,error='COMPLETE',None
    try:
        with partial.open('xb') as stream,tarfile.open(fileobj=stream,mode='w:gz',compresslevel=6,dereference=False) as archive:
            for row in manifest['entries']:
                source=OUT/row['path']
                if metadata(source)!=row:raise RuntimeError(f'source changed after inventory: {row["path"]}')
                archive.add(source,arcname=row['path'],recursive=False)
            content=json.dumps(manifest,indent=2,ensure_ascii=False).encode()
            header=tarfile.TarInfo('PACKAGE_METADATA.json');header.size=len(content);header.mtime=int(time.time())
            archive.addfile(header,io.BytesIO(content))
        partial.replace(destination)
    except BaseException as exc:status,error='FAILED_PARTIAL_PRESERVED',repr(exc);raise
    finally:
        dump(output/f'STAGE_{stage}_PACK_RECEIPT.json',dict(stage=stage,status=status,error=error,
            wall_seconds=time.perf_counter()-wall,CPU_seconds=time.process_time()-cpu,
            archive=str(destination),partial=str(partial),SHA_checks=False,new_physics_RHS=0,
            source_files_modified=False,inventory=str(manifest_path)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['A','B','C','D'],required=True)
    parser.add_argument('--mode',choices=['plan','inventory','pack'],default='plan')
    parser.add_argument('--output',type=Path,default=HERE/'output')
    args=parser.parse_args()
    if not args.output.resolve().is_relative_to(HERE.resolve()):raise ValueError('packaging outputs must stay under packaging/')
    if args.mode=='plan':
        print(json.dumps(dict(status='PLAN_ONLY_NO_SCAN_NO_COMPRESSION',stage=args.stage,
            future_primary_root=f'track_{args.stage.lower()}',include='same-stage raw/jobs, common sources/config/ledger; C/D corresponding parent raw observations/truth only after closure',
            output=str(args.output),new_RHS=0,SHA_checks=False),indent=2));return
    args.output.mkdir(parents=True,exist_ok=True)
    if args.mode=='inventory':
        result=inventory(args.stage,args.output);print(json.dumps({k:v for k,v in result.items() if k!='entries'},indent=2))
    else:pack(args.stage,args.output)


if __name__=='__main__':main()

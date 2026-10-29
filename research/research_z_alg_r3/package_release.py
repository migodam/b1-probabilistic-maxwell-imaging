#!/usr/bin/env python3
"""Build the full R3 reproducibility archive and a curated reading kit.

No model, forward solver, derivative, scorer, network request, or hash check
is executed. Default invocation is a read-only packaging plan; pass --build
only after the final reports and rendering review are complete.
"""
import argparse
import json
from pathlib import Path, PurePosixPath
import sys
import zipfile


ROOT = Path(__file__).resolve().parent
PREFIX = 'research_z_alg_r3'
FULL_NAME = 'FULL_REPRODUCIBILITY.zip'
READING_NAME = 'Research_Z_ALG_R3_Reading_Kit.zip'
MANIFEST_NAMES = {'BUNDLE_MANIFEST.json','READING_KIT_MANIFEST.json'}
MAX_FILE_BYTES = 100_000_000
FULL_DIRECTORIES = {'src','experiments','public_data','results','audit','theory',
                    'figures','truth_after_freeze','frozen_predictions','reports',
                    'snapshots','assets','docs','site'}
FORBIDDEN_PARTS = {'private_truth','private','.git','__pycache__','.pytest_cache',
                   '.mypy_cache','.venv','node_modules'}
FORBIDDEN_JSON_KEYS = {'seed_decimal','private_seed','generation_seed',
                       'private_generation_seed','generator_secret'}
READING_AUDIT_STEMS = {'FINAL_SCIENTIFIC_REVIEW','PAPER_NOVELTY',
                       'MECHANISM_AND_IMPLEMENTATION','HISTORY_AUDIT_CLOSURE',
                       'FINAL_REVIEW','EXPERIMENT_FINAL_VERIFICATION'}


def eligible(path, root=ROOT):
    relative = path.relative_to(root)
    if path.is_symlink():
        return False
    if any(part in FORBIDDEN_PARTS for part in relative.parts):
        return False
    lower = path.name.lower()
    if lower.endswith(('.zip','.pyc','.pyo','.tmp','.partial','.part')):
        return False
    if '.zip.' in lower or path.name in MANIFEST_NAMES:
        return False
    if lower in {'.env','.ds_store'} or lower.startswith(('generator_secret','sealed_targets')):
        return False
    if 'audit_worktrees_tmp' in lower:
        return False
    if len(relative.parts)>1 and relative.parts[0] not in FULL_DIRECTORIES:
        return False
    return True


def full_files(root=ROOT):
    result = []
    for path in sorted(root.rglob('*')):
        if path.is_file() and eligible(path,root):
            result.append((path,f'{PREFIX}/{path.relative_to(root).as_posix()}'))
    sibling = root.parent/'research_z_alg_r2'/'src'/'physics.py'
    if not sibling.is_file() or sibling.is_symlink():
        raise FileNotFoundError('Required minimal immutable generation dependency: '+str(sibling))
    result.append((sibling,'research_z_alg_r2/src/physics.py'))
    return sorted(result,key=lambda item:item[1])


def reading_files(full):
    selected = []
    for source, member in full:
        relative = PurePosixPath(member)
        if relative.parts[0] != PREFIX:
            continue
        parts = relative.parts[1:]
        name = parts[-1]
        suffix = Path(name).suffix.lower()
        include = False
        if len(parts)==1:
            include = suffix in {'.md','.html','.json','.csv'} or name in {
                'run_solver.py','run_portfolio.py','package_release.py','requirements.txt'}
        elif parts[0] in {'src','theory','figures'}:
            include = True
        elif parts[0]=='audit':
            include = len(parts)==2 and (Path(name).stem in READING_AUDIT_STEMS or
                       ('FINAL' in name.upper() and 'REVIEW' in name.upper()))
        elif parts[0]=='experiments':
            include = ('CONTRACT' in name.upper() and suffix=='.json') or name=='requirements.txt'
            include = include or ('configs' in '/'.join(parts) and suffix=='.json')
            include = include or (name.lower().startswith(('reproduce','replay')) and suffix in {'.py','.md','.txt','.sh'})
        elif parts[0]=='results':
            include = ('SUMMARY' in name.upper() and suffix in {'.md','.json','.csv'}) or name.upper().endswith('METRICS.CSV')
        if include:
            selected.append((source,member))
    return selected


def inspect_private_fields(value, filename):
    if isinstance(value,dict):
        for key,item in value.items():
            if key.lower() in FORBIDDEN_JSON_KEYS:
                raise ValueError(f'Private generation field {key!r} in export candidate {filename}')
            inspect_private_fields(item,filename)
    elif isinstance(value,list):
        for item in value:
            inspect_private_fields(item,filename)


def validate_payload(files):
    names = set()
    for source,member in files:
        if member in names:
            raise ValueError('Duplicate archive member: '+member)
        names.add(member)
        if source.stat().st_size>MAX_FILE_BYTES:
            raise ValueError('Export member exceeds 100 MB: '+member)
        parts = PurePosixPath(member).parts
        if any(part in FORBIDDEN_PARTS for part in parts) or '..' in parts:
            raise ValueError('Forbidden archive path: '+member)
        if source.suffix.lower()=='.npz':
            with zipfile.ZipFile(source) as archive:
                arrays = set(archive.namelist())
            if 'truth_after_freeze' in parts and arrays!={'m.npy'}:
                raise ValueError('Disclosed truth must contain exactly m, without seeds: '+member)
            if any('seed' in Path(array).stem.lower() or 'secret' in Path(array).stem.lower() for array in arrays):
                raise ValueError('Private seed/secret array in NPZ: '+member)
        elif source.suffix.lower()=='.json':
            inspect_private_fields(json.loads(source.read_text()),member)


def manifest(files, kind, archive_name, own_manifest):
    entries = [dict(path=member,bytes=source.stat().st_size) for source,member in files]
    return dict(schema='R3_RELEASE_MANIFEST_V1',kind=kind,archive=archive_name,
                archive_root=PREFIX,inventory=entries,payload_file_count=len(entries),
                payload_bytes=sum(item['bytes'] for item in entries),
                inventory_uses='Paths and byte sizes only; no cryptographic hashes are computed.',
                manifest_member=f'{PREFIX}/{own_manifest}',
                manifest_inventory_exclusion='This manifest is included as an auxiliary archive member but is omitted from its own payload inventory. Both root manifest files are excluded during payload discovery, preventing self/cross recursion.',
                compression=dict(method='ZIP_DEFLATED',level=9,canonical_member_timestamp='1980-01-01 00:00:00'),
                maximum_export_file_bytes=MAX_FILE_BYTES,
                privacy=dict(private_truth_exported=False,private_generation_seeds_exported=False,
                             disclosed_truth_policy='truth_after_freeze NPZ members must be exactly {m.npy}'),
                scope=('Full R3 sources, experiment protocols/drivers/configs, cached public observations/operators, complete results, frozen predictions, audit snapshots, reports, figures and disclosed m-only truth.' if kind=='full' else
                       'Curated explanatory reports, root summaries/cost CSV, theory, selected final reviews, figures, core source/CLIs, contracts/configs and replay commands. Intentionally excludes raw solver results, public NPZ data, truth NPZs, frozen predictions and source snapshots; use the full archive for complete replay.'),
                generation_dependency=(dict(archive_member='research_z_alg_r2/src/physics.py',
                                             purpose='experiments/prepare_data.py imports the immutable sibling Model for fresh target generation',
                                             available_in_this_archive=True) if kind=='full' else
                                       dict(available_in_this_archive=False,full_archive='FULL_REPRODUCIBILITY.zip',
                                            purpose='Fresh generation and complete replay require the full archive.')))


def write_json(path,value):
    temporary = path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    temporary.replace(path)


def write_archive(path,files,manifest_path):
    temporary = path.with_name(path.name+'.tmp')
    members = files+[(manifest_path,f'{PREFIX}/{manifest_path.name}')]
    try:
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
            for source,member in sorted(members,key=lambda item:item[1]):
                info = zipfile.ZipInfo(member,date_time=(1980,1,1,0,0,0))
                info.external_attr = 0o100644<<16
                archive.writestr(info,source.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
        if temporary.stat().st_size>MAX_FILE_BYTES:
            raise ValueError('Compressed archive exceeds 100 MB: '+path.name)
        # Readback checks names, sizes and compression metadata; no hash pass.
        with zipfile.ZipFile(temporary) as archive:
            infos = archive.infolist()
            expected = {name:source.stat().st_size for source,name in members}
            actual = {info.filename:info.file_size for info in infos}
            if actual!=expected:
                raise ValueError('ZIP path/byte-size readback mismatch: '+path.name)
            receipt = dict(name=path.name,archive_bytes=temporary.stat().st_size,
                           archive_members=len(infos),uncompressed_bytes=sum(i.file_size for i in infos),
                           compressed_member_bytes=sum(i.compress_size for i in infos),
                           every_member_at_most_100MB=all(i.file_size<=MAX_FILE_BYTES for i in infos),
                           archive_at_most_100MB=True,
                           compression_ratio=temporary.stat().st_size/max(sum(i.file_size for i in infos),1))
        temporary.replace(path)
        return receipt
    finally:
        if temporary.exists():
            temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',action='store_true',help='Write manifests and archives after reports/render review are final')
    parser.add_argument('--replace',action='store_true',help='Explicitly replace existing owned manifests/ZIPs')
    args = parser.parse_args()
    full = full_files()
    reading = reading_files(full)
    validate_payload(full)
    summaries = [manifest(full,'full',FULL_NAME,'BUNDLE_MANIFEST.json'),
                 manifest(reading,'reading',READING_NAME,'READING_KIT_MANIFEST.json')]
    if not args.build:
        print(json.dumps(dict(mode='PLAN_ONLY_NO_WRITES',bundles=[{key:value for key,value in summary.items() if key in {'kind','archive','payload_file_count','payload_bytes','scope','generation_dependency'}} for summary in summaries]),indent=2))
        return
    outputs = [ROOT/'BUNDLE_MANIFEST.json',ROOT/'READING_KIT_MANIFEST.json',ROOT/FULL_NAME,ROOT/READING_NAME]
    existing = [str(path) for path in outputs if path.exists()]
    if existing and not args.replace:
        raise FileExistsError('Existing release outputs preserved; use --replace for an authorized rebuild: '+', '.join(existing))
    for summary,name in zip(summaries,('BUNDLE_MANIFEST.json','READING_KIT_MANIFEST.json')):
        write_json(ROOT/name,summary)
    receipts = [write_archive(ROOT/FULL_NAME,full,ROOT/'BUNDLE_MANIFEST.json'),
                write_archive(ROOT/READING_NAME,reading,ROOT/'READING_KIT_MANIFEST.json')]
    print(json.dumps(dict(mode='BUILT',archives=receipts,
                          manifests=['BUNDLE_MANIFEST.json','READING_KIT_MANIFEST.json'],
                          extraction='Unzip into one directory; cd research_z_alg_r3. The full ZIP also places its minimal immutable R2 generation dependency beside R3.'),indent=2))


if __name__=='__main__':
    main()

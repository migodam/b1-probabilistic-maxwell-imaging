"""Publish-directory preparation and exact stored-array package; no hashing."""
from pathlib import Path
import argparse, json, shutil, zipfile
root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser();parser.add_argument('--destination',required=True);parser.add_argument('--archive',required=True)
a=parser.parse_args();dest=Path(a.destination).resolve();archive=Path(a.archive).resolve()
if dest==root or root in dest.parents:raise ValueError('Destination must be outside source project')
dest.mkdir(parents=True,exist_ok=True);archive.parent.mkdir(parents=True,exist_ok=True)
def excluded(p):return any(x in ('.deps','.venv','__pycache__','.git','mpl','.DS_Store') for x in p.relative_to(root).parts)
arrays=sorted(p for p in root.rglob('*.npz') if not excluded(p))
manifest=dict(release_tag='research-w-continuation-20261023',asset=archive.name,
              download='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/research-w-continuation-20261023/'+archive.name,
              unpack_relative_to='research/research_w_continuation/',
              files=[dict(path=p.relative_to(root).as_posix(),bytes=p.stat().st_size) for p in arrays],
              array_files=len(arrays),total_array_bytes=sum(p.stat().st_size for p in arrays),
              scope='Actual saved arrays, including earlier failed/repaired variants; scientific eligibility is defined by reports/ledgers. No fabricated missing states.')
(root/'ARRAY_MANIFEST.json').write_text(json.dumps(manifest,indent=2))
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_STORED,allowZip64=True) as z:
    for p in arrays:z.write(p,p.relative_to(root).as_posix())
    z.write(root/'ARRAY_MANIFEST.json','ARRAY_MANIFEST.json')
copied=[]
for p in sorted(root.rglob('*')):
    if not p.is_file() or excluded(p) or p.suffix=='.npz':continue
    if p.suffix not in ('.py','.md','.json','.jsonl','.txt','.csv','.png','.svg','.html','.log'):continue
    rel=p.relative_to(root);target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);copied.append(dict(path=rel.as_posix(),bytes=p.stat().st_size))
audit=dict(copied_files=len(copied),files=copied,large_arrays_in_release=manifest['array_files'],
           third_party_installs_redistributed=False,original_P0_modified=False,
           no_locked_data=True,no_credentials_required_in_artifacts=True,
           publication_scope='New continuation directory only; add separate reading entry in research/README.md')
(root/'PUBLICATION_FILES.json').write_text(json.dumps(audit,indent=2));shutil.copy2(root/'PUBLICATION_FILES.json',dest/'PUBLICATION_FILES.json')
print(json.dumps(dict(public_files=len(copied),array_files=len(arrays),array_bytes=manifest['total_array_bytes'],archive_bytes=archive.stat().st_size,destination=str(dest))))

"""R13 byte-manifest wrapper around the existing R12 restoration implementation.

Downloads selected logical assets into a NEW/empty destination. Never runs
experiments, checks SHA, installs dependencies or writes into R12 directories.
"""
import argparse,json,runpy,stat,sys,zipfile,shutil
from pathlib import Path,PurePosixPath


def main():
    root=Path(__file__).resolve().parents[1];p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,default=root/'DOWNLOADS.json');p.add_argument('--asset',action='append',required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--extract-dir',type=Path);p.add_argument('--legacy-restorer',type=Path,default=root.parent/'research_max_struct_r12/restore_release_assets.py');a=p.parse_args()
    manifest=json.loads(a.manifest.read_text())
    if manifest.get('published') is not True:raise RuntimeError('R13 restore requires a published DOWNLOADS.json with actually verified public part URLs, not the private local transport plan.')
    assets={x['name']:x for x in manifest['logical_assets']}
    if len(set(a.asset))!=len(a.asset):raise ValueError('Select each logical archive once')
    for name in a.asset:
        if name not in assets:raise KeyError('Unknown logical asset '+name)
        if PurePosixPath(name).name!=name:raise ValueError('Logical archive name must be a basename')
        asset=assets[name]
        if sum(part['bytes'] for part in asset['parts'])!=asset['bytes']:raise ValueError('Manifest byte sum mismatch '+name)
        if [part.get('order',i+1) for i,part in enumerate(asset['parts'])]!=list(range(1,len(asset['parts'])+1)):raise ValueError('Transport part ordering is not sequential')
        if any(not part.get('url','').startswith('https://') for part in asset['parts']):raise ValueError('Manifest has no published HTTPS URL for every part')
    dest=a.output_dir.resolve()
    legacy_root=(root.parent/'research_max_struct_r12').resolve()
    if dest.is_relative_to(legacy_root) or a.extract_dir and a.extract_dir.resolve().is_relative_to(legacy_root):raise ValueError('R13 restoration never writes into R12 directories')
    if dest.exists() and any(dest.iterdir()):raise FileExistsError('Choose a new/empty external restore directory; existing files are not overwritten.')
    if a.extract_dir and a.extract_dir.exists():raise FileExistsError('--extract-dir must be a NEW directory')
    if not a.legacy_restorer.is_file():raise FileNotFoundError('Restore existing project R12 helper or supply --legacy-restorer; no dependency install is performed.')
    argv=[str(a.legacy_restorer),'--manifest',str(a.manifest.resolve()),'--output-dir',str(dest)]
    for name in a.asset:argv+=['--asset',name]
    previous=sys.argv
    try:sys.argv=argv;runpy.run_path(str(a.legacy_restorer),run_name='__main__')
    finally:sys.argv=previous
    if a.extract_dir:
        a.extract_dir.mkdir(parents=True,exist_ok=False);target=a.extract_dir.resolve()
        for name in a.asset:
            if not name.endswith('.zip'):raise ValueError('R13 extraction mode expects plainZIP; logical transport parts have already been concatenated.')
            asset_target=target/Path(name).stem;asset_target.mkdir(exist_ok=False)
            with zipfile.ZipFile(dest/name) as z:
                for item in z.infolist():
                    rel=PurePosixPath(item.filename.replace('\\','/'))
                    if rel.is_absolute() or '..' in rel.parts or stat.S_ISLNK(item.external_attr>>16):raise ValueError('Unsafe ZIP path/member '+item.filename)
                    if not(asset_target/Path(*rel.parts)).resolve().is_relative_to(asset_target):raise ValueError('ZIP path leaves new destination')
                for item in z.infolist():
                    rel=PurePosixPath(item.filename.replace('\\','/'));destination=asset_target/Path(*rel.parts)
                    if item.is_dir():destination.mkdir(parents=True,exist_ok=True);continue
                    destination.parent.mkdir(parents=True,exist_ok=True)
                    with z.open(item) as source,destination.open('xb') as output:shutil.copyfileobj(source,output,8*1024*1024)
        print(target,flush=True)
    print('Restored files only. Original absolute paths/mtime seals are not portable: create a new prediction seal at the restored root before scoring; parent authorizes it separately.',flush=True)
if __name__=='__main__':main()

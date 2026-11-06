"""Package this study and the already-public direct dependency archive; no science."""
from pathlib import Path
import argparse,zipfile,json
P=Path(__file__).resolve().parents[1];REPO=P.parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--dependency-zip',required=True);p.add_argument('--output',required=True);a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
    files=[f for f in sorted(P.rglob('*'))if f.is_file()and '__pycache__'not in f.parts and 'artifacts'not in f.parts and not f.name.endswith('.pyc')]
    with zipfile.ZipFile(out/'reading-pack.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6)as z:
        for f in files:
            if f.suffix!='.npz':z.write(f,str(f.relative_to(REPO)))
    with zipfile.ZipFile(a.dependency_zip)as dep,zipfile.ZipFile(out/'full-reproduction.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=1)as z:
        for name in dep.namelist():z.writestr(name,dep.read(name))
        for f in files:z.write(f,str(f.relative_to(REPO)))
        z.writestr('GN_GREEN_R2_PACKAGE_MANIFEST.json',json.dumps(dict(study='research/research_gn_search_space_r2',new_study_files=len(files),dependency='gn-green-error-geometry-20261106 published full-reproduction.zip',array_policy='All saved NPZ arrays included; GT already disclosed; no private seed',hash_policy='No hash checks requested'),indent=2))
    print(json.dumps({f.name:f.stat().st_size for f in out.glob('*.zip')},indent=2))
if __name__=='__main__':main()

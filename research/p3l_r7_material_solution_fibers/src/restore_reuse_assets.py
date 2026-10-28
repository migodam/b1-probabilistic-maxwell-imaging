"""Restore only declared frozen reuse assets, never overwrite differing bytes."""
from pathlib import Path
import argparse,hashlib,json,os,shutil
P=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--gaussian-root',type=Path,required=True);ap.add_argument('--apply',action='store_true');a=ap.parse_args()
    rows=json.loads((P/'ASSET_MANIFEST.json').read_text())['files'];restored=[]
    for row in rows:
        rel=Path(row['restore_relative_to_Gaussian'])
        if rel.is_absolute() or '..' in rel.parts:raise ValueError('unsafe asset path')
        source=P/row['package_path'];b=source.read_bytes()
        if hashlib.blake2b(b).hexdigest()!=row['blake2b']:raise ValueError('changed packaged input')
        dest=a.gaussian_root.resolve()/rel
        if dest.exists() and dest.read_bytes()!=b:raise FileExistsError('differing file; restore to an isolated Gaussian root: '+str(dest))
        if a.apply and not dest.exists():
            dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
            os.utime(dest,ns=(row['mtime_ns'],row['mtime_ns']))
            if dest.read_bytes()!=b:raise ValueError('asset byte readback')
        restored.append(str(dest))
    print(json.dumps(dict(mode='APPLY' if a.apply else 'PREVIEW_ONLY',files=len(restored),targets=restored),indent=2))
if __name__=='__main__':main()

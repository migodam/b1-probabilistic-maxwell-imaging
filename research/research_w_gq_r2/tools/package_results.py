"""One local archive; no Git mutation, network publication or hash validation."""
from pathlib import Path
import json,time,zipfile
from build_evidence import cost
ROOT=Path(__file__).resolve().parents[1]
def main():
    archive=ROOT/'Research_W_GQ_R2_Results.zip'
    if archive.exists():raise RuntimeError('archive already exists; do not silently replace')
    cost()
    files=[]
    for p in ROOT.rglob('*'):
        rel=p.relative_to(ROOT)
        if not p.is_file() or p==archive or any(x in rel.parts for x in ['__pycache__','.mpl']):continue
        if p.suffix in ['.pyc','.lock']:continue
        files.append(p)
    manifest=dict(created_unix=time.time(),archive=archive.name,source_root='Research_W_GQ_R2',files=len(files),uncompressed_bytes=sum(p.stat().st_size for p in files),includes_original_inputs=True,includes_exposed_old_bank=True,includes_locked=False,automatic_publication=False,hash_checks=False,accounting='staged snapshot; final completed packaging receipt/accounting appended under delivery/')
    (ROOT/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2));files.append(ROOT/'PACKAGE_MANIFEST.json')
    with zipfile.ZipFile(archive,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=3,allowZip64=True) as z:
        for p in sorted(files):
            compressed=p.suffix.lower() in ['.npz','.pdf','.png','.jpg','.jpeg']
            z.write(p,'Research_W_GQ_R2/'+p.relative_to(ROOT).as_posix(),compress_type=zipfile.ZIP_STORED if compressed else zipfile.ZIP_DEFLATED)
    with zipfile.ZipFile(archive) as z:
        names=z.namelist();assert 'Research_W_GQ_R2/REPORT_ZH.md' in names
        start=z.read('Research_W_GQ_R2/START_HERE.md').decode();assert '1.13137' in start
    receipt=dict(archive=archive.name,bytes=archive.stat().st_size,entries=len(names),ordinary_saved_report_readback=True,full_archive_hash_or_sweep=False,publication=False)
    (ROOT/'PACKAGE_RECEIPT.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt))
if __name__=='__main__':main()

"""Light reading distribution; full scientific arrays are separate release shards."""
from pathlib import Path
import zipfile,json,time,hashlib
R=Path(__file__).resolve().parents[3];P=R/'research/p3l_r9_stabilized_subspace';O=P/'packaging';A=R/'research/maxwell_formula_atlas';dst=O/'Research_P3LR9S_Reading_Kit.zip'
if dst.exists():raise RuntimeError('refuse overwrite')
files=[f for f in P.iterdir() if f.is_file() and f.suffix in ['.md','.html','.json','.csv']]
for d in ['src','config','formula_collection','figures','GPT_PRO_HANDOFF']:
 files.extend(f for f in (P/d).rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix not in ['.npz','.npy','.pyc'])
files.extend(f for f in A.rglob('*') if f.is_file())
for f in (P/'audit').rglob('*'):
 if f.is_file() and f.suffix in ['.md','.json','.csv','.txt'] and 'reading' not in f.parts and f.stat().st_size<2_000_000:files.append(f)
files=sorted(set(files));start=time.monotonic()
with zipfile.ZipFile(dst,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for f in files:z.write(f,str(f.relative_to(R)))
h=hashlib.blake2b(digest_size=32)
with dst.open('rb') as f:
 for b in iter(lambda:f.read(8<<20),b''):h.update(b)
rec=dict(file=dst.name,files=len(files),source_bytes=sum(f.stat().st_size for f in files),archive_bytes=dst.stat().st_size,blake2b_256=h.hexdigest(),wall_seconds=time.monotonic()-start,role='READING_KIT_NOT_FULL_RAW_EVIDENCE',physics_calls=0)
(O/(dst.name+'.receipt.json')).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))

"""Byte-only completed continuation harvest, with source seal verification."""
import time,json,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];out=ROOT/'RESULTS/vector_continuation_v04';checks=[]
for arm in ('sourceCIE','GlobalDual','AplusGlobal'):
 p=out/arm/'SEALED.json';seal=json.loads(p.read_text())
 for row in seal['files']:
  f=REPO/row['path'];s=f.stat()
  if s.st_size!=row['bytes'] or s.st_mtime_ns!=row['mtime_ns']:raise RuntimeError('Source seal changed')
  checks.append(row)
receipt=dict(status='ALL3_SOURCE_SEALS_VERIFIED',file_count=len(checks),GT_read=False)
(out/'SOURCE_SEAL_VERIFICATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
paths=list(out.rglob('*'))+list((ROOT/'delegated/vector/windows_v04_continuation_queue').glob('*.json'))
archive=REPO/'VECTOR_CONTINUATION_FINAL.zip'
with zipfile.ZipFile(archive,'w',zipfile.ZIP_STORED) as z:
 for p in paths:
  if p.is_file():z.write(p,p.relative_to(REPO).as_posix())
 p=REPO/'CONTINUATION_EXTRACTION_COST.json'
 if p.exists():z.write(p,(out/'REMOTE_EXTRACTION_COST.json').relative_to(REPO).as_posix())
cost=dict(platform='XINAN',CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,scope='Own source verification and byte packaging only; no scientific solves or GT',physical_RHS=0)
p=out/'HARVEST_COST.json';p.write_text(json.dumps(cost,indent=2)+'\n')
with zipfile.ZipFile(archive,'a',zipfile.ZIP_STORED) as z:z.write(p,p.relative_to(REPO).as_posix())
print(json.dumps(dict(**cost,archive_bytes=archive.stat().st_size)))

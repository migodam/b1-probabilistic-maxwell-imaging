"""Lossless final evidence packaging; archives are never physics receipts."""
from pathlib import Path
import argparse,csv,json,hashlib,subprocess,time
R=Path(__file__).resolve().parents[3];P=R/'research/p3l_r9_stabilized_subspace';O=P/'packaging';INV=R.parents[2]
def pack(name,root,rels,role):
 out=O/name
 if out.exists():raise RuntimeError('Archive already exists: '+str(out))
 rels=sorted(set(rels));listing=O/(name+'.files.txt');listing.write_text('\n'.join(rels)+'\n');start=time.monotonic()
 with (O/(name+'.tar_stderr')).open('wb') as err:
  a=subprocess.Popen(['/usr/bin/tar','-cf','-','-T',str(listing)],cwd=root,stdout=subprocess.PIPE,stderr=err)
  b=subprocess.run(['/opt/homebrew/bin/zstd','-T2','-3','-q','-o',str(out)],stdin=a.stdout,capture_output=True);a.stdout.close();code=a.wait()
 if code or b.returncode:raise RuntimeError((code,b.returncode,b.stderr.decode()))
 h=hashlib.blake2b(digest_size=32)
 with out.open('rb') as f:
  for block in iter(lambda:f.read(8<<20),b''):h.update(block)
 rec=dict(file=name,files=len(rels),source_bytes=sum((root/q).stat().st_size for q in rels),archive_bytes=out.stat().st_size,blake2b_256=h.hexdigest(),wall_seconds=time.monotonic()-start,role=role,physics_calls=0)
 (O/(name+'.receipt.json')).write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec))
def dependencies():
 rows=list(csv.DictReader((P/'audit/reproducibility/EXTERNAL_DEPENDENCY_INVENTORY.csv').open()))
 for row in rows:
  assert (INV/row['archive_path']).is_file(),row['archive_path']
 pack('R9_historical_dependencies.tar.zst',INV,[x['archive_path'] for x in rows],'GAUSSIAN_RELATIVE_HISTORICAL_REUSE_NOT_NEW_PHYSICS')
def metadata():
 archived=set()
 for listing in O.glob('R9_fresh_*.tar.zst.files.txt'):archived.update(listing.read_text().splitlines())
 archived.update((O/'R9_development_and_probability_arrays.tar.zst.files.txt').read_text().splitlines())
 files=[]
 for f in P.rglob('*'):
  if not f.is_file():continue
  rel=f.relative_to(P)
  if any(x in rel.parts for x in ['packaging','inputs','__pycache__']):continue
  if rel.parts[:3]==('audit','theory','reading'):continue
  if f.suffix in ['.lock','.pyc']:continue
  rr=str(f.relative_to(R))
  if rr not in archived:files.append(rr)
 for f in (R/'research/maxwell_formula_atlas').rglob('*'):
  if f.is_file():files.append(str(f.relative_to(R)))
 pack('R9_complete_metadata_scoring_and_sources.tar.zst',R,files,'REPO_RELATIVE_METADATA_WITH_OFFLINE_ONLY_SYNTHETIC_TRUTH')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('part',choices=['dependencies','metadata']);v=a.parse_args();globals()[v.part]()

"""Lossless packaging only: no solves, statistics, input mutation, or deletion."""
from pathlib import Path
import subprocess,json,time,hashlib
R=Path(__file__).resolve().parents[3];P=R/'research/p3l_r9_stabilized_subspace';V=P/'results/possibility_volume';O=P/'packaging';O.mkdir(exist_ok=True)
assert (V/'offline_statistics/SCORING_COMPLETE.json').exists()
groups={}
for start in [0,8,16,24]:
 fs=[]
 for i in range(start,start+8):
  name=f'R9V_FRESH_{i:03d}'
  for d in [V/'parents'/name,V/'jobs'/name]:fs.extend(x for x in d.rglob('*') if x.is_file())
 groups[f'R9_fresh_{start:03d}_{start+7:03d}.tar.zst']=fs
fs=[]
for branch in ['exploration','exploration_rev1','stable_roi','probability']:
 d=P/'results'/branch
 fs.extend(x for x in d.rglob('*') if x.is_file() and x.suffix in ['.npz','.npy'])
fs.extend(x for x in (P/'audit').rglob('*.npz'))
groups['R9_development_and_probability_arrays.tar.zst']=fs
receipts=[]
for name,fs in groups.items():
 dest=O/name
 if dest.exists():raise RuntimeError('refuse overwrite '+str(dest))
 rels=sorted(set(str(x.relative_to(R)) for x in fs));fl=O/(name+'.files.txt');fl.write_text('\n'.join(rels)+'\n')
 t=time.monotonic();t0=time.process_time()
 with (O/(name+'.tar_stderr')).open('wb') as err:
  a=subprocess.Popen(['/usr/bin/tar','-cf','-','-T',str(fl)],cwd=R,stdout=subprocess.PIPE,stderr=err)
  b=subprocess.run(['/opt/homebrew/bin/zstd','-T2','-3','-q','-o',str(dest)],stdin=a.stdout,capture_output=True);a.stdout.close();ar=a.wait()
 if ar or b.returncode:raise RuntimeError((ar,b.returncode,b.stderr.decode()))
 h=hashlib.blake2b(digest_size=32)
 with dest.open('rb') as f:
  for block in iter(lambda:f.read(8<<20),b''):h.update(block)
 rec={'file':name,'files':len(rels),'source_bytes':sum((R/f).stat().st_size for f in rels),'archive_bytes':dest.stat().st_size,'blake2b_256':h.hexdigest(),'wall_seconds':time.monotonic()-t,'role':'LOSSLESS_ARCHIVE_NO_NEW_PHYSICS'}
 receipts.append(rec);(O/'RAW_SHARDS.json').write_text(json.dumps(receipts,indent=2)+'\n');print(json.dumps(rec),flush=True)

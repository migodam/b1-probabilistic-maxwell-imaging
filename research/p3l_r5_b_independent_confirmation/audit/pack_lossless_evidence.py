"""Local immutable evidence delivery. No upload, deletion, hash or science."""
from pathlib import Path
import os,json,tarfile,subprocess,time,resource,mimetypes,collections
OUT=Path(__file__).resolve().parent
BASE=OUT.parents[2]/'Gaussian'
REPO=BASE/'B1/b1-p3l-r5-b-confirm'
S=REPO/'research/p3l_r5_b_independent_confirmation'
ZSTD='/opt/homebrew/bin/zstd'
LEDGER=OUT/'DELIVERY_LEDGER.jsonl'
def event(**kw):
 payload=dict(utc=time.time(),new_physics_RHS=0,new_training=0);payload.update(kw)
 with LEDGER.open('a') as f:f.write(json.dumps(payload)+'\n');f.flush();os.fsync(f.fileno())
def files(p):
 return [x for x in sorted(p.rglob('*')) if x.is_file() and not x.is_symlink() and '__pycache__' not in x.parts and x.suffix!='.pyc'] if p.is_dir() else [p]
def cpu():
 a=resource.getrusage(resource.RUSAGE_SELF);b=resource.getrusage(resource.RUSAGE_CHILDREN)
 return a.ru_utime+a.ru_stime+b.ru_utime+b.ru_stime
def pack(name,paths):
 paths=sorted(set(paths));dest=OUT/(name+'.tar.zst');partial=OUT/(name+'.tar.zst.partial');receipt=OUT/(name+'.json')
 if dest.exists() and receipt.exists() and not partial.exists():
  r=json.loads(receipt.read_text());assert r['status']=='FULL_BYTE_ROUNDTRIP_PASS' and dest.stat().st_size==r['archive_bytes']
  for m in r['members']:
   st=Path(m['source']).stat();assert (st.st_size,st.st_mtime_ns)==(m['bytes'],m['mtime_ns'])
  event(event='RESUME_VERIFIED_EXISTING_ARCHIVE',name=name,files=r['files'],archive_bytes=r['archive_bytes'],reason='Original full byte comparison finished; bookkeeping exception corrected without repacking or source changes')
  return {k:v for k,v in r.items() if k not in ['members','mime_counts']}
 if dest.exists() or partial.exists() or receipt.exists():raise FileExistsError('No overwrite: '+name)
 start=time.perf_counter();c0=cpu();manifest=[];mime=collections.Counter()
 for p in paths:
  if not p.is_file() or p.is_symlink():raise ValueError('Regular source required '+str(p))
  arc=p.relative_to(BASE).as_posix()
  if '..' in Path(arc).parts or Path(arc).is_absolute():raise ValueError('Unsafe member '+arc)
  st=p.stat();kind=mimetypes.guess_type(p.name)[0] or 'application/octet-stream';mime[kind]+=1
  manifest.append(dict(path=arc,source=str(p),bytes=st.st_size,mtime_ns=st.st_mtime_ns,mode=st.st_mode&0o777,mime_by_extension=kind))
 event(event='ARCHIVE_BEGIN',name=name,files=len(paths),source_bytes=sum(x['bytes'] for x in manifest))
 try:
  with partial.open('xb') as output:
   proc=subprocess.Popen([ZSTD,'-3','-T1','-q','-c'],stdin=subprocess.PIPE,stdout=output)
   with tarfile.open(fileobj=proc.stdin,mode='w|',format=tarfile.PAX_FORMAT) as tar:
    for p,m in zip(paths,manifest):
     st=p.stat()
     if (st.st_size,st.st_mtime_ns)!=(m['bytes'],m['mtime_ns']):raise RuntimeError('Source changed before pack '+str(p))
     tar.add(p,arcname=m['path'],recursive=False)
     st=p.stat()
     if (st.st_size,st.st_mtime_ns)!=(m['bytes'],m['mtime_ns']):raise RuntimeError('Source changed during pack '+str(p))
   proc.stdin.close()
   if proc.wait()!=0:raise RuntimeError('zstd compression failed')
   output.flush();os.fsync(output.fileno())
  if partial.stat().st_size>=2_000_000_000:raise RuntimeError('Archive exceeds2GB cap')
  event(event='ROUNDTRIP_BEGIN',name=name,archive_bytes=partial.stat().st_size)
  proc=subprocess.Popen([ZSTD,'-q','-d','-c',str(partial)],stdout=subprocess.PIPE)
  expected={m['path']:m for m in manifest};seen=set();compared=0
  with tarfile.open(fileobj=proc.stdout,mode='r|') as tar:
   for member in tar:
    if not member.isfile() or member.name not in expected or member.name in seen or '..' in Path(member.name).parts or Path(member.name).is_absolute():raise RuntimeError('Unexpected unsafe member '+member.name)
    m=expected[member.name]
    if member.size!=m['bytes']:raise RuntimeError('Size mismatch '+member.name)
    with tar.extractfile(member) as extracted,Path(m['source']).open('rb') as original:
     while True:
      block=extracted.read(2*1024*1024)
      if block!=original.read(len(block)):raise RuntimeError('Byte mismatch '+member.name)
      compared+=len(block)
      if not block:
       if original.read(1):raise RuntimeError('Source grew '+member.name)
       break
    st=Path(m['source']).stat()
    if (st.st_size,st.st_mtime_ns)!=(m['bytes'],m['mtime_ns']):raise RuntimeError('Source mutated '+member.name)
    seen.add(member.name)
  while proc.stdout.read(1024*1024):pass
  if proc.wait()!=0 or seen!=set(expected):raise RuntimeError('Incomplete decompression/member set')
  partial.rename(dest)
  result=dict(status='FULL_BYTE_ROUNDTRIP_PASS',archive=str(dest),files=len(paths),source_bytes=sum(m['bytes'] for m in manifest),archive_bytes=dest.stat().st_size,roundtrip_compared_bytes=compared,wall_seconds=time.perf_counter()-start,CPU_seconds_self_plus_children=cpu()-c0,compression='tar PAX +zstd3 single thread; decompressed members compared directly to originals',paths_relative_to='Gaussian',mime_scope='extension classification, not semantic content validation',mime_counts=dict(mime),new_physics_RHS=0,new_training=0,source_deletion=False,hashes_computed=False,SHA256=False,members=manifest)
  receipt.write_text(json.dumps(result,indent=2)+'\n');event(event='ARCHIVE_COMPLETE',**{k:v for k,v in result.items() if k not in ['members','mime_counts']})
  print(json.dumps({k:v for k,v in result.items() if k not in ['members','mime_counts']}),flush=True)
  return {k:v for k,v in result.items() if k not in ['members','mime_counts']}
 except BaseException as e:
  event(event='FAILURE_PARTIAL_RETAINED',name=name,error=repr(e),wall_seconds=time.perf_counter()-start,CPU_seconds=cpu()-c0);raise
def main():
 start=time.perf_counter();c0=cpu();groups={'test-a-raw':[],'test-b-raw':[],'shared-probe-scores-ledger':[]};raw=[]
 for p in sorted((S/'raw').rglob('*')):
  if p.is_symlink():raise ValueError('Raw symlink requires explicit handling '+str(p))
  if not p.is_file():continue
  rel=p.relative_to(S/'raw');r=str(rel);key='test-a-raw' if rel.parts[:2]==('predictions','A') or '_A_256_' in r else ('test-b-raw' if rel.parts[:2]==('predictions','B') or '_B_256_' in r else 'shared-probe-scores-ledger');groups[key].append(p);raw.append(p)
 groups['shared-probe-scores-ledger'] += [S/'RESOURCE_LEDGER.jsonl',S/'RESOURCE_LEDGER_SNAPSHOT.json']
 checkpoints=json.loads((S/'audit/model_replay/CHECKPOINT_IDENTITIES_BLAKE2B.json').read_text())['checkpoints'];assert len(checkpoints)==16
 hist=[Path(r['path']) for r in checkpoints]
 R4=BASE/'B1/b1-p3l-r4-conditional-transport';R5=BASE/'B1/b1-p3l-r5-parallel-ab/research/p3l_r5_parallel_ab';A30=BASE/'A30/a30-maxwell-material-manifold'
 hist += files(R4/'src')+[R4/'research/p3l_r4_conditional_transport/configs/TARGET.json',R4/'research/p3l_r4_conditional_transport/configs/PROTOCOL.json',R4/'research/p3l_r4_conditional_transport/inputs/background/weak.npz',R4/'research/p1_material_target_20261018/inputs/geometry.json']
 hist += files(A30/'src')+files(A30/'vendor/frozen_a29/src')
 hist += [p for p in files(R5/'stage_b_material_risk') if p.suffix=='.py']+[R5/'stage_b_material_risk/ASSET_INVENTORY.json']
 hist += files(S/'src')+files(S/'configs')+[S/'audit/model_replay/CHECKPOINT_IDENTITIES_BLAKE2B.json']
 # Frozen target metadata and exact runtime feature compatibility clarification.
 for p in [R5/'shared/P1_TARGET_FROZEN.json',R4/'research/p3l_r4_conditional_transport/AUDIT/fusion/providers/CLOSED_K_METADATA_CLARIFICATION.json']:
  if p.is_file():hist.append(p)
 # Packaging dependencies are metadata only; do not ship virtual environments.
 for base in [R4,A30,R5.parents[1]]:
  for name in ['pyproject.toml','requirements.txt','requirements-nn.txt','uv.lock']:
   p=base/name
   if p.is_file() and p.stat().st_size<2_000_000:hist.append(p)
 groups['essential-reused-source-checkpoints']=hist
 flattened=[p for paths in groups.values() for p in paths]
 if len(set(raw))!=len(raw):raise ValueError('Raw duplicates')
 planned=dict(groups={k:dict(files=len(set(v)),source_bytes=sum(p.stat().st_size for p in set(v))) for k,v in groups.items()},all_raw_files=len(raw),all_raw_bytes=sum(p.stat().st_size for p in raw),checkpoint_count=16,archive_base=str(BASE),no_upload=True,no_source_deletion=True)
 (OUT/'PLAN.json').write_text(json.dumps(planned,indent=2)+'\n');event(event='DELIVERY_BEGIN',**planned)
 results=[pack(k,v) for k,v in groups.items()]
 receipt=dict(status='ALL_LOSSLESS_BYTE_ROUNDTRIPS_PASS',archives=results,raw_coverage=planned,wall_seconds=time.perf_counter()-start,CPU_seconds_self_plus_children=cpu()-c0,new_physics_RHS=0,new_training=0,no_upload=True,originals_preserved=True,exclusions=['Other historical TRAIN32768 or other giant historical data arrays','Historical campaign ledgers except new independent confirmation ledger','Virtual environments and installed dependency binaries','In-progress report renders/reading ZIP; parent creates later','Source __pycache__ and pyc excluded; all raw files included without such exclusions'],limitations=['Archive+source byte equality proves local lossless packaging only, not scientific acceptance','No source ownership/license inference','No external upload or anonymous download verification'],paths_relative_to='Gaussian')
 (OUT/'DELIVERY_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');event(event='DELIVERY_COMPLETE',status=receipt['status'],wall_seconds=receipt['wall_seconds'],CPU_seconds=receipt['CPU_seconds_self_plus_children']);print(json.dumps(receipt,indent=2),flush=True)
if __name__=='__main__':main()

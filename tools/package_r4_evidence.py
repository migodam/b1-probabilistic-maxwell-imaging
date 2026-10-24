#!/usr/bin/env python3
"""R4 standard-library streaming evidence inventory and confirmed split packager.

Default inventory is metadata-only, one os.scandir traversal per selected tree.
Full packing NEVER runs without a separate parent confirmation JSON. No hashes,
physics/training/inference, GT analysis, external symlink following, or upload.
Every invocation writes a NEW directory below OUT/packaging; prior runs survive.
"""
from pathlib import Path
from contextlib import ExitStack
from collections import defaultdict
import argparse,json,os,stat,sys,time,subprocess,fcntl,gzip,tarfile,shutil,math,resource,heapq
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport';PACK=OUT/'packaging'
A30=ROOT.parents[1]/'A30/a30-maxwell-material-manifold';R2=ROOT.parent/'b1-p3l-r2-risk-posterior'
LIMIT=1_500_000_000;VERSION='R4_STREAM_EVIDENCE_V1';MIB=1024**2

def dump(path,value):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():raise FileExistsError(path)
 temp=Path(str(path)+'.writing')
 with temp.open('x') as f:json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
 os.link(temp,path);temp.unlink()
 fd=os.open(path.parent,os.O_RDONLY)
 try:os.fsync(fd)
 finally:os.close(fd)

def identity(path):
 p=Path(path);s=p.lstat();return dict(path=str(p.absolute()),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
def git(repo,*args):
 p=subprocess.run(['git','-C',str(repo),*args],capture_output=True,text=True)
 return p.stdout.strip() if p.returncode==0 else 'UNAVAILABLE: '+p.stderr.strip()
def provenance():
 return dict(format=VERSION,workspace=str(ROOT),workspace_commit=git(ROOT,'rev-parse','HEAD'),runtime_source='actual working-tree bytes, including uncommitted R4 implementation; commit is provenance, not clean-tree claim',
  pinned_R2_declared_commit='6eadb4690c962aaa95821fc4056b677f0f924a6b',R2_worktree_commit=git(R2,'rev-parse','HEAD'),A30_declared_pin='717d3d26a290248bd2b9aeb55f3d03694b407b84',A30_current_commit=git(A30,'rev-parse','HEAD'),
  A30_backend_vendored_origin='vendor/frozen_a29 including nested vendor; A30 SOURCE_COMMIT5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4',
  bundled_dependencies=['workspace/src/p1_material (actual unchanged P1 runtime)','workspace/src/p3r2 and p3r2_mechanism (actual baseline runtime)','dependencies/a30/src plus entire vendor tree and provenance/config/requirements','dependencies/r2_baseline frozen n512 local_geometry_cnn_s33 directory','workspace/research/p1_material_target_20261018/inputs/geometry.json'],
  external_dependencies=['Python/NumPy/SciPy/PyTorch/MPS runtimes and OS are not bundled; use archived ENVIRONMENT and requirements','the entire historical R2 cache is not recursively copied; all real R4-local inputs are included, external symlinks are provenance-only','frozen manifests contain original absolute paths/mtime identities; restore original mount layout or document an explicit relocation mapping, never silently rewrite closed evidence'],
  hashes_used=False,external_symlinks_followed=False)
def roots():
 selected=[(OUT,'workspace/'+str(OUT.relative_to(ROOT)),'R4'),(ROOT/'src','workspace/src','runtime_source'),(ROOT/'tools','workspace/tools','runtime_source'),(ROOT/'tests','workspace/tests','runtime_source'),(A30/'src','dependencies/a30/src','pinned_A30_source'),(A30/'vendor','dependencies/a30/vendor','pinned_A30_vendor'),
 (R2/'research/p3l_r2_risk_posterior_exploration/A/residual_training/n512/local_geometry_cnn_s33','dependencies/r2_baseline/local_geometry_cnn_s33','reused_R2_baseline'),(ROOT/'research/p1_material_target_20261018/inputs/geometry.json','workspace/research/p1_material_target_20261018/inputs/geometry.json','reused_P1_geometry')]
 for base,prefix in [(ROOT,'workspace'),(A30,'dependencies/a30')]:
  with os.scandir(base) as entries:
   for e in entries:
    if e.name.startswith('.') and e.name!='.gitignore':continue
    if e.is_file(follow_symlinks=False) and (e.name.endswith(('.md','.txt','.toml','.yaml','.yml','.html','.css','.js','.svg','.pdf','.png','.jpg','.webp','.json')) or e.name in ['LICENSE','a30_config.json']):selected.append((Path(e.path),prefix+'/'+e.name,'source_metadata'))
 if PACK.exists():
  with os.scandir(PACK) as entries:
   for e in entries:
    if e.is_file(follow_symlinks=False) and e.name.endswith(('.md','.json','.py','.txt')):selected.append((Path(e.path),'workspace/'+str(Path(e.path).relative_to(ROOT)),'packaging_metadata'))
 # Roots are deliberately disjoint: no duplicate traversal of the same file.
 existing=[];missing=[]
 for p,arc,origin in selected:
  if not os.path.lexists(p):missing.append(dict(source=str(p),archive_path=arc));continue
  existing.append((p,arc,origin))
 return existing,missing

def category(source,origin):
 if origin!='R4':return 'compact_core',origin
 r=source.relative_to(OUT);p=r.parts;name=source.name
 if len(p)>1 and p[0]=='inputs':
  if p[1] in ['joint','background']:return 'reused_inputs','reused_historical_R2_or_paid_background; not new independent R4 parents'
  if p[1] in ['online','offline_only']:
   if name.startswith('R2_'):return 'reused_inputs','reused original R2 TRAIN/VAL observation/label or metadata export; not new R4 independent parents'
   if name.startswith(('R4DEV_','R4L1_','R4L2_')):return 'new_R4_inputs','new R4 simulated parents/observations/truth; exact seeds in original manifests'
   return 'derived_training_inputs','observation/truth export without new R4 cohort prefix; use original manifest provenance, not new-independent-parent claim'
  return 'derived_training_inputs','new_R4 features/contexts, often on reused historical R2 TRAIN/VAL parents; not new independent parent claim'
 if source.suffix in ['.npz','.npy','.bin','.dat']:
  if p and p[0]=='results' and any(s in p for s in ['risk45','density','samples','features','online']):return 'raw_neural_and_scores','new_R4 neural samples/actions/scores; category is storage only, not posterior qualification'
  return 'raw_physics_and_arrays','new_R4 paid physical/cache/attempt arrays and numerical arrays; exact experiment provenance retained'
 return 'compact_core','R4 reports/models/checkpoints/logs/ledgers/source snapshots/figures; includes failures and unqualified results'

def enumerate_once(destination,deadline=None):
 """One scandir pass. No data arrays/checkpoints are opened for inventory."""
 start=time.perf_counter();selected,missing=roots();handles={};stats=defaultdict(lambda:dict(files=0,logical_bytes=0,stored_content_bytes=0,allocated_bytes=0,hardlink_aliases=0,directories=0));links={};dircount=0;special=[];largest=[];seen_roots=[];entries=0
 destination.mkdir(parents=True,exist_ok=False);(destination/'entries').mkdir();dep=provenance();dump(destination/'DEPENDENCIES.json',dep)
 def writer(cat):
  if cat not in handles:handles[cat]=(destination/'entries'/(cat+'.jsonl')).open('x')
  return handles[cat]
 def emit(cat,row):writer(cat).write(json.dumps(row,separators=(',',':'),ensure_ascii=False)+'\n')
 def walk(path,arc,origin):
  nonlocal entries,dircount
  if deadline:deadline.check()
  info=path.lstat();entries+=1;mode=stat.S_IMODE(info.st_mode);row=dict(source=str(path.absolute()),archive_path=arc,size=info.st_size,mtime_ns=info.st_mtime_ns,mode=mode,origin=origin)
  if stat.S_ISLNK(info.st_mode):
   special.append(dict(row,type='symlink',archived=False,link_target=os.readlink(path),reason='provenance only; never followed/restored'));return
  if path==PACK or path.name in ['.git','__pycache__','.venv','node_modules'] or path.name=='.DS_Store' or path.suffix=='.pyc':
   special.append(dict(row,type='excluded',archived=False,reason='packaging self/output or Git/runtime housekeeping; no scientific arrays excluded by extension/ignore rules'));return
  if stat.S_ISDIR(info.st_mode):
   row.update(type='directory',category='compact_core',provenance='directory metadata');emit('compact_core',row);dircount+=1;stats['compact_core']['directories']+=1
   with os.scandir(path) as it:
    for child in it:walk(Path(child.path),arc+'/'+child.name,origin)
  elif stat.S_ISREG(info.st_mode):
   cat,history=category(path,origin);requested=cat;row.update(type='file',provenance=history,allocated_bytes=getattr(info,'st_blocks',0)*512)
   key=(info.st_dev,info.st_ino)
   if info.st_nlink>1 and key in links:
    target,cat=links[key];row.update(hardlink_to=target,requested_category=requested);stats[cat]['hardlink_aliases']+=1;content=0
   else:
    if info.st_nlink>1:links[key]=(arc,cat)
    content=info.st_size
   row['category']=cat;emit(cat,row);z=stats[cat];z['files']+=1;z['logical_bytes']+=info.st_size;z['stored_content_bytes']+=content;z['allocated_bytes']+=row['allocated_bytes']
   item=(info.st_size,arc)
   if len(largest)<25:heapq.heappush(largest,item)
   elif item>largest[0]:heapq.heapreplace(largest,item)
  else:special.append(dict(row,type='special',archived=False,reason='unsupported nonregular OS object; not silently archived'))
 try:
  for path,arc,origin in selected:seen_roots.append(dict(source=str(path),archive_path=arc,origin=origin));walk(path,arc,origin)
 finally:
  for h in handles.values():h.close()
 dump(destination/'EXCLUSIONS_AND_SYMLINKS.json',special)
 metadata_bytes=sum(p.stat().st_size for p in (destination/'entries').iterdir())+(destination/'DEPENDENCIES.json').stat().st_size+(destination/'EXCLUSIONS_AND_SYMLINKS.json').stat().st_size
 usage=shutil.disk_usage(PACK);content=sum(v['stored_content_bytes'] for v in stats.values());files=sum(v['files'] for v in stats.values());logical=sum(v['logical_bytes'] for v in stats.values());tar_overhead=(files+dircount)*3072+10*MIB;upper=content+tar_overhead+metadata_bytes+max(64*MIB,int(content*.005));reserve=max(5*1024**3,int(upper*.10));elapsed=time.perf_counter()-start
 summary=dict(format=VERSION,status='METADATA_INVENTORY_ONLY_NO_ARCHIVE',created_unix=time.time(),scan_wall_seconds=elapsed,selected_roots=seen_roots,missing_roots=missing,categories=dict(stats),logical_file_bytes=logical,unique_hardlink_content_bytes=content,files=files,directories=dircount,enumerated_entries=entries,symlinks=sum(r['type']=='symlink' for r in special),excluded_objects=len(special),largest_files=[dict(bytes=n,archive_path=p) for n,p in sorted(largest,reverse=True)],
  journal_files={cat:'entries/'+cat+'.jsonl' for cat in handles},part_limit_bytes=LIMIT,estimated_uncompressed_tar_upper_bytes=upper,inventory_metadata_bytes=metadata_bytes,estimated_max_parts_no_compression=sum(math.ceil((v['stored_content_bytes']+(v['files']+v['directories'])*3072+(metadata_bytes if cat=='compact_core' else 0)+max(64*MIB,int(v['stored_content_bytes']*.005)))/LIMIT) for cat,v in stats.items()),disk=dict(total_bytes=usage.total,free_bytes=usage.free,required_bundle_upper_plus_reserve=upper+reserve,sufficient_for_default_local_bundle=usage.free>=upper+reserve),
  time_estimate=dict(assumption='not a benchmark; raw tar mostly I/O-bound, small-file overhead and core gzip add uncertainty',at_200_MiB_per_second_minutes=content/(200*MIB)/60,at_50_MiB_per_second_minutes=content/(50*MIB)/60,extra_metadata_overhead_minutes_range=[max(2,elapsed/60*2),max(5,elapsed/60*6)]),
  deduplication='disjoint scopes, each path one category; hardlinked file content once with aliases in same category; different historical file versions are not content-hash-deduplicated',inventory_consistency='provisional metadata snapshot while final reports may be produced; confirmed pack performs ONE fresh final scan and holds quiet locks; no repeated giant du',hashes_used=False,archive_started=False,upload_started=False)
 dump(destination/'SUMMARY.json',summary);return summary

def rows(path):
 with path.open() as f:
  for line in f:yield json.loads(line)
def verify_row(row):
 p=Path(row['source']);s=p.lstat()
 if stat.S_ISLNK(s.st_mode) or (s.st_size,s.st_mtime_ns,stat.S_IMODE(s.st_mode))!=(row['size'],row['mtime_ns'],row['mode']):raise RuntimeError('source changed: '+str(p))
 return p
class Deadline:
 def __init__(self,seconds):self.start=time.perf_counter();self.end=self.start+seconds
 def check(self):
  if time.perf_counter()>self.end:raise TimeoutError('delivery max-wall exceeded; preserve completed parts and in-progress .writing file')
def event(destination,name,row):
 row=dict(unix_ns=time.time_ns(),**row)
 with (destination/name).open('a') as f:f.write(json.dumps(row,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno())
 print(json.dumps(row),flush=True)
 return row
class SplitWriter:
 def __init__(self,dest,prefix,limit,*,category=None,codec=None,deadline=None):
  self.dest=dest;self.prefix=prefix;self.limit=limit;self.file=None;self.current=0;self.total=0;self.parts=[];self.category=category;self.codec=codec;self.deadline=deadline
 def write(self,data):
  if self.deadline:self.deadline.check()
  view=memoryview(data);length=len(view)
  while view:
   if self.deadline:self.deadline.check()
   if self.file is None:
    self.final=self.dest/(self.prefix+f'.part{len(self.parts)+1:04d}');self.pending=Path(str(self.final)+'.writing');self.file=self.pending.open('xb');self.current=0;self.parts.append(dict(file=self.final.name,bytes=0,status='WRITING_NOT_UPLOADABLE'))
   n=min(len(view),self.limit-self.current);self.file.write(view[:n]);self.current+=n;self.total+=n;self.parts[-1]['bytes']=self.current;view=view[n:]
   if self.current==self.limit:self.close_part()
  return length
 def flush(self):
  if self.file:self.file.flush()
 def close_part(self):
  if self.file:
   self.file.flush();os.fsync(self.file.fileno());self.file.close();self.file=None
   os.rename(self.pending,self.final)
   # Durably expose the immutable final filename before publishing its marker.
   fd=os.open(self.dest,os.O_RDONLY)
   try:os.fsync(fd)
   finally:os.close(fd)
   row=self.parts[-1];row.update(status='CLOSED_IMMUTABLE',closed_unix_ns=time.time_ns())
   if self.category:
    # Main uploader consumes only atomically published closed_parts JSON.
    marker=dict(category=self.category,codec=self.codec,absolute_path=str(self.final),**row)
    dump(self.dest/'closed_parts'/(self.final.name+'.json'),marker)
    event(self.dest,'PARTS_CLOSED.jsonl',dict(event='PART_CLOSED',**marker))
 def abort(self):
  if self.file:
   self.file.flush();os.fsync(self.file.fileno());self.file.close();self.file=None
   self.parts[-1].update(status='INCOMPLETE_NOT_UPLOADABLE',incomplete_file=self.pending.name)

def add(archive,row):
 p=verify_row(row);t=tarfile.TarInfo(row['archive_path']);t.mode=row['mode'];t.uid=t.gid=0;t.uname=t.gname='';t.mtime=row['mtime_ns']//10**9;t.pax_headers={'mtime':str(t.mtime)+'.'+str(row['mtime_ns']%10**9).zfill(9)}
 if row['type']=='directory':t.type=tarfile.DIRTYPE;archive.addfile(t)
 elif row.get('hardlink_to'):t.type=tarfile.LNKTYPE;t.linkname=row['hardlink_to'];archive.addfile(t)
 else:
  t.size=row['size'];fd=os.open(p,os.O_RDONLY|os.O_NOFOLLOW)
  with os.fdopen(fd,'rb') as f:
   s=os.fstat(f.fileno())
   if (s.st_size,s.st_mtime_ns)!=(row['size'],row['mtime_ns']):raise RuntimeError('source changed at open')
   archive.addfile(t,f)
  verify_row(row)
 archive.members.clear()
def last_ledger():
 p=OUT/'RESOURCE_LEDGER.jsonl'
 with p.open('rb') as f:f.seek(0,2);n=f.tell();f.seek(max(0,n-2*MIB));lines=f.read().splitlines()
 return json.loads(lines[-1])
def hold_quiet(stack):
 for name in ['physics.lock','training.lock','training_cpu.lock','budget.lock']:
  p=OUT/name
  if p.is_symlink():raise ValueError('symlink lock refused')
  h=stack.enter_context(p.open('a+'))
  try:fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
  except BlockingIOError:raise RuntimeError('campaign lock busy: '+name)
 row=last_ledger()
 if row['after']['active_wall']:raise RuntimeError('active campaign jobs: '+repr(row['after']['active_wall']))
 return row

def pack(a,dest):
 confirmation=json.loads(a.confirmation.read_text())
 if confirmation.get('released') is not True or confirmation.get('experiment')!='r4_evidence_stream_v1' or confirmation.get('final_reports_complete') is not True or confirmation.get('no_active_writers') is not True:raise ValueError('main explicit final-report/quiet packaging confirmation required')
 if confirmation.get('workspace')!=str(ROOT):raise ValueError('confirmation workspace mismatch')
 if confirmation.get('source_commit') and git(ROOT,'rev-parse','HEAD')!=confirmation['source_commit']:raise ValueError('confirmed final source commit differs')
 started=time.perf_counter();cpu=time.process_time();deadline=Deadline(a.max_wall);receipt=dict(format=VERSION,status='FAILED',parts=[],started_unix=time.time(),max_wall_seconds=a.max_wall,part_limit_bytes=a.part_bytes,hashes_used=False,upload_started=False)
 try:
  with ExitStack() as stack:
   ledger=hold_quiet(stack);dest.mkdir(parents=True,exist_ok=False)
   receipt['prepack_final_research_ledger_sequence']=ledger['sequence']
   try:
    summary=enumerate_once(dest/'inventory',deadline)
    if summary['missing_roots']:raise RuntimeError('required source/dependency roots missing')
    if not summary['disk']['sufficient_for_default_local_bundle']:raise RuntimeError('insufficient free disk for conservative raw-tar upper bound plus reserve')
    for cat,journal in summary['journal_files'].items():
     deadline.check();compressed=cat=='compact_core' or a.raw_compression=='gzip';prefix=cat+('.tar.gz' if compressed else '.tar');codec_name='gzip1' if compressed else 'tar';writer=SplitWriter(dest,prefix,a.part_bytes,category=cat,codec=codec_name,deadline=deadline);cat_start=time.perf_counter();done=content=0;last_progress=cat_start
     event(dest,'CATEGORY_EVENTS.jsonl',dict(event='CATEGORY_START',category=cat,codec=codec_name,files=summary['categories'][cat]['files'],source_content_bytes=summary['categories'][cat]['stored_content_bytes']))
     try:
      with ExitStack() as codec:
       output=codec.enter_context(gzip.GzipFile(filename='',fileobj=writer,mode='wb',compresslevel=1,mtime=0)) if compressed else writer
       # Matching stream/copy buffers avoid tarfile's default10KiB repeated
       # tail-slicing of1MiB chunks. Memory is bounded independently of file size.
       with tarfile.open(fileobj=output,mode='w|',bufsize=MIB,format=tarfile.PAX_FORMAT,dereference=False) as tar:
        tar.copybufsize=MIB
        for row in rows(dest/'inventory'/journal):
         deadline.check();add(tar,row)
         if row['type']=='file':done+=1;content+=0 if row.get('hardlink_to') else row['size']
         if time.perf_counter()-last_progress>=10:
          event(dest,'CATEGORY_EVENTS.jsonl',dict(event='CATEGORY_PROGRESS',category=cat,files_done=done,input_content_bytes=content,closed_parts=sum(v['status']=='CLOSED_IMMUTABLE' for v in writer.parts),wall_seconds=time.perf_counter()-cat_start));last_progress=time.perf_counter()
        if cat=='compact_core':
         metadata=[dest/'inventory'/n for n in ['SUMMARY.json','DEPENDENCIES.json','EXCLUSIONS_AND_SYMLINKS.json']]+[dest/'inventory'/v for v in summary['journal_files'].values()]
         for p in metadata:
          deadline.check();ss=p.stat();add(tar,dict(source=str(p),archive_path='bundle_metadata/'+str(p.relative_to(dest/'inventory')),size=ss.st_size,mtime_ns=ss.st_mtime_ns,mode=0o644,type='file'))
      writer.close_part();event(dest,'CATEGORY_EVENTS.jsonl',dict(event='CATEGORY_COMPLETE',category=cat,files_done=done,output_bytes=writer.total,parts=len(writer.parts),wall_seconds=time.perf_counter()-cat_start))
     except BaseException:
      writer.abort();event(dest,'CATEGORY_EVENTS.jsonl',dict(event='CATEGORY_FAILED_NOT_COMPLETE',category=cat,files_done=done,output_bytes=writer.total));raise
     finally:receipt['parts'].extend([dict(category=cat,codec=codec_name,**p) for p in writer.parts])
    event(dest,'CATEGORY_EVENTS.jsonl',dict(event='FINAL_METADATA_RECHECK_START'))
    # No second tree traversal. Existing file and directory times detect late
    # edits/additions/removals. No arrays reopened for analysis or hashes.
    for journal in summary['journal_files'].values():
     for row in rows(dest/'inventory'/journal):deadline.check();verify_row(row)
    if last_ledger()['sequence']!=ledger['sequence']:raise RuntimeError('ledger changed despite quiet lock')
    receipt.update(status='COMPLETE',inventory_summary=summary,ledger_sequence=ledger['sequence'],final_prepack_research_ledger_included=True)
   except BaseException as exc:receipt['error']=repr(exc);raise
   finally:
    receipt.update(wall_seconds=time.perf_counter()-started,completed_unix=time.time(),peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024));dump(dest/'RECEIPT.json',receipt)
 finally:
  # ExitStack releases every campaign lock BEFORE delivery accounting. Do not
  # run inside Book, do not append to the immutable researchclosed ledger.
  if dest.exists():
   delivery=dict(kind='DELIVERY_NOT_RESEARCH',status=receipt['status'],bundle=str(dest),wall_seconds=time.perf_counter()-started,CPU_seconds=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),campaign_locks_released=True,research_ledger_modified=False,research_RHS=0,inference=0,training=0,hashes_used=False,completed_unix=time.time())
   dump(dest/'DELIVERY_RECEIPT.json',delivery)
   with (PACK/'DELIVERY_ACCOUNTING.jsonl').open('a') as f:
    fcntl.flock(f,fcntl.LOCK_EX);f.write(json.dumps(delivery,separators=(',',':'))+'\n');f.flush();os.fsync(f.fileno());fcntl.flock(f,fcntl.LOCK_UN)
 return receipt

def verify_restored(a):
 manifest=a.inventory/'SUMMARY.json';summary=json.loads(manifest.read_text());restored=a.restored.absolute();checked=0
 if not restored.is_dir():raise ValueError('existing NEW extraction directory required')
 for journal in summary['journal_files'].values():
  for row in rows(a.inventory/journal):
   rel=Path(row['archive_path'])
   if rel.is_absolute() or '..' in rel.parts:raise ValueError('unsafe archive path')
   target=restored/rel;info=target.lstat()
   if stat.S_ISLNK(info.st_mode) or (row['type']=='file' and (not stat.S_ISREG(info.st_mode) or info.st_size!=row['size'])):raise RuntimeError('restored type/size mismatch: '+str(target))
   if a.restore_times:os.utime(target,ns=(row['mtime_ns'],row['mtime_ns']),follow_symlinks=False)
   if target.stat().st_mtime_ns!=row['mtime_ns']:raise RuntimeError('mtime precision mismatch; use explicit --restore-times on extraction only')
   checked+=1
 return dict(status='RESTORED_SIZES_AND_NANOSECOND_TIMES_MATCH',entries_checked=checked,hashes_used=False,content_digest_verification=False)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['inventory','pack','verify-restored']);p.add_argument('--out',type=Path);p.add_argument('--inventory',type=Path);p.add_argument('--restored',type=Path);p.add_argument('--restore-times',action='store_true');p.add_argument('--confirmation',type=Path);p.add_argument('--part-bytes',type=int,default=LIMIT);p.add_argument('--raw-compression',choices=['none','gzip'],default='none');p.add_argument('--max-wall',type=float,default=14400);a=p.parse_args()
 if a.command=='verify-restored':
  if not a.inventory or not a.restored:p.error('verification requires --inventory and --restored')
  print(json.dumps(verify_restored(a)));return
 if a.out is None:p.error('--out required')
 dest=a.out.expanduser().absolute();PACK.mkdir(parents=True,exist_ok=True)
 if PACK not in dest.parents or dest.exists():raise ValueError('new output directory strictly below OUT/packaging required')
 if not 0<a.max_wall<=14400:raise ValueError('delivery maxwall must be >0 and <=4hours')
 if not 1<=a.part_bytes<=LIMIT:raise ValueError('every split asset must be <=1,500,000,000 bytes')
 if a.command=='pack' and not a.confirmation:p.error('pack requires --confirmation; no default archive')
 result=enumerate_once(dest) if a.command=='inventory' else pack(a,dest)
 print(json.dumps(dict(status=result['status'],out=str(dest),files=result.get('files'),bytes=result.get('logical_file_bytes'),wall=result.get('scan_wall_seconds',result.get('wall_seconds')))))
if __name__=='__main__':main()

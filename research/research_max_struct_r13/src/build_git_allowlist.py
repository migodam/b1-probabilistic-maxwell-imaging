"""Explicit metadata-first Git publication inventory; never stages or deletes.

Run again after final docs/site generation. Raw studies/arrays are Release-only;
private operational/chat/fulltext files excluded from Git and release assets.
"""
import argparse,json,re,time
from pathlib import Path
from datetime import datetime,timezone
PRIVATE_PARTS={'.git','.venv','venv','__pycache__','.cache','node_modules','private_operational','scorer_fixture','raw_retrieval','fulltext','fulltexts','third_party_papers'}
RAW_PARTS={'RESULTS','sealed_truth','data_public','archives'}
RAW_NAMES=('CANDIDATE','HISTORY','PAID_PREFIX','ALL_ESTIMATOR_ROWS','ALL_PAID_AND_OVERRUN','SCORES','COST_INVENTORY','ORACLE_POOL_DIAGNOSTIC','MEASURED_PREDICTIONS','GEOMETRY')

def classify(path,root):
 rel=path.relative_to(root);parts=rel.parts;name=path.name;lower=name.lower();suffix=path.suffix.lower();size=path.stat().st_size
 if path.is_symlink():return 'exclude_all','symlink'
 if any(x in PRIVATE_PARTS or 'fixture' in x.lower() or 'fake' in x.lower() for x in parts):return 'exclude_all','private operation/environment/cache/fixture/fulltext'
 if any(parts[i:i+2] in [('audit','transfers'),('audit','remote_metadata')] for i in range(len(parts)-1)) or name=='REMOTE_TRANSFER_PROGRESS.json':return 'exclude_all','private operational transfer/controller metadata'
 if any(x in lower for x in ['chat_readback','private_chat','conversation_readback','thread_readback']) or suffix in ['.pem','.key','.pfx','.exp'] or lower in ['.env','credentials.json','token.json','id_rsa','id_ed25519','.ds_store'] or name.endswith('.lock'):return 'exclude_all','private payload/key/raw third-party input'
 if suffix=='.html' and 'audit' in parts and 'readback' in lower:return 'release_only','Archived prior website readback, not an active local page'
 if name=='QUEUE_ATTEMPT_HISTORY.json' and size<2_000_000:return 'git_allowlist','Compact failed-attempt cost/provenance receipt, not numerical candidate history'
 if 'sources' in parts and not('pro' in parts or 'user' in parts or name in ['SOURCE_INDEX.json','SOURCE_ACCESS.json']):return 'exclude_all','third-party retrieved source/fulltext; preserve access metadata elsewhere'
 if suffix in ['.py','.sh','.ps1'] and size<=2_000_000 and name!='build_git_allowlist.py':
  text=path.read_text(errors='replace')
  if re.search(r'''["'](?:ssh|scp|sftp)["']|IdentityFile|BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY''',text):return 'exclude_all','operational SSH/key argument source, retained locally'
 if any(x in RAW_PARTS for x in parts) or 'predictions' in parts and suffix!='.json':return 'release_only','raw result/input/truth/archive/state'
 if suffix in ['.npz','.npy','.pt','.pth','.zip','.gz','.7z','.tar','.exp','.log','.jsonl'] or '.part' in name or name.startswith('runtime_'):return 'release_only','raw arrays/models/archives/logs'
 if suffix=='.pdf' and not('FIGURES' in parts or 'figures' in parts or name=='MAX_STRUCT_Research_20261112.pdf'):return 'exclude_all','nonfigure third-party/full-paper PDF'
 if size>50_000_000:return 'release_only','single file above50MB review limit'
 if suffix in ['.png','.jpg','.jpeg','.svg','.webp','.pdf']:return 'git_allowlist','actual figure or own Pro report'
 if suffix in ['.py','.sh','.ps1','.md','.tex','.bib','.html','.css','.js'] or name in ['.gitignore','LICENSE','LICENSE.txt']:return 'git_allowlist','own source/configuration/docs/site/attribution'
 if suffix in ['.json','.csv','.tsv']:
  if size>=2_000_000 or any(x in name.upper() for x in RAW_NAMES) and 'review_tables' not in parts:return 'release_only','bulk history/candidate/prefix/raw table or compact limit'
  if 'queues' in parts and name not in ['RECEIPT.json','SOURCE_SNAPSHOT.json','RESOURCE_PROBE.json','initial_RECEIPT.json']:return 'release_only','noncanonical queue snapshots'
  if name=='DOWNLOADS.json' and json.loads(path.read_text()).get('published') is not True:return 'exclude_all','unpublished download plan has no public manifest status'
  return 'git_allowlist','compact scientific table/config/audit/selection index'
 return 'parent_review','unclassified format'

def main():
 c,w=time.process_time(),time.perf_counter();p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();root=a.root.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);rows=[]
 for f in sorted(root.rglob('*')):
  if not f.is_file() or f.is_relative_to(out):continue
  dest,reason=classify(f,root);rows.append(dict(path=f.relative_to(root).as_posix(),repository_path='research/research_max_struct_r13/'+f.relative_to(root).as_posix(),bytes=f.stat().st_size,destination=dest,reason=reason))
 allowed=[r for r in rows if r['destination']=='git_allowlist'];summary={dest:dict(files=sum(r['destination']==dest for r in rows),bytes=sum(r['bytes'] for r in rows if r['destination']==dest)) for dest in ['git_allowlist','release_only','exclude_all','parent_review']}
 d=dict(status='EXPLICIT_PROSPECTIVE_ALLOWLIST_NO_GIT_MUTATION',created_at_utc=datetime.now(timezone.utc).isoformat(),root=str(root),summary=summary,git_allowlist=[r['path'] for r in allowed],repository_git_allowlist=[r['repository_path'] for r in allowed],files=rows,files_over50MB=[r for r in rows if r['bytes']>50_000_000],CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,host='Apple_local',NPZ_payloads_read=False,source_files_preserved=True,staging_or_publication_performed=False,refresh='Rerun into a NEW output after root scientific docs/site changes; inventory output itself excluded to avoid self-reference')
 (out/'GIT_ALLOWLIST.json').write_text(json.dumps(d,indent=2)+'\n');(out/'GIT_ALLOWLIST.txt').write_text('\n'.join(r['repository_path'] for r in allowed)+'\n');print(json.dumps(dict(summary=summary,CPU_seconds=d['CPU_seconds'],output=str(out))))
if __name__=='__main__':main()

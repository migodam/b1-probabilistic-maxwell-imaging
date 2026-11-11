from pathlib import Path
import json,time,zipfile,sys
folder=Path(sys.argv[1]);plan=json.loads((folder/'PLAN.json').read_text());rows=[];cpu=time.process_time();wall=time.perf_counter()
for asset in plan['archives']:
 out=folder/asset['name'];index=out.with_suffix('.index.json')
 if out.exists()and index.exists():
  old=json.loads(index.read_text())
  if old['archive_bytes']!=out.stat().st_size or old['entries']!=asset['entries']:raise RuntimeError('Existing archive does not match frozen plan')
  rows.append(old);print(json.dumps(dict(asset=out.name,resumed=True,bytes=out.stat().st_size)),flush=True);continue
 if out.exists():raise RuntimeError('Partial archive preserved; inspect and choose explicit new name, never overwrite')
 c,t=time.process_time(),time.perf_counter()
 with zipfile.ZipFile(out,'x',compression=zipfile.ZIP_STORED,allowZip64=True)as z:
  for e in asset['entries']:
   p=Path(e['source_root'])/e['source_relative']
   if not p.is_file()or p.stat().st_size!=e['bytes']:raise RuntimeError('Source size changed: '+str(p))
   z.write(p,e['archive_path'])
 if out.stat().st_size>=plan['max_asset_bytes']:raise RuntimeError('Asset cap exceeded')
 with zipfile.ZipFile(out)as z:
  actual=[(x.filename,x.file_size)for x in z.infolist()]
  expected=[(e['archive_path'],e['bytes'])for e in asset['entries']]
  if actual!=expected:raise RuntimeError('Central directory index mismatch')
 row=dict(name=out.name,archive_bytes=out.stat().st_size,file_count=len(asset['entries']),source_bytes=asset['source_bytes'],entries=asset['entries'],CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,hashes='NOT_COMPUTED',verification='ZIP central-directory names/counts/byte sizes')
 index.write_text(json.dumps(row,indent=2)+'\n');rows.append(row);print(json.dumps(dict(asset=out.name,bytes=row['archive_bytes'],files=row['file_count'])),flush=True)
receipt=dict(status='COMPLETE',archives=[{k:v for k,v in x.items()if k!='entries'}for x in rows],archive_phase_CPU_seconds_sum=sum(x['CPU_seconds']for x in rows),CPU_seconds=time.process_time()-cpu,process_CPU_seconds_near_exit=time.process_time(),wall_seconds=time.perf_counter()-wall,scientific_forward_calls=0,source_files_deleted=0,hashes='NOT_COMPUTED')
if(folder/'PACKAGING_RECEIPT.json').exists():
 with(folder/'PACKAGING_RESUMES.jsonl').open('a')as f:f.write(json.dumps(receipt)+'\n')
else:(folder/'PACKAGING_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')

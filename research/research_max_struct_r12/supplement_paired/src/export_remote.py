"""Archive only completed paired outputs; no source deletion or scientific run."""
from pathlib import Path
import json,time,zipfile,argparse
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();root=a.root
receipt=json.loads((root/'queues/main/RECEIPT.json').read_text())
if receipt['status'] not in ('COMPLETED','COMPLETED_WITH_ERRORS'):raise RuntimeError('Finite queue still active')
start,wall=time.process_time(),time.perf_counter();target=root/'paired_results_full.zip'
with zipfile.ZipFile(target,'x',zipfile.ZIP_DEFLATED,compresslevel=3) as z:
 for folder in ('RESULTS','queues'):
  for path in sorted((root/folder).rglob('*')):
   if path.is_file():z.write(path,path.relative_to(root))
 info=[dict(name=x.filename,bytes=x.file_size,compressed_bytes=x.compress_size) for x in z.infolist()]
record=dict(archive=target.name,bytes=target.stat().st_size,files=len(info),entries=info,CPU_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall,source_deleted=False,new_physics_runs=0)
(root/'REMOTE_ARCHIVE_RECEIPT.json').write_text(json.dumps(record,indent=2));print(json.dumps({k:v for k,v in record.items() if k!='entries'}))

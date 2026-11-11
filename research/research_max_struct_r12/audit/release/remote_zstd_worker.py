from pathlib import Path
import json,time,sys,zstandard
folder=Path(sys.argv[1]);plan=json.loads((folder/'PLAN.json').read_text());receipt=folder/'ZSTD_RECEIPT.json';old=json.loads(receipt.read_text())if receipt.exists()else dict(archives=[]);prior={x['zip_name']:x for x in old['archives']};rows=[]
params=zstandard.ZstdCompressionParameters.from_level(6,window_log=27,enable_ldm=True,threads=2,write_checksum=0)
for asset in plan['archives']:
 src=folder/asset['name'];out=folder/(asset['name']+'.zst')
 if out.exists()and asset['name']in prior:
  row=prior[asset['name']]
  if out.stat().st_size!=row['bytes']:raise RuntimeError('Existing zstd size mismatch')
  rows.append(row);continue
 if out.exists():raise RuntimeError('Partial zstd preserved; explicit inspection required')
 c,t=time.process_time(),time.perf_counter()
 with src.open('rb')as f,out.open('xb')as target:zstandard.ZstdCompressor(compression_params=params).copy_stream(f,target)
 row=dict(zip_name=src.name,name=out.name,bytes=out.stat().st_size,zip_bytes=src.stat().st_size,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,version=zstandard.__version__,window_log=27,enable_ldm=True,threads=2,frame_checksum=False)
 rows.append(row);receipt.write_text(json.dumps(dict(status='COMPLETE'if len(rows)==len(plan['archives'])else'RUNNING',archives=rows,scientific_forward_calls=0,source_files_deleted=0),indent=2)+'\n');print(json.dumps(row),flush=True)
receipt.write_text(json.dumps(dict(status='COMPLETE',archives=rows,scientific_forward_calls=0,source_files_deleted=0),indent=2)+'\n')

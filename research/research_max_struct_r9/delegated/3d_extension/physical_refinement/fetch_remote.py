"""Remote only: verify exact sealed files then package completed public outputs."""
import json,time,zipfile
from pathlib import Path
T=time.perf_counter();C=time.process_time();REPO=Path('D:/AI/MAX_STRUCT_R9_PHYSREF_20261106');ROOT=REPO/'research/research_max_struct_r9';OUT=ROOT/'RESULTS/original_vector3d/physical_refinement'
seal=json.loads((OUT/'SEALED_ALL_THREE.json').read_text());assert len(seal['entries'])==3
checks=[]
for e in seal['entries']:
    for f in e['files']:
        p=REPO/f['path'];s=p.stat();assert s.st_size==f['bytes'] and s.st_mtime_ns==f['mtime_ns'];checks.append(f['path'])
(OUT/'FETCH_SEAL_VERIFICATION.json').write_text(json.dumps(dict(status='ALL_REMOTE_BYTES_AND_MTIMES_MATCH_SEAL',files=len(checks),GT_access=False)))
with zipfile.ZipFile(REPO/'COMPLETED_PUBLIC_OUTPUTS.zip','w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob('*')):
        if p.is_file():z.write(p,str(p.relative_to(REPO)))
    z.write(REPO/'EXTRACTION_COST.json','EXTRACTION_COST.json')
receipt=dict(platform='XINAN',CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,operation='VERIFY_ALL_THREE_SEAL_AND_ZIP_PUBLIC_OUTPUTS',scope='Own remote helper only; GPU child process CPU excluded',additive=True,included_in=None,counts=dict(new_physical_RHS=0))
(OUT/'REMOTE_FETCH_PROCESS_COST.json').write_text(json.dumps(receipt,indent=2))
with zipfile.ZipFile(REPO/'COMPLETED_PUBLIC_OUTPUTS.zip','a',zipfile.ZIP_DEFLATED) as z:z.write(OUT/'REMOTE_FETCH_PROCESS_COST.json',str((OUT/'REMOTE_FETCH_PROCESS_COST.json').relative_to(REPO)))
print(json.dumps(receipt))

"""Own remote artifacts byte harvest. No NPZ loads, GT scores or solves."""
import time
T=time.perf_counter();C=time.process_time()
from pathlib import Path
import json,zipfile
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]


def main():
    tasks=['v08_'+a for a in ('sourceCIE','GlobalDual','AplusGlobal')]
    checks=[]
    for name in tasks:
        folder=ROOT/'RESULTS/vector'/name
        if not (folder/'SEALED.json').exists():raise RuntimeError('PRIMARY_NOT_READY '+name)
        seal=json.loads((folder/'SEALED.json').read_text())
        for row in seal['files']:
            p=REPO/row['path'];s=p.stat()
            if s.st_size!=row['bytes'] or s.st_mtime_ns!=row['mtime_ns']:raise RuntimeError('Source seal changed '+str(p))
            checks.append(row)
    verification=dict(status='ALL_V08_PRIMARY_SOURCE_BYTES_AND_MTIMES_MATCH_SEALS',tasks=tasks,files=len(checks),GT_arrays_read=False)
    (ROOT/'delegated/vector/PRIMARY_REMOTE_SEAL_VERIFICATION.json').write_text(json.dumps(verification,indent=2)+'\n')
    paths=[]
    for name in tasks:paths.extend(p for p in (ROOT/'RESULTS/vector'/name).rglob('*') if p.is_file())
    # Own producer metadata, including stopped attempts; no unrelated app files.
    for parent in [ROOT/'RESULTS/vector_engineering_attempts',ROOT/'delegated/vector']:
        if parent.exists():
            for p in parent.rglob('*'):
                if p.is_file() and (p.suffix in ('.json','.jsonl') or parent.name=='vector_engineering_attempts'):
                    paths.append(p)
    for name in ('gpu_backend_probe','gpu_dense_backend_probe','generation_windows_v04_attempt1_missing_dependency','generation_windows_v04','generation_gpu_v08','case_v08_cuda_dense_process','queue_probe'):
        folder=ROOT/'RESULTS/vector'/name
        if folder.exists():paths.extend(p for p in folder.rglob('*') if p.is_file() and p.suffix in ('.json','.jsonl'))
    paths=sorted(set(paths))
    archive=REPO/'VECTOR_FINAL_PUBLIC_AND_ENGINEERING.zip'
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_STORED) as z:
        for p in paths:z.write(p,p.relative_to(REPO).as_posix())
        for p in REPO.glob('*COST*.json'):
            arc='research/research_max_struct_r10/delegated/vector/remote_repository_costs/'+p.name
            z.write(p,arc)
    receipt=dict(status='REMOTE_ARTIFACTS_BYTE_HARVEST_NO_GT_LOAD',CPU_seconds=time.process_time()-C,
        wall_seconds=time.perf_counter()-T,platform='XINAN',archive_bytes=archive.stat().st_size,
        archived_files=len(paths),source_seal_verification=verification,
        scope='Own byte packaging and metadata source verification; no arrays interpreted, no score',physical_RHS=0)
    (ROOT/'delegated/vector/FINAL_HARVEST_COST.json').write_text(json.dumps(receipt,indent=2)+'\n')
    with zipfile.ZipFile(archive,'a',zipfile.ZIP_STORED) as z:z.write(ROOT/'delegated/vector/FINAL_HARVEST_COST.json',(ROOT/'delegated/vector/FINAL_HARVEST_COST.json').relative_to(REPO).as_posix())
    print(json.dumps(receipt),flush=True)


if __name__=='__main__':main()

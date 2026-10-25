"""Lossless independent R5 evidence archives; never delete scientific sources.

Archive paths are relative to Gaussian/. Streaming byte verification preserves
original file bytes. There are no SHA-256/data hash checks and no new physics.
"""
from pathlib import Path
import sys,json,tarfile,subprocess,time,os

P=Path(__file__).resolve().parents[1];REPO=P.parents[1];GAUSSIAN=REPO.parents[1]
OUT=GAUSSIAN.parent/'research/delegated/r5_delivery_20261026'
LEDGER=P/'DELIVERY_LEDGER.jsonl'

def event(x):
    with LEDGER.open('a') as f:f.write(json.dumps(dict(time=time.time(),physics_RHS=0,**x))+'\n')
def files(directory):
    return [p for p in sorted(directory.rglob('*')) if p.is_file() and not p.is_symlink()
            and '__pycache__' not in p.parts and p.suffix!='.pyc']
def pack(name,paths):
    dest=OUT/(name+'.tar.zst');receipt=OUT/(name+'.json')
    if dest.exists():
        if receipt.exists():return json.loads(receipt.read_text())
        raise RuntimeError('Existing unreceipted archive; do not overwrite')
    paths=sorted(set(paths));start=time.perf_counter();cpu=time.process_time()
    manifest=[dict(path=str(p.relative_to(GAUSSIAN)),bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns) for p in paths]
    partial=dest.with_suffix(dest.suffix+'.partial')
    event(dict(event='PACK_BEGIN',name=name,files=len(paths),source_bytes=sum(x['bytes'] for x in manifest)))
    try:
        with partial.open('xb') as out:
            proc=subprocess.Popen(['/opt/homebrew/bin/zstd','-9','-T1','-q','-c'],stdin=subprocess.PIPE,stdout=out)
            with tarfile.open(fileobj=proc.stdin,mode='w|') as tar:
                for p,meta in zip(paths,manifest):
                    st=p.stat()
                    tar.add(p,arcname=meta['path'],recursive=False)
                    after=p.stat()
                    if (st.st_size,st.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise RuntimeError('Source changed during evidence pack')
            proc.stdin.close()
            if proc.wait()!=0:raise RuntimeError('zstd pack failed')
            out.flush();os.fsync(out.fileno())
        # Read every archive member and compare with the existing source bytes.
        proc=subprocess.Popen(['/opt/homebrew/bin/zstd','-q','-d','-c',str(partial)],stdout=subprocess.PIPE)
        expected={x['path']:x for x in manifest};seen=set()
        with tarfile.open(fileobj=proc.stdout,mode='r|') as tar:
            for member in tar:
                if not member.isfile() or member.name not in expected or member.name in seen:raise RuntimeError('Unexpected archive member')
                meta=expected[member.name]
                if member.size!=meta['bytes']:raise RuntimeError('Member size differs')
                archived=tar.extractfile(member)
                with (GAUSSIAN/member.name).open('rb') as original:
                    while True:
                        block=archived.read(1024*1024)
                        if block!=original.read(len(block)):raise RuntimeError('Archive roundtrip byte mismatch')
                        if not block:
                            if original.read(1):raise RuntimeError('Source length differs')
                            break
                seen.add(member.name)
        while proc.stdout.read(1024*1024):pass
        if proc.wait()!=0 or seen!=set(expected):raise RuntimeError('Incomplete zstd archive verification')
        partial.rename(dest)
        result=dict(name=dest.name,status='FULL_BYTE_ROUNDTRIP_PASS',files=len(paths),
                    original_bytes=sum(x['bytes'] for x in manifest),archive_bytes=dest.stat().st_size,
                    wall_seconds=time.perf_counter()-start,parent_CPU_seconds=time.process_time()-cpu,
                    compression='tar + zstd9 single thread',paths_relative_to='Gaussian',
                    scientific_source_deletion=False,physics_RHS=0,members=manifest)
        receipt.write_text(json.dumps(result,indent=2)+'\n')
        event(dict(event='PACK_COMPLETE',name=name,files=len(paths),archive_bytes=result['archive_bytes'],wall_seconds=result['wall_seconds']))
        print(json.dumps({k:v for k,v in result.items() if k!='members'}),flush=True)
        return result
    except BaseException as exc:
        event(dict(event='PACK_FAILED_RETAINED',name=name,error=repr(exc),wall_seconds=time.perf_counter()-start))
        raise

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    groups=[('r5-track-a',files(P/'stage_a_response_coordinates')),
            ('r5-track-b',files(P/'stage_b_material_risk')),
            ('r5-shared-physics',files(P/'shared')),
            ('r5-parent-probes',files(P/'audit/curvature_bounds')+files(P/'audit/seam_maxwell'))]
    refs=json.loads((P/'audit/EVIDENCE_INPUT_REFERENCES.json').read_text())['indispensable_historical_references']
    hist=[]
    for item in refs:
        p=Path(item['path'])
        if p.is_relative_to(REPO):continue
        hist+=files(p) if p.is_dir() else [p]
    # Source helpers used by the original checkpoint formats, not old datasets.
    hist+=files(REPO.parent/'b1-p3l-r4-conditional-transport/src')
    groups.append(('r5-essential-reused-inputs',hist))
    results=[pack(name,paths) for name,paths in groups]
    compact=dict(status='ALL_BYTE_ROUNDTRIPS_PASS',archives=[{k:v for k,v in r.items() if k!='members'} for r in results],
                 source_arrays_preserved=True,physics_RHS=0)
    (P/'EVIDENCE_ARCHIVES.json').write_text(json.dumps(compact,indent=2)+'\n')
    print(json.dumps(dict(status=compact['status'],archives=len(results),archive_bytes=sum(r['archive_bytes'] for r in results))),flush=True)
if __name__=='__main__':main()

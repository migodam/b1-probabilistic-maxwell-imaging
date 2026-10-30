"""Archive the already-frozen physical arrays once. No scientific replay."""
from pathlib import Path
import subprocess,json,time

P=Path(__file__).resolve().parents[1];REPO=P.parents[1]
def main():
    dest=P/'packaging';dest.mkdir(exist_ok=True)
    batches=[('R12_raw_P1.tar.zst',[P/'results/P1']),('R12_raw_N1.tar.zst',[P/'results/N'])]
    for j in range(3):batches.append((f'R12_raw_N2_part{j+1}.tar.zst',[P/'results/N2'/f'N_{i:03d}' for i in range(3*j,3*j+3)]))
    receipt=[]
    for name,paths in batches:
        out=dest/name
        if out.exists():raise RuntimeError('archive exists; no overwrite')
        rel=[str(p.relative_to(REPO)) for p in paths]
        st=time.perf_counter()
        tar=subprocess.Popen(['/usr/bin/tar','-cf','-','-C',str(REPO),*rel],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        z=subprocess.run(['/opt/homebrew/bin/zstd','-T2','-3','-o',str(out)],stdin=tar.stdout,capture_output=True)
        tar.stdout.close();err=tar.stderr.read().decode();rc=tar.wait()
        row=dict(name=name,bytes=out.stat().st_size,source_paths=rel,tar_exit=rc,zstd_exit=z.returncode,wall_seconds=time.perf_counter()-st,stderr=err+z.stderr.decode())
        receipt.append(row);print(name,row['bytes'],round(row['wall_seconds'],2),flush=True)
        (dest/'RAW_ARCHIVE_CREATION.json').write_text(json.dumps(receipt,indent=2)+'\n')
        if rc or z.returncode:raise RuntimeError('archive creation failed; retained')
if __name__=='__main__':main()

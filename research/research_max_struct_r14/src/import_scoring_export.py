"""Restore explicit material-only transport, retaining remote full originals."""
import argparse,json,time,zipfile
from pathlib import Path,PurePosixPath


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--archive',type=Path,required=True);p.add_argument('--campaign',required=True);a=p.parse_args()
    r=a.root.resolve();cpu=time.process_time();wall=time.perf_counter();count=0;same=[]
    metadata=r/'audit/transports';metadata.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.archive) as z:
        for info in z.infolist():
            if info.is_dir():continue
            rel=PurePosixPath(info.filename.replace('\\','/'))
            if rel.is_absolute() or '..' in rel.parts:raise ValueError(info.filename)
            dest=metadata/(a.campaign+'_OFFLINE_EXPORT_RECEIPT.json') if rel.as_posix()=='OFFLINE_EXPORT_RECEIPT.json' else r.joinpath(*rel.parts)
            dest.parent.mkdir(parents=True,exist_ok=True);data=z.read(info)
            if dest.exists():
                old=dest.read_bytes()
                if old!=data:
                    if dest.suffix=='.json' and json.loads(old)==json.loads(data):same.append(str(rel))
                    else:raise FileExistsError(str(dest))
            else:dest.write_bytes(data)
            count+=1
    value=dict(status='MATERIAL_ONLY_EXPORT_IMPORTED',campaign=a.campaign,files=count,
        archive_bytes=a.archive.stat().st_size,CPU_seconds=time.process_time()-cpu,
        wall_seconds=time.perf_counter()-wall,semantically_identical_JSON_kept=same,
        material_transform='NONE',truth_access='NONE',new_physics=0,
        full_fields='Full remote original archives required; these selected-view NPZ contain only unchanged m')
    with (metadata/(a.campaign+'_IMPORT_RECEIPT.json')).open('x') as f:json.dump(value,f,indent=2);f.write('\n')
    print(json.dumps(value))


if __name__=='__main__':main()

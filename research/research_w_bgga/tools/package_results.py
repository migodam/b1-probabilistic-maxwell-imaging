"""Package saved, public stage-relative evidence, without digests or physics."""
from pathlib import Path
import json,zipfile
ROOT=Path(__file__).resolve().parents[1]

def main():
    target=ROOT/'release_assets';target.mkdir(exist_ok=True)
    groups={
        'g0-physics':['results/preflight','results/g0_qualification','results/g0_derivatives',
            'results/g0_operators','results/g0_stress','results/g0_fixed_posterior'],
        'g1-weights-and-intervals':['results/g1_replay_001','results/g1_normalization_repair','results/g1_offline_extra'],
        'receipts-and-figures':['results/jobs','figures','research/delegated/g1_figures'],
    }
    manifest=[]
    for group,paths in groups.items():
        files=sorted({p for rel in paths for p in (ROOT/rel).rglob('*') if p.is_file() and p.suffix not in ('.lock','.pyc')})
        dest=target/f'research-w-bgga-{group}.zip'
        if dest.exists():raise FileExistsError(f'Unique publication asset required: {dest}')
        with zipfile.ZipFile(dest,'x',zipfile.ZIP_DEFLATED,compresslevel=6,allowZip64=True) as archive:
            for p in files: archive.write(p,p.relative_to(ROOT).as_posix())
        with zipfile.ZipFile(dest) as archive:
            members=archive.infolist()
            assert len(members)==len(files)
            assert all(m.file_size==p.stat().st_size for m,p in zip(members,files))
        manifest.append(dict(asset=dest.name,bytes=dest.stat().st_size,uncompressed_bytes=sum(p.stat().st_size for p in files),
            members=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size) for p in files],
            extract_into='research/research_w_bgga',scope='saved actual arrays/traces, no new physics'))
    (ROOT/'PUBLICATION_ASSETS.json').write_text(json.dumps(dict(tag='research-w-bgga-20261028',assets=manifest,
        file_digest_checks=False,old_bank='separate unchanged research-w-slam-20261025 controlled_E2E3 asset',
        accounting='closing accounting asset added after anonymous readback; science arrays in these immutable stage archives'),indent=2)+'\n')
    print(json.dumps([dict(asset=r['asset'],bytes=r['bytes'],files=len(r['members'])) for r in manifest]))
if __name__=='__main__':main()

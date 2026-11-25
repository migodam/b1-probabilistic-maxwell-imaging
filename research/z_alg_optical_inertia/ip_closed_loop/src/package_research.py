"""Package the closed-loop study at its original repository-relative paths."""
import argparse
from collections import Counter
import json
import os
from pathlib import Path
import re
from tempfile import NamedTemporaryFile
from time import perf_counter, process_time
from zipfile import ZipFile, ZIP_DEFLATED


PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[2]
STUDY = PAPER.parent.relative_to(ROOT)
R8 = Path('research/delegated/max_struct_r1/r8_worktree/research')
MANIFEST = 'RESEARCH_ARCHIVE_MANIFEST.json'
DEFAULT = PAPER / 'closed_loop_reproduction.zip'

# Keep borrowed data and imports at precisely the paths used by each program.
STUDY_FILES = {
    'inherited_map_operator': ('fixed_operator_paper/experiments/B/public/native_geometry.npz',),
    'inherited_reduction_inputs': (
        'synthetic/public/pw360_10x20.npz',
        'ip_guided_v4/boxes/results/shiftbox32_01__pw360_10x20__n02.json',
        'ip_guided_v4/boxes/results/shiftbox32_01__pw360_10x20__n02.npz',
        'ip_guided_v4/boxes/results/shiftbox32_01__pw360_10x20__n02__offline.npz',
        'ip_guided_v4/boxes/results/disk_radius.json',
        'ip_final_v5/strong_control/strong_ring_arrays.npz',
        'ip_final_v5/strong_control/strong_ring_verified.json',
    ),
    'laboratory_code': ('validation/manitoba/physical.py', 'fresnel/physical.py'),
}
R8_FILES = (
    'research_max_struct_r8/src/fresnel_adapter.py',
    'research_max_struct_r8/src/paper_physics.py',
    'research_max_struct_r6/src/fft_tsom_ops.py',
    'research_max_struct_r4/src/cie_onehop.py',
    'research_max_struct_r2/src/cie_global.py',
    'research_fbg_r3_experiments/src/__init__.py',
    'research_fbg_r3_experiments/src/cie_tsom_baseline.py',
    'research_fbg_r3_experiments/src/physics.py',
    'research_z_alg_r4_diversity/src/__init__.py',
    'research_z_alg_r4_diversity/src/physics.py',
    'research_z_alg_r4_diversity/src/solver.py',
)
MANDATORY = (
    'REPRODUCE.md', 'requirements.txt', 'manuscript_ip.tex', 'manuscript_ip.pdf',
    'maps/public/acquisition.json', 'maps/public/geometry_32.npz', 'maps/public/geometry_96.npz',
    'maps/results/ipcl-fov-f1__M96-R01__weak__7tx9rx__n03__b0.03125__broad_sector.npz',
    'maps/results/arb_M96-R01_broad.json', 'beams/dependencies.txt', 'beams/paired_replay_power_96.json',
    'pair/power_96/proposal_arrays.npz', 'reduction/result_32.json', 'reduction/result_96.json',
    'reduction/psd_maps/dark_replay_full_32.json', 'lab/measured.py', 'lab/held_ports.csv',
)
STUDY_DIRECTORIES = {'adaptive','beams','figures','lab','maps','pair','reduction','results','reviews','src','theory','sources'}
FORMATS = {'.cls','.css','.csv','.html','.json','.md','.npz','.pdf','.png','.py','.svg','.tex','.txt'}
LAB_RELEASE = {'DESIGN.md','held_ports.csv','measured.py','results.json','loaded_headroom.json'}
REVIEW_RELEASE = {'RESPONSE_ZH.md','neon-01.md','neon-02.md','neon-03.md','neon-04.md'}
HIDDEN_TEXT = re.compile(r'(?:file://|/(?:Users|Volumes)/|/private/tmp/|/home/[^/\s]+/)')


def public_study_file(path):
    relative = path.relative_to(PAPER)
    parts = relative.parts
    if path.suffix.lower() not in FORMATS or any(
        part.startswith('.') or part in {'__pycache__','tmp','cache','private_prompts'} for part in parts
    ):
        return False
    if path.name in {'ORIGINAL_REQUEST.txt','TASK_RESUME.md','COVER_LETTER_IP.md'} or '.raw.' in path.name:
        return False
    if len(parts) > 1 and parts[0] not in STUDY_DIRECTORIES:
        return False
    if parts[0] == 'lab':
        return len(parts) == 2 and path.name in LAB_RELEASE
    if parts[0] == 'reviews':
        return len(parts) == 2 and path.name in REVIEW_RELEASE
    return True


def inventory(output):
    records = {}
    def add(path, role):
        source = ROOT / path
        if not source.is_file() or source.is_symlink() or source.resolve() != source:
            raise FileNotFoundError(f'Expected local regular file at {path.as_posix()}')
        if (source != Path(__file__).resolve() and source.suffix in FORMATS - {'.npz','.pdf','.png'}
                and HIDDEN_TEXT.search(source.read_text(encoding='utf-8'))):
            raise ValueError(f'Remove local machine path from release material {path.as_posix()}')
        records[path.as_posix()] = {'path':path.as_posix(),'bytes':source.stat().st_size,'role':role}

    for path in PAPER.rglob('*'):
        if path.is_file() and path.resolve() != output and public_study_file(path):
            add(path.relative_to(ROOT),'closed_loop_study')
    for name in MANDATORY:
        path = (PAPER / name).relative_to(ROOT).as_posix()
        if path not in records:
            raise FileNotFoundError(f'Missing study input {path}')
    lines = (PAPER/'beams/dependencies.txt').read_text().splitlines()
    inherited = [Path(line.strip()) for line in lines if line.strip()]
    if len(inherited) != 5:
        raise ValueError('Expected the five declared inherited beam dependencies')
    for path in inherited:
        add(path,'inherited_beam_operator')
    for role, entries in STUDY_FILES.items():
        for name in entries:
            add(STUDY/name,role)
    for name in R8_FILES:
        add(R8/name,'laboratory_loader_dependency')
    return [records[path] for path in sorted(records)]


def package(output, records, manifest):
    output.parent.mkdir(parents=True,exist_ok=True)
    with NamedTemporaryFile(prefix='.closed-loop-',suffix='.part',dir=output.parent,delete=False) as stage:
        temporary = Path(stage.name)
    try:
        with ZipFile(temporary,'w',compression=ZIP_DEFLATED,compresslevel=4,allowZip64=True) as archive:
            for row in records:
                archive.write(ROOT/row['path'],row['path'])
                if archive.getinfo(row['path']).file_size!=row['bytes']:
                    raise RuntimeError(f"Source changed size during packaging: {row['path']}")
            archive.writestr(MANIFEST,json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
        temporary.chmod(0o644)
        os.replace(temporary,output)
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=DEFAULT)
    parser.add_argument('--inventory',action='store_true')
    args = parser.parse_args()
    output = args.output.resolve()
    start,cpu = perf_counter(),process_time()
    records = inventory(output)
    manifest = dict(study=PAPER.relative_to(ROOT).as_posix(),python='3.13',
                    source_files=len(records),input_bytes=sum(row['bytes'] for row in records),
                    files_by_role=dict(sorted(Counter(row['role'] for row in records).items())),
                    laboratory_data='Source and derived scalar results; raw instruments and calibrated port caches supplied separately by readers',
                    members=records)
    if not args.inventory:
        package(output,records,manifest)
        manifest['archive_bytes'] = output.stat().st_size
    manifest.update(wall_seconds=perf_counter()-start,cpu_seconds=process_time()-cpu)
    print(json.dumps({k:v for k,v in manifest.items() if k!='members'},ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()

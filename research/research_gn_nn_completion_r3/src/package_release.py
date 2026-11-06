"""Package frozen results without running science or reading array contents."""
from pathlib import Path
import argparse
import json
import zipfile

STUDY = Path(__file__).resolve().parents[1]
REPO = STUDY.parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--part-bytes', type=int, default=1_300_000_000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    files = [p for p in sorted(STUDY.rglob('*')) if p.is_file()
             and '__pycache__' not in p.parts and p.suffix not in ('.zip', '.pyc')
             and not p.is_relative_to(args.output.resolve())]
    dependencies = [REPO/'research/research_z_alg_r4_diversity/src/physics.py',
                    REPO/'research/research_z_alg_r4_diversity/public_data/operators.npz',
                    REPO/'research/research_gn_search_space_r2/src/search.py',
                    REPO/'research/research_gn_search_space_r2/src/score_search.py']
    files.extend(dependencies)
    groups = {'source-and-reports': [], 'data-and-anchors': [],
              'trained-models': [], 'predictions-and-evaluation': []}
    for p in files:
        rel = p.relative_to(STUDY) if p.is_relative_to(STUDY) else None
        group = 'source-and-reports'
        if rel and rel.parts[0] in ('data', 'anchors', 'cache'):
            group = 'data-and-anchors'
        elif rel and (rel.parts[0] == 'models' or p.suffix == '.pt'):
            group = 'trained-models'
        elif rel and rel.parts[0] in ('predictions', 'evaluation'):
            group = 'predictions-and-evaluation'
        groups[group].append(p)
    inventory = []
    for name, paths in groups.items():
        parts, chunk, size = [], [], 0
        for path in paths:
            n = path.stat().st_size
            if chunk and size + n > args.part_bytes:
                parts.append(chunk)
                chunk, size = [], 0
            chunk.append(path)
            size += n
        if chunk:
            parts.append(chunk)
        for i, part in enumerate(parts):
            target = args.output/f'{name}-{i+1:02d}.zip'
            with zipfile.ZipFile(target, 'w', zipfile.ZIP_DEFLATED, compresslevel=1) as z:
                for path in part:
                    z.write(path, str(path.relative_to(REPO)))
            inventory.append(dict(file=target.name, bytes=target.stat().st_size,
                                  files=len(part), group=name))
    manifest = dict(study=str(STUDY.relative_to(REPO)), artifacts=inventory,
                    layout='Extract every part into the same directory, preserving research/ paths.',
                    scope='Saved material arrays, model weights, candidates, costs, code and direct frozen physics dependencies. No recomputation or hash checks.',
                    geometry='Use included data/operators.npz for re-generation; platform fonts may differ, so saved material arrays are authoritative.')
    (args.output/'RELEASE_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps(manifest, indent=2))


if __name__ == '__main__':
    main()

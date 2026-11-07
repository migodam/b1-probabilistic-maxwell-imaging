"""Package frozen R5 bytes and the exact TRAIN/VAL dependencies; no numerical run."""
from pathlib import Path
import argparse
import json
import subprocess
import time
import zipfile
from datetime import datetime, timezone

STUDY = Path(__file__).resolve().parents[1]
REPO = STUDY.parents[1]
R4 = STUDY.parent/'research_gn_nn_preservation_r4'


def members(folder):
    return sorted(p for p in folder.rglob('*') if p.is_file()
                  and '__pycache__' not in p.parts
                  and p.suffix != '.pyc' and p.name != '.DS_Store')


def write_zip(output, paths):
    started = time.perf_counter()
    entries = []
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED,
                         compresslevel=1, allowZip64=True) as archive:
        for path in sorted(set(paths)):
            relative = str(path.relative_to(REPO))
            archive.write(path, relative)
            entries.append({'path': relative, 'bytes': path.stat().st_size})
    return {'asset': output.name, 'bytes': output.stat().st_size,
            'files': len(entries), 'members': entries,
            'packaging_wall_seconds': time.perf_counter()-started}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    records = json.loads((R4/'data/manifests/evaluation.json').read_text())['records']
    selected = set(json.loads((R4/'audit/OPTIONAL_FT_CONTRACT.json').read_text())['parent_ids'])
    records = [r for r in records if r['split'] == 'val'
               or (r['split'] == 'train' and r['parent_id'] in selected)]
    models = members(STUDY/'models') + members(STUDY/'models_noiseaware')
    own = [p for p in members(STUDY)
           if p.suffix != '.zip' and 'artifacts' not in p.relative_to(STUDY).parts]
    frozen = [R4/'data/operators.npz', R4/'data/reference_jacobian.npz',
              R4/'data/reference_jacobian.json', R4/'data/manifests/evaluation.json',
              R4/'audit/OPTIONAL_FT_CONTRACT.json', R4/'audit/FINAL_ENVIRONMENT.json',
              STUDY.parent/'research_gn_green_error_20261106/results/fixed_spectra.npz']
    for record in records:
        frozen += [R4/'data'/record['public_path'], R4/'data'/record['truth_path'],
                   R4/'anchors'/(record['obs_id']+'.npz')]
    # Preserve original Python packages in full, including relative imports and
    # package initializers; do not regenerate operators, projectors or physics.
    packages = ['research_gn_nn_preservation_r4', 'research_gn_nn_completion_r3',
                'research_gn_search_space_r2', 'research_z_alg_r4_diversity',
                'research_z_alg_r4_basin', 'research_z_alg_r3_digits',
                'research_z_alg_r3', 'research_gn_green_error_20261106']
    for package in packages:
        frozen += [p for p in members(STUDY.parent/package) if p.suffix == '.py']
    frozen += members(R4/'cache')
    assets = [write_zip(args.output/'gn-nn-r5-models.zip', models),
              write_zip(args.output/'gn-nn-r5-reproducibility.zip', own),
              write_zip(args.output/'gn-nn-r5-frozen-inputs.zip', frozen)]
    index = {'schema': 'gn-nn-r5-frozen-release-v1',
             'created_utc': datetime.now(timezone.utc).isoformat(),
             'source_commit': subprocess.check_output(['git','rev-parse','HEAD'], cwd=REPO, text=True).strip(),
             'scope': 'DEVELOPMENT TRAIN128 VAL78,11 models,no new TEST',
             'restore': 'Extract reproducibility and frozen-inputs ZIPs into the same repository root, preserving research/... paths. Models ZIP is a small convenience copy.',
             'prior_results_unchanged': True, 'numerical_runs_during_packaging': False,
             'assets': assets}
    (args.output/'gn-nn-r5-release-index.json').write_text(json.dumps(index, indent=2)+'\n')
    print(json.dumps({'assets': [{k:v for k,v in a.items() if k != 'members'} for a in assets],
                      'index': str(args.output/'gn-nn-r5-release-index.json')}, indent=2))


if __name__ == '__main__':
    main()

"""Create an isolated public-only Windows package after generation; never launch.

Explicit source allowlist, no recursive repository/data copy, no sealed_truth or
scoring code. Remote execution_authorized remains false for parent review.
"""
from pathlib import Path, PurePosixPath
import argparse
import json
import shutil
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
PROJECT = 'research/research_max_struct_r4'
# Existing remote runner/vendor live in publication repository, outside worktree.
INV_SLAM = ROOT.parents[5]
REMOTE_SOURCE = INV_SLAM/'research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r4_diversity'
SOURCE_FILES = [
    PROJECT+'/src/confirmation_solver.py',
    PROJECT+'/src/cie_block_absorption.py',
    PROJECT+'/src/cie_onehop.py',
    PROJECT+'/src/cooperative_geometry.py',
    PROJECT+'/experiments/cooperative_paths.py',
    'research/research_max_struct_r2/src/cie_global.py',
    'research/research_max_struct_r2/src/profile_current.py',
    'research/research_max_struct_r2/experiments/run_blind.py',
    'research/research_fbg_r3_experiments/src/physics.py',
    'research/research_fbg_r3_experiments/src/cie_tsom_baseline.py',
    'research/research_fbg_r3_experiments/src/full_maxwell_corrector.py',
    'research/research_z_alg_r4_diversity/src/physics.py',
    'research/research_z_alg_r4_basin/src/rational_rays.py',
]
OPERATOR_KEYS = {'G','E0','S','coords','sources','receivers','k','radius','shape','grid_shape'}
OBSERVATION_KEYS = {'Y','sigma','noise_radius','host','amplitude','material_bounds',
                    'smooth_weight','tv_weight','tv_epsilon','prior_width'}


def check_npz(path, keys):
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != keys:
            raise ValueError('Refuse unexpected or missing public keys: '+str(path))


def package(destination):
    cfg_path = ROOT/'configs/confirmation_frozen.json'
    cfg = json.loads(cfg_path.read_text())
    if not (cfg['status'] == 'FROZEN_BEFORE_GENERATION' and cfg['algorithm_finalized'] is True):
        raise ValueError('Parent-finalized frozen config required')
    generation_cfg = json.loads((ROOT/'audit/confirmation/FROZEN_CONFIG.json').read_text())
    if cfg != generation_cfg:
        raise ValueError('Configuration changed since generation; preserve realization and resolve provenance')
    manifest = json.loads((ROOT/'data_public/confirmation/manifest.json').read_text())
    records = manifest['records']
    if cfg['primary_methods'] != ['onehop','global','native_extended']:
        raise ValueError('Frozen18-task primary method list mismatch')
    if manifest['status'] != 'FRESH_UNSCORED' or len(records) != 6 or len({r['observation_id'] for r in records}) != 6:
        raise ValueError('Six unique fresh public observations required')
    if destination.exists():
        raise FileExistsError('Preserve existing package; use a new destination')
    copies = [(REPO/name, PurePosixPath(name)) for name in SOURCE_FILES]
    copies += [(REMOTE_SOURCE/'remote_runner.py', PurePosixPath('remote_runner.py')),
               (REMOTE_SOURCE/'vendor/threadpoolctl.py', PurePosixPath('vendor/threadpoolctl.py'))]
    snapshot = ROOT/'audit/confirmation/executed_source_snapshot'
    for source,target in copies:
        frozen = snapshot/str(target)
        if not frozen.is_file() or source.read_bytes() != frozen.read_bytes():
            raise ValueError('Runtime source differs from pre-generation snapshot: '+str(target))
    # Preserve vendor license/metadata without copying caches or arbitrary directories.
    metadata = REMOTE_SOURCE/'vendor/threadpoolctl-3.6.0.dist-info'
    copies += [(path, PurePosixPath('vendor')/metadata.name/path.relative_to(metadata).as_posix())
               for path in sorted(metadata.rglob('*')) if path.is_file() and '__pycache__' not in path.parts]
    public_root = ROOT/'data_public/confirmation'
    check_npz(public_root/'operators.npz', OPERATOR_KEYS)
    copies += [(public_root/'operators.npz', PurePosixPath('data_public/confirmation/operators.npz')),
               (public_root/'manifest.json', PurePosixPath('data_public/confirmation/manifest.json'))]
    for row in records:
        ident = row['observation_id']
        if not ident.startswith('obs_') or len(ident) != 20 or any(c not in '0123456789abcdef' for c in ident[4:]):
            raise ValueError('Invalid opaque observation ID')
        if row['path'] != ident+'.npz':
            raise ValueError('Nonlocal public observation path')
        check_npz(public_root/row['path'], OBSERVATION_KEYS)
        copies.append((public_root/row['path'], PurePosixPath('data_public/confirmation')/row['path']))
    for source, target in copies:
        if not source.is_file() or any(word in target.parts for word in ('sealed_truth','scoring','predictions','results')):
            raise ValueError('Refuse source/path: '+str(source))
    destination.mkdir(parents=True)
    included = []
    for source, target in copies:
        dest = destination/str(target)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
        included.append(dict(path=target.as_posix(), bytes=dest.stat().st_size, source=str(source)))
    # Per-method configs contain no family, seed, label, material truth or clean Y.
    config_dir = destination/'configs'
    config_dir.mkdir()
    public_cfg = {key:value for key,value in cfg.items() if key != 'private_shape_families'}
    (config_dir/'confirmation_frozen_PUBLIC.json').write_text(json.dumps(public_cfg,indent=2)+'\n')
    remote = PurePosixPath(cfg['remote_root'])
    tasks = []
    for row in records:
        for method in cfg['primary_methods']:
            ident = row['observation_id']+'__'+method
            tasks.append(dict(task_id=ident, backend='cpu', method=method, config=cfg['methods'][method],
                operators_path=(remote/'data_public/confirmation/operators.npz').as_posix(),
                public_path=(remote/'data_public/confirmation'/row['path']).as_posix(),
                trajectory_path=(remote/'predictions'/row['observation_id']/(method+'.npz')).as_posix(),
                hard_timeout_s=900))
    execution = dict(status='PACKAGED_PENDING_PARENT_EXECUTION_REVIEW',execution_authorized=False,
        solver_module='research.research_max_struct_r4.src.confirmation_solver',
        output_dir=(remote/'results/confirmation').as_posix(),tasks=tasks,
        expected_tasks=18,requested_cpu_workers=4,requested_cuda_workers=0,scientific_wallcap=None,
        timeout_semantics='900s emergency process cap only; timeout is incomplete scientific disposition, no automatic retry/budget expansion')
    (config_dir/'confirmation_manifest.json').write_text(json.dumps(execution,indent=2)+'\n')
    for name in ('configs/confirmation_frozen_PUBLIC.json','configs/confirmation_manifest.json'):
        included.append(dict(path=name,bytes=(destination/name).stat().st_size,source='derived public protocol'))
    (destination/'PACKAGE_FILELIST.json').write_text(json.dumps(dict(status='PUBLIC_ONLY_NO_LAUNCH',files=included,
        excluded=['sealed_truth','scoring code','generator code','old observations','old results','credentials','environment installs'],
        truth_access='NONE',remote_root=str(remote),source_validation='Exact public NPZ key allowlists; explicit source filelist; no hashes'),indent=2)+'\n')
    print(json.dumps(dict(status='PUBLIC_ONLY_PACKAGE_PREPARED_NOT_LAUNCHED',path=str(destination),tasks=18,execution_authorized=False)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    package(args.output.resolve())

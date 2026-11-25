#!/usr/bin/env python3
"""Package the experiment outputs and the exact project-owned replay inputs."""

import argparse
from collections import Counter
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZIP_STORED, ZipFile

import numpy as np

HERE = Path(__file__).resolve().parent
GUIDED = HERE.parent
ROOT = HERE.parents[3]
STUDY = GUIDED.parent.relative_to(ROOT)
OUTPUT = ROOT / '.artifacts/z_alg_guided_v4/release/Z_ALG_Guided_Evidence.zip'
FLOW = Path('Flow Matching/FlowMultiDecode_Local')
FLOW_RUN = FLOW / 'results/phase_b/runs/main/seed_11/material'
FLOW_TEST = FLOW / 'results/phase_b/data/main/test'
FLOW_INDEX = (0, 37, 111)
WORKTREES = Path('research/delegated/max_struct_r1')
SRAP = Path('research/srap_cie_array_study_20261121')
GENERAL6 = Path('research/delegated/z_alg_r2_publication/repository/research/research_z_alg_r2')
HIGH = Path('research/z_alg_secant_tasks/highdim')
OUTCOME_SUFFIXES = {'.csv', '.json', '.md', '.npy', '.npz', '.pdf', '.png', '.py', '.tex', '.tsv'}
ALREADY_COMPRESSED = {'.npz', '.pdf', '.png', '.pt'}
CHAMBER_FIELDS = {'material', 'complex_born_prediction', 'train_mask',
                  'held_transmitter_mask', 'held_receiver_mask'}

HISTORICAL_LINKS = {
    'unified': 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/optical-inertia-unified',
    'v3': 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/optical-inertia-ip-v3',
    'validation': 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-optical-inertia-validation-v1',
    'cie': 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-xudong-cie-20261124',
    'general6': 'https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_z_alg_r2/',
}
CHAMBER_ORIGINS = {
    'Institut Fresnel': 'https://doi.org/10.1088/0266-5611/21/6/S09',
    'University of Manitoba': 'https://home.cc.umanitoba.ca/~lovetrij/EMILab/mwt_data.html',
}

TASK_FILES = ('REPRODUCE.md', 'RESEARCH_REPORT_ZH.md')
TASK_DIRS = ('benchmark', 'boxes', 'families', 'figures', 'pool', 'routing', 'theory', 'threshold')
CORE_RECORDS = {
    'theory/check_forms.py', 'theory/check_forms.json',
    'families/nested.py', 'families/nested.json', 'families/same_site.py', 'families/same_site.json',
    'families/low_amplitude_site_01.npz', 'families/low_amplitude_site_02.npz',
    'pool/cie_pair/paired.py', 'pool/cie_pair/direct.py', 'pool/cie_pair/direct_diagnostic.py',
    'pool/cie_pair/summary.json', 'pool/cie_pair/per_object.csv',
    *(f'pool/cie_pair/{arm}/{artifact}'
      for arm in ('results', 'corrected_results', 'direct_results')
      for artifact in ('online.json', 'offline.json', 'predictions.npz', 'offline_ports.npz')),
    'threshold/confirmation/predictions.npz', 'threshold/confirmation/offline.json',
    'threshold/confirmation/public/cases.json',
    'threshold/same_stop.py', 'threshold/same_stop.json', 'threshold/same_stop.md',
}
STUDY_INPUTS = {
    'core': 'certificate.py inertia.py',
    'synthetic/public': 'fine_geometry.npz coarse_green.npz pw360_10x20.npz',
    'submission/curve': 'search.py curve_inputs.py noise_curve.csv',
    'physical_revision/array': 'screen.py array_scan.csv',
    'physical_revision/born': 'extinction.py run.py fixed_port_weights.npz',
    'physical_revision/born/controls': 'arrays.npz results.json',
    'physical_revision/born/ports': 'results.json',
    'ip_unified/numerics/feasibility/strict_endpoints': 'control_disk075__pw360_32x32__pair2pct.json',
    'ip_unified/numerics/loss': '''
        run.py finite_0_inputs.npz finite_1_inputs.npz
        finite_0_state_verified.json finite_1_state_verified.json
        finite_0_verified.json finite_1_verified.json finite_0_directional.json
    ''',
    'ip_unified/numerics/native': 'u100_best.json required_currents.json',
    'validation/srap': 'public_observations.npz online.json cases.csv',
    'validation/max_struct': 'public_inputs.npz paired_results.csv summary.json public_provenance.json',
    'fresnel': 'physical.py',
    'results/fresnel': 'SENSITIVITY.csv',
    'validation/manitoba': 'physical.py',
    'validation/transfer': 'symmetric_margins.csv',
}
R8_CODE = '''
    research_fbg_r3_experiments/src/__init__.py
    research_fbg_r3_experiments/src/cie_tsom_baseline.py
    research_fbg_r3_experiments/src/physics.py
    research_max_struct_r2/src/cie_global.py
    research_max_struct_r4/src/cie_onehop.py
    research_max_struct_r6/src/fft_tsom_ops.py
    research_max_struct_r8/src/fresnel_adapter.py
    research_max_struct_r8/src/paper_physics.py
    research_z_alg_r4_diversity/src/__init__.py
    research_z_alg_r4_diversity/src/physics.py
    research_z_alg_r4_diversity/src/solver.py
'''


def task_files():
    files = {GUIDED / name for name in TASK_FILES if (GUIDED / name).is_file()}
    for name in TASK_DIRS:
        for path in (GUIDED / name).rglob('*'):
            rel = path.relative_to(GUIDED)
            if path.is_file() and path.suffix in OUTCOME_SUFFIXES and all(
                not piece.startswith('.') and piece not in ('__pycache__', 'private')
                for piece in rel.parts
            ):
                if rel.parts[:2] == ('routing', 'laboratory') and path.suffix == '.npz':
                    with np.load(path, allow_pickle=False) as chamber:
                        if set(chamber.files) != CHAMBER_FIELDS:
                            raise ValueError(f'Unexpected arrays in measured-output artifact {rel}')
                files.add(path)
    return {path.relative_to(ROOT) for path in files}


def historical_files():
    paths = {STUDY / parent / name for parent, names in STUDY_INPUTS.items() for name in names.split()}
    for directory in ('submission/curve/observations', 'physical_revision/array/observations'):
        p = STUDY / directory
        paths.update(f.relative_to(ROOT) for f in (ROOT / p).glob('*') if f.suffix in {'.json', '.npz'})
    paths.update(f.relative_to(ROOT) for f in (ROOT / STUDY / 'synthetic/evaluation').glob('optfresh_[0-9][0-9].npz'))
    paths.update({
        HIGH / 'public_vie_inputs.npz', HIGH / 'observation_radii.json',
        GENERAL6 / 'src/physics.py', GENERAL6 / 'experiments/GENERAL6_CONTRACT.json',
        GENERAL6 / 'experiments/general6_public/GENERAL6_00.npz',
        GENERAL6 / 'results/general6_GENERAL6_00_general_quotient20.json',
        Path('research/limited_aperture_routed_dr_20261120/run_la_study.py'),
        *(SRAP / name for name in ('engine.py', 'legacy_la.py', 'data2d.py')),
        *(Path('research/z_alg_xudong_cie/minpeak_native') / name for name in ('search.py', 'cie_seed.py')),
        *(WORKTREES / 'r8_worktree/research' / name for name in R8_CODE.split()),
        WORKTREES / 'r13_worktree/research/research_z_alg_r4_diversity/src/physics.py',
        WORKTREES / 'r13_worktree/research/research_z_alg_r4_diversity/src/__init__.py',
        WORKTREES / 'r13_worktree/research/research_z_alg_r4_diversity/src/solver.py',
        WORKTREES / 'r13_worktree/research/research_max_struct_r10/src/lowrank_denominator.py',
        WORKTREES / 'r14_worktree/research/research_fbg_r3_experiments/src/cie_tsom_baseline.py',
        FLOW / 'README.md', FLOW / 'requirements-local.txt', FLOW / 'configs/protocol.json',
        FLOW / 'configs/physics_protocol.json',
        *(FLOW / 'src' / name for name in ('common.py', 'physics.py', 'physics_learning.py',
                                           'physics_dc.py', 'models.py')),
        *(FLOW_RUN / name for name in ('flow/checkpoint.pt', 'flow/done.json',
                                      'probe/checkpoint.pt', 'latent_stats.npz')),
        FLOW_TEST / 'manifest.json',
    })
    with np.load(ROOT / STUDY / 'validation/srap/public_observations.npz', allow_pickle=False) as data:
        paths.update(SRAP / 'data/truth' / (str(name) + '.npz') for name in data['ids'])
    return paths


def masked_flow_array(name):
    original = np.load(ROOT / FLOW_TEST / name, allow_pickle=False)
    subset = np.full_like(original, np.nan)
    subset[list(FLOW_INDEX)] = original[list(FLOW_INDEX)]
    stream = BytesIO()
    np.save(stream, subset, allow_pickle=False)
    return stream.getvalue()


def selected_npz(name, keys):
    # A .npz is a collection of .npy members. Keep their exact represented
    # coefficients without transforming the array contents.
    stream = BytesIO()
    with ZipFile(ROOT / name) as original, ZipFile(stream, 'w', compression=ZIP_DEFLATED, compresslevel=3) as out:
        for key in keys:
            item = key + '.npy'
            out.writestr(item, original.read(item))
    return stream.getvalue()


def package(output):
    new = task_files()
    new_names = {path.relative_to(GUIDED.relative_to(ROOT)).as_posix() for path in new}
    missing_records = sorted(CORE_RECORDS - new_names)
    if missing_records:
        raise FileNotFoundError('Missing experiment records:\n' + '\n'.join(missing_records))
    history = historical_files()
    chosen = {
        STUDY / 'ip_unified/numerics/native/u100_best.npz': ('negative_vec', 'negative_eig', 'center'),
        STUDY / 'ip_unified/numerics/feasibility/strict_endpoints/control_disk075__pw360_32x32__pair2pct.npz':
            ('Y', 'material'),
        HIGH / 'results/r10_twod_0001.npz': ('Y_archive',),
    }
    with np.load(ROOT / STUDY / 'validation/max_struct/public_inputs.npz', allow_pickle=False) as data:
        for name in data['cases']:
            chosen[WORKTREES / 'r14_worktree/research/research_max_struct_r14/sealed_truth/independent24' /
                   (str(name) + '.npz')] = ('material_forward', 'Y_clean')
    generated = {name: selected_npz(name, keys) for name, keys in chosen.items()}
    for name in ('y.npy', 'chi.npy'):
        generated[FLOW_TEST / name] = masked_flow_array(name)
    manifest = {
        'entry_point': str(GUIDED.relative_to(ROOT) / 'REPRODUCE.md'),
        'historical_release_links': HISTORICAL_LINKS,
        'measured_provider_links': CHAMBER_ORIGINS,
        'new_artifacts': sorted(v.as_posix() for v in new),
        'core_experiment_records': sorted(CORE_RECORDS),
        'archived_inputs': sorted(v.as_posix() for v in history),
        'repacked_npz_arrays': {v.as_posix(): keys for v, keys in chosen.items()},
        'flow_test_index_subset': FLOW_INDEX,
        'flow_omitted_rows_representation': 'NaN rows at original omitted indices. Required indices retain exact archived values.',
        'laboratory_outputs': sorted(CHAMBER_FIELDS),
    }
    generated[GUIDED.relative_to(ROOT) / 'benchmark/dependency_manifest.json'] = (
        json.dumps(manifest, indent=2) + '\n'
    ).encode()
    missing = sorted(str(p) for p in new | history if not (ROOT / p).is_file())
    if missing:
        raise FileNotFoundError('\n'.join(missing))
    output = output.resolve()
    output.parent.mkdir(exist_ok=True, parents=True)
    partial = output.with_suffix('.zip.part')
    paths = new | history
    with ZipFile(partial, 'w', compression=ZIP_DEFLATED, compresslevel=3) as archive:
        for name in sorted(paths, key=Path.as_posix):
            archive.write(ROOT / name, str(name),
                          compress_type=ZIP_STORED if name.suffix in ALREADY_COMPRESSED else ZIP_DEFLATED)
        for name, content in sorted(generated.items(), key=lambda item: item[0].as_posix()):
            archive.writestr(str(name), content,
                             compress_type=ZIP_STORED if name.suffix in ALREADY_COMPRESSED else ZIP_DEFLATED)
    partial.replace(output)
    types = Counter(path.suffix for path in paths | generated.keys())
    print(f'{output} | {sum(types.values())} members | {types[".npz"]} array archives'
          f' | {output.stat().st_size:,} bytes', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    package(parser.parse_args().output)

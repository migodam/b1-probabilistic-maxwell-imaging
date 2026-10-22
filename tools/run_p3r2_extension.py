"""Resume independent data-scale learning and exact learned-proposal controls.

The preceding physics queue owns all solves until its explicit closure marker.
No automatic paid replay is performed after an error or an existing output.
"""
from pathlib import Path
import os, sys, time, json, subprocess

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from p3r2.core import OUT, write


def run(name, arguments):
    with (OUT/(name+'.log')).open('x') as log:
        child = subprocess.run([sys.executable, *arguments], cwd=ROOT,
            stdout=log, stderr=subprocess.STDOUT)
    if child.returncode:
        write(OUT/'EXTENSION_ERROR.json', dict(task=name, exit_code=child.returncode,
            no_automatic_replay=True, raw_outputs_preserved=True))
        raise RuntimeError('Extension stopped at '+name)


def main():
    marker = OUT/'ROUND2_PHYSICS_COMPLETE.json'
    while not marker.exists():
        if (OUT/'ROUND2_ERROR.json').exists():
            raise RuntimeError('Previous physics queue stopped; inspect its error')
        time.sleep(15)
    if (OUT/'EXTENSION_COMPLETE.json').exists():
        raise FileExistsError('No duplicate extension campaign')
    # Low-cost direct-risk scale study is independent of posterior qualification.
    run('A_training32768', ['tools/run_p3r2_accounted.py', '--stage', 'A',
        '--job', 'supervised_material_risk_scale2', '--max-wall', '5400', '--',
        sys.executable, 'tools/train_p3r2_risk.py',
        '--train', str(OUT/'A/data/train32768.npz'),
        '--validation', str(OUT/'A/data/validation.npz'),
        '--target', str(OUT/'configs/TARGET.json'),
        '--output', str(OUT/'A/training32768'), '--sizes', '32768',
        '--seeds', '11', '22', '33', '--models', 'mean', 'ridge',
        'mlp_small', 'mlp_medium', '--epochs', '100', '--max-seconds', '360'])
    # Same FULL posterior, equal paid transition counts, two distinct frozen frames.
    for frame, selected in (
        ('informative', OUT/'A/density_training/SELECTED.json'),
        ('riskweak', OUT/'A/weak_density_training/SELECTED.json')):
        run('D2_corrected_'+frame, ['tools/run_p3r2_learned_correction.py',
            '--job', 'D2_corrected_'+frame, '--selected', str(selected),
            '--chains', '2', '--warmup', '128', '--production', '256',
            '--rhs-cap', '131072', '--max-wall', '3600'])
    write(OUT/'EXTENSION_COMPLETE.json', dict(status='COMPLETE',
        dataset_parent_count=32768, corrected_proposals='two frozen 10D frames, FULL66 MH target',
        fresh_locked_truth_read=False))


if __name__ == '__main__':
    main()

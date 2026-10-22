"""Execute the fixed fresh-parent study once, with truth opened only at the end.

Method/model freeze is prepared separately. An interrupted or failed child is
preserved and stops the sequence, rather than silently replaying paid states.
"""
from pathlib import Path
import json, subprocess, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from p3r2.core import OUT, write


def run(label, arguments):
    log = OUT/('fresh_'+label+'.log')
    write(OUT/'FRESH_CHECKPOINT.json', dict(status='RUNNING', task=label,
        command=arguments, utc_unix=time.time()))
    with log.open('x') as handle:
        result = subprocess.run([sys.executable, *arguments], cwd=ROOT,
            stdout=handle, stderr=subprocess.STDOUT)
    if result.returncode:
        write(OUT/'FRESH_CHECKPOINT.json', dict(status='STOPPED_ERROR',
            task=label, exit_code=result.returncode, log=str(log),
            paid_outputs_preserved=True, automatic_replay=False))
        raise RuntimeError('Fresh sequence stopped at '+label)
    print(label+' COMPLETE', flush=True)


def analysis(label, script):
    run(label, ['tools/run_p3r2_accounted.py', '--stage', 'L', '--job', label,
        '--max-wall', '3600', '--', sys.executable, script])


def main():
    manifest = json.loads((OUT/'METHODS_FROZEN_FOR_LOCKED.json').read_text())
    if manifest.get('status') != 'FROZEN' or not manifest.get('include_matched_point176'):
        raise ValueError('Complete frozen protocol, including the ordinary-cost control, required')
    if (OUT/'FRESH_CHECKPOINT.json').exists():
        raise FileExistsError('Review existing fresh sequence before any manual continuation')
    run('generate', ['tools/generate_p3r2_data.py', 'locked', '--count', '256'])
    run('anchors', ['tools/run_p3r2_anchors.py', 'locked', '--count', '256',
        '--native', '96', '--probit', '64', '--matched', '16'])
    analysis('locked_inference', 'tools/infer_p3r2_locked.py')
    for entry in manifest['trajectory_jobs']:
        args = ['tools/run_p3r2_geometry.py', '--job', entry['job'],
            '--parents', *entry['parents'], '--arms', *entry['arms'],
            '--anchor-directory', str(OUT/'L/contexts/cache'),
            '--split', 'fresh_eval', '--steps', '8', '--max-wall', '1200',
            '--seed', str(entry['seed'])]
        if 'dynamic_bandit' in entry['arms']:
            args += ['--policy', str(OUT/'C/policy/FROZEN_POLICY.json')]
        if 'dynamic_mbrl' in entry['arms']:
            args += ['--model-policy', str(OUT/'C/policy/FROZEN_MODEL_POLICY.json')]
        run(entry['job'], args)
    analysis('locked_trajectory_integration', 'tools/integrate_p3r2_fresh_trajectories.py')
    run('locked_decision_forward', ['tools/evaluate_p3r2_decision_forward.py'])
    # Both online families have completed before either offline evaluator starts.
    analysis('locked_risk_evaluation', 'tools/evaluate_p3r2_locked.py')
    analysis('locked_trajectory_evaluation', 'tools/evaluate_p3r2_fresh_trajectories.py')
    analysis('locked_offline_geometry', 'tools/analyze_p3r2_locked_geometry.py')
    write(OUT/'FRESH_CHECKPOINT.json', dict(status='COMPLETE', parents=256,
        independent_locked_cohort=True, no_online_truth_access=True))


if __name__ == '__main__':
    main()

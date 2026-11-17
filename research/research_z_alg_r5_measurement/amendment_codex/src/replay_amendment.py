"""Replay exposed Codex amendment only; preserve every published output."""
from complete_quotient import *
import argparse,shutil,subprocess
ap=argparse.ArgumentParser();ap.add_argument('--output-name',required=True);args=ap.parse_args()
if Path(args.output_name).name!=args.output_name:raise ValueError('A new sibling directory name is required')
out=BASE.parent/args.output_name;out.mkdir(exist_ok=False)
for name in ['src','inputs','private']:
 if name=='inputs':
  (out/name).mkdir()
  for f in (BASE/name).glob('T*.npz'):shutil.copy2(f,out/name/f.name)
 else:shutil.copytree(BASE/name,out/name,ignore=shutil.ignore_patterns('__pycache__'))
# Original frozen inference is input to exposed amendment designs; never relabel it new work.
shutil.copytree(BASE/'results',out/'results');shutil.copytree(OUT/'src',out/'amendment_codex/src',ignore=shutil.ignore_patterns('__pycache__'));(out/'amendment_codex/results').mkdir();(out/'amendment_codex/figures').mkdir()
steps=['run_complete_batch.py','verify_saved_certificates.py','test_amendment.py','noise_coverage.py','structured_noise.py','fold_diagnostic.py','nearest_fold.py','minimax_fold.py','certify_fold.py','constrained_source.py','known_host_probe.py','certify_source_gap.py','present_results.py','plot_host_revision.py']
for script in steps:subprocess.run([sys.executable,str(out/'amendment_codex/src'/script)],check=True)
save(out/'amendment_codex/results/REPLAY_SCOPE.json',dict(status='EXPOSED_REPLAY_NOT_NEW_BLIND',copied_original_results_are_not_new_experiments=True,source=str(BASE)))
print(out)

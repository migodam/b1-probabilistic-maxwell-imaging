"""Replay exposed stored targets in a new sibling directory; never overwrite."""
from common import *
import argparse,shutil,subprocess
ap=argparse.ArgumentParser();ap.add_argument('--output-name',required=True);args=ap.parse_args()
if Path(args.output_name).name!=args.output_name:raise ValueError('Use a sibling directory name, not a path')
out=ROOT.parent/args.output_name
out.mkdir(exist_ok=False)
shutil.copytree(ROOT/'src',out/'src',ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy(ROOT/'EXPERIMENT_PROTOCOL.md',out/'EXPERIMENT_PROTOCOL.md')
(out/'inputs').mkdir();(out/'results').mkdir();(out/'figures').mkdir()
shutil.copytree(ROOT/'private',out/'private')
for f in (ROOT/'inputs').glob('T*.npz'):shutil.copy(f,out/'inputs'/f.name)
shutil.copy(ROOT/'inputs/cases.json',out/'inputs/cases.json')
shutil.copy(ROOT/'results/generation_receipt.json',out/'results/generation_receipt.json')
save(out/'inputs/EXPOSED_REPLAY.json',dict(status='EXPOSED_REPLAY_NOT_NEW_BLIND',source=str(ROOT),generation_receipt='Historical generation cost copied for provenance; no fresh generation charged in actual replay process'))
for script in ['debug.py','test_formulae.py','freeze.py','run_gate1.py','qualify_predictions.py','freeze_predictions.py','score_v2.py','account_and_plot.py']:
 subprocess.run([sys.executable,str(out/'src'/script)],check=True)
print('Exposed replay complete:',out)

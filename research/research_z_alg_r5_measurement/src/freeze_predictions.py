from common import *
import datetime
for case in ['T0','T1','T2','T3','T4','T5']:assert (ROOT/'results'/case/'summary.json').exists()
files=list((ROOT/'results').rglob('*.json'))+list((ROOT/'results').rglob('*.npz'))+list((ROOT/'inputs').glob('*.npz'))+[ROOT/'src/score.py',ROOT/'EXPERIMENT_PROTOCOL.md']
original=json.loads((ROOT/'inputs/solver_freeze.json').read_text());changes={f:digest(ROOT/f) for f,h in original['files'].items() if digest(ROOT/f)!=h};assert not changes,changes
save(ROOT/'inputs/prediction_freeze.json',dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),all_inference_complete=True,GT_materials_not_read_by_inference=True,scoring_rule='RMS .1; raw per-root errors are coverage only; online estimates selected without GT',files={str(f.relative_to(ROOT)):digest(f) for f in files},solver_changed_files=changes))
print('All predictions, retention decisions and scoring code frozen. GT may now be read by scorer.')

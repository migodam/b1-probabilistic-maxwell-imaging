"""After final scoring, export material-only replay targets without secrets."""
import datetime
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def main():
    required=['results/dev_initial_scores.json','results/holdout_scores.json','results/confirm_v7_scores.json',
              'audit/HOLDOUT_V6_FINAL_FREEZE.json','audit/CONFIRM_V7_FINAL_FREEZE.json']
    for path in required:
        if not (ROOT/path).exists():raise RuntimeError('Required completed score/freeze absent: '+path)
    confirm=json.loads((ROOT/'results'/'confirm_v7_scores.json').read_text())
    if not confirm.get('truth_read_after_prediction_freeze'):raise RuntimeError('Confirm not scored after freeze')
    destination=ROOT/'truth_after_freeze'
    if destination.exists():raise FileExistsError(destination)
    rows=[]
    for split in ('dev','holdout','confirm'):
        for source in sorted((ROOT/'private_truth'/split).glob('*.npz')):
            # Deliberately select only m; never read the seed field.
            with np.load(source,allow_pickle=False) as private:material=private['m'].copy()
            target=destination/split/source.name;target.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(target,m=material)
            with np.load(target,allow_pickle=False) as check:
                if check.files!=['m'] or not np.array_equal(check['m'],material):raise RuntimeError('Material-only export failed')
            target.chmod(0o444)
            rows.append(dict(split=split,case_id=source.stem,file=str(target.relative_to(ROOT)),fields=['m'],bytes=target.stat().st_size))
    text='''# Disclosed material targets for replay

These targets were disclosed only after the corresponding frozen predictions
were scored. They are for reproducibility and are not hidden or new blind
targets. DEV is exploratory; V6 holdouts became development evidence for the
separately frozen V7 revision; V7 used four fresh NEWCONFIRM12 targets.

Every NPZ contains exactly one array, `m`. Target generation seeds remain in
the ignored local `private_truth/` directory and are not included here. Public
optimization seeds in method configs are separate and remain public.

Replay from the R3 directory, writing a new output path:

```sh
python experiments/score_predictions.py --split confirm --freeze-manifest audit/CONFIRM_V7_FINAL_FREEZE.json --truth-root truth_after_freeze --disclosed-truth-replay --output /tmp/r3_v7_replay_scores.json
```

The replay flag permits changed file modification times after a Git clone and
checks the recorded prediction byte lengths. It explicitly labels the result
as disclosed-truth replay; it cannot establish fresh blindness. It does not
change frozen predictions or run a physics solver.
'''
    (destination/'README.md').write_text(text)
    receipt=dict(disclosed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                 publication_role='Material-only disclosed-truth replay, not new blind evidence',
                 original_private_truth_retained_local=True,target_generation_seed_values_read=False,
                 target_generation_seed_values_disclosed=False,files=rows,
                 final_confirm_prediction_freeze='audit/CONFIRM_V7_FINAL_FREEZE.json',
                 final_confirm_scores='results/confirm_v7_scores.json',new_physics_rhs=0)
    (ROOT/'audit'/'DISCLOSED_TRUTH_PACKAGE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(material_only_files=len(rows),target_generation_seeds_disclosed=False,new_physics_rhs=0)))


if __name__=='__main__':main()

"""Explicit same-holdout adapter for the pre-GT budget-fill supplement.

The sixteen primary frozen files are untouched. The aliases reuse the exact
same nine measurements and private truth; no new target is generated.
"""
import argparse,json
from pathlib import Path
import score_predictions as frozen_scorer
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','score']);a=ap.parse_args()
 if not (ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_HOLDOUT.json').exists():raise SystemExit('Original 27 must freeze before budget-filled prediction freeze')
 for base in ['public_data','private_truth']:
  alias=ROOT/base/'budget_filled';target=ROOT/base/'holdout'
  if alias.exists() and alias.resolve()!=target.resolve():raise RuntimeError('Unexpected existing split alias')
  if not alias.exists():alias.symlink_to('holdout',target_is_directory=True)
 receipt=dict(split='budget_filled',base_split='holdout',same_nine_measurements=True,new_target_generation=False,source_freeze_timing='supplement frozen after main Y generation, before all GT scoring',primary='original holdout predictions remain primary and unchanged',budget_rule='only remaining 90 seconds minus original wall time, original online incumbent; all methods eligible',scorer='unchanged primary frozen experiments/score_predictions.py invoked with explicit same-data aliases')
 (ROOT/'audit/EXPERIMENT_BUDGET_FILLED_DATA_ALIAS.json').write_text(json.dumps(receipt,indent=2)+'\n')
 if a.action=='freeze':frozen_scorer.freeze('budget_filled')
 else:
  if not (ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_BUDGET_FILLED.json').exists():raise SystemExit('All 27 filled predictions must freeze before any scoring')
  frozen_scorer.score('budget_filled')
if __name__=='__main__':main()

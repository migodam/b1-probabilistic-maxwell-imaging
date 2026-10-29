"""Single-method adapter for the pre-GT ordinary data-adjoint control.

Primary frozen sources remain untouched; the control reuses the same nine Y.
"""
import argparse,json
from pathlib import Path
import score_predictions as frozen_scorer
ROOT=Path(__file__).resolve().parents[1]
def main():
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['freeze','score']);a=ap.parse_args()
 for marker in ['audit/EXPERIMENT_PREDICTIONS_FREEZE_HOLDOUT.json','audit/EXPERIMENT_PREDICTIONS_FREEZE_BUDGET_FILLED.json']:
  if not (ROOT/marker).exists():raise SystemExit('Original and budget-filled predictions must freeze first')
 for base in ['public_data','private_truth']:
  alias=ROOT/base/'data_adjoint';target=ROOT/base/'holdout'
  if alias.exists() and alias.resolve()!=target.resolve():raise RuntimeError('Unexpected split alias')
  if not alias.exists():alias.symlink_to('holdout',target_is_directory=True)
 frozen_scorer.METHODS=['data_adjoint'];frozen_scorer.LABELS['data_adjoint']='Data adjoint L-BFGS + Full-J'
 receipt=dict(split='data_adjoint',base_split='holdout',same_nine_measurements=True,new_target_generation=False,supplement='ordinary data-adjoint matched optimizer attribution control',source_freeze_timing='after Y generation and before any holdout GT scoring',frozen_primary_scorer_unchanged=True,methods_override='in-memory single-method adapter only; metrics unchanged')
 (ROOT/'audit/EXPERIMENT_DATA_ADJOINT_DATA_ALIAS.json').write_text(json.dumps(receipt,indent=2)+'\n')
 if a.action=='freeze':frozen_scorer.freeze('data_adjoint')
 else:
  if not (ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_DATA_ADJOINT.json').exists():raise SystemExit('All9 control predictions must freeze before any GT scoring')
  frozen_scorer.score('data_adjoint')
if __name__=='__main__':main()

"""Zero-physics adversarial IO/qualification invariants, not Maxwell evidence."""
import json
import tempfile
from pathlib import Path
import numpy as np
from offline_score import load_predictions, public_npz, score, cloud_score, read, write, DEFAULT_CONFIG


def main():
    checks={};cfg=read(DEFAULT_CONFIG)
    with tempfile.TemporaryDirectory(prefix='r12_c_io_qa_') as td:
        p=Path(td);cloud=p/'CLOUD.npz';np.savez(cloud,materials=np.zeros((1,576),complex),phi=np.zeros(1),qualified=np.ones(1,bool))
        group=dict(dataset='N',parent='anonymous',method='anonymous',cloud_path=str(cloud),roi=list(range(256)))
        manifest=p/'manifest.json';manifest.write_text(json.dumps(dict(predictions_sealed=False,truth_read=False,groups=[group])))
        try:load_predictions(manifest)
        except RuntimeError:checks['unsealed_predictions_refused_before_truth_IO']=True
        contaminated=p/'contaminated.npz';np.savez(contaminated,material=np.zeros(576),xgt=np.ones(66))
        try:public_npz(contaminated)
        except ValueError:checks['GT_bearing_prediction_rejected_by_keys']=True
        manifest.write_text(json.dumps(dict(predictions_sealed=True,truth_read=False,groups=[group])))
        try:score(manifest,p/'NONEXISTENT_PRIVATE_TRUTH_MANIFEST.json',p/'score')
        except RuntimeError as error:checks['digit_score_without_main_blind_read_refused_before_truth_mapping']=('blind' in str(error))
        manifest.write_text(json.dumps(dict(predictions_sealed=True,truth_read=False,expected_parents={'N':['anonymous','missing']},groups=[group])))
        try:load_predictions(manifest)
        except ValueError:checks['missing_parent_not_silently_excluded']=True
        empty=cloud_score([],np.zeros(576),dict(roi=np.arange(256)),cfg)
        checks['empty_cloud_counts_as_coverage_failure']=not empty['coverage']['0.075'] and empty['nearest_material_MSE'] is None
        checks['empty_cloud_has_no_invented_nearest_penalty']=empty['nearest_error_right_censored'] is True
    if len(checks)!=6 or not all(checks.values()):raise RuntimeError('adversarial IO invariant failed: '+str(checks))
    receipt=dict(status='PASS',checks=checks,new_physics_calls=0,real_GT_files_opened=0,
         scope='synthetic IO/empty-cloud invariants only; no physical/scientific acceptance inference')
    write(Path(__file__).resolve().parent/'GT_ISOLATION_ADVERSARIAL_QA.json',receipt)
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()

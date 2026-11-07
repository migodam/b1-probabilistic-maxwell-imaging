"""R4 public evaluation reuses the unchanged R3 full-VIE corrector."""
from research.research_gn_nn_completion_r3.src.evaluate import (
    evaluate_arm, evaluate_control, evaluate_full45_baseline, save_result,
)

def evaluate_independent_baseline(public,common_anchor,baseline_anchor,baseline_receipt,method,noise):
    receipt,arrays=evaluate_full45_baseline(public,common_anchor,baseline_anchor,baseline_receipt,noise=noise)
    receipt.update(arm=method,control=method,cost_scope='INDEPENDENT_COLD45_PLUS_PAID_VERIFICATION_IN55_ENVELOPE',
                   budget={'wall_seconds':55.,'max_LU':600,'max_complex_RHS':36000},
                   online_selection='own independent baseline endpoint if physically feasible; common strict GN is comparison only',
                   control_stopping='independent45s method; no extra optimization once its method has stopped')
    receipt['wall_overshoot_seconds']=max(0.,receipt['wall_seconds']-55.)
    return receipt,arrays


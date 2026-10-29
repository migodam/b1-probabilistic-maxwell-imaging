"""One actual support-birth shape verifies the inversion/cost-ledger interface."""
from __future__ import annotations
import os
for n in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[n]='1'
import sys
from pathlib import Path
import numpy as np
D=Path(__file__).resolve().parent
sys.path.insert(0,str(D.parent/'shared'));sys.path.insert(0,str(D))
import core
from controlled_solver import GlobalControlledSolver


def main():
    s=core.load_setup()
    z=np.load(core.R10/'followup_shape_reuse/R10BI_FRESH_003_t2.npz')
    a=z['fraction'][s['ROI']];delta=complex(z['delta'][0])
    y=np.load(core.R10/'results/candidates_fresh10/inputs/online/R10BI_FRESH_003.npz')['observation']
    ledger=core.Ledger(D/'API_COST_LEDGER.jsonl',context=dict(task='real support-birth API compatibility'),cap=80)
    model=GlobalControlledSolver(s['C'],s['O'],s['E'],s['Fh'],ledger=ledger)
    result=model.evaluate(a,delta,need_jac=True,y=y)
    reference=core.ExactSolver(s['C'],s['O'],s['E'],s['Fh'],ledger=ledger).evaluate(a,delta,need_jac=True,y=y)
    core.write(D/'API_VERIFICATION.json',dict(status='COMPLETE',shape='R10BI_FRESH_003_t2',
        source_major_J_occ_shape=list(result['J_occ'].shape),J_delta_shape=list(result['J_delta'].shape),
        F_shape=list(result['F'].shape),finite=all(np.isfinite(result[k]).all() for k in ('F','J_occ','J_delta')),
        F_difference_Fnorm=float(np.linalg.norm(result['F']-reference['F'])),
        J_occ_relative_difference=float(np.linalg.norm(result['J_occ']-reference['J_occ'])/np.linalg.norm(reference['J_occ'])),
        J_delta_relative_difference=float(np.linalg.norm(result['J_delta']-reference['J_delta'])/np.linalg.norm(reference['J_delta'])),
        own_conservative_RHS=model.ledger.summary()['RHS_conservative_total'],
        external_conservative_RHS=ledger.rhs,external_LU=ledger.lu,
        own_diagnostics=result['diagnostics'],derivative_status='EXACT_DERIVATIVE_NOT_APPROXIMATE_DERIVATIVE',
        interpretation='interface and source-major compatibility audit, not scientific acceptance',GT_access=False,RHS576=0))
    ledger.close()


if __name__=='__main__':main()

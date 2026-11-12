"""Public-only finite-pool minimum-change retention, before truth scoring.

No candidate generation, optimization or new forward. Main min-Phi remains
unchanged. This alternative explicitly prefers the raw anchor among already-paid
discrepancy-feasible materials, with primary fallback if none are available.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,time
from pathlib import Path
from seal_results import read,validate_seal


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--seal',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    c,w=time.process_time(),time.perf_counter();root=a.root.resolve();seal=read(a.seal);validate_seal(seal,require_authorization=False)
    if Path(seal['root']).resolve()!=root:raise ValueError('Wrong root')
    import numpy as np
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);cfg=seal['manifest_snapshot'];campaign=cfg['campaign'];rows=[]
    for task in cfg['tasks']:
        cpu,wall=time.process_time(),time.perf_counter();folder=root/'RESULTS'/campaign/task['task_id'];reason=None;selected=None;distance=None;fallback=True;n=0
        try:
            r=read(folder/'RECEIPT.json');history=read(folder/'HISTORY.json');cap=float(task.get('CPU_cap',cfg['protocol']['CPU_cap']))
            with np.load(folder/'CANDIDATES.npz',allow_pickle=False) as z:materials=z['m'].copy()
            raw=r.get('raw_candidate');primary=r.get('primary_candidate');selected=primary
            if raw is None:reason='RAW_UNAVAILABLE_PRIMARY_FALLBACK'
            else:
                anchor=materials[raw]
                qualified=[h for h in history if h.get('legal') and h.get('eligible') and h.get('noise_feasible') and h['CPU_seconds']<=cap]
                n=len(qualified)
                if qualified:
                    def key(h):
                        difference=materials[h['candidate']]-anchor
                        return (float(np.mean(difference*difference)),h['Phi'],h['candidate'])
                    chosen=min(qualified,key=key);selected=chosen['candidate'];distance=float(np.sqrt(key(chosen)[0]));fallback=False
                    reason='MIN_RMS_CHANGE_TO_RAW_IN_LEGAL_PAID_BALL_POOL'
                else:reason='NO_LEGAL_PAID_BALL_CANDIDATE_PRIMARY_FALLBACK'
        except Exception as error:reason='UNAVAILABLE:'+repr(error)
        rows.append(dict(task_id=task['task_id'],case_id=task['case_id'],method=task['method'],selection_candidate=selected,fallback_to_primary=fallback,raw_anchor_RMS_distance=distance,qualified_candidates=n,reason=reason,CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,truth_access='NONE'))
    result=dict(status='PUBLIC_SELECTION_NO_GT_NO_NEW_SOLVES',campaign=campaign,prediction_seal=str(a.seal.resolve()),rule='minimum RMS material change from raw among legal paid discrepancy-feasible candidates; Phi then candidate index break exact ties; primary fallback',prior_semantics='Explicit raw-anchor preference, different from minPhi or minR. This is selection/retention, not improved same-objective optimization.',cost_policy='All original proposal/corrector costs retained; selection scan is extra CPU. No measured skipped-corrector or full-wave saving.',tasks=rows,CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-w,truth_access='NONE',new_Maxwell_LU=0,new_optimization=0)
    (out/'SELECTIONS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(campaign=campaign,tasks=len(rows),CPU_seconds=result['CPU_seconds'],truth_access='NONE',new_Maxwell_LU=0)))
if __name__=='__main__':main()

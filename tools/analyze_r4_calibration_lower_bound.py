"""Post-hoc standard data-processing bounds, using closed fixed-draw counts.

No posterior teacher, new forward, truth array, training or model selection.
The bound is conditional on the frozen discrete P1 law and ideal IID draws;
it is not a bit-level or continuous-PDE certificate.
"""
from pathlib import Path
import argparse
import json
import math
import time

from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'research/p3l_r4_conditional_transport'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scores',type=Path,default=S/'evidence_summary/FULL_PROBABILITY.json')
    p.add_argument('--roster',type=Path,default=S/'configs/FROZEN_DENSITY_BOTH_COHORTS.json')
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    start=time.perf_counter();cpu=time.process_time()
    scores=json.loads(a.scores.read_text());roster=json.loads(a.roster.read_text())
    assert scores['status']=='OFFLINE_LIKELIHOOD_SCORE_COMPLETE' and scores['no_residual_selection']
    assert roster['fixed_sample_index']==0 and len(roster['methods'])==26
    family_size=2*len(roster['methods']);tail=.025;delta=.05/family_size
    rows=[]
    for name,r in scores['Gamma256'].items():
        n=r['planned_parents'];k=r['PIT_above_975']
        assert r['fixed_sample_index']==0 and r['Gamma_shape']==256 and r['Gamma_scale']==1
        assert n in (512,256) and n==r['finite_parents'] and not r['failed_parents']
        assert type(k) is int and 0<=k<=n
        lower=float(beta.ppf(delta,k,n-k+1)) if k else 0.
        tv=max(0.,lower-tail)
        kl=(tail*math.log(tail/lower)+(1-tail)*(math.log1p(-tail)-math.log1p(-lower))) if lower>tail else 0.
        rows.append(dict(group=name,n=n,tail_event_count=k,model_event_probability_lower=lower,
                         average_conditional_material_TV_lower=tv,
                         average_conditional_material_forward_KL_lower_nats=kl))
    assert len(rows)==8
    def identity(path):
        v=path.stat();return dict(path=str(path.resolve()),bytes=v.st_size,mtime_ns=v.st_mtime_ns)
    result=dict(status='POST_HOC_INTERPRETATION_ONLY',confidence_simultaneous=.95,
                Bonferroni_family_size=family_size,measured_model_cohort_rows=len(rows),
                target_event_probability=tail,fixed_threshold=.975,fixed_sample_index=0,
                conditional_target='original FULL material posterior pushforward, averaged over prior-predictive Y',
                assumptions=['fixed discrete P1 forward/noise law','frozen original26 model roster',
                             'ideal IID parent and generator random streams','negligible classification error from numerical forward'],
                density_required=False,material_Lebesgue_density_assumed=False,
                CFM_scope='actual frozen NFE32 sample law, not continuous ODE density',
                roster_CFM_NFE=roster['CFM_NFE'],rows=rows,
                new_Maxwell_RHS=0,new_GT_arrays_opened=0,model_changes=False,
                original_scores_modified=False,original_gate_status_changed=False,
                identity=identity(a.scores),frozen_roster_identity=identity(a.roster),
                wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('x') as stream:json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print(json.dumps(dict(status=result['status'],rows=len(rows),family_size=family_size,
                         wall_seconds=result['wall_seconds'],RHS=0)))


if __name__=='__main__':main()

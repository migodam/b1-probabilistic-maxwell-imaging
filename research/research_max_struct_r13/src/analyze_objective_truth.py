"""Post-seal GT-objective diagnostic using generation-paid coarse predictions.

No optimization, new Maxwell forward, parameter selection or GT-selected roots.
This distinguishes an objective that prefers a different legal material from
failure to find a known lower-objective reference. Cross-grid reference is the
coarse cell average, not an exactly representable continuous true material.
Noisy least-squares truth is ordinarily not an exact minimizer: this fact alone
is not novel or evidence of a serious failure. Interpret the material deviation
and objective margin jointly; arbitrarily small noise fitting is not a basin
breakthrough or a harmful-selection certificate.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
import argparse,json,sys,time
from pathlib import Path
from seal_results import read,validate_seal


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--seal',type=Path,required=True);p.add_argument('--campaign',action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    cpu,wall=time.process_time(),time.perf_counter();root=a.root.resolve();seal=read(a.seal)
    validate_seal(seal,True)
    if Path(seal['root']).resolve()!=root or seal['campaign']!='confirmation_core48':raise ValueError('Authorized complete core48 seal required')
    import numpy as np
    from score_results import table
    sys.path.insert(0,str(root.parents[1]))
    from research.research_fbg_r3_experiments.src.physics import Regularizer
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);truth_rows=[];lookup={};rows=[]
    tasks=seal['manifest_snapshot']['tasks'];case_inputs={t['case_id']:t['input_file'] for t in tasks}
    for cid,filename in sorted(case_inputs.items()):
        with np.load(root/filename.replace('\\','/'),allow_pickle=False) as z:
            Y=z['Y'].copy();delta=float(z['delta']);lower=float(z['lower']);upper=float(z['upper']);shape=tuple(z['grid_shape']);crossgrid=int(z['forward_side'])!=int(z['inverse_side'])
        with np.load(root/'sealed_truth/fresh48'/f'{cid}.npz',allow_pickle=False) as z:
            m=z['material_coarse'].copy();prediction=z['Y_coarse_truth'].copy()
        r=prediction-Y;D=.5*float(np.vdot(r,r).real/np.vdot(Y,Y).real);R=float(Regularizer(shape).evaluate(m)[0]);ratio=float(np.linalg.norm(r)/delta)
        row=dict(case_id=cid,matched_discrete_model=not crossgrid,GT_reference='coarse cell-average material; generation-paid inverse-model prediction',GT_legal=bool(np.min(m)>=lower and np.max(m)<=upper),GT_D=D,GT_R=R,GT_Phi=D+R,GT_residual_over_delta=ratio,GT_discrepancy_feasible=ratio<=1,qualification='OFFLINE_DIAGNOSTIC; not available to solver/selection. Cross-grid reference does not exactly parameterize the fine generative target.')
        truth_rows.append(row);lookup[cid]=row
    for campaign in a.campaign:
        score_root=root/'analysis'/campaign;scoring=read(score_root/'SCORING_RECEIPT.json');validate_seal(read(scoring['authorized_seal']),True)
        for row in read(score_root/'SCORES.json'):
            cid=row['case_id']
            if cid not in lookup or row['view']=='initial':continue
            gt=lookup[cid];available=row.get('available') and row.get('Phi') is not None
            item=dict(campaign=campaign,case_id=cid,method=row['method'],view=row['view'],available=bool(available),matched_discrete_model=gt['matched_discrete_model'],GT_Phi=gt['GT_Phi'],GT_D=gt['GT_D'],GT_R=gt['GT_R'],GT_discrepancy_feasible=gt['GT_discrepancy_feasible'],candidate_Phi=row.get('Phi'),candidate_D=row.get('D'),candidate_R=row.get('R'),candidate_NRMSE=row.get('NRMSE'),candidate_boxrange_NRMSE=row.get('boxrange_NRMSE'),candidate_discrepancy_feasible=row.get('discrepancy_feasible'),candidate_Phi_minus_GT=row['Phi']-gt['GT_Phi'] if available else None,candidate_Phi_below_GT_reference=bool(row['Phi']<gt['GT_Phi']) if available else None,interpretation='A legal material with numerically lower Phi disproves GT being the minimizer of this finite noisy objective at these points; it is not exact-data nonuniqueness or a global-optimum proof. A higher candidate Phi leaves a known better reference; it does not alone prove a strict local minimum.')
            rows.append(item)
    table(out,'TRUTH_OBJECTIVE_REFERENCES',truth_rows);table(out,'CANDIDATE_MINUS_TRUTH_REFERENCE',rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f=out/'FIGURES';f.mkdir()
    for campaign in a.campaign:
        group=[r for r in rows if r['campaign']==campaign and r['view']=='primary' and r.get('available')]
        if not group:continue
        methods=sorted({r['method'] for r in group});fig,axes=plt.subplots(1,2,figsize=(13,5))
        for ax,matched in zip(axes,(True,False)):
            for method in methods:
                q=[r for r in group if r['method']==method and r['matched_discrete_model']==matched]
                ax.scatter([r['candidate_Phi_minus_GT'] for r in q],[r['candidate_NRMSE'] for r in q],s=20,label=method,alpha=.75)
            ax.axvline(0,color='gray',ls=':');ax.set_xscale('symlog',linthresh=1e-5);ax.set_xlabel('Paid primary Phi minus coarse truth-reference Phi');ax.set_ylabel('Amplitude NRMSE');ax.set_title('Matched discrete model' if matched else 'Cross-grid model mismatch');ax.grid(alpha=.2)
        axes[0].legend(fontsize=7,ncol=2);fig.suptitle(campaign+': objective preference and material recovery are separate');fig.tight_layout();fig.savefig(f/(campaign+'_objective_truth.png'),dpi=160);fig.savefig(f/(campaign+'_objective_truth.pdf'));plt.close(fig)
    (out/'RECEIPT.json').write_text(json.dumps(dict(status='POST_SEAL_OFFLINE_DIAGNOSTIC',truth_access='AUTHORIZED_OFFLINE_ONLY',new_Maxwell_LU=0,new_adjoint_RHS=0,new_optimization=0,source_predictions='generation-paid Y_coarse_truth',CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,cases=len(truth_rows),comparison_rows=len(rows)),indent=2)+'\n')
    print(json.dumps(dict(cases=len(truth_rows),comparison_rows=len(rows),new_Maxwell_LU=0)))
if __name__=='__main__':main()

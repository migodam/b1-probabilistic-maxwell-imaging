"""Post-freeze fresh scorer. This module is never imported by inverse workers."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import json,csv,sys
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]))
from research.research_z_alg_r3_digits.experiments.score_predictions import metrics
METHODS=['trf','data_lbfgs','globaldual','uniform_bridge','uniform_best_only','uniform_multistart']
COLORS=['#0072b2','#009e73','#d55e00','#882255','#44aabb','#332288']

def main():
    folder=ROOT/'results/fresh_v3'
    manifest=json.loads((folder/'MANIFEST.json').read_text())
    mapping=json.loads((ROOT/'truth_after_freeze/fresh/mapping.json').read_text())
    metadata={r['case_id']:r for r in mapping['records']}
    rows=[];images={};truths={};curves={}
    for job in manifest['jobs']:
        case,method=job['case_id'],job['method']
        result=json.loads((folder/(case+'_'+method+'.json')).read_text())
        with np.load(ROOT/'truth_after_freeze/fresh'/(case+'.npz'),allow_pickle=False) as f:
            truth,mask=f['m'].ravel(),f['mask'].ravel()
        truths[case]=truth
        top=result['top1'];material=np.asarray(top['m'])
        images[(case,method)]=material
        scores=metrics(material,truth,mask)
        candidate_scores=[dict(index=i,origin=c['origin'],objective=c['objective'],relative_residual=c['relative_residual'],
           **metrics(np.asarray(c['m']),truth,mask)) for i,c in enumerate(result['candidates'])]
        data_gate=top['relative_residual']<=.001
        material_gate=scores['RMSE']<=1.5 and scores['IoU']>=.7
        row=dict(case_id=case,digit=metadata[case]['digit'],method=method,status=result['status'],
            scope='FRESH_AFTER_V3_FREEZE',**scores,relative_data_residual=top['relative_residual'],
            objective=top['objective'],K=len(candidate_scores),candidate_scores=candidate_scores,
            oracle_bestK_RMSE=min((c['RMSE'] for c in candidate_scores),default=scores['RMSE']),
            online_selection='minimum original full data loss + common smooth/TV; no GT',
            oracle_metric='posthoc GT ranking, not online recovery',
            data_gate_pass=data_gate,material_gate_pass=material_gate,joint_gate_pass=data_gate and material_gate,
            all_576_unknowns_released=result['all_unknowns_released'],
            wall_seconds=result['wall_seconds'],cpu_seconds=result['cpu_seconds'],
            costcounts=result['costcounts'],gate=result.get('gate'),
            lu_budget_overshoot=result.get('lu_budget_overshoot'),rhs_budget_overshoot=result.get('complex_rhs_budget_overshoot'),
            raw_proposal_count=len(result['raw_proposals']),corrected_proposal_count=sum(bool(p.get('corrected')) for p in result['raw_proposals']))
        rows.append(row)
        trace=result['trajectory']
        seconds=np.asarray([p['seconds'] for p in trace]);objective=np.asarray([p['objective'] for p in trace]);residual=np.asarray([p['relative_residual'] for p in trace])
        rhs=np.asarray([p['counts']['actual_complex_RHS'] for p in trace]);rmse=np.asarray([np.sqrt(np.mean((np.asarray(p['m'])-truth)**2)) for p in trace])
        # At each resource budget select by the online objective among completed trials.
        best=np.minimum.accumulate(objective);indices=[];which=0
        for i in range(len(objective)):
            if objective[i]<=objective[which]:which=i
            indices.append(which)
        selected=np.asarray(indices)
        curves[(case,method)]=dict(seconds=seconds,rhs=rhs,RMSE=rmse,online_RMSE=rmse[selected],online_residual=residual[selected],objective=objective)
        np.savez_compressed(folder/(case+'_'+method+'.quality_curve.npz'),**curves[(case,method)])
        del result,trace
    summary=dict(schema='r4-fresh-scores-v3',scored_utc=datetime.now(timezone.utc).isoformat(),
       predictions_closed_utc=manifest['completed_utc'],truth_first_read_after_predictions_closed=True,
       scope='3 fresh targets, 6 frozen methods, unchanged 576-cell VIE',rows=rows,
       joint_success_by_method={m:sum(r['joint_gate_pass'] for r in rows if r['method']==m) for m in METHODS},
       scoring_physics_cost={'LU':0,'complex_RHS':0},
       known_limits=['3 targets only','finite candidate archive, not branch or posterior coverage','same-grid synthetic discrete model','data and material gates separate'])
    (folder/'SCORES.json').write_text(json.dumps(summary,indent=2)+'\n')
    scalar=[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in rows]
    with (folder/'SCORES.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=sorted(set().union(*(r.keys() for r in scalar))));writer.writeheader();writer.writerows(scalar)
    cases=sorted(metadata,key=lambda c:metadata[c]['digit'])
    fig,axes=plt.subplots(3,7,figsize=(18,9),layout='constrained')
    for i,case in enumerate(cases):
        for j,key in enumerate(['truth']+METHODS):
            ax=axes[i,j];m=truths[case] if key=='truth' else images[(case,key)]
            image=ax.imshow(m.reshape(24,24),vmin=.5,vmax=20,cmap='viridis',interpolation='nearest')
            ax.set_xticks([]);ax.set_yticks([])
            if i==0:ax.set_title('Ground truth' if key=='truth' else key,fontsize=10)
            if j==0:ax.set_ylabel('Fresh digit '+str(metadata[case]['digit']))
            else:
                r=next(r for r in rows if r['case_id']==case and r['method']==key)
                ax.set_xlabel('RMSE %.2f | IoU %.2f\nr %.2g | K %d'%(r['RMSE'],r['IoU'],r['relative_data_residual'],r['K']),fontsize=8)
    fig.colorbar(image,ax=axes.ravel().tolist(),shrink=.5,label='Real material contrast; common [0.5,20]')
    fig.suptitle('R4 frozen V3 fresh confirmation: all three targets and all six methods\nOriginal-physics top1; GT used only after all eighteen predictions closed',fontsize=13)
    for ext in ('png','pdf'):fig.savefig(ROOT/'figures'/('fresh_reconstruction.'+ext),dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(15,10),layout='constrained')
    for i,case in enumerate(cases):
        for method,color in zip(METHODS,COLORS):
            c=curves[(case,method)]
            axes[i,0].plot(c['seconds'],c['online_RMSE'],color=color,label=method,lw=1)
            axes[i,1].plot(c['rhs'],c['online_RMSE'],color=color,lw=1)
            axes[i,2].semilogy(c['seconds'],c['online_residual'],color=color,lw=1)
        axes[i,0].set_ylabel('Fresh digit %s\nGT material RMSE (posthoc)'%metadata[case]['digit'])
        axes[i,0].set_xlabel('Paid solver wall seconds');axes[i,1].set_xlabel('Paid complex RHS');axes[i,1].set_ylabel('GT material RMSE (posthoc)')
        axes[i,2].set_xlabel('Paid solver wall seconds');axes[i,2].set_ylabel('Relative data residual')
        for ax in axes[i]:ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('Resource-matched prefixes: material of the lowest online objective seen so far\nThese curves use GT only to score a selection already made from measurements',fontsize=12)
    for ext in ('png','pdf'):fig.savefig(ROOT/'figures'/('fresh_quality_vs_budget.'+ext),dpi=180)
    plt.close(fig)
    print(json.dumps({'joint_success_by_method':summary['joint_success_by_method'],'rows':[{k:r[k] for k in ('digit','method','RMSE','IoU','relative_data_residual','K','oracle_bestK_RMSE')} for r in rows]},indent=2))

if __name__=='__main__':main()

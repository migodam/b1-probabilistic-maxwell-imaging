"""Authorized postseal finite-noise/raw-shadow diagnostics and figure layout QA.
No forward/projection/optimizer; main scores and frozen sources are retained.
"""
from pathlib import Path
import json,csv,time,shutil
import numpy as np
from scipy.stats import chi2
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'RESULTS/offline_scoring_v1'
FIG=ROOT/'FIGURES/offline_scoring_v1'
METHODS=('cie_single','cie_dual','global_dual','a_plus','fullG19','dr','multistart','nn_direct')
LABELS=dict(zip(METHODS,('CIE','Dual','Global','A+','G19','DR','Multi','NN')))
def read(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).open('x').write(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(name,rows):
    dump(OUT/(name+'.json'),rows)
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (OUT/(name+'.csv')).open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def main():
    start=time.process_time();rows=read(OUT/'per_object.json');raw=read(OUT/'raw_DR_fixed_selections.json')
    # No GT material or observation was regenerated. Check the actual saved
    # simulator Yclean against the exact same public observation/radius.
    gt=[];gmap={}
    for stage in ('dev6','confirm24'):
        cases=sorted({r['case_id'] for r in rows if r['stage']==stage})
        for case in cases:
            with np.load(ROOT/'sealed_truth'/stage/(case+'.npz'),allow_pickle=False) as z:
                t=z['m'].copy();clean=z['Y_clean'].copy() if 'Y_clean' in z.files else None
            r=dict(stage=stage,case_id=case,saved_Yclean_available=clean is not None,GT_exact_box_legal=bool(t.min()>=.02 and t.max()<=20),new_FW=0)
            if clean is not None:
                with np.load(ROOT/'data_public'/stage/(case+'.npz'),allow_pickle=False) as z:Y=z['Y'];sigma=float(z['sigma'])
                radius=float(sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2));norm=float(np.linalg.norm(Y-clean))
                r.update(GT_raw_noise_norm=norm,delta_CN99=radius,GT_noise_over_bound=norm/radius,
                    GT_noise_feasible=bool(norm<=radius),GT_relative_data_residual=norm/float(np.linalg.norm(Y)),
                    qualification='VERIFIED_FROM_SAVED_YCLEAN_NO_FW')
            else:r.update(GT_noise_feasible=None,qualification='UNVERIFIED_NO_SAVED_YCLEAN_NO_EXTRA_SOLVE')
            gt.append(r);gmap[case]=r
    table('GT_saved_noise_qualification',gt)
    joined=[]
    for r in rows:
        g=gmap[r['case_id']];a=dict(r,GT_noise_qualification=g['qualification'],GT_noise_feasible=g['GT_noise_feasible'],GT_exact_box_legal=g['GT_exact_box_legal'])
        a['finite_noise_two_distinct_legal_candidates']=bool(r.get('available') and r.get('legal_exact') and r.get('noise_feasible') and g['GT_noise_feasible'] and g['GT_exact_box_legal'] and r.get('RMSE',0)>0)
        a['ambiguity_scope']='At most finite-noise observational ambiguity of two candidates; not exactnonuniqueness or impossibility for everyprior'
        joined.append(a)
    table('per_object_noise_qualification',joined)
    diag=[];comparisons=[];rawcontrols=[]
    for stage in ('dev6','confirm24'):
        ids=sorted({r['case_id'] for r in rows if r['stage']==stage})
        for case in ids:
            p=next(r for r in joined if r['stage']==stage and r['case_id']==case and r['method']=='dr' and r['selector']=='primary')
            s=next(r for r in joined if r['stage']==stage and r['case_id']==case and r['method']=='dr' and r['selector']=='secondary')
            rawcase={r['view']:r for r in raw if r['stage']==stage and r['case_id']==case}
            views=dict(primary=p,secondary=s,**rawcase)
            for name,r in views.items():
                a=dict(stage=stage,case_id=case,view=name,available=r.get('available',False),GT_noise_feasible=gmap[case]['GT_noise_feasible'])
                for metric in ('NRMSE','RMSE','relativeL2','IoU','Phi','R','relative_data_residual','noise_feasible','descriptive_success','physical_material_class'):
                    a[metric]=r.get(metric)
                if r.get('available'):
                    a.update(inherited_NRMSE_pass=r['NRMSE']<=.15,inherited_IoU_pass=r['IoU']>=.7,
                        failure_components='PASS' if r['descriptive_success'] else '+'.join(x for x,test in [('NRMSE',r['NRMSE']>.15),('IoU',r['IoU']<.7),('box',not r['legal_inherited'])] if test),
                        finite_noise_two_distinct_legal_candidates=bool(r['legal_exact'] and r['noise_feasible'] and gmap[case]['GT_noise_feasible'] and gmap[case]['GT_exact_box_legal'] and r['RMSE']>0))
                diag.append(a)
            c=dict(stage=stage,case_id=case,secondary_available=s['available'])
            for name,r in views.items():
                for metric in ('NRMSE','relativeL2','IoU','Phi','R','relative_data_residual','noise_feasible'):
                    c[name+'_'+metric]=r.get(metric)
            for name,r in rawcase.items():
                for metric in ('NRMSE','relativeL2','IoU','Phi','R','relative_data_residual'):
                    c['primary_minus_'+name+'_'+metric]=p[metric]-r[metric]
                c[name+'_to_primary_noise_transition']=str(r['noise_feasible'])+'->'+str(p['noise_feasible'])
            comparisons.append(c)
            for view,r in views.items():
                for method in ('cie_dual','fullG19'):
                    # Raw diagnostic shadows compare to finalprimary controls;
                    # this is explicitly not an equal-spend newexperiment.
                    selector=view if view in ('primary','secondary') else 'primary'
                    control=next(x for x in rows if x['stage']==stage and x['case_id']==case and x['method']==method and x['selector']==selector)
                    d=dict(stage=stage,case_id=case,DR_view=view,control_method=method,control_selector=selector,
                        DR_available=r.get('available',False),control_available=control['available'],scope='Postseal object diagnostic; rawproposal andfullcontrol are different paid budgets')
                    for metric in ('NRMSE','relativeL2','IoU','Phi','R','relative_data_residual'):
                        d['DR_minus_control_'+metric]=r.get(metric)-control.get(metric) if r.get('available') and control['available'] else None
                    rawcontrols.append(d)
    table('DR_failure_continuous_diagnosis',diag);table('DR_raw_primary_secondary_comparison',comparisons);table('DR_all_views_vs_Dual_G19',rawcontrols)
    # Figure QA: move pair availability into ylabels instead of covering CIs.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42})
    pairs=read(OUT/'pairs.json');colors={m:plt.get_cmap('tab10')(i) for i,m in enumerate(METHODS)}
    history=FIG/'layout_history';history.mkdir(exist_ok=True)
    for suffix in ('png','pdf'):
        old=FIG/('confirm24_DR_paired_bootstrap_forest.'+suffix)
        shutil.copy2(old,history/(old.stem+'_original_layout.'+suffix))
    def forest(name,controls):
        f,axes=plt.subplots(1,2,figsize=(12,4.6))
        for ax,sel in zip(axes,('primary','secondary')):
            labels=[]
            for y,control in enumerate(controls):
                p=next(x for x in pairs if x['stage']=='confirm24' and x['selector']==sel and x['view']=='endpoint' and {x['left'],x['right']}=={'dr',control})
                vals=(p['NRMSE_difference_mean'],p['NRMSE_difference_mean_low'],p['NRMSE_difference_mean_high'])
                mid,lo,hi=vals if p['left']=='dr' else (-vals[0],-vals[2],-vals[1])
                ax.errorbar(mid,y,xerr=[[max(0,mid-lo)],[max(0,hi-mid)]],fmt='o',color=colors[control],capsize=3)
                labels.append(LABELS[control]+f"  ({p['paired_available']}/24)")
            ax.set_yticks(range(len(controls)),labels);ax.axvline(0,color='gray',ls='--');ax.grid(alpha=.2);ax.set_title(sel)
            ax.set_xlabel('Mean NRMSE difference: DR minus control\nnegative favors DR; paired n in labels')
        f.suptitle('confirm24:10000 material-object bootstrap, seed20261112,95% percentile CI\n'+('All controls; secondary conditional on jointlyavailable' if 'multistart' in controls else 'Focused axis excludes Multi; see full forest for all7controls'))
        f.tight_layout()
        for ext in ('png','pdf'):f.savefig(FIG/(name+'.'+ext),dpi=180,bbox_inches='tight')
        plt.close(f)
    forest('confirm24_DR_paired_bootstrap_forest',[m for m in METHODS if m!='dr'])
    forest('confirm24_DR_paired_bootstrap_forest_focused',[m for m in METHODS if m not in ('dr','multistart')])
    # Continuous per-object values behind the old binary descriptive labels.
    f,axes=plt.subplots(3,1,figsize=(12,8),sharex=True)
    for view,label,color in [('raw_proposal_minPhi','raw proposal minPhi','tab:blue'),('primary','primary','tab:red'),('secondary','secondary','tab:green')]:
        rs=[r for r in diag if r['stage']=='confirm24' and r['view']==view]
        by={r['case_id']:r for r in rs};ids=sorted(by)
        for ax,key in zip(axes,('NRMSE','relativeL2','IoU')):
            ax.plot(range(24),[by[c].get(key) if by[c]['available'] else np.nan for c in ids],'-o',ms=3,label=label,color=color)
            ax.set_ylabel(key);ax.grid(alpha=.2)
    axes[0].legend(ncol=3);axes[2].set_xticks(range(24),[f'{i:02d}' for i in range(24)]);axes[2].set_xlabel('Every confirmation objectID suffix; gaps=unavailable, not removed')
    f.suptitle('DR full24 continuous errors: fixed rawproposal selector vs primary/secondary');f.tight_layout()
    for ext in ('png','pdf'):f.savefig(FIG/('confirm24_DR_continuous_all_objects.'+ext),dpi=180,bbox_inches='tight')
    plt.close(f)
    # Append the explicit qualification; do not rewrite frozen scores.
    confirm_gt=[r for r in gt if r['stage']=='confirm24'];inside=[r for r in diag if r['stage']=='confirm24' and r['view']=='primary' and r.get('noise_feasible') and r.get('descriptive_success') is False]
    pairs_ok=sum(r.get('finite_noise_two_distinct_legal_candidates',False) for r in inside)
    cmp=[r for r in comparisons if r['stage']=='confirm24'];transitions={}
    for r in cmp:
        key=r.get('raw_proposal_minPhi_to_primary_noise_transition');transitions[key]=transitions.get(key,0)+1
    delta=[r['primary_minus_raw_proposal_minPhi_NRMSE'] for r in cmp if 'primary_minus_raw_proposal_minPhi_NRMSE' in r]
    text=['','## 追加：保存的GT噪声资格与连续材料误差（无新FW）','',
        f"确认24的Yclean均已保存；直接核对Y−Yclean与同一CN99半径：{sum(r['GT_noise_feasible'] is True for r in confirm_gt)}/24个GT本次也在球内，其余保留资格差异。逐对象表见GT_saved_noise_qualification。",'',
        f"DR主规则11个球内旧描述未过对象中，NRMSE>.15为{sum(r['NRMSE']>.15 for r in inside)}个、IoU<.7为{sum(r['IoU']<.7 for r in inside)}个；均值NRMSE={np.mean([r['NRMSE'] for r in inside]):.6f}，均值relativeL2={np.mean([r['relativeL2'] for r in inside]):.6f}，relativeL2范围{min(r['relativeL2'] for r in inside):.6f}–{max(r['relativeL2'] for r in inside):.6f}。这主要是超水平集/细节精度未过旧标准，不能把它们统称完全没有成像。",'',
        f"其中{pairs_ok}/{len(inside)}个对象同时满足：GT合法且用已保存Yclean核验在球内，DR主材料也合法、在同球内、与GT不同。这只支持该次有限噪声半径下两个候选都与观测相容的诊断，不是exact nonuniqueness、连续模型不可辨识定理，也不证明所有先验都无法选择GT。",'',
        f"固定raw proposal minPhi→primary的noise状态转移计数：{transitions}；NRMSE下降{sum(d<0 for d in delta)}、上升{sum(d>0 for d in delta)}，最大下降{min(delta):+.6f}、最大上升{max(delta):+.6f}。没有按GT挑某一步rawshadow。",'',
        '全部DR连续诊断/旧标准分项在DR_failure_continuous_diagnosis；raw/primary/secondary并排数值在DR_raw_primary_secondary_comparison；raw与Dual/G19的全部诊断差异在DR_all_views_vs_Dual_G19，注明rawproposal与完整control的预算不同，不能偷换为新等成本比较。']
    with (ROOT/'CONFIRMATION_RESULTS_ZH.md').open('a') as f:f.write('\n'.join(text)+'\n')
    dump(OUT/'ADDITIONAL_DIAGNOSTICS_RECEIPT.json',dict(status='COMPLETE_POSTSEAL_DIAGNOSTICS',new_FW=0,new_algorithms=0,GT='EXPLICIT_ROOT_OFFLINE_AUTHORIZATION',
        GT_confirm_noise_feasible=sum(r['GT_noise_feasible'] is True for r in confirm_gt),DR_inside_old_failure_two_candidate_qualified=pairs_ok,
        process_CPU_seconds_near_write=time.process_time()-start,priority_scope='Singleactualmidrunsnapshotnotwholehistory',old_figure_layout_preserved=True))
    print(json.dumps(dict(status='COMPLETE_GT_SAVED_NOISE_QUALIFICATION',GT_confirm_CN99=sum(r['GT_noise_feasible'] is True for r in confirm_gt),DR_inside_failure_qualified=pairs_ok,new_FW=0)),flush=True)
if __name__=='__main__':main()

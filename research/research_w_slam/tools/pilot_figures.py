#!/usr/bin/env python3
"""Saved-array pilot figures only; no solver import. Run via stage wrapper.

All quantitative summaries average repeats within parent first. Paired bootstrap
is descriptive over twelve parents, seed2026102510,2000 draws; no new gates.
"""
from pathlib import Path
import argparse,csv,json,math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
SEED=2026102510


def jsonread(path):return json.loads(path.read_text()) if path.exists() else None


def number(value):
    try:
        result=float(value)
        return result if np.isfinite(result) else np.nan
    except (ValueError,TypeError):return np.nan


def grouped(records,category,metrics):
    """Keep repeat counts; missing/failed scores remain missing, never zero."""
    groups=defaultdict(list)
    for r in records:groups[(int(r['parent']),str(r[category]))].append(r)
    out=[]
    for (parent,method),rr in sorted(groups.items()):
        row=dict(parent=parent,**{category:method},saved_rows=len(rr),repeat_count=len(set(r.get('repeat',0) for r in rr)))
        for metric in metrics:
            values=np.array([number(r.get(metric)) for r in rr]);valid=np.isfinite(values)
            repeats=defaultdict(list)
            for item,value in zip(rr,values):
                if np.isfinite(value):repeats[item.get('repeat',0)].append(value)
            repeat_means=[np.mean(value) for value in repeats.values()]
            row[metric]=float(np.mean(repeat_means)) if repeat_means else np.nan
            row[metric+'_valid_rows']=int(valid.sum())
            row[metric+'_valid_repeats']=len(repeats)
        out.append(row)
    return out


def paired_bootstrap(rows,category,metric,baseline):
    by=defaultdict(dict)
    for r in rows:
        if r['repeat_count']==2 and r.get(metric+'_valid_repeats',2)==2 and np.isfinite(number(r.get(metric))):by[r[category]][r['parent']]=r[metric]
    result=[]
    for method in sorted(by):
        ids=sorted(set(by[method])&set(by.get(baseline,{})))
        differences=np.array([by[method][i]-by[baseline][i] for i in ids])
        row=dict(method=method,baseline=baseline,metric=metric,paired_parents=len(ids),
                 estimate=float(differences.mean()) if len(ids) else np.nan,
                 descriptive_only=True,bootstrap_draws=2000,seed=SEED)
        if len(ids)==12:
            rng=np.random.default_rng(SEED);draws=differences[rng.integers(0,12,size=(2000,12))].mean(axis=1)
            row.update(low=float(np.quantile(draws,.025)),high=float(np.quantile(draws,.975)),status='DESCRIPTIVE_12_PARENT_PAIRED_BOOTSTRAP')
        else:row.update(low=np.nan,high=np.nan,status='NOT_RUN_INCOMPLETE_12_PARENT_REPEAT_PAIRS')
        result.append(row)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=ROOT)
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args();root=args.root.resolve()
    out=(args.output or root/'figures/pilot_actual').resolve()
    out.mkdir(parents=True,exist_ok=False)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    manifest=[];gaps=[]
    def table(name,rows):
        if not rows:return
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (out/(name+'.csv')).open('w',newline='') as f:
            wr=csv.DictWriter(f,keys);wr.writeheader();wr.writerows(rows)
    def save(fig,name,en,zh,sources):
        if not getattr(fig,'_pilot_manual_layout',False):fig.tight_layout()
        fig.savefig(out/(name+'.png'),dpi=180);fig.savefig(out/(name+'.pdf'));plt.close(fig)
        manifest.append(dict(name=name,caption_en=en,caption_zh=zh,sources=sources))
    def curves(rows,category,metric,ax,label=None):
        for method in sorted(set(r[category] for r in rows)):
            rr=sorted([r for r in rows if r[category]==method],key=lambda r:r['parent'])
            ax.plot([r['parent'] for r in rr],[r.get(metric,np.nan) for r in rr],'o-',ms=3,label=method)
        ax.set_xlabel('Independent parent');ax.set_ylabel(label or metric)
    # A0 is continuous-material optimization; kept separate from finite beliefs.
    a0=jsonread(root/'results/a0/scoring/RESULTS.json')
    if a0:
        records=[]
        for r in a0['rows']:
            records.append(dict(r,parent=r['parent_id'],repeat=r['noise_repeat']))
        metrics=['training_noise_envelope_compatible','holdout_noise_envelope_compatible',
                 'both_noise_envelopes_compatible','pose_error_lambda0','support_iou']
        a0group=grouped(records,'method',metrics);table('A0_PARENT_MEANS',a0group)
        table('A0_RAW_ROWS',a0['rows'])
        table('A0_PAIRED_BOOTSTRAP',[b for m in metrics for b in paired_bootstrap(a0group,'method',m,'FULL_J_TR')])
        fig,axes=plt.subplots(2,2,figsize=(11,7))
        for ax,metric in zip(axes.ravel(),[metrics[0],metrics[1],metrics[3],metrics[4]]):curves(a0group,'method',metric,ax)
        axes[0,0].legend(fontsize=6);save(fig,'A0_parent_outcomes',
            'Continuous A0: repeats averaged within parent. Train/holdout envelope compatibility, pose error and support IoU remain distinct; missing/failed scores are not imputed.',
            '连续A0：先在parent内平均噪声副本；训练/留出兼容、位姿误差和支撑IoU分开。缺失或失败得分不填零。',
            ['results/a0/scoring/RESULTS.json','results/a0/scoring/FULL_TABLE.csv'])
        truth=jsonread(root/'inputs/offline/CONTINUOUS_TRUTHS.json')
        if truth:
            methods=sorted(set(r['method'] for r in records));images=[];gridrows=[]
            def material(support,eps,ell):
                s=np.asarray(support,bool);e=np.broadcast_to(np.asarray(eps,float),s.shape);l=np.broadcast_to(np.asarray(ell,float),s.shape)
                repeat=(6//s.shape[0],24//s.shape[1]);mask=np.repeat(np.repeat(s,repeat[0],0),repeat[1],1)
                er=np.repeat(np.repeat(e,repeat[0],0),repeat[1],1);loss=np.repeat(np.repeat(l,repeat[0],0),repeat[1],1)
                return np.where(mask,er-1,0.),np.where(mask,loss,0.),mask
            fig,axes=plt.subplots(4,1+len(methods),figsize=(3*(1+len(methods)),7),squeeze=False)
            archive={}
            for ip,parent in enumerate((0,1,2,3)):
                label=next((r for r in truth if r['parent_id']==parent),None)
                choices=[('TRUTH',label)]+[(method,next((r for r in records if r['parent']==parent and r['repeat']==0 and r['method']==method and 'selected_theta' in r),None)) for method in methods]
                for im,(name,row) in enumerate(choices):
                    ax=axes[ip,im];key=f'parent{parent}_{name}'
                    if row is None:
                        ax.text(.5,.5,'NOT_SAVED',ha='center',transform=ax.transAxes);ax.set_axis_off();continue
                    if name=='TRUTH':contrast,loss,mask=material(row['support'],row['epsinf'],row['ell'])
                    else:
                        theta=row['selected_theta'];contrast,loss,mask=material(row['selected_support'],theta[0],theta[1])
                    archive[key+'_contrast']=contrast;archive[key+'_ell']=loss;archive[key+'_support']=mask
                    for iy,ix in np.ndindex(contrast.shape):gridrows.append(dict(parent=parent,method=name,repeat='truth' if name=='TRUTH' else 0,
                        row=iy,col=ix,x=-1.5+(ix+.5)*.125,y=.5+(iy+.5)*.125,contrast=contrast[iy,ix],ell=loss[iy,ix],support=int(mask[iy,ix])))
                    m=ax.imshow(contrast,origin='lower',extent=(-1.5,1.5,.5,1.25),vmin=0,vmax=5,aspect='auto',interpolation='nearest')
                    ax.set_title(f'P{parent} {name}',fontsize=7);ax.set_xlabel('x');ax.set_ylabel('y')
            np.savez_compressed(out/'A0_ACTUAL_MATERIAL_GRIDS.npz',**archive)
            table('A0_ACTUAL_MATERIAL_GRIDS',gridrows)
            save(fig,'A0_fixed_materials','Fixed parents0–3, repeat0: actual selected and true wall coefficients on the public24×6 grid; no cartoon geometry or invented wave fields. Selection failures remain blank.',
                 '固定parent0–3、副本0：公共24×6网格上的实际真值与已选墙系数；非示意几何，不补造波场。未保存候选留空。',
                 ['results/a0/scoring/RESULTS.json','inputs/offline/CONTINUOUS_TRUTHS.json'])
    else:gaps.append('A0 scoring RESULTS.json absent; no A0 outcomes/materials invented')
    # Original analysis/v1 had all72 mean layouts invalid. It is retained as an
    # implementation failure, never interpreted as zero physical candidates.
    repaired_path=root/'results/a0/analysis_repaired/RESULTS.json'
    repaired=jsonread(repaired_path)
    valid_repaired=False
    if repaired:
        summary=repaired.get('summary',{});rr=repaired.get('rows',[])
        evaluable=sum(np.isfinite(number(r.get('noiseless_joint_mean_white_distance'))) for r in rr)
        valid_repaired=(len(rr)==72 and evaluable==72 and
                        summary.get('noiseless_pair_evaluable_count')==72 and
                        summary.get('noiseless_pair_diagnosis_complete') is True and
                        summary.get('status')=='OFFLINE_SAVED_ARRAY_ANALYSIS_COMPLETE')
        if valid_repaired:
            table('A0_REPAIRED_ALL72_OFFLINE_MEANS',rr)
            table('A0_REPAIRED_METHOD_STATS',summary.get('method_stats',[]))
            table('A0_REPAIRED_PARENT_MEANS',repaired.get('parent_method_means',[]))
            table('A0_REPAIRED_PAIRED_BOOTSTRAP',repaired.get('paired_bootstrap',[]))
    if not valid_repaired:
        gaps.append('A0 repaired noiseless analysis not complete/evaluable72; original all72-layout-invalid v1 is not a negative physical finding')
    # 3x3 extension: use paid primary NPZ candidate coefficients and offline
    # declared truths. Both epsilon panels share one fixed physical color norm.
    wide=root/'results/a0/wide_extension'
    wide_truth=jsonread(wide/'offline/TRUTHS.json')
    wide_terminal=jsonread(wide/'ONLINE_RESULTS.json')
    if wide_truth and wide_terminal and wide_terminal.get('status')=='ALL_ONLINE_TERMINAL_DURABLE':
        import matplotlib.colors as colors
        from matplotlib.cm import ScalarMappable
        fixed=sorted(wide_truth,key=lambda r:r['parent_id'])
        fig,axes=plt.subplots(len(fixed),2,figsize=(9,2.2*len(fixed)),squeeze=False)
        norm=colors.Normalize(vmin=1.,vmax=6.);archive={};wide_rows=[]
        def physical_grid(support,eps,ell):
            support=np.asarray(support,bool)
            if support.shape not in ((2,2),(3,3),(6,24)):raise ValueError('Unexpected wide material grid')
            e=np.broadcast_to(np.asarray(eps,float),support.shape);l=np.broadcast_to(np.asarray(ell,float),support.shape)
            ry,rx=6//support.shape[0],24//support.shape[1]
            def expand(x):return np.repeat(np.repeat(x,ry,0),rx,1)
            mask=expand(support)
            return np.where(mask,expand(e),1.),np.where(mask,expand(l),0.),mask
        for ip,label in enumerate(fixed):
            parent=int(label['parent_id']);path=wide/f'parent_{parent:02d}_primary.npz'
            actual=[('truth',physical_grid(label['support'],label['epsinf'],label['ell']))]
            if path.exists():
                with np.load(path,allow_pickle=False) as z:theta=z['theta'].copy();support=z['support'].copy()
                actual.append(('selected',physical_grid(support,theta[0],theta[1])))
            else:actual.append(('selected',None));gaps.append(f'Wide parent{parent} paid primary NPZ absent; selected material left blank')
            for col,(name,fields) in enumerate(actual):
                ax=axes[ip,col]
                if fields is None:ax.text(.5,.5,'NOT_SAVED',transform=ax.transAxes,ha='center');ax.set_axis_off();continue
                eps,ell,mask=fields;prefix=f'parent{parent}_{name}'
                archive[prefix+'_epsinf']=eps;archive[prefix+'_ell']=ell;archive[prefix+'_support']=mask
                ax.imshow(eps,origin='lower',extent=(-1.5,1.5,.5,1.25),cmap='viridis',norm=norm,aspect='auto',interpolation='nearest')
                ax.set_title(f'3x3 extension P{parent}: {name} epsilon_inf');ax.set_xlabel('x');ax.set_ylabel('y')
                for iy,ix in np.ndindex(eps.shape):wide_rows.append(dict(parent=parent,kind=name,row=iy,col=ix,
                    x=-1.5+(ix+.5)*.125,y=.5+(iy+.5)*.125,epsinf=eps[iy,ix],ell=ell[iy,ix],support=int(mask[iy,ix])))
        # Reserve an explicit shared colorbar axis so both comparison panels
        # retain equal physical scales and equal plotting areas.
        fig.subplots_adjust(right=.87,hspace=.8,wspace=.3)
        fig._pilot_manual_layout=True
        cax=fig.add_axes([.9,.15,.02,.7]);fig.colorbar(ScalarMappable(norm=norm,cmap='viridis'),cax=cax,label='epsilon_inf (same 1–6 scale)')
        table('A0_WIDE_MATCHED_MATERIAL_GRIDS',wide_rows)
        np.savez_compressed(out/'A0_WIDE_MATCHED_MATERIAL_GRIDS.npz',**archive)
        save(fig,'A0_wide_matched_materials','Saved3x3-extension offline truth and paid selected wall coefficients on the same public24×6 grid. Every panel uses the identical epsilon_inf1–6 color scale; selection/nonconvergence status is unchanged. No new wavefield or PDE.',
             '已保存3x3扩展的离线真值与已计费候选墙系数，位于同一公共24×6网格。全部面板共用epsilon_inf1–6色标；不改写候选/不收敛状态，不生成新波场或PDE。',
             ['results/a0/wide_extension/offline/TRUTHS.json','results/a0/wide_extension/parent_*_primary.npz'])
    else:gaps.append('Optional3x3 extension terminal/truth inputs absent; no matched material panel invented')
    # Finite atom coordinates and joint weights; no wslam module import.
    infer=root/'results/e2e3/inference';gridpath=root/'results/e2e3/finite_bank/online/FINITE_GRID.npz'
    methods=jsonread(infer/'METHOD_RESULTS.json');scores=jsonread(infer/'OFFLINE_SCORES.json')
    postpath=infer/'POSTERIORS.npz'
    if postpath.exists() and gridpath.exists():
        with np.load(gridpath,allow_pickle=False) as z:grid={k:z[k].copy() for k in ('pose','support','epsinf','ell','support_id')}
        with np.load(postpath,allow_pickle=False) as z:
            refs=z['weights'].copy();mw=z['method_weights'].copy();methodnames=z['methods'].astype(str).tolist()
            completed=int(z['completed_parent'])
        # Shape/normalization checks do not silently normalize partial zero rows.
        if refs.ndim!=4 or mw.ndim!=5 or refs.shape[-1]!=len(grid['pose']):raise ValueError('Finite posterior schema mismatch')
        selected=[p for p in (0,3,6,9) if p<=completed and p<refs.shape[0]]
        support=np.asarray(grid['support']).reshape(-1,4);pose=grid['pose'];eps=grid['epsinf'];ell=grid['ell']
        coords=np.column_stack((support*(eps[:,None]-1),support*ell[:,None],pose))
        names=['wall_c0','wall_c1','wall_c2','wall_c3','wall_ell0','wall_ell1','wall_ell2','wall_ell3','target_x','target_y']
        def repeat_mean_covariance(parent,t):
            covariances=[]
            for w in refs[parent,:,t]:
                mu=w@coords;center=coords-mu
                covariances.append((center*w[:,None]).T@center)
            return np.mean(covariances,axis=0)
        beliefrows=[];covrows=[]
        fig,axes=plt.subplots(len(selected),6,figsize=(14,3*len(selected)),squeeze=False) if selected else (None,None)
        for ip,parent in enumerate(selected):
            for t in range(3):
                w=refs[parent,:,t].mean(axis=0)
                if not np.isclose(w.sum(),1.,rtol=0,atol=1e-10):raise ValueError('Incomplete or nonnormalized saved reference')
                occupancy=w@support;mean=w@pose
                for cell,value in enumerate(occupancy):beliefrows.append(dict(parent=parent,prefix=t+1,quantity='occupancy',index=cell,value=value))
                for dim in range(2):beliefrows.append(dict(parent=parent,prefix=t+1,quantity='pose_mean',index=dim,value=mean[dim]))
                pmass=np.zeros((len(np.unique(pose[:,1])),len(np.unique(pose[:,0]))))
                xv,yv=np.unique(pose[:,0]),np.unique(pose[:,1])
                for index,mass in enumerate(w):pmass[np.searchsorted(yv,pose[index,1]),np.searchsorted(xv,pose[index,0])]+=mass
                for iy,ix in np.ndindex(pmass.shape):beliefrows.append(dict(parent=parent,prefix=t+1,quantity='pose_atom_mass',index=iy*len(xv)+ix,value=pmass[iy,ix],x=xv[ix],y=yv[iy]))
                axes[ip,2*t].imshow(occupancy.reshape(2,2),vmin=0,vmax=1,origin='lower',interpolation='nearest')
                axes[ip,2*t].set_title(f'P{parent} prefix{t+1} occupancy')
                # Unequal CC-node spacing: display atoms at their actual coordinates.
                xx,yy=np.meshgrid(xv,yv);m=axes[ip,2*t+1].scatter(xx.ravel(),yy.ravel(),c=pmass.ravel(),s=130,cmap='magma',vmin=0)
                axes[ip,2*t+1].set_title(f'P{parent} prefix{t+1} pose masses');axes[ip,2*t+1].set_xlim(-1.1,1.1);axes[ip,2*t+1].set_ylim(1.5,2.7)
                cov=repeat_mean_covariance(parent,t)
                for i in range(10):
                    for j in range(10):covrows.append(dict(parent=parent,prefix=t+1,row=names[i],column=names[j],covariance=cov[i,j]))
        if fig is not None:
            table('E2_FIXED_REFERENCE_BELIEFS',beliefrows);table('E2_WALL_TARGET_COVARIANCE',covrows)
            save(fig,'E2_fixed_beliefs','Explicit finite-prior reference, fixed parents0/3/6/9: two-repeat mean distributions at prefixes1–3. Pose colors are atom masses at actual CC nodes, not continuous density; no representative selected by outcome.',
                 '显式有限先验参考，固定parent0/3/6/9：前缀1–3的两副本平均分布。位姿颜色为真实CC节点上的原子质量，非连续密度；不按结果选例。',
                 ['results/e2e3/inference/POSTERIORS.npz','results/e2e3/finite_bank/online/FINITE_GRID.npz'])
            fig,axes=plt.subplots(1,len(selected),figsize=(4*len(selected),4),squeeze=False)
            for ip,parent in enumerate(selected):
                cov=repeat_mean_covariance(parent,-1)
                image=axes[0,ip].imshow(cov[:8,8:],aspect='auto',cmap='RdBu_r');axes[0,ip].set_title(f'P{parent}, final prefix');axes[0,ip].set_xticks([0,1],['target x','target y']);axes[0,ip].set_yticks(range(8),names[:8]);fig.colorbar(image,ax=axes[0,ip],shrink=.8)
            save(fig,'E2_wall_target_crosscov','Finite joint posterior wall–target cross covariance, averaged across the two within-parent covariances; units differ between contrast and loss rows. This is not a local-Gaussian or continuous-posterior certificate.',
                 '先计算两份副本各自的有限联合后验墙–目标交叉协方差，再在parent内平均；对比度与损耗行单位不同。非局部高斯或连续后验证书。',
                 ['results/e2e3/inference/POSTERIORS.npz','results/e2e3/finite_bank/online/FINITE_GRID.npz'])
    else:gaps.append('Finite grid or POSTERIORS absent; no E2 belief/covariance figure fabricated')
    if methods:
        flat=[]
        for r in methods:
            m=r['posterior_metrics'];flat.append(dict(parent=r['parent'],repeat=r['repeat'],time=r['time'],method=r['method'],
                action_regret=m['action_regret'],target_cdf_error=max(m['coordinate_cdf_sup'][8:]),event_error=m['event_px_ge_0_error'],
                mode_error=m['mode_mass_max_error_ref_ge_005'],cpu_seconds=r['cpu_seconds'],
                standalone_cpu_seconds=number(r.get('standalone_cpu_seconds')),
                standalone_cpu_source='METHOD_RESULTS.standalone_cpu_seconds' if np.isfinite(number(r.get('standalone_cpu_seconds'))) else 'FALLBACK_OLD_DERIVED_SUM_NOT_CANONICAL',
                shared_likelihood_cpu_seconds=r['shared_likelihood_cpu_seconds'],shared_batch_reference_cpu_seconds=r['shared_batch_reference_cpu_seconds'],
                setup_rhs_attribution=r['setup_rhs_attribution'],actual_new_maxwell_rhs=r['actual_new_maxwell_rhs']))
        table('E3_RAW_METHOD_METRICS',flat)
        final=[r for r in flat if r['time']==3];g=grouped(final,'method',['action_regret','target_cdf_error','event_error','mode_error'])
        table('E3_FINAL_PARENT_METRICS',g)
        table('E3_PAIRED_BOOTSTRAP',[r for metric in ['action_regret','target_cdf_error'] for r in paired_bootstrap(g,'method',metric,'INCREMENTAL_JOINT_SMOOTHER')])
        # Aggregate each method's three updates once per observation. Common
        # likelihood/batch costs remain explicit standalone attribution, not new FE.
        by=defaultdict(list)
        for r in flat:by[(r['parent'],r['repeat'],r['method'])].append(r)
        costs=[]
        bank_final=jsonread(root/'results/e2e3/finite_bank/online/JOB_COST_FINAL.json') or {}
        bank_cost_complete=bank_final.get('complete') is True
        bank_cpu=number(bank_final.get('measured_process_cpu_seconds')) if bank_cost_complete else np.nan
        bank_wall=number(bank_final.get('summed_job_wall_seconds')) if bank_cost_complete else np.nan
        if not bank_cost_complete:gaps.append('Finite-bank JOB_COST_FINAL.json complete=true absent; setupCPU is unavailable, never taken from in-flight COST_RECEIPT/COMPLETE snapshots')
        for (p,r,m),items in sorted(by.items()):
            canonical=all(np.isfinite(x['standalone_cpu_seconds']) for x in items)
            online_cpu=sum(x['standalone_cpu_seconds'] if np.isfinite(x['standalone_cpu_seconds']) else
                           x['cpu_seconds']+x['shared_likelihood_cpu_seconds']+x['shared_batch_reference_cpu_seconds'] for x in items)
            costs.append(dict(parent=p,repeat=r,method=m,method_update_cpu=sum(x['cpu_seconds'] for x in items),
                shared_likelihood_cpu=sum(x['shared_likelihood_cpu_seconds'] for x in items),
                shared_reference_cpu=sum(x['shared_batch_reference_cpu_seconds'] for x in items),
                standalone_online_cpu=online_cpu,
                standalone_cpu_source='CANONICAL_METHOD_RESULTS_FIELD' if canonical else 'FALLBACK_OLD_DERIVED_SUM_NOT_CANONICAL',
                common_shared_bank_setup_rhs=11800,common_shared_bank_setup_factorizations=2950,
                actual_new_maxwell_rhs=sum(x['actual_new_maxwell_rhs'] for x in items),
                common_bank_setup_cpu_seconds=bank_cpu,
                common_bank_setup_summed_job_wall_seconds=bank_wall,
                setup_cpu_source='online/JOB_COST_FINAL.json measured wrapper jobs only' if bank_cost_complete else 'NOT_AVAILABLE_INCOMPLETE_FINAL_JOB_RECEIPT',
                standalone_total_cpu_with_shared_setup=bank_cpu+online_cpu if np.isfinite(bank_cpu) else np.nan,
                attribution='shared bank charged once; standalone alternatives are not additive campaign costs'))
        table('E3_STANDALONE_COST_ATTRIBUTION',costs)
        cg=grouped(costs,'method',['method_update_cpu','shared_likelihood_cpu','shared_reference_cpu','standalone_online_cpu'])
        table('E3_PARENT_COSTS',cg)
        costmean={m:np.mean([x['standalone_online_cpu'] for x in cg if x['method']==m]) for m in set(x['method'] for x in cg)}
        fig,axes=plt.subplots(1,2,figsize=(11,4))
        for method in sorted(set(x['method'] for x in g)):
            rr=[x for x in g if x['method']==method];cpu=costmean[method]
            axes[0].scatter(cpu,np.mean([x['action_regret'] for x in rr]));axes[0].annotate(method,(cpu,np.mean([x['action_regret'] for x in rr])),fontsize=5)
            axes[1].scatter(cpu,np.mean([x['target_cdf_error'] for x in rr]));axes[1].annotate(method,(cpu,np.mean([x['target_cdf_error'] for x in rr])),fontsize=5)
        fallback_used=any(r['standalone_cpu_source']!='CANONICAL_METHOD_RESULTS_FIELD' for r in costs)
        for ax in axes:ax.set_xlabel('Standalone CPU: METHOD_RESULTS field' if not fallback_used else 'Standalone CPU: includes labeled legacy fallback')
        axes[0].set_ylabel('Finite-reference action regret');axes[1].set_ylabel('Target CDF sup error')
        fig.suptitle('All methods also require shared setup: 11800 RHS / 2950 LU',fontsize=10)
        save(fig,'E3_error_vs_CPU','Finite-prior final errors versus canonical METHOD_RESULTS standalone_cpu_seconds, including each method’s actual likelihood/reference dependence; any old derived fallback is explicitly labeled. Common11800RHS/2950LU setup is separate. SetupCPU comes only from complete JOB_COST_FINAL wrapper receipts, never in-flight snapshots. No online acceleration/equal-risk claim.',
             '有限先验最终误差对照METHOD_RESULTS规范standalone_cpu_seconds，包含各方法真实似然/参考依赖；旧派生回退会明确标注。共同11800RHS/2950LU准备另列；准备CPU只读完成的JOB_COST_FINAL wrapper收据，不使用运行中快照。不宣称在线加速或等风险优势。',
             ['results/e2e3/inference/METHOD_RESULTS.json','results/e2e3/finite_bank/online/JOB_COST_FINAL.json'])
        if scores:
            sg=grouped([r for r in scores if r['time']==3],'method',['pose_error','occupancy_brier'])
            table('E3_FINAL_OFFLINE_PARENT_ERRORS',sg)
            table('E3_OFFLINE_PAIRED_BOOTSTRAP',paired_bootstrap(sg,'method','pose_error','INCREMENTAL_JOINT_SMOOTHER'))
            fig,axes=plt.subplots(1,2,figsize=(10,4))
            chosen=[r for r in sg if r['method'] in ('HARD_FROZEN_WALL','INCREMENTAL_JOINT_SMOOTHER')]
            curves(chosen,'method','pose_error',axes[0]);curves(g,'method','action_regret',axes[1]);axes[0].legend(fontsize=7)
            save(fig,'E3_parent_hard_vs_smoother','Final finite-prior offline pose errors and reference-action regret, parent first. Hard freezing and ordinary smoother are explicit controls; all twelve independent parents, not24 independent noise cases.',
                 '有限先验最终离线位姿误差及相对参考动作风险，先按parent聚合。硬冻结与普通smoother作为显式对照；独立单位是12个parent，不是24份噪声。',
                 ['results/e2e3/inference/OFFLINE_SCORES.json','results/e2e3/inference/METHOD_RESULTS.json'])
    else:gaps.append('METHOD_RESULTS absent; E3 metrics/cost figures skipped')
    e4path=root/'results/e4/comparison/FULL_TABLE.csv'
    if e4path.exists():
        with e4path.open() as f:rows=list(csv.DictReader(f))
        table('E4_RAW_COMPARISON',rows)
        rec=[]
        for r in rows:rec.append(dict(r,parent=int(r['parent_id']),repeat=int(r['repeat_id']),method=str(r['method'])+' n='+str(r['particles'])))
        pg=grouped(rec,'method',['target_cdf_max','event_error','important_mode_mass_error','action_regret','process_cpu_seconds'])
        for row in pg:
            row['common_shared_bank_setup_rhs']=11800
            row['common_shared_bank_setup_factorizations']=2950
        table('E4_PARENT_METRICS',pg)
        baselines=sorted(set(r['method'] for r in pg if r['method'].startswith('MAP_POINT n=')))
        if baselines:table('E4_PAIRED_BOOTSTRAP',[r for metric in ['target_cdf_max','action_regret'] for r in paired_bootstrap(pg,'method',metric,baselines[0])])
        fig,axes=plt.subplots(1,2,figsize=(11,4));curves(pg,'method','target_cdf_max',axes[0]);curves(pg,'method','action_regret',axes[1]);axes[0].legend(fontsize=5)
        save(fig,'E4_saved_comparison','Saved E4 finite-prior comparison; repeats and MC seeds averaged within parent. Offline projections using exact posterior moments are not deployable algorithms; no continuous convergence claim.',
             '已保存E4有限先验比较，副本及MC种子先在parent内平均。使用精确后验矩的离线投影不是可部署算法，不声称连续收敛。',
             ['results/e4/comparison/FULL_TABLE.csv'])
    else:gaps.append('Optional E4 raw CSV absent')
    e5=root/'results/e5/finite_sensing';online=jsonread(e5/'ONLINE_RESULTS.json')
    gt=jsonread(e5/'offline/PARENT_GROUPED_GT_SCORES.json')
    if online:
        pg=grouped(online,'policy',['target_bayes_risk','initial_target_bayes_risk','selection_cpu_seconds_standalone','inference_cpu_seconds'])
        for row in pg:
            row['common_shared_bank_setup_rhs']=11800
            row['common_shared_bank_setup_factorizations']=2950
        table('E5_PARENT_ONLINE_METRICS',pg);counts=[]
        table('E5_RISK_PAIRED_BOOTSTRAP',paired_bootstrap(pg,'policy','target_bayes_risk','Uniform'))
        for policy in sorted(set(r['policy'] for r in online)):
            for parent in sorted(set(r['parent'] for r in online)):
                rr=[r for r in online if r['policy']==policy and r['parent']==parent]
                for action in range(4):counts.append(dict(parent=parent,policy=policy,action=action,selection_fraction=np.mean([r['action']==action for r in rr]),saved_repeats=len(rr)))
        table('E5_PARENT_SELECTION_COUNTS',counts)
        policies=sorted(set(r['policy'] for r in online));fig,axes=plt.subplots(1,2,figsize=(11,4))
        bottom=np.zeros(len(policies))
        for action in range(4):
            heights=np.array([np.mean([r['selection_fraction'] for r in counts if r['policy']==p and r['action']==action]) for p in policies])
            axes[0].bar(np.arange(len(policies)),heights,bottom=bottom,label=f'Action{action}');bottom+=heights
        axes[0].set_xticks(range(len(policies)),policies,rotation=40,ha='right');axes[0].set_ylabel('Parent-mean selection fraction');axes[0].legend(fontsize=6)
        curves(pg,'policy','target_bayes_risk',axes[1]);axes[1].legend(fontsize=6)
        save(fig,'E5_policy_choices_risk','Conditional finite-prior active sensing: parent-grouped choices and posterior target risk. Independent future noise is preserved; all policies share the costly bank. These are conditional pilot outcomes, not general sensing superiority.',
             '条件有限先验主动感知：按parent聚合的选点与目标后验风险。保留独立未来噪声，各策略共用已付费bank；仅为条件pilot，不宣称一般选点优势。',
             ['results/e5/finite_sensing/ONLINE_RESULTS.json'])
        if gt:
            rows=[dict(r,repeat=0,repeat_count=2,saved_rows=2) for r in gt]
            table('E5_PARENT_GT_ERRORS',rows)
            table('E5_PAIRED_BOOTSTRAP',paired_bootstrap(rows,'policy','mean_target_squared_error','Uniform'))
            fig,ax=plt.subplots(figsize=(8,4));curves(rows,'policy','mean_target_squared_error',ax);ax.legend(fontsize=6)
            save(fig,'E5_parent_GT_errors','Offline finite-prior target squared error, two repeats averaged within parent. Descriptive paired bootstrap compares twelve shared parents; no new acceptance threshold.',
                 '有限先验离线目标平方误差，两副本先在parent内平均。描述性配对bootstrap以12个共同parent比较，不引入新接受门槛。',
                 ['results/e5/finite_sensing/offline/PARENT_GROUPED_GT_SCORES.json'])
    else:gaps.append('Optional E5 online results absent')
    text=['# Pilot figures / 条件试验图表','',
        'Continuous A0 and explicit finite-prior E2–E5 are different statistical experiments. All figures read saved evidence only.','',
        '量化图先在parent内平均副本；2000次配对bootstrap、seed2026102510，仅12个完整parent时显示描述性区间；不作新科学门判断。','']
    for item in manifest:text += ['## '+item['name'],'',item['caption_zh'],'',item['caption_en'],'',f"![{item['name']}]({item['name']}.png)",'','Sources: '+', '.join(item['sources']),'']
    (out/'CAPTIONS_ZH_EN.md').write_text('\n'.join(text)+'\n')
    (out/'MANIFEST.json').write_text(json.dumps(dict(figures=manifest,missing_inputs=gaps,
        execution='SAVED_ARRAY_RENDER_ONLY_NO_PDE',paired_bootstrap=dict(seed=SEED,draws=2000,statistical_unit='12 independent parent means',descriptive_only=True),
        A0_noiseless_analysis_source='results/a0/analysis_repaired/RESULTS.json' if valid_repaired else 'UNAVAILABLE_NOT_COMPLETE_EVALUABLE72',
        A0_original_analysis_status='PRESERVED_LAYOUT_FAILURE_NOT_NEGATIVE_PHYSICAL_FINDING',
        finite_setup_CPU_policy='Only online/JOB_COST_FINAL.json with complete=true; never in-flight COST_RECEIPT or COMPLETE snapshot',
        no_scientific_gate_assigned=True),ensure_ascii=False,indent=2))
    print(json.dumps(dict(output=str(out),figure_count=len(manifest),gaps=gaps)))


if __name__=='__main__':main()

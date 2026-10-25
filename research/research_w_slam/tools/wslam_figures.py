#!/usr/bin/env python3
"""Render actual saved E0/E1 arrays only: no solver import or physical action."""
from pathlib import Path
import csv,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'figures/actual_e0_e1'


def read(path):
    return json.loads(path.read_text()) if path.exists() else None


def write_csv(name,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/(name+'.csv')).open('w',newline='') as file:
        writer=csv.DictWriter(file,keys);writer.writeheader();writer.writerows(rows)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.bbox':'tight','pdf.fonttype':42})
    figures=[];missing=[]
    def save(fig,name,en,zh,sources):
        fig.savefig(OUT/(name+'.png'),dpi=180);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
        figures.append(dict(name=name,caption_en=en,caption_zh=zh,sources=sources))
    original=read(ROOT/'results/e0/qualification/RESULTS.json')
    refined=read(ROOT/'results/e0/refinement01_executable/RESULTS.json')
    if original and refined:
        old={r['case']:r for r in original['rows'] if 'point_grid_white' in r}
        new=[r for r in refined['rows'] if 'point_grid_white' in r]
        rows=[dict(case=r['case'],initial_white=old[r['case']]['point_grid_white'],
                   refined_white=r['point_grid_white'],threshold=.1) for r in new]
        write_csv('01_precision_convergence',rows)
        fig,ax=plt.subplots(figsize=(10,3.8));x=np.arange(len(rows))
        ax.semilogy(x,[r['initial_white'] for r in rows],'o-',label='Initial mesh pair')
        ax.semilogy(x,[r['refined_white'] for r in rows],'s-',label='Refined mesh pair')
        ax.axhline(.1,color='r',ls='--',label='Frozen 0.1 threshold')
        ax.set_xticks(x,[r['case'].replace('support_','S') for r in rows],rotation=65,ha='right',fontsize=7)
        ax.set_ylabel('Full point-data difference / sigma');ax.legend();fig.tight_layout()
        save(fig,'01_precision_convergence','Fixed E0 data convergence; initial failures retained. These finite discrete cases are not a uniform continuum certificate.',
             '固定 E0 数据精度对照，保留最初失败。有限离散案例不等于统一连续 PDE 认证。',
             ['results/e0/qualification/RESULTS.json','results/e0/refinement01_executable/RESULTS.json'])
    else:missing.append('E0 original/refined precision results')
    blocks=read(ROOT/'results/e0/four_blocks/RESULTS.json')
    if blocks:
        rows=[{k:r[k] for k in ['target_contrast','relative_full_identity','target_induced_relative_identity','Born_target_induced_relative_error','setup_rhs']} for r in blocks]
        write_csv('02_four_blocks_born',rows)
        fig,ax=plt.subplots(figsize=(6.8,3.6))
        for key,label in [('relative_full_identity','Full four-block identity'),('target_induced_relative_identity','Target-induced identity'),('Born_target_induced_relative_error','Target-Born error')]:
            ax.semilogy([r['target_contrast'] for r in rows],[r[key] for r in rows],'o-',label=label)
        ax.set_xlabel('Intrinsic target contrast');ax.set_ylabel('Relative error');ax.legend();fig.tight_layout()
        save(fig,'02_four_blocks_born','Same-discretization Maxwell four-block identity versus target-Born error. Exact algebra uses paid port setup; it does not establish continuum accuracy or speedup.',
             '同一离散 Maxwell 算子的四块恒等式与目标 Born 误差。精确代数使用已计费端口准备，不证明连续精度或加速。',
             ['results/e0/four_blocks/RESULTS.json'])
    e1=read(ROOT/'results/e1/dynamic/RESULTS.json')
    erows=e1['rows'] if e1 else []
    confusion=[r for r in erows if r['case']=='confusion_model']
    if confusion:
        rows=[]
        for r in confusion:
            for i,value in enumerate(r['rho']):
                rows.append(dict(model=r['label'],direction=i,retention=value,
                    information_lambda_min=float(np.linalg.eigvalsh(r['Keff'])[0]),
                    bare_angle_rad=r['angles'][i],statistical_prior_assisted=bool(np.any(r['statistical_nuisance_precision']))))
        write_csv('03_confusion_models',rows)
        fig,axes=plt.subplots(1,2,figsize=(9,3.6));x=np.arange(len(confusion))
        for i in range(2):axes[0].bar(x+(i-.5)*.32,[r['rho'][i] for r in confusion],width=.32,label=f'rho {i+1}')
        axes[1].bar(x,[np.linalg.eigvalsh(r['Keff'])[0] for r in confusion],color='#487c89')
        labels=['Pose','Pose + intrinsic','+ gain + anchored nav'][:len(confusion)]
        for ax in axes:ax.set_xticks(x,labels,rotation=15,ha='right')
        axes[0].set_ylabel('Generalized information retention');axes[0].legend();axes[1].set_ylabel('Minimum absolute information eigenvalue')
        fig.tight_layout()
        save(fig,'03_confusion_models','Full16RX local linear information: pose, intrinsic target, then gain and anchored navigation. Last model includes declared statistical priors; bare principal angles alone do not equal prior-assisted retention.',
             '完整 16 接收器的局部线性信息：位姿、目标内禀参数、增益和已锚定导航。最后一项包含声明的统计先验，裸主角不等于先验辅助保留率。',
             ['results/e1/dynamic/RESULTS.json'])
    information=[r for r in erows if r['case']=='information_structure']
    if information:
        rows=[]
        for r in information:
            for mode,key in [('shared','rho_shared'),('independent','rho_independent')]:
                for i,value in enumerate(r[key]):rows.append(dict(panel=r['label'],mode=mode,direction=i,rho=value,scope=r['scope']))
        write_csv('04_aperture_retention',rows)
        fig,ax=plt.subplots(figsize=(10,3.8));x=np.arange(len(information))
        for mode,key,offset,color in [('Shared','rho_shared',-.15,'#237a8a'),('Independent','rho_independent',.15,'#b36047')]:
            minima=[min(r[key]) for r in information];maxima=[max(r[key]) for r in information]
            ax.vlines(x+offset,minima,maxima,color=color,lw=3);ax.scatter(x+offset,minima,color=color,label=mode);ax.scatter(x+offset,maxima,color=color)
        ax.set_xticks(x,[r['label'].replace('full_','16RX ').replace('minimal_','RX8 ') for r in information],rotation=30,ha='right',fontsize=8)
        ax.set_ylabel('Retention eigenvalue range');ax.legend(loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2);fig.tight_layout()
        save(fig,'04_aperture_retention','Local retention depends on aperture and nuisance sharing. RX8 is a single-complex-port mechanism control; it does not demonstrate full-aperture finite ambiguity. Actual-motion and hypothetical equality panels retain their original labels.',
             '局部信息保留依赖孔径与干扰参数共享方式。RX8 为单复端口机制对照，不证明完整孔径的有限歧义；实际运动与假设相等约束沿用原标签。',
             ['results/e1/dynamic/RESULTS.json'])
    profiles=[r for r in erows if r['case']=='finite_compensation_profile']
    if profiles:
        rows=[dict(mode=r['mode'],wall_step=r['wall_step'],endpoint_white=r['endpoint_white'],
            local_predicted_white=r['local_predicted_white'],actual_minus_local_white=r['actual_minus_local_white'],
            success=r['success'],nfev=r['nfev'],independent_grid_white=r.get('independent_grid_white','')) for r in profiles]
        write_csv('05_finite_corrector',rows)
        fig,axes=plt.subplots(1,2,figsize=(9,3.7))
        for mode in dict.fromkeys(r['mode'] for r in profiles):
            group=sorted([r for r in profiles if r['mode']==mode],key=lambda r:r['wall_step'])
            axes[0].semilogy([r['wall_step'] for r in group],[max(r['endpoint_white'],1e-16) for r in group],'o-',label=mode)
            axes[1].semilogy([r['wall_step'] for r in group],[max(r['actual_minus_local_white'],1e-16) for r in group],'o-')
            bad=[r for r in group if not r['success']]
            axes[0].scatter([r['wall_step'] for r in bad],[max(r['endpoint_white'],1e-16) for r in bad],marker='x',s=90,color='r')
        axes[0].set_ylabel('Finite endpoint residual / sigma');axes[1].set_ylabel('Actual minus local prediction / sigma')
        for ax in axes:ax.set_xlabel('Frozen wall-direction step')
        axes[0].legend(fontsize=7);fig.tight_layout()
        save(fig,'05_finite_corrector','Paid finite correctors at RX8, compared with the fixed local prediction. Red crosses retain optimizer nonconvergence. Local and finite residuals are separate evidence; these profiles are not A0 inference or full-aperture ambiguity.',
             'RX8 已计费有限补偿剖面与固定局部预测对照。红叉保留优化不收敛；局部与有限残差分开记录，不属于 A0 推断或完整孔径歧义。',
             ['results/e1/dynamic/RESULTS.json'])
    fieldpath=ROOT/'results/e0/qualification/derivative_base.npz'
    if fieldpath.exists():
        with np.load(fieldpath,allow_pickle=False) as z:
            points=z['points'];epsilon=z['epsilon'];field=z['field'][:,0]
        # Deterministic display subsample only; no interpolated or invented field.
        idx=np.arange(0,len(points),max(1,math.ceil(len(points)/12000)))
        pts=points[idx];eps=epsilon[idx];values=field[idx]
        write_csv('06_sampled_fields',[dict(x=p[0],y=p[1],eps_real=e.real,eps_imag=e.imag,E_real=v.real,E_imag=v.imag) for p,e,v in zip(pts,eps,values)])
        fig,axes=plt.subplots(1,3,figsize=(10,4))
        for ax,value,title,cmap in zip(axes,[eps.real,values.real,values.imag],['Re epsilon','Re total E (TX0)','Im total E (TX0)'],['viridis','RdBu_r','RdBu_r']):
            m=ax.scatter(pts[:,0],pts[:,1],c=value,s=2,cmap=cmap,rasterized=True)
            ax.set_aspect('equal');ax.set_title(title);ax.set_xlabel('x');ax.set_ylabel('y');fig.colorbar(m,ax=ax,shrink=.7)
        fig.tight_layout()
        save(fig,'06_sampled_fields','Actual saved E0 representative material and complex total field at quadrature points (initial derivative grid). Deterministic display subsampling only; this image is not an independent error check.',
             '实际保存的 E0 代表材料与复总场，取自最初导数网格的积分点。仅确定性抽样显示，该图不是独立误差审核。',
             ['results/e0/qualification/derivative_base.npz'])
    else:missing.append('representative quadrature field archive')
    # Future evidence is inventoried, never fabricated or interpreted via a guessed schema.
    optional=[]
    for folder in ('a0','probability'):
        directory=ROOT/'results'/folder
        if directory.exists():
            for path in sorted(directory.rglob('*.json')):
                value=read(path)
                optional.append(dict(path=str(path.relative_to(ROOT)),top_level_keys=list(value) if isinstance(value,dict) else 'list'))
    audit=dict(E1_snapshot_exists=e1 is not None,profile_count=len(profiles),
        information_panel_count=len(information),confusion_model_count=len(confusion),
        profile_failures=[r for r in erows if r['case']=='finite_compensation_profile_failure'],
        missing_explicit_raw_vectors=['Confusion NPZ archives do not embed nuisance-column names; ordering comes from the saved driver source.',
          'Profile archives contain final residual/Jacobian, not every optimizer iterate or trial field.',
          'No standalone E1 full-field archive; representative fields belong to E0.'],
        present_raw=['Frequency/pose data,J,names,TX,RX','NAV Jacobian and +/- observations',
         'A/B, C/R/DC, retention and rank arrays in per-panel NPZ','Profile endpoint poses/material/data/residual/local prediction/J',
         'Direction and fixed steps in RESULTS.json','Proper local Gaussian covariance and declared prior precision'],
        inference_policy='Missing arrays are not reconstructed as observed measurements or posterior evidence')
    (OUT/'SCHEMA_AND_MISSING.json').write_text(json.dumps(dict(audit=audit,optional_evidence=optional,missing_plot_inputs=missing),indent=2))
    (OUT/'FIGURES.json').write_text(json.dumps(figures,ensure_ascii=False,indent=2))
    text=['# 实际 E0/E1 图表 / Saved evidence figures','',
          '只读取已有数据；未执行 Maxwell、LU、RHS 或新的实验。 / Saved-array rendering only; no new physical actions.','']
    for f in figures:
        text += [f"## {f['name']}",'',f['caption_zh'],'',f['caption_en'],'',
                 f"![{f['name']}]({f['name']}.png)",'',
                 'Sources: '+', '.join('`'+s+'`' for s in f['sources']),'']
    text += ['## 原始量保存范围 / Raw scope','',*['- '+s for s in audit['missing_explicit_raw_vectors']],
             '', 'Optional future A0/probability JSON inventory is in SCHEMA_AND_MISSING.json; no missing data were invented.']
    (OUT/'CAPTIONS_ZH_EN.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(dict(figure_count=len(figures),output=str(OUT),missing=missing)))


if __name__=='__main__':main()

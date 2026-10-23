#!/usr/bin/env python3
"""Fixed observation000 OFFLINE_RENDER gallery; zero Maxwell/Book/optimization.

Import is inert. main must pass four completion gates BEFORE any GT access.
Only parameters.npy row0 is streamed from OFFLINE_TRUTH.npz; true F is never
read. No selection, fitting, sampling, inference or generation plan is run.
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
import types
import zipfile
from pathlib import Path

import numpy as np

def _stage_root():
    # Supports the repository layout and the adjacent stage/history Release layout.
    name='research_w_p0r_to_p1'
    for parent in Path(__file__).resolve().parents:
        for candidate in (parent,parent/'research'/name):
            if candidate.name==name and (candidate/'configs/AUTHORIZATION.json').is_file():
                return candidate
    raise RuntimeError('Research W evidence root was not found')

ROOT = _stage_root()
FIG = ROOT/'figures/controlled_gallery'
TAB = ROOT/'tables/controlled_gallery'
OBSERVATION = 0
FREQUENCY = 1.
RSS_LIMIT = 1_000_000_000
GATES = {
    'results/P1B/offline_scores/COMPLETE.json': 'OFFLINE_SCORES_READY_FOR_PARENT_REVIEW',
    'results/P1B/inference/COMPLETE.json': 'FULL_NODE_INFERENCE_COMPLETE',
    'results/P1B/point_validation/COMPLETE.json': 'FULL_POINT_VALIDATION_COMPLETE',
    'results/P1B/continuous_point_control/COMPLETE.json': 'FULL_CONTINUOUS_POINT_CONTROL_COMPLETE',
}
SCOPE = ('固定预注册 observation 000；受控公开两墙 prior；模型外质量 = null（未量化）。\n'
         '16 个独立参数 parents / 32 噪声观测；GT 仅离线；非一般 unknown-wall 在线成功；原 minimum_feature 未解决。')


def read_json(path):
    return json.loads(Path(path).read_text())


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT.resolve()))


def inside(path):
    p=(ROOT/Path(path)).resolve()
    if not p.is_relative_to(ROOT.resolve()):
        raise ValueError('Input path escapes this controlled study')
    return p


def write_json(path,obj):
    def encode(v):
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,np.ndarray):return v.tolist()
        raise TypeError(type(v).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=encode,allow_nan=False)+'\n')


def write_csv(path,rows):
    with Path(path).open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


def peak_rss():
    value=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform=='darwin' else value*1024)


def memory_guard():
    peak=peak_rss()
    if peak>RSS_LIMIT:
        raise MemoryError(f'Peak RSS {peak} exceeds explicit {RSS_LIMIT} byte limit')
    return peak


def require_gates(root=ROOT):
    """Fail closed, before reading either offline score GT or truth parameters."""
    answer={}
    for rel,status in GATES.items():
        path=Path(root)/rel
        if not path.exists():raise RuntimeError('OFFLINE_RENDER blocked: missing '+rel)
        value=read_json(path)
        if value.get('status')!=status:
            raise RuntimeError('OFFLINE_RENDER blocked: unfinished '+rel)
        if value.get('observations')!=32:
            raise RuntimeError('OFFLINE_RENDER blocked: completion population changed')
        if rel=='results/P1B/offline_scores/COMPLETE.json' and value.get('parent_count')!=16:
            raise RuntimeError('OFFLINE_RENDER blocked: independent parameter-parent count changed')
        if rel=='results/P1B/inference/COMPLETE.json' and value.get('independent_parents')!=16:
            raise RuntimeError('OFFLINE_RENDER blocked: inference parent count changed')
        answer[rel]=value
    return answer


def read_truth0(root=ROOT):
    """Read ONLY the first C-order five-float parameter row, never true F."""
    gates=require_gates(root)  # Helper cannot bypass the GT gate either.
    path=Path(root)/'results/P1B/data/OFFLINE_TRUTH.npz'
    with zipfile.ZipFile(path) as archive:
        with archive.open('parameters.npy') as stream:
            version=np.lib.format.read_magic(stream)
            if version==(1,0):shape,fortran,dtype=np.lib.format.read_array_header_1_0(stream)
            elif version==(2,0):shape,fortran,dtype=np.lib.format.read_array_header_2_0(stream)
            else:raise ValueError('Unsupported GT parameter NPY header; refusing broader GT read')
            if shape!=(32,5) or fortran or dtype.hasobject or dtype.kind!='f' or dtype.itemsize!=8:
                raise ValueError('GT parameters must be C-order float64 (32,5)')
            raw=stream.read(5*dtype.itemsize)
            if len(raw)!=5*dtype.itemsize:raise ValueError('Incomplete GT row0')
            parameters=np.frombuffer(raw,dtype=dtype).astype(float,copy=True)
    if not np.isfinite(parameters).all() or parameters[0] not in (0.,1.):
        raise ValueError('Invalid fixed GT row0')
    return parameters,dict(role='OFFLINE_RENDER_ONLY',observation=0,
        path='results/P1B/data/OFFLINE_TRUTH.npz',zip_member='parameters.npy',
        rows_read=[0],parameter_order=['branch','u','contrast','amplitude','phase'],
        data_payload_bytes_read=40,truth_F_read=False,generation_plan_read=False,
        gates_checked_before_read=gates,selection='observation000 fixed before outcomes',
        inference_or_method_parameters_changed=False)


def npz(path,keys=None):
    with np.load(path,allow_pickle=False) as z:
        result={k:z[k].copy() for k in (z.files if keys is None else keys)}
    if sum(a.nbytes for a in result.values())>128_000_000:
        raise MemoryError('Unexpectedly large input; field-history loading is prohibited')
    memory_guard()
    return result


def load_public():
    """Completed fixed case inputs only; no truth, plans, sample generation."""
    base=ROOT/'results/P1B'
    example_path=base/'handoff/observation_000/EXAMPLE.json'
    example=read_json(example_path)
    if example.get('observation')!=0 or example.get('family_id')!='W_FIXED_SUPPORT_JOINT_V1':
        raise ValueError('Only the predeclared controlled observation000 may be rendered')
    rows=read_json(base/'inference/RESULTS.json')
    if not rows or rows[0].get('observation')!=0:raise ValueError('Fixed observation000 identity missing')
    row=rows[0]
    factor_path=inside(example['downstream']['finite_joint_weights']['path'])
    if factor_path.parent.parent!=base/'inference/observation_000':
        raise ValueError('Finite reference does not belong to fixed observation000')
    material_path=inside(example['downstream']['public_material_response']['path'])
    if material_path!=base/'inference/PUBLIC_MATERIAL_GRID.npz':raise ValueError('Unexpected public response grid')
    material=npz(material_path,['F_full','branches','u','contrast','x','shape','response_kind'])
    finite=npz(factor_path,['posterior_weights','grid_shape','material_branches','material_u',
        'material_contrast','material_x','amplitude_axis','phase_axis','response_kind'])
    if (tuple(finite['grid_shape'][:3])!=(2,65,5) or tuple(material['shape'])!=(2,65,5)
        or str(material['response_kind'])!='FULL_CACHED_NODE' or str(finite['response_kind'])!='FULL_CACHED_NODE'):
        raise ValueError('Gallery requires the actual650 full-node finite target')
    for left,right in [('material_branches','branches'),('material_u','u'),
                       ('material_contrast','contrast'),('material_x','x')]:
        if not np.array_equal(finite[left],material[right]):raise ValueError('Finite factors/grid coordinates mismatch')
    shape=tuple(int(i) for i in finite['grid_shape'])
    w=np.asarray(finite['posterior_weights'],float).reshape(650,-1)
    if w.shape[1]!=shape[3]*shape[4] or np.any(w<0) or not np.isfinite(w).all() or abs(float(w.sum())-1)>1e-10:
        raise ValueError('Invalid joint posterior weights; no silent normalization')
    predictions_path=base/'handoff/observation_000/PREDICTIONS_AND_MATERIALS.npz'
    predictions=npz(predictions_path)
    observed=npz(factor_path.parent/'OBSERVATION.npz',['y'])['y']
    if observed.shape!=(5,4,12) or not np.iscomplexobj(observed) or not np.isfinite(observed).all():
        raise ValueError('Invalid actual noisy observation')
    point_path=base/'point_validation/observation_000/RESULT.json'
    control_path=base/'continuous_point_control/observation_000/RESULT.json'
    point=read_json(point_path);control=read_json(control_path)
    scores_path=base/'offline_scores/SCORES.json'
    scores=[r for r in read_json(scores_path) if r['observation']==0]
    if not scores or point.get('observation')!=0 or control.get('observation')!=0:
        raise ValueError('Missing fixed-case score/control identity')
    sources=[relative(p) for p in [example_path,factor_path,material_path,predictions_path,
        factor_path.parent/'OBSERVATION.npz',point_path,control_path,scores_path,base/'inference/RESULTS.json']]
    return dict(example=example,row=row,material=material,finite=finite,joint_weights=w,
        predictions=predictions,observed=observed,point=point,control=control,scores=scores,sources=sources)


def finite_quantiles(x,weights,probabilities=(.025,.5,.975)):
    x=np.asarray(x,float);weights=np.asarray(weights,float)
    if x.shape!=weights.shape or np.any(weights<0) or abs(float(weights.sum())-1)>1e-10:
        raise ValueError('Expected normalized finite atom weights')
    order=np.argsort(x,kind='stable');cdf=np.cumsum(weights[order]);cdf[-1]=1.
    return x[order][np.searchsorted(cdf,probabilities,side='left')]


def branch_stream_statistics(family,material,masses,grid):
    """Two wall decodes + streaming650 exact frozen target overrides.

    No650 x ny x nx field history. The decoder's target is an overwrite,
    so moments use differences of first/second powers on the exact disk.
    """
    from p0.geometry import wall_values
    from wr.endpoint_geometry import linear_material_values
    decoder=types.SimpleNamespace(wall_values=wall_values)
    seed=family.config['wall_support']['seed']
    flat=grid.reshape(-1,2)
    walls=[]
    for branch in (0,1):
        er,sigma=linear_material_values(seed,family.walls[branch],flat,decoder,probe=None)
        walls.append(np.stack((er,sigma/FREQUENCY)))
    walls=np.asarray(walls)
    branch_mass=np.array([float(masses[material['branches']==b].sum()) for b in (0,1)])
    mean=np.einsum('b,bdp->dp',branch_mass,walls)
    second=np.einsum('b,bdp->dp',branch_mass,walls*walls)
    radius=float(family.config['target']['radius'])
    for i,mass in enumerate(masses):
        b=int(material['branches'][i]);u=float(material['u'][i]);c=float(material['contrast'][i])
        node=family.coordinates(b,u,c)
        active=np.linalg.norm(flat-np.asarray(node['center']),axis=1)<radius
        values=np.array([1+c,float(family.config['target']['conductivity'])/FREQUENCY])[:,None]
        baseline=walls[b][:,active]
        mean[:,active]+=mass*(values-baseline)
        second[:,active]+=mass*(values*values-baseline*baseline)
        if i in (0,162,324,325,487,649):
            reconstructed=walls[b].copy();reconstructed[:,active]=values
            actual=family.material(b,u,c,grid,f=FREQUENCY)
            if not np.allclose(reconstructed[0]+1j*reconstructed[1],actual.ravel(),rtol=0,atol=1e-12):
                raise ValueError('Streamed target overwrite differs from frozen decoder')
        if i%50==0:memory_guard()
    variance=second-mean*mean
    minimum=float(variance.min())
    if minimum < -1e-10:raise FloatingPointError('Material second moment lost numerical consistency')
    variance=np.maximum(variance,0.)
    return dict(mean=mean.reshape(2,*grid.shape[:-1]),variance=variance.reshape(2,*grid.shape[:-1]),
        branch_mass=branch_mass,minimum_variance_before_roundoff_clip=minimum,
        method='Two public frozen wall decodes plus streamed650 exact disk overwrites; no field history',
        actual_decoder_spotcheck_node_indices=[0,162,324,325,487,649])


def correlation_preserving_prediction(public):
    f=public['finite'];w=public['joint_weights'];F=public['material']['F_full']
    a,p=np.meshgrid(f['amplitude_axis'],f['phase_axis'],indexing='ij')
    gamma=(a*np.exp(1j*p)).ravel()
    factors=(w@gamma.real)+1j*(w@gamma.imag)
    expected=np.tensordot(factors,F,axes=(0,0))
    saved=public['predictions']['posterior_predictive_approx']
    if not np.allclose(expected,saved,rtol=1e-12,atol=1e-13):raise ValueError('Saved E[gammaF] disagrees with actual joint factors')
    masses=w.sum(axis=1)
    if not np.allclose(masses,public['predictions']['material_posterior_masses'],rtol=0,atol=1e-12):
        raise ValueError('Saved material masses differ from joint factors')
    wrong=(w.sum(axis=0)@gamma)*np.tensordot(masses,F,axes=(0,0))
    return expected,dict(factorization_not_used=True,
        factorized_difference_raw_norm=float(np.linalg.norm(expected-wrong)),
        factorized_difference_white_norm=float(np.linalg.norm(expected-wrong)/.01),
        formula='sum_m F_full[m]*sum_gamma(gamma*w[m,gamma]); not E[gamma]*E[F]')


def setup_plotting():
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    names={f.name for f in font_manager.fontManager.ttflist}
    font=next((f for f in ('PingFang SC','Heiti SC','Arial Unicode MS') if f in names),None)
    if font is None:raise RuntimeError('No CJK font available for the Chinese gallery captions')
    matplotlib.rcParams.update({'font.family':[font,'DejaVu Sans'],'axes.unicode_minus':False,
        'font.size':10,'figure.facecolor':'#faf9f5','axes.facecolor':'#faf9f5',
        'savefig.facecolor':'#faf9f5','svg.fonttype':'path','axes.spines.top':False,'axes.spines.right':False})
    import matplotlib.pyplot as plt
    return plt,font


def save_figure(plt,fig,name,caption,sources):
    fig.text(.045,.025,caption+'\n'+SCOPE,fontsize=8.5,va='bottom',ha='left',linespacing=1.5)
    fig.savefig(FIG/f'{name}.png',dpi=160)
    fig.savefig(FIG/f'{name}.svg')
    plt.close(fig);memory_guard()
    (FIG/f'{name}.caption.md').write_text(caption+'\n\n'+SCOPE+'\n')
    write_json(FIG/f'{name}.sources.json',dict(observation=0,sources=sources,new_RHS=0,new_LU=0,
        truth_role='OFFLINE_RENDER_ONLY',outside_family_mass=None,scope=SCOPE))


def material_figure(plt,public,truth,representatives,statistics,x,y):
    fields=[np.stack((truth.real,truth.imag))]+[np.stack((z.real,z.imag)) for z in representatives]+[statistics['mean']]
    events=public['example']['event_partition']['events']
    titles=['GT（仅离线）']
    for b in (0,1):
        mass=events[b]['posterior_mass']
        titles.append(f'候选世界 {"A" if b==0 else "B"}\n事件质量 '+('UNQUALIFIED' if mass is None else f'{mass:.3%}'))
    titles.append('有限 posterior 材料均值\n不是合法的单一世界')
    fig,axes=plt.subplots(2,4,figsize=(16,9.2),sharex=True,sharey=True)
    fig.subplots_adjust(left=.055,right=.92,top=.82,bottom=.22,wspace=.12,hspace=.18)
    fig.suptitle('同一观测的两种候选世界与材料均值：图像来自真实保存参数',x=.045,ha='left',fontsize=19,y=.97)
    fig.text(.045,.90,'two fixed wall fields（公开 q0/q1）  |  16 independent parameter parents  |  posterior mean 不是合法世界',fontsize=12)
    extent=[x[0],x[-1],y[0],y[-1]]
    for row,label in enumerate(('Re ε（相对介电参数）','Im ε = σ/f，f=1')):
        vmin=min(float(z[row].min()) for z in fields);vmax=max(float(z[row].max()) for z in fields)
        for col,z in enumerate(fields):
            image=axes[row,col].imshow(z[row],extent=extent,origin='lower',cmap='viridis',
                vmin=vmin,vmax=vmax,interpolation='nearest',aspect='equal')
            if row==0:axes[row,col].set_title(titles[col],fontsize=12,pad=12)
            if col==0:axes[row,col].set_ylabel(label+'\ny（原生坐标）')
            if row==1:axes[row,col].set_xlabel('x（原生坐标）')
        cax=fig.add_axes([.94,.56 if row==0 else .25,.012,.23])
        fig.colorbar(image,cax=cax)
    save_figure(plt,fig,'01_material_worlds',
        '每行四图共享色标；ε 图为固定400×300节点采样，未做平滑。候选列质量来自独立连续参考事件积分，非单点概率。\n'
        '右列按有限联合后验权重平均材料，墙与目标均平均；可出现两个弱目标，不能解释为第三个合法世界。Im ε 为本研究 σ/f，并非耗散功率。',public['sources']+['results/P1B/data/OFFLINE_TRUTH.npz:parameters[0]'])


def task_figure(plt,public,true_x,masses):
    material=public['material'];row=public['row']
    bounds=[e['definition']['target_x_interval'] for e in public['example']['event_partition']['events']]
    finite_q=finite_quantiles(material['x'],masses)
    point=public['point']['methods']['direct joint point']
    control=public['control']['best']
    actions=[dict(method='GT（仅离线）',x=true_x,interval=None,kind='GT'),
        dict(method='Direct-grid 点',x=float(point['x']),interval=None,kind='point'),
        dict(method='强化连续 control 点',x=float(control['x']),interval=None,kind='point')]
    for schedule,label in [('direct','Direct SMC 均值'),('progressive','Progressive SMC 均值')]:
        e=row['ensemble'][schedule];interval=e['target_interval95']
        actions.append(dict(method=label,x=float(e['average_x_mean']),
            interval=[float(interval['lower']),float(interval['upper'])],kind='finite_Bayes_action'))
    actions.extend([dict(method='有限求积参考均值',x=float(masses@material['x']),
        interval=[float(finite_q[0]),float(finite_q[2])],kind='finite_Bayes_action'),
        dict(method='连续参考求积均值',x=float(row['continuous_ref']['x_mean']),interval=None,kind='continuous_quadrature_Bayes_action')])
    fig,(ax,cax)=plt.subplots(1,2,figsize=(15.5,8.8),gridspec_kw={'width_ratios':[1.2,1]})
    fig.subplots_adjust(left=.16,right=.96,top=.81,bottom=.26,wspace=.24)
    fig.suptitle('目标 x：合法世界位置与平方损失 Bayes action 的区别',x=.045,ha='left',y=.96,fontsize=19)
    for a in (ax,cax):
        for lo,hi in bounds:a.axvspan(lo,hi,color='#6e9e83',alpha=.17)
        a.axvspan(bounds[0][1],bounds[1][0],color='#bfa98d',alpha=.13,hatch='//')
        a.set_xlim(bounds[0][0]-.007,bounds[1][1]+.007)
        a.set_xlabel('目标 x（原生位置坐标）')
    colors=['#20252b','#af5f2c','#9653a4','#2c6eac','#357f5b','#313d4c','#a93440']
    for i,(record,color) in enumerate(zip(actions,colors)):
        in_gap=bounds[0][1]<record['x']<bounds[1][0]
        record.update(in_unsupported_prior_gap=bool(in_gap),
            lower95=None if record['interval'] is None else record['interval'][0],
            upper95=None if record['interval'] is None else record['interval'][1],
            interval_semantics='finite weighted atoms, no bin spreading' if record['interval'] else 'NOT_ASSIGNED')
        if record['interval'] is not None:ax.plot(record['interval'],[i,i],color=color,lw=2,alpha=.8)
        ax.scatter(record['x'],i,color=color,marker='D' if i==0 else 'o',s=45,zorder=5)
        if in_gap:ax.annotate('仅 Bayes action',xy=(record['x'],i),xytext=(7,7),textcoords='offset points',fontsize=8,color=color)
    ax.set_yticks(range(len(actions)),[a['method'] for a in actions]);ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    ax.set_title('阴影为两段 prior 位置支撑；斜线区无合法 target center',fontsize=10,pad=15)
    order=np.argsort(material['x'],kind='stable');sx=material['x'][order];sw=masses[order]
    unique,starts=np.unique(sx,return_index=True);atom_mass=np.add.reduceat(sw,starts)
    cax.step(np.r_[bounds[0][0]-.007,unique,bounds[1][1]+.007],
        np.r_[0,np.cumsum(atom_mass),1],where='post',color='#324c72',lw=1.4,label='有限原子 CDF（step）')
    continuous=row['continuous_ref'];cdfrows=[]
    for cut,values in continuous['CDF'].items():
        t=float(cut);atomic=float(masses[material['x']<=t].sum())
        cdfrows.append(dict(cut=t,finite_atomic_CDF=atomic,continuous_order5=values['order5'],
            continuous_order9=values['order9'],order5_vs9_change=values['error']))
    cax.scatter([r['cut'] for r in cdfrows],[r['continuous_order5'] for r in cdfrows],marker='x',
        color='#8c8d8f',s=45,label='独立连续积分 GL5 检查点',zorder=6)
    cax.scatter([r['cut'] for r in cdfrows],[r['continuous_order9'] for r in cdfrows],
        facecolors='none',edgecolors='#a94832',s=65,label='独立连续积分 GL9 检查点',zorder=7)
    cax.set(ylim=(-.025,1.025),ylabel='累积概率（无量纲）');cax.grid(alpha=.15)
    cax.legend(fontsize=8,frameon=False,loc='upper left');cax.set_title('有限 CDF 不展宽；连续参考只画已积分检查点',fontsize=10,pad=15)
    write_csv(TAB/'TASK_ACTIONS.csv',[{k:v for k,v in r.items() if k!='interval'} for r in actions])
    write_csv(TAB/'CDF_CHECKPOINTS.csv',cdfrows)
    save_figure(plt,fig,'02_target_actions_and_cdf',
        '区间为各有限加权分布的等尾95%原子分位区间；不保证连续覆盖率。SMC 按4个既有seed的等权混合汇总，无新增采样。\n'
        '均值若落在 unsupported gap，只是平方损失下的 action，不是第三世界。强化control仅保存的Full候选最优，未获 actual-Full KKT/global最优认证。',public['sources']+['results/P1B/data/OFFLINE_TRUTH.npz:parameters[0]'])
    return actions,cdfrows


def prediction_figure(plt,public,predictive):
    observed=public['observed'];pred=public['predictions']
    curves=[('噪声观测 y',observed,'#22272d'),('候选 A：γA × paid Full',pred['gammaF'][0],'#bb6830'),
        ('候选 B：γB × paid Full',pred['gammaF'][1],'#376dad'),('有限 E[γF]',predictive,'#34805e')]
    fig,axes=plt.subplots(2,2,figsize=(15.5,8.8),sharex=True)
    fig.subplots_adjust(left=.07,right=.96,top=.78,bottom=.25,wspace=.19,hspace=.23)
    fig.suptitle('完整复响应：噪声观测、实际 paid Full 候选和相关性保留的 E[γF]',x=.045,ha='left',y=.96,fontsize=18)
    for col,component in enumerate(('real','imag')):
        for label,F,color in curves:
            axes[0,col].plot(np.arange(240),getattr(F,component).ravel(),color=color,
                lw=1.1 if F is observed else 1.4,alpha=.75,label=label)
            if F is not observed:
                axes[1,col].plot(np.arange(240),getattr((F-observed)/.01,component).ravel(),
                    color=color,lw=1.1,alpha=.85)
        axes[0,col].set_title('实部' if col==0 else '虚部')
        axes[0,col].set_ylabel('复响应原始数值（本研究归一化）')
        axes[1,col].set_ylabel('残差 / σ，σ=0.01（无量纲）')
        axes[1,col].axhline(0,color='#808080',lw=.6)
        for a in axes[:,col]:
            for boundary in (47.5,95.5,143.5,191.5):a.axvline(boundary,color='#cccccc',lw=.6)
            a.grid(axis='y',alpha=.14)
            a.set_xticks([23.5+48*i for i in range(5)],
                [f'f={f:g}' for f in public['example']['acquisition']['frequencies']])
        axes[1,col].set_xlabel('频率块内依次排列4源×12接收器；共240复通道')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.51,.88),ncol=4,frameon=False,fontsize=9)
    rows=[]
    for label,F,_ in curves[1:]:
        distance=float(np.linalg.norm(F-observed)/.01)
        rows.append(dict(prediction=label,white_noisy_observation_residual_norm=distance,
            Phi=.5*distance**2,new_RHS=0,meaning='Noisy fit diagnostic, not truth error or a two-mean witness test'))
    write_csv(TAB/'RESPONSE_RESIDUALS.csv',rows)
    save_figure(plt,fig,'03_observation_and_paid_predictions',
        '候选预测直接取已付费 full9 原始复 F，再乘其保存代表 γ；不是重新求解。E[γF] 使用650个实际Full节点的联合权重，保留 γ—材料相关性。\n'
        '原始值未标为SI场强；频率为本研究归一化频率。残差是对噪声观测的fit，不是真值误差；后验预测均值也不要求对应单一物理世界。',public['sources']+['results/P1A/rules/full_9.npz'])
    return rows


def variance_figure(plt,statistics,x,y,sources):
    fig,axes=plt.subplots(1,2,figsize=(11.5,6.7))
    fig.subplots_adjust(left=.07,right=.94,top=.8,bottom=.26,wspace=.25)
    fig.suptitle('有限 posterior 材料方差：墙分支与目标连续参数共同贡献',x=.045,ha='left',fontsize=16,y=.96)
    for i,title in enumerate(('Var(Re ε)','Var(Im ε)，f=1')):
        image=axes[i].imshow(statistics['variance'][i],origin='lower',extent=[x[0],x[-1],y[0],y[-1]],
            vmin=0,cmap='magma',interpolation='nearest',aspect='equal')
        axes[i].set(title=title,xlabel='x（原生坐标）',ylabel='y（原生坐标）')
        fig.colorbar(image,ax=axes[i],fraction=.045,pad=.03)
    save_figure(plt,fig,'04_finite_material_variance',
        '逐像素有限后验 Var=E[field²]−E[field]²；单位分别为相对介电实部²与(σ/f)²。非噪声、网格或模型误差图。\n'
        '仅保存两墙场和流式650个目标覆盖的矩统计，不保存650份场历史；像素为展示采样，不宣称连续边界证书。',sources)


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    gates=require_gates()  # MUST be the first read; no GT/offline scores before this.
    memory_guard()
    public=load_public()
    sys.path.insert(0,str(ROOT/'src'))
    sys.path.insert(0,str(ROOT.parent/'research_w_continuation/vendor'))
    from wr.probability_family import Family
    family=Family.load()
    p=public['predictions'];x=p['grid_x'];y=p['grid_y']
    if x.shape!=(400,) or y.shape!=(300,) or float(p['frequency'])!=FREQUENCY:
        raise ValueError('Unexpected fixed400x300, f1 handoff material grid')
    xx,yy=np.meshgrid(x,y);grid=np.stack((xx,yy),axis=-1)
    events=public['example']['event_partition']['events']
    if len(events)!=2 or [e['representative']['branch'] for e in events]!=[0,1]:
        raise ValueError('Expected both predeclared branch events in order')
    representatives=[]
    paid=npz(ROOT/'results/P1A/rules/full_9.npz',['F','branches','u','contrast'])
    for i,event in enumerate(events):
        rep=event['representative'];source=rep['full9_prediction_source'];idx=int(source['node_index'])
        if source['rule_file']!='results/P1A/rules/full_9.npz':raise ValueError('Representative is not the declared paid full9 node')
        if (int(paid['branches'][idx])!=rep['branch'] or float(paid['u'][idx])!=rep['u']
            or float(paid['contrast'][idx])!=rep['contrast'] or not np.array_equal(paid['F'][idx],p['Fraw'][i])):
            raise ValueError('Representative raw prediction lost paid-node provenance')
        gained=rep['amplitude']*np.exp(1j*rep['phase'])*p['Fraw'][i]
        if not np.allclose(gained,p['gammaF'][i],rtol=0,atol=1e-13):raise ValueError('Representative gain prediction mismatch')
        z=family.material(rep['branch'],rep['u'],rep['contrast'],grid,f=FREQUENCY)
        if not np.allclose(z.real,p['epsilon_real'][i],rtol=0,atol=1e-12) or not np.allclose(z.imag,p['epsilon_imag_loss'][i],rtol=0,atol=1e-12):
            raise ValueError('Saved representative material differs from frozen decoder')
        representatives.append(z)
    predictive,correlation=correlation_preserving_prediction(public)
    masses=public['joint_weights'].sum(axis=1)
    statistics=branch_stream_statistics(family,public['material'],masses,grid)
    # Only now read GT, after checking all public/frozen material inputs.
    truth_parameters,access=read_truth0()
    true_branch=int(truth_parameters[0]);true_u=float(truth_parameters[1]);true_c=float(truth_parameters[2])
    truth=family.material(true_branch,true_u,true_c,grid,f=FREQUENCY)
    true_x=float(family.coordinates(true_branch,true_u,true_c)['center'][0])
    for score in public['scores']:
        if score['true_branch']!=true_branch or abs(float(score['true_target_x'])-true_x)>1e-12:
            raise ValueError('Offline fixed-case scores disagree with GT row0 identity')
    FIG.mkdir(parents=True,exist_ok=True);TAB.mkdir(parents=True,exist_ok=True)
    write_json(TAB/'GT_ACCESS.json',access)
    plt,font=setup_plotting()
    material_figure(plt,public,truth,representatives,statistics,x,y)
    actions,cdf=task_figure(plt,public,true_x,masses)
    residuals=prediction_figure(plt,public,predictive)
    variance_figure(plt,statistics,x,y,public['sources'])
    np.savez_compressed(TAB/'FIELD_STATISTICS.npz',grid_x=x,grid_y=y,frequency=FREQUENCY,
        finite_material_mean=statistics['mean'],finite_material_variance=statistics['variance'],
        finite_branch_mass=statistics['branch_mass'],material_node_masses=masses,
        posterior_predictive= predictive)
    numbers=[dict(quantity='GT_target_x_offline',value=true_x,unit='native_position'),
        dict(quantity='factorized_prediction_white_difference_not_used',value=correlation['factorized_difference_white_norm'],unit='dimensionless_white_response')]
    for b,event in enumerate(events):
        numbers.extend([dict(quantity=f'branch{b}_continuous_event_mass',value=event['posterior_mass'],unit='probability_or_null'),
            dict(quantity=f'branch{b}_finite_event_mass',value=float(statistics['branch_mass'][b]),unit='probability')])
    write_csv(TAB/'MAJOR_NUMBERS.csv',numbers)
    write_json(TAB/'SOURCES.json',dict(observation=0,selection='predeclared fixed000, no outcome selection',
        two_fixed_wall_fields=True,independent_parameter_parents=16,noise_observations=32,
        posterior_mean_is_legal_world=False,
        gates=gates,inputs=public['sources']+['results/P1A/rules/full_9.npz','configs/W_CONTROLLED_PRIOR.json',
        'inputs/W_FAMILY_WALLS.npz','results/P1B/data/OFFLINE_TRUTH.npz:parameters[0]'],
        decoder_sources=['src/wr/probability_family.py:Family.material','src/wr/endpoint_geometry.py:linear_material_values',
        '../research_w_continuation/vendor/p0/geometry.py:wall_values'],
        prediction_correlation=correlation,material_statistics={k:v for k,v in statistics.items() if k not in ('mean','variance')},
        source_inference_status=public['row']['numerical_probability_status'],
        outside_family_mass=None,original_minimum_feature='UNRESOLVED',
        new_Maxwell=0,new_RHS=0,new_LU=0,Book_constructed=False,sampling=False,optimization=False))
    (TAB/'README.md').write_text('# 固定 observation 000 离线成像图库\n\n'
        '四组图为实际保存参数/完整响应与冻结 decoder 的离线展示。完成门槛先于 GT 读取；只读 GT parameters 第0行，不读真值F或生成计划。\n\n'
        '候选世界列质量属于两个公开墙分支与目标位置区间事件；不是单个材料点概率。材料均值/方差使用有限联合posterior，不能视为合法第三世界。\n\n'
        '目标图的95%区间来自有限原子权重，不宣称连续覆盖率；CDF保留原子step，连续参考只显示独立积分检查点。预测均值使用相关性保留的E[gamma F]。\n\n'
        '650材料节点仅流式累计目标覆盖的场矩，公开墙各解码一次，没有场历史堆叠。模型外质量为null，一般unknown-wall在线成功不在图的主张范围，原minimum_feature语义保持未解决。\n\n'
        '具体来源与主要数字见SOURCES.json、MAJOR_NUMBERS.csv及其余表；GT访问范围见GT_ACCESS.json，CPU/墙钟/峰值RSS见RUNTIME.json。\n')
    runtime=dict(status='OFFLINE_RENDER_COMPLETE_FOR_PARENT_REVIEW',observation=0,font=font,
        cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        peak_RSS_bytes=memory_guard(),RSS_limit_bytes=RSS_LIMIT,new_Maxwell=0,new_RHS=0,new_LU=0,
        GT_role='OFFLINE_RENDER_ONLY',truth_F_read=False,generation_plan_read=False,
        sampling=False,optimization=False,Book_constructed=False,scientific_GO='NOT_EVALUATED')
    write_json(TAB/'RUNTIME.json',runtime)
    print(json.dumps(runtime,ensure_ascii=False),flush=True)


if __name__=='__main__':
    main()

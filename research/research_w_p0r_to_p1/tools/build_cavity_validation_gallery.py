#!/usr/bin/env python3
"""Decoder/paid-array cavity gallery. Import inert; zero PDE/Book/backend.

Only the fixed cal-only nomination and its completed B2 independent validation
are read. No P1 truth, observation or generation plan is accessed. A partial
24-task sweep or unfinished argmax check produces NOT_RUN, never a final plot.
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
import types
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
FIG=ROOT/'figures/cavity_validation'
TAB=ROOT/'tables/cavity_validation'
NOMINATION=ROOT/'configs/SECOND_CAVITY_NOMINATION.json'
IDENTIFIER='cavity_seed2_narrow_bestcal'
VALIDATION=ROOT/'results/B2/validation'/IDENTIFIER
CENTER=(-.85,.85)
RADIUS=.06
SIGMA_REAL=.01
FREQUENCY=1.


def load_json(path):
    return json.loads(Path(path).read_text())


def write_json(path,value):
    def encode(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False,default=encode)+'\n')


def inside(rel):
    p=(ROOT/rel).resolve()
    if not p.is_relative_to(ROOT.resolve()):raise ValueError('Candidate path escapes this stage')
    return p


def rss_bytes():
    n=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(n if sys.platform=='darwin' else n*1024)


def completion_state():
    required=[VALIDATION/'RESULTS.json',VALIDATION/'argmax_validation/TASK_ARGMAX_GRID_CHECK.json',
              VALIDATION/'argmax_validation/COMPLETE.json',VALIDATION/'ALL_TASKS.json',
              VALIDATION/'argmax_validation/TASK_ARGMAX_GRID_CHECK.npz']
    missing=[str(p.relative_to(ROOT)) for p in required if not p.exists()]
    if missing:return dict(status='NOT_RUN',reason='Independent validation still incomplete',missing=missing)
    result=load_json(required[0]);argmax=load_json(required[1]);tasks=load_json(required[3])
    if argmax.get('status') not in ('PASS_EMPIRICAL_ARGMAX_GRID_CHECK','FAIL_EMPIRICAL_ARGMAX_GRID_CHECK'):
        return dict(status='NOT_RUN',reason='Argmax numerical check unfinished/incomplete',source_status=argmax.get('status'))
    if (len(tasks)!=24 or {r.get('probe') for r in tasks}!=set(range(24))
        or not all((VALIDATION/f'probe{i}_PAIR.npz').exists() for i in range(24))):
        return dict(status='NOT_RUN',reason='All24 actual paid task pairs required')
    if not isinstance(result.get('numerical_pass'),bool) or 'd_T_sup_24' not in result:
        return dict(status='NOT_RUN',reason='Base validation result incomplete')
    return dict(status='READY',result=result,argmax=argmax,tasks=tasks)


def load_pair(path,keys):
    with np.load(path,allow_pickle=False) as z:
        arrays=[z[k].copy() for k in keys]
    for a in arrays:
        if a.shape!=(5,4,12) or not np.iscomplexobj(a) or not np.isfinite(a).all():
            raise ValueError('Expected complete actual5x4x12 complex response: '+str(path))
    return arrays


def norm_difference(pair):
    return float(np.linalg.norm(pair[1]-pair[0])/SIGMA_REAL)


def checked_sources(state):
    cfg=load_json(NOMINATION)
    if (cfg.get('id')!=IDENTIFIER or cfg.get('seed')!=2026102202 or cfg.get('kind')!='cavity'
            or cfg.get('band')!='narrow' or cfg.get('chart')!='linear'
            or cfg.get('selection_uses_only_calibration') is not True
            or cfg.get('cavities')!=[[list(CENTER),RADIUS]]):
        raise ValueError('Frozen nominee/strict-air cavity contract changed')
    if cfg!=load_json(VALIDATION/'MANIFEST.json'):
        raise ValueError('Independent validation is not the exact frozen nominee')
    with np.load(inside(cfg['x_file']),allow_pickle=False) as z:q=z['x'].copy()
    if q.shape!=(108,) or np.iscomplexobj(q) or not np.isfinite(q).all() or np.any(abs(q)>1):
        raise ValueError('Invalid completed nominee affine material coordinates')
    original=load_json(ROOT/'configs/ORIGINAL_P0.json')
    if original['noise']['each_real_quadrature_std']!=SIGMA_REAL or len(original['frequencies']['narrow'])!=5:
        raise ValueError('Original noise/acquisition differs; no silent rescaling')
    for world,expected in ((0,np.zeros_like(q)),(1,q)):
        p=VALIDATION/f'fields/cal_world{world}_fine/endpoint_0000/PARAMETERS.json'
        metadata=load_json(p)
        if not np.array_equal(np.asarray(metadata['x']),expected):
            raise ValueError('Decoded wall parameters differ from paid cal endpoint')
        recorded=metadata.get('metadata',{})
        expected_cavities=cfg['cavities'] if world else []
        if recorded.get('cavities')!=expected_cavities or recorded.get('seed')!=cfg['seed']:
            raise ValueError('Paid cal endpoint cavity/support differs from nominee')
    pair=load_pair(VALIDATION/'cal_PAIR.npz',('fine0','fine1'))
    dcal=norm_difference(pair)
    if not np.isclose(dcal,state['result']['d_cal'],rtol=1e-12,atol=1e-12):
        raise ValueError('Actual saved calibration pair disagrees with result')
    return cfg,q,dcal,pair,original


def task_arrays(state,public_probes):
    rows=[];arrays=[]
    source_by_id={r['probe']:r for r in state['tasks']}
    for index in range(24):
        source=source_by_id[index];probe=public_probes[index]
        if (not np.array_equal(source['center'],probe[0]) or source['radius']!=probe[1]
                or source['contrast']!=probe[2]):
            raise ValueError('The fixed original24 task class changed')
        pair=load_pair(VALIDATION/f'probe{index}_PAIR.npz',('data0','data1'))
        distance=norm_difference(pair)
        if not np.isclose(distance,source['distance'],rtol=1e-12,atol=1e-12):
            raise ValueError('Task distance is not reproduced by actual paid response')
        rows.append(dict(probe=index,center_x=probe[0][0],center_y=probe[0][1],radius=probe[1],
            contrast=probe[2],distance_white=distance,threshold=3.,meets_task_threshold=bool(distance>=3),
            source=f'results/B2/validation/{IDENTIFIER}/probe{index}_PAIR.npz'))
        arrays.append(pair)
    best=max(range(24),key=lambda i:rows[i]['distance_white'])
    if best!=state['result']['task_argmax'] or best!=state['argmax']['task_argmax']:
        raise ValueError('Frozen24-task argmax disagrees across source checks')
    for recorded in (state['result']['d_T_sup_24'],state['argmax']['fine_distance']):
        if not np.isclose(rows[best]['distance_white'],recorded,rtol=1e-12,atol=1e-12):
            raise ValueError('Actual task supremum/argmax response inconsistent')
    return rows,np.asarray(arrays),best


def check_argmax_arrays(state,q,responses,argmax):
    path=VALIDATION/'argmax_validation/TASK_ARGMAX_GRID_CHECK.npz'
    with np.load(path,allow_pickle=False) as z:
        if int(z['probe_index'])!=argmax or not np.array_equal(z['original_q'],q):
            raise ValueError('Argmax check was not performed at the frozen nominee')
    coarse=load_pair(path,('coarse0','coarse1'))
    fine=load_pair(path,('fine0','fine1'))
    if not np.array_equal(np.asarray(fine),responses[argmax]):
        raise ValueError('Argmax grid check fine responses differ from actual fixed task')
    changes=[float(np.linalg.norm(a-b)/SIGMA_REAL) for a,b in zip(coarse,fine)]
    result=state['argmax']
    for actual,saved in [(norm_difference(coarse),result['coarse_distance']),
        (changes[0],result['world0_grid_distance']),(changes[1],result['world1_grid_distance'])]:
        if not np.isclose(actual,saved,rtol=1e-12,atol=1e-12):
            raise ValueError('Argmax grid metadata is not reproduced by actual arrays')
    if bool(max(changes)<=.1)!=result['task_argmax_numerical_pass']:
        raise ValueError('Argmax grid qualification contradicts saved actual arrays')
    return np.asarray(coarse)


def plot(cfg,fields,grid_x,grid_y,rows,dcal,state,profile):
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import font_manager
    names={f.name for f in font_manager.fontManager.ttflist}
    chosen=next((x for x in ('PingFang SC','Heiti SC','Arial Unicode MS') if x in names),None)
    if chosen is None:raise RuntimeError('A Chinese font is required for this gallery')
    matplotlib.rcParams.update({'font.family':[chosen,'DejaVu Sans'],'svg.fonttype':'path',
        'axes.unicode_minus':False,'font.size':10,'figure.facecolor':'#faf9f5',
        'axes.facecolor':'#faf9f5','savefig.facecolor':'#faf9f5'})
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle
    fig=plt.figure(figsize=(16.5,10.))
    gs=fig.add_gridspec(2,3,left=.055,right=.96,top=.82,bottom=.23,
        width_ratios=[1,1,1.65],wspace=.3,hspace=.20)
    axes=np.array([[fig.add_subplot(gs[i,j]) for j in range(2)] for i in range(2)])
    task=fig.add_subplot(gs[:,2])
    fig.suptitle('严格空气孔洞不能被材料调制抹除：两个校准墙世界与全部24任务',
        x=.045,y=.97,ha='left',fontsize=19)
    fig.text(.045,.90,'ORACLE fixed-support only  |  原 minimum_feature = 0.1 语义 UNRESOLVED  |  nomination 仅按 cal 选定',fontsize=12)
    extent=[grid_x[0],grid_x[-1],grid_y[0],grid_y[-1]]
    for component,label in enumerate(('Re ε','Im ε = σ/f，f=1')):
        arrays=fields.real if component==0 else fields.imag
        lo,hi=float(arrays.min()),float(arrays.max())
        for world in (0,1):
            ax=axes[component,world]
            image=ax.imshow(arrays[world],origin='lower',extent=extent,interpolation='nearest',
                aspect='equal',cmap='viridis',vmin=lo,vmax=hi)
            if world==1:
                ax.add_patch(Circle(CENTER,RADIUS,fill=False,color='#e95131',lw=1.6))
                ax.annotate('严格空气孔',xy=(CENTER[0],CENTER[1]),xytext=(-1.75,1.65),
                    arrowprops={'arrowstyle':'->','color':'#e95131'},fontsize=9,color='#a43220')
            ax.set_xlabel('x（原生坐标）');ax.set_ylabel(label+'\ny（原生坐标）')
            if component==0:ax.set_title('世界0：原墙，无新增孔' if world==0 else '世界1：固定 q + 严格空气孔',fontsize=12,pad=12)
        fig.colorbar(image,ax=list(axes[component,:]),fraction=.025,pad=.025)
    distances=np.array([r['distance_white'] for r in rows]);ids=np.arange(24)
    task.bar(ids,distances,color=np.where(distances>=3,'#387f62','#386ea5'),width=.77)
    task.axhline(3,color='#a43427',ls='--',lw=1.6,label='任务门槛 d ≥ 3')
    best=int(np.argmax(distances))
    task.scatter(best,distances[best],color='#a43427',marker='D',s=40,zorder=5)
    task.annotate(f'argmax probe {best}\ndT = {distances[best]:.5f}',
        xy=(best,distances[best]),xytext=(0.05,.83),textcoords='axes fraction',
        arrowprops={'arrowstyle':'->','color':'#555555'},fontsize=11)
    task.set(xlabel='原公开任务 probe 索引（固定顺序）',ylabel='完整复响应白化距离（无量纲）',
        xticks=ids[::2],ylim=(0,max(3.6,float(distances.max())*1.3)))
    task.set_title(f'全部24个已付费实际任务响应\ncal d = {dcal:.5f} '+('≤ 1' if dcal<=1 else '> 1'),fontsize=13,pad=14)
    task.grid(axis='y',alpha=.16);task.legend(loc='upper right',frameon=False,fontsize=9)
    status=('PASS' if state['argmax']['task_argmax_numerical_pass'] else 'FAIL')
    fig.text(.045,.11,
        f'同一600×450节点网格、每行共同色标；世界1孔内严格 Re ε=1、Im ε=0，中心(-0.85,0.85)，r=0.06；世界0不新增该孔。\n'
        f'距离 = ||F1−F0|| / 0.01，5频率×4源×12接收器，全Re/Im范数；固定gamma=1，实际总场响应，无surrogate、无拼接24任务。\n'
        f'argmax两网格经验检查：{status}；不构成连续PDE误差证书。仅oracle已知支撑构造，非一般未知墙在线成功；模型概率未定义。\n'
        '本次只解码/重算已付数组并绘图：新增RHS=0、LU=0；既有候选搜索与独立验证物理费用已付，并非免费。',
        fontsize=9,ha='left',va='bottom',linespacing=1.65)
    fig.savefig(FIG/'CAVITY_WORLDS_AND_24_TASKS.png',dpi=170)
    fig.savefig(FIG/'CAVITY_WORLDS_AND_24_TASKS.svg')
    plt.close(fig)
    profile['font']=chosen


def main():
    cpu,wall=time.process_time(),time.perf_counter()
    state=completion_state()
    if state['status']!='READY':
        TAB.mkdir(parents=True,exist_ok=True)
        pending=dict(**state,new_RHS=0,new_LU=0,PDE_run=False,Book_constructed=False,
            cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,peak_RSS_bytes=rss_bytes())
        write_json(TAB/'PROFILE.json',pending)
        print(json.dumps(pending),flush=True);return
    cfg,q,dcal,cal_pair,original=checked_sources(state)
    # Only pure public geometry and endpoint decoder are imported after gates.
    sys.path.insert(0,str(ROOT/'src'))
    sys.path.insert(0,str(ROOT.parent/'research_w_continuation/vendor'))
    from p0.geometry import wall_values,probes
    from wr.endpoint_geometry import linear_material_values
    if any(n in sys.modules for n in ('wr.physics','wr.native_world','wr.physical_witness','wr.accounting','wq.fem')):
        raise RuntimeError('Pure renderer imported a forbidden physical/accounting backend')
    task_rows,responses,argmax=task_arrays(state,probes())
    argmax_coarse=check_argmax_arrays(state,q,responses,argmax)
    gx=np.linspace(-2.,2.,600);gy=np.linspace(0.,3.,450)
    xx,yy=np.meshgrid(gx,gy);points=np.c_[xx.ravel(),yy.ravel()]
    decoder=types.SimpleNamespace(wall_values=wall_values)
    fields=[]
    for which,parameters in ((0,np.zeros_like(q)),(1,q)):
        er,sigma=linear_material_values(cfg['seed'],parameters,points,decoder,
            cavities=((CENTER,RADIUS),) if which else (),probe=None)
        fields.append((er+1j*sigma/FREQUENCY).reshape(len(gy),len(gx)))
    fields=np.asarray(fields)
    hole=(np.linalg.norm(points-np.asarray(CENTER),axis=1)<RADIUS).reshape(len(gy),len(gx))
    if not hole.any() or not np.all(fields[1][hole]==1.+0j):
        raise ValueError('Strict-air cavity was erased or lost in the decoder')
    if not np.all(fields[0][hole].real>1):raise ValueError('Original reference is unexpectedly already air at the new cavity')
    profile=dict(status='BUILT_DECODER_AND_ACTUAL_PAID_ARRAYS_ONLY',new_RHS=0,new_LU=0,
        PDE_run=False,Book_constructed=False,backend_imported=False,P1_truth_or_plan_access=False)
    FIG.mkdir(parents=True,exist_ok=True);TAB.mkdir(parents=True,exist_ok=True)
    plot(cfg,fields,gx,gy,task_rows,dcal,state,profile)
    np.savez_compressed(FIG/'MATERIALS_AND_ACTUAL_RESPONSES.npz',grid_x=gx,grid_y=gy,
        frequency=FREQUENCY,epsilon_real=fields.real,epsilon_imag_loss=fields.imag,
        conductivity_sigma=fields.imag*FREQUENCY,strict_cavity_mask=hole,
        reference_parameters=np.zeros_like(q),candidate_parameters=q,
        calibration_actual_pair=np.asarray(cal_pair),task_actual_pairs=responses,
        argmax_actual_coarse_pair=argmax_coarse,
        task_distances_white=np.array([r['distance_white'] for r in task_rows]))
    with (TAB/'ALL24_TASK_DISTANCES.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(task_rows[0]));writer.writeheader();writer.writerows(task_rows)
    summary=dict(d_cal=dcal,cal_threshold_le1=bool(dcal<=1),d_T_sup_24=task_rows[argmax]['distance_white'],
        task_argmax=argmax,task_threshold_ge3=bool(task_rows[argmax]['distance_white']>=3),
        combined_cal_task_threshold=bool(dcal<=1 and task_rows[argmax]['distance_white']>=3),
        base_validation_numerical_pass=state['result']['numerical_pass'],
        argmax_grid_numerical_pass=state['argmax']['task_argmax_numerical_pass'],
        original_minimum_feature='UNRESOLVED',strict_original_prior_GO=False,
        strict_air_cavity=dict(center=CENTER,radius=RADIUS,retained=True,decoded_interior_pixels=int(hole.sum()),
            exact_sample_values=dict(epsilon_real=1.,conductivity=0.,loss=0.)),
        scope='ORACLE fixed-support only; cal-only nominated fixed pair; no online unknown-wall or posterior claim',
        science_decision='PARENT_REVIEW_ONLY')
    write_json(TAB/'NUMBERS.json',summary)
    sources=['configs/SECOND_CAVITY_NOMINATION.json',cfg['x_file'],'configs/ORIGINAL_P0.json',
        f'results/B2/validation/{IDENTIFIER}/MANIFEST.json',
        f'results/B2/validation/{IDENTIFIER}/RESULTS.json',
        f'results/B2/validation/{IDENTIFIER}/ALL_TASKS.json',
        f'results/B2/validation/{IDENTIFIER}/cal_PAIR.npz',
        f'results/B2/validation/{IDENTIFIER}/argmax_validation/TASK_ARGMAX_GRID_CHECK.json',
        f'results/B2/validation/{IDENTIFIER}/argmax_validation/TASK_ARGMAX_GRID_CHECK.npz']
    sources += [r['source'] for r in task_rows]
    write_json(TAB/'SOURCES.json',dict(inputs=sources,
        decoder=['src/wr/endpoint_geometry.py:linear_material_values',
        '../research_w_continuation/vendor/p0/geometry.py:wall_values/probes'],
        display='Same600x450 nodal grid; no interpolation smoothing; row-shared colorbars; f=1',
        noise='sigma=.01 for each real component, no channel normalization; original covariance unchanged',
        unit='Material relative epsilon/loss in original normalized convention; spatial native coordinates, not silently metres',
        task_selection='All24 fixed public tasks, no omission or surrogate; argmax checked after cal-only nomination',
        source_physics_cost=dict(original_validation_B2_cumulative_snapshot=state['result'].get('totals',{}).get('B2'),
            snapshot_scope='Stage aggregate, not attributed as this renderer/candidate marginal cost',
            argmax_previously_paid_actual_cost=state['argmax'].get('actual_cost'),
            source_search_and_validation_were_already_paid=True),
        renderer_incremental_physics=dict(rhs=0,LU=0),summary=summary))
    profile.update(cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
        peak_RSS_bytes=rss_bytes(),decoded_complex_field_array_bytes=int(fields.nbytes),
        saved_task_response_array_bytes=int(responses.nbytes),output=str(FIG))
    write_json(TAB/'PROFILE.json',profile)
    print(json.dumps(dict(**summary,profile=profile),ensure_ascii=False),flush=True)


if __name__=='__main__':
    main()

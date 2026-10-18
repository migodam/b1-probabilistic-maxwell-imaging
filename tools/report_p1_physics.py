"""Analyze persisted P1 qualification arrays and publish them; no solver calls."""
from pathlib import Path
import json,sys,time,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p1_material.prior import MaterialPrior,topology
from p1_material.target import pack,write_json
from publish_lxz import STYLE
from publish_research_f import rendered

OUT=ROOT/'research/p1_material_target_20261018'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p1-physics-qualification-20261018'


def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(OUT/f'figures/{name}.{ext}',dpi=175,bbox_inches='tight')
    plt.close(fig)


def main():
    start=time.perf_counter()
    if not (OUT/'PHYSICS_CLOSED.json').exists():raise RuntimeError('only report a completed physical job')
    cfg=json.loads((OUT/'configs/ONLINE_TARGET.json').read_text());points=json.loads((OUT/'inputs/geometry.json').read_text())['points'];P=MaterialPrior(cfg,points);sigma=cfg['noise']['sigma']
    rows=[json.loads(x) for x in (OUT/'PHYSICS_LEDGER.jsonl').read_text().splitlines()];job=[x for x in rows if x['event']=='JOB_END'][-1]
    endpoints=[]
    for path in sorted((OUT/'results/physical/endpoints').glob('*.npz')):
        with np.load(path) as z:state={k:z[k].copy() for k in ['q','chi','prediction','current','field','backward_residual']}
        index=int(path.stem.split('_')[-1]);state['path']=str(path.relative_to(ROOT))
        with np.load(OUT/f'results/physical/attempts/input_{index:04}.npz') as z:state['purpose']=str(z['purpose'])
        state['shared_material_constitutive_relative_error']=float(np.linalg.norm(state['current']-state['chi'][:,None]*state['field'])/max(np.linalg.norm(state['current']),1e-300));endpoints.append(state)
    jacobians=[]
    for path in sorted((OUT/'results/physical/jacobians').glob('*.npz')):
        with np.load(path) as z:jacobians.append(dict(q=z['q'].copy(),A=z['A'].copy(),path=str(path.relative_to(ROOT))))
    # Offline scoring starts only after the inverse program and its ledger closed.
    with np.load(OUT/'inputs/offline_materials/dev.npz') as z:qtruth=z['latent'][:2].copy();truth=z['chi'][:2].copy()
    q0=np.r_[np.zeros(64),.5,.5];initial=P.decode(q0)[0];allmetrics=[];plotdata=[]
    compact={}
    for i,name in enumerate(['P1_obs_A','P1_obs_B']):
        r=json.loads((OUT/f'results/physical/{name}_qualification.json').read_text());history=json.loads((OUT/f'results/physical/{name}_history.json').read_text());q=np.array(r['q'])
        state=next(x for x in endpoints if np.array_equal(x['q'],q) and x['purpose'].startswith(name));truestate=next(x for x in endpoints if x['purpose']=='dataset_truth_'+name)
        with np.load(OUT/f'inputs/online/{name}.npz') as z:y=z['observation'].copy()
        chi,a=P.decode(q);_,atruth=P.decode(qtruth[i]);J=next(x for x in jacobians if np.array_equal(x['q'],q));sv=np.linalg.svd(J['A'],compute_uv=False)
        invstates=[x for x in endpoints if x['purpose'].startswith(name)];steps=[];secants=[]
        for old,new in zip(invstates[:-1],invstates[1:]):
            step=new['q']-old['q'];steps.append(float(np.linalg.norm(step[:64])))
            previous=next(x for x in jacobians if np.array_equal(x['q'],old['q']));actual=pack(new['prediction']-old['prediction'])*np.sqrt(2)/sigma;linear=previous['A']@step
            secants.append(dict(actual_white_data_step=float(np.linalg.norm(actual)),linear_white_data_step=float(np.linalg.norm(linear)),relative_linear_secant_error=float(np.linalg.norm(actual-linear)/max(np.linalg.norm(actual),1e-300))))
        r.update(dict(parent_id=['P1_dev_000','P1_dev_001'][i],offline_GT_scoring_only=True,full_complex_MAE=float(np.mean(abs(chi-truth[i]))),full_complex_RMS=P.material_rms(chi,truth[i]),Re_MAE=float(np.mean(abs(chi.real-truth[i].real))),Im_MAE=float(np.mean(abs(chi.imag-truth[i].imag))),initial_complex_MAE=float(np.mean(abs(initial-truth[i]))),initial_white_residual=float(np.sqrt(2*history[0]['Phi'])),last_projected_gradient_inf=history[-1]['projected_gradient_inf'],last_material_rank=history[-1]['material_rank'],white_J_singular_values=sv.tolist(),white_J_rank_relative_1e8=int(sum(sv>sv[0]*1e-8)),active_property_bounds=dict(lower=np.flatnonzero(q[64:]<=1e-10).tolist(),upper=np.flatnonzero(q[64:]>=1-1e-10).tolist()),latent_Gaussian_step_lengths=steps,linear_secant_diagnostics=secants,noise_at_generating_material_white_norm=float(np.linalg.norm(pack(truestate['prediction']-y))*np.sqrt(2)/sigma),total_scattered_signal_RMS=float(np.sqrt(np.mean(abs(truestate['prediction'])**2))),actual_total_field_SNR_dB=float(20*np.log10(np.sqrt(np.mean(abs(truestate['prediction'])**2))/sigma)),topology_reconstruction=topology(a['fraction']),topology_truth_OFFLINE=topology(atruth['fraction']),qualified_state_array=state['path'],truth_state_array_OFFLINE=truestate['path'],failure_reason=None if r['qualified'] else '12-state per-parent allocation exhausted while feasible projected gradient remains large; no absence-of-solution claim',new_discovered_alternative_worlds=0,posterior_mass=None))
        allmetrics.append(r);plotdata.append((name,r,history,chi,truth[i],state['prediction'],truestate['prediction'],y))
        compact.update({name+'_q':q,name+'_chi':chi,name+'_prediction':state['prediction'],name+'_observation':y,name+'_truth_chi_OFFLINE':truth[i],name+'_truth_prediction_OFFLINE':truestate['prediction']})
    np.savez_compressed(OUT/'results/physical/QUALIFICATION_COMPACT_OFFLINE.npz',**compact)
    metrics=dict(target_id=cfg['target_id'],executed_source_commit='117359b',authorization='configs/AUTHORIZED_PHYSICS_20261018.json',parents=allmetrics,qualified_anchors=sum(x['qualified'] for x in allmetrics),new_discovered_alternative_worlds=0,physical_JOB=job,physical_endpoints_persisted=len(endpoints),paid_receiver_Jacobians_persisted=len(jacobians),max_backward_residual=max(float(x['backward_residual']) for x in endpoints),max_shared_material_constitutive_error=max(x['shared_material_constitutive_relative_error'] for x in endpoints),GT_used_online=False,posthoc_analysis_new_solver_calls=0,NN_GPU_Mie_new3D=0,cold_warm_method_comparison='NOT_RUN',A17_X_competition='NOT_RUN',posterior_reference='NOT_RUN',stage_closed=True,reserve_unspent=dict(states=32-job['totals']['states'],RHS=1024-job['totals']['RHS']))
    write_json(OUT/'PHYSICS_METRICS.json',metrics)
    write_json(OUT/'PHYSICS_DECISION.json',dict(material_target='ACTUALLY_RUN',qualified_anchors=metrics['qualified_anchors'],A='QUALIFIED_FEASIBLE_BASE_NOT_STATIONARY_MAP',B='UNQUALIFIED_NONCONVERGED_AT_FROZEN_LIMIT',new_alternative_worlds=0,A17_X_increment='NOT_RUN',next='retain A as a real anchor; diagnose B optimization without treating it as physical impossibility; no new paid comparison authorized',new_physics_states=job['totals']['states'],new_complex_RHS=job['totals']['RHS'],no_extra_iterations=True))
    write_json(OUT/'evidence/PHYSICS_ENDPOINT_INVENTORY.json',dict(rows=[dict(path=x['path'],purpose=x['purpose'],backward_residual=float(x['backward_residual']),shared_material_constitutive_relative_error=x['shared_material_constitutive_relative_error']) for x in endpoints],all_sources_share_one_chi=True,offline_analysis_only=True))
    for name,r,h,x,t,p,pt,y in plotdata:
        fig,ax=plt.subplots(2,3,figsize=(11,7.5),constrained_layout=True)
        for col,(a,label) in enumerate([(initial,'Independent prior-mean init'),(x,'Full-J bounded LM'),(t,'Generating material / OFFLINE')]):
            for row,(value,lim,cmap) in enumerate([(1+a.real,(1.1,2.9),'viridis'),(a.imag,(.04,.22),'magma')]):
                im=ax[row,col].imshow(value.reshape(24,24).T,origin='lower',extent=[-.5,.5,-.5,.5],vmin=lim[0],vmax=lim[1],cmap=cmap);ax[row,col].set_xlabel('x [m]')
                if col==0:ax[row,col].set_ylabel(('Re epsilon_r' if row==0 else 'Im epsilon_r (loss)')+'\ny [m]')
                if row==0:ax[row,col].set_title(label)
                if col==2:fig.colorbar(im,ax=ax[row,:],shrink=.85)
        fig.suptitle(f"{name}: qualified={r['qualified']}; white residual={r['white_residual']:.3f} / tau={r['tau']:.3f}\nSingle base reconstruction, not multi-world discovery; truth never enters inversion",fontsize=12)
        save(fig,'04_physics_'+name)
    fig,ax=plt.subplots(1,2,figsize=(11,4.2),constrained_layout=True)
    for name,r,h,*_ in plotdata:
        ax[0].semilogy([x['states_used'] for x in h],[x['white_residual'] for x in h],'o-',label=name)
        ax[1].semilogy([x['states_used'] for x in h],[x['projected_gradient_inf'] for x in h],'o-',label=name)
    ax[0].axhline(allmetrics[0]['tau'],color='k',ls='--',label='Frozen noise-ball tau');ax[0].set(xlabel='Per-parent full states used',ylabel='Exact full-Maxwell whitened residual');ax[0].legend()
    ax[1].set(xlabel='Per-parent full states used',ylabel='Projected posterior-objective gradient (inf norm)');ax[1].legend();fig.suptitle('Frozen 12-state cap; no extra iteration to rescue case B');save(fig,'05_physics_convergence')
    fig,ax=plt.subplots(2,2,figsize=(12,7),constrained_layout=True)
    for row,(name,r,h,x,t,p,pt,y) in enumerate(plotdata):
        ax[row,0].plot(y[:,0].real,'k.',label='Noisy Re, source 0');ax[row,0].plot(p[:,0].real,label='Reconstructed Re');ax[row,0].plot(pt[:,0].real,'--',label='Generating Re / offline');ax[row,0].set_title(name+' — representative source only');ax[row,0].legend(fontsize=8)
        d=abs((p-y).T)*np.sqrt(2)/sigma;im=ax[row,1].imshow(d,origin='lower',aspect='auto');fig.colorbar(im,ax=ax[row,1]);ax[row,1].set(xlabel='Receiver',ylabel='Source',title='All-source |complex residual| × sqrt(2)/sigma')
    save(fig,'06_physics_predictions')
    a,b=allmetrics
    text=f'''# P1真实Maxwell资格：取得1个合格基础解，另1个优化未完成

2026-10-18。本轮已按新增明确授权真实运行，不再只是材料生成。模型、prior、sigma、gain、两个parent与12states/parent限制均沿用已发布协议。未重跑旧P0/F3，未训练NN。

## 实际结果

| 条件 | 初始白化残差 | 最终白化残差 | noise-ball上限 | 资格 | 最后投影梯度∞范数 |
|---|---:|---:|---:|---|---:|
| A / dev000 | {a['initial_white_residual']:.4f} | {a['white_residual']:.4f} | {a['tau']:.4f} | PASS：合法且数据相容 | {a['last_projected_gradient_inf']:.6g} |
| B / dev001 | {b['initial_white_residual']:.4f} | {b['white_residual']:.4f} | {b['tau']:.4f} | FAIL：未相容 | {b['last_projected_gradient_inf']:.6g} |

两者均从z=0、u=.5的共同先验均值开始。在线反演只读取y、sigma和target identity，不读真实latent、GT、N3或结构标签。两份观测按事先保存的noise生成，无挑选seed或放大噪声。

**目前有1个合格P1 physics-only基础材料、0个新发现的替代世界。** A不是认证MAP，noise-ball也不是posterior支持边界；没有posterior概率或多模式结论。

![A实际重建](figures/04_physics_P1_obs_A.png)

![B失败保留](figures/04_physics_P1_obs_B.png)

## B为什么失败，当前能排除什么？

B的Φ由152020.76降至47984.08，投影梯度仍约2.24e5；每一步都降低目标，全部11次Gaussian更新长度均触及.5上限，但12states限额先到，不能视为收敛。固定Gaussian步长限幅与复杂材料形状的优化是当前瓶颈候选，尚未隔离各自影响；这也不保证追加迭代必然成功。

生成材料与反演使用同一P1、同一full-forward、mesh、source、receiver和固定gain/noise。闭账后仅用已保存数组检查：真实生成材料的白化噪声范数为A={a['noise_at_generating_material_white_norm']:.6f}、B={b['noise_at_generating_material_white_norm']:.6f}，均低于固定阈值。因此B的本轮失败不支持“先验里没有可行解”或“物理唯一性阻止反演”的结论。没有用这些GT信息修改在线初始化或重跑。

sigma仍为.00119211798514519，没有为新物体重新定标。按总散射响应RMS/sigma定义的实际SNR分别{a['actual_total_field_SNR_dB']:.3f}/{b['actual_total_field_SNR_dB']:.3f}dB，不是每物体固定20dB；总散射包含已知宿主响应。

![真实收敛轨迹](figures/05_physics_convergence.png)

## 物理与材料质量

26个端点全已立即保存，最大状态backward residual为{metrics['max_backward_residual']:.3g}；多源共享材料本构误差最大{metrics['max_shared_material_constitutive_error']:.3g}。完整白化J来自32个receiver共享adjoint，共24份，原数组均保留。

离线复材料MAE：A由{a['initial_complex_MAE']:.6f}降至{a['full_complex_MAE']:.6f}，B由{b['initial_complex_MAE']:.6f}降至{b['full_complex_MAE']:.6f}。这只用于评分，不用于选择endpoint。A与GT的材料RMS仍为{a['full_complex_RMS']:.6f}，说明数据相容不等于精确材料恢复；不能因此将真实生成材料和逆解自动称两个已发现的重要posterior世界。

A损耗参数处于允许上界，投影梯度与active bounds在PHYSICS_METRICS中记录；未认证驻点。两例末端Dg数值rank均66，不把这些材料维数或小奇异值当exact-fiber存在证据。

![实际散射拟合](figures/06_physics_predictions.png)

## 全部成本与停止

新增26/32个full-material states、976/1024个complex full-solver RHS：208 forward RHS + 768 receiver-adjoint RHS；26次LU。受控job wall={job['wall_seconds']:.6f}s，含两truth、初始化、全部反演、导数、检查与端点持久化。不是算法独立cold/warm计时，也不是完整研究/上传host时间。旧paid背景导入单列，历史setup不当作免费。

没有数值异常或line-search拒绝；B属于科学/优化资格失败，仍计费。剩余6states/48RHS是原失败预留，不转用来追加迭代或寻找正例，阶段封账。物理之后的分析只读已有数组，没有新LU/RHS/forward。

## 科学决定

继续保留A作为真实P1物理anchor；B标记NONCONVERGED_AT_FROZEN_LIMIT。P1比旧双Gaussian chart更丰富，但现在不足以证明多解机制或NN价值。

F3-native、A17-X finite feedback、exact T-matrix、普通Full-J和goal-oriented ROM的共同机制比较仍NOT_RUN：本次授权只覆盖观测与独立物理资格，不以预留额度启动另一campaign。下一有价值的实验是在已合格A上冻结同材料位移/同成本对照；B应先做零新增物理的步长、局部模型与已有Jacobian分析，不能无限增加迭代。没有A17增量或强TAP算法贡献结论。

[机器可读结果](PHYSICS_METRICS.json) · [裁决](PHYSICS_DECISION.json) · [逐项账本](PHYSICS_LEDGER.jsonl) · [原始数组]({RELEASE})
'''
    (OUT/'PHYSICS_QUALIFICATION_REPORT_ZH.md').write_text(text)
    nav='<nav><a href="../../">B1首页</a><a href="index.html">P1入口</a><a href="physics_report.html">真实物理报告</a><a href="physics_gallery.html">实际重建图</a><a href="report.html">材料模型报告</a><a href="'+RELEASE+'">物理证据</a></nav>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P1 · 真实Maxwell资格实验</title><style>'+STYLE+'</style><main>'+nav
    (OUT/'physics_report.html').write_text(top+rendered(text)+'</main></html>')
    gallery=top+'<h1>P1实际Maxwell重建：1/2合格，无新多解发现</h1><p>GT仅在物理作业闭账后的离线评分使用。图像共用材料色标，不把未合格B隐藏，不将noise-ball解释为posterior质量。</p>'
    for name,title in [('04_physics_P1_obs_A','A：合法、数据相容的独立物理基础解'),('04_physics_P1_obs_B','B：同限制下未完成的反演'),('05_physics_convergence','实际残差与投影梯度'),('06_physics_predictions','散射预测与全部8源残差')]:
        gallery+='<section><h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+title+'"><a href="figures/'+name+'.pdf">PDF</a></section>'
    (OUT/'physics_gallery.html').write_text(gallery+'</main></html>')
    banner='<section><h2>更新：真实Maxwell资格已运行</h2><div class="metrics"><span><b>1 / 2</b>合格基础材料</span><span><b>26</b>新states</span><span><b>976</b>complex RHS</span><span><b>5.48 s</b>受控job wall</span></div><p>A残差22.29通过阈值23.79；B残差309.79未合格。未追加迭代，无新替代world；A17-X机制比较尚未运行。</p><p><a href="physics_report.html">真实物理中文报告 →</a> · <a href="physics_gallery.html">实际重建图 →</a></p></section>'
    path=OUT/'index.html';old=path.read_text()
    if '更新：真实Maxwell资格已运行' not in old:old=old.replace('<main>','<main>'+banner,1)
    old=old.replace('<h1>P1已经实现：材料支持就绪，真实Maxwell资格尚未运行</h1>','<h1>P1：材料模型与第一轮真实Maxwell资格</h1>').replace('在线观测与独立物理anchor仍为NOT_RUN。','下方是原材料构造阶段记录；最新真实物理结果见页首更新。').replace('当前物理预算为0，调用在导入物理backend前终止。','默认预算仍为0；本次明确授权通过独立overlay执行，现已封账。').replace('第一轮资格建议：32个states、1024个complex RHS、600秒CPU。只在新的明确授权后执行；历史预算不借用。','第一轮资格已按独立授权完成：26states/976RHS；历史预算未借用。');path.write_text(old)
    notice='## 2026-10-18更新：真实P1资格已运行\n\n[真实物理报告](PHYSICS_QUALIFICATION_REPORT_ZH.md) · [网页](physics_report.html) · [实际重建图](physics_gallery.html)\n\n26states / 976RHS / 5.477s封账：两个固定对象取得1个合法且数据相容的physics-only基础解，另1个未收敛。没有新发现替代world；A17-X比较未运行。下面保留原零物理材料构造记录，不能把历史状态当最新结果。\n\n'
    path=OUT/'START_HERE.md'
    if not path.read_text().startswith('## 2026-10-18更新：真实P1资格已运行'):path.write_text(notice+path.read_text())
    original=OUT/'evidence/MATERIAL_CONSTRUCTION_DECISION.json'
    if not original.exists():original.write_text((OUT/'DECISION.json').read_text())
    decision=json.loads(original.read_text());decision.update(dict(historical_material_decision='evidence/MATERIAL_CONSTRUCTION_DECISION.json',actual_observations='COMPLETED_TWO_FIXED_PARENTS',independent_physics_anchors='ONE_OF_TWO_QUALIFIED',ambiguity_witness='NOT_RUN',A17_increment='NOT_RUN',new_physics_authorization='AUTHORIZED_AND_CLOSED',new_physical_states=26,new_complex_RHS=976,next_required='no automatically authorized continuation; scientific interpretation in PHYSICS_DECISION.json'));write_json(OUT/'DECISION.json',decision)
    binding=json.loads((OUT/'A17X_P1_BINDING.json').read_text());binding.update(dict(actual_observations='COMPLETED',physics_only_anchors=dict(qualified=['P1_obs_A'],unqualified=['P1_obs_B']),qualification_cost=dict(states=26,RHS=976,controlled_wall_seconds=job['wall_seconds']),future_A17_competition_status='NOT_RUN_NOT_AUTHORIZED',source='PHYSICS_METRICS.json'));write_json(OUT/'A17X_P1_PHYSICS_STATUS.json',binding)
    path=OUT/'PHYSICS_NEXT_RUN.md'
    if not path.read_text().startswith('> 历史冻结建议'):path.write_text('> 历史冻结建议：用户随后明确授权，本资格阶段现已完成并封账。实际结果为26states/976RHS、1/2合格；见[真实物理报告](PHYSICS_QUALIFICATION_REPORT_ZH.md)。以下原预算建议保留，不代表当前仍未授权，也不授予后续比较额度。\n\n'+path.read_text())
    write_json(OUT/'P1_RESOURCE_SUMMARY.json',dict(material_stage='RESOURCE_LEDGER.jsonl; zero Maxwell; successful analysis .976s, first failure wall unmeasured',physical_stage='PHYSICS_LEDGER.jsonl',physical_states=26,physical_LU=26,forward_RHS=208,receiver_adjoint_RHS=768,total_complex_RHS=976,controlled_physics_job_wall=job['wall_seconds'],independent_cold_warm_times='NOT_RUN',posthoc_scoring_new_solves=0,NN_GPU_Mie_new3D=0,legacy_setup='known paid background; cost not remeasured',host_research_and_upload_wall='NOT_MEASURED',stage_closed=True))
    for name in ['README.md','START_HERE.md']:
        path=ROOT/name;s=path.read_text();s=s.replace('192份先验材料、24个材料侧运动端点已保存；真实P1观测和物理anchor尚未运行。新增Maxwell states/RHS均0，不训练NN。','材料构造已完成；新增真实物理资格26states/976RHS：1/2基础材料合格，B未收敛。无新替代world，不训练NN。[真实物理报告](research/p1_material_target_20261018/PHYSICS_QUALIFICATION_REPORT_ZH.md)。');path.write_text(s)
    for name in ['index.html','research/index.html']:
        path=ROOT/name;s=path.read_text().replace('192份实际先验材料，0新物理；观测与anchor待独立授权。','真实资格26states/976RHS，1/2基础材料合格；未发现新多解。').replace('模型、实际材料数组和付费接口已实现；物理资格尚未运行。','真实P1已运行：1/2合格基础材料，26states/976RHS。');path.write_text(s)
    archive=ROOT/'publication_local/P1_Physics_Qualification_20261018_Evidence.zip'
    excluded={'physics.lock','PUBLICATION_RECEIPT.json','WEB_READBACK.json','EVIDENCE_RELEASE.json','PHYSICS_PUBLICATION_RECEIPT.json'}
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=4) as z:
        for path in sorted(OUT.rglob('*')):
            if path.is_file() and path.name not in excluded and '__pycache__' not in path.parts:z.write(path,path.relative_to(ROOT))
        for path in (ROOT/'src/p1_material').glob('*.py'):z.write(path,path.relative_to(ROOT))
        for name in ['run_p1_physics.py','report_p1_physics.py','prepare_p1_material_target.py','publish_p1_material.py','publish_lxz.py','publish_research_f.py']:z.write(ROOT/'tools'/name,'tools/'+name)
    with zipfile.ZipFile(archive) as z:count=len(z.namelist())
    write_json(OUT/'PHYSICS_EVIDENCE_RELEASE.json',dict(asset=archive.name,bytes=archive.stat().st_size,files=count,full_endpoints=26,receiver_Jacobians=24,experiment_source_commit='117359b',release=RELEASE,postprocessing_no_new_solver_calls=True))
    print(json.dumps(dict(qualified=metrics['qualified_anchors'],totals=job['totals'],wall=job['wall_seconds'],archive=str(archive),archive_bytes=archive.stat().st_size,reporting_wall_seconds=time.perf_counter()-start),ensure_ascii=False))


if __name__=='__main__':main()

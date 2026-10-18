"""Render completed zero-physics P1 audits and frozen protocols; no solver."""
from pathlib import Path
import sys,json,platform,zipfile,html
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from publish_lxz import STYLE
from publish_research_f import rendered
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p1_zero_physics_audit_20261018';REL='research/p1_zero_physics_audit_20261018/'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p1-zero-physics-audit-20261018'


def write(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(OUT/f'figures/{name}.{ext}',dpi=175,bbox_inches='tight')
    plt.close(fig)


def main():
    a=json.loads((OUT/'A_ONLINE_GEOMETRY.json').read_text());t=json.loads((OUT/'A_DECODER_ONLY_TRIALS.json').read_text());off=json.loads((OUT/'OFFLINE_GT_PAIR_EVALUATION.json').read_text());b=json.loads((OUT/'B_OPTIMIZER_CACHED_AUDIT.json').read_text())['rows']
    (OUT/'figures').mkdir(exist_ok=True)
    fig,ax=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    ax[0].semilogy(a['latent_response_singular_values'],'o-',ms=3);ax[0].set(title='Whitened data / unit latent',xlabel='Ordered direction',ylabel='Singular value (latent metric)')
    ax[1].semilogy(a['material_normalized_response_singular_values'],'o-',ms=3,label='Data only');ax[1].semilogy(a['prior_aware_material_singular_values'],'o-',ms=3,label='Data + latent-prior quadratic');ax[1].axhline(50,color='k',ls='--',label='.02 RMS gives unit local norm');ax[1].set(title='Unit actual material RMS',xlabel='Ordered direction',ylabel='Local quadratic norm / material RMS');ax[1].legend(fontsize=8)
    ax[2].semilogy([x['unit_latent_material_RMS'] for x in a['weak_latent_directions']],'o-');ax[2].axhline(.02,color='k',ls='--');ax[2].set(title='Weak latent directions barely move material',xlabel='Weak latent direction order',ylabel='Material RMS / unit latent')
    fig.suptitle('A: coordinate conditioning is not Maxwell uncertainty; no exact-fiber certificate');save(fig,'01_material_geometry')
    fig,ax=plt.subplots(1,2,figsize=(12,4),constrained_layout=True);colors={'raw_latent_weak':'#b96e3f','material_metric_weak':'#8064a2','prior_aware_material_weak':'#246b81'}
    for family,color in colors.items():
        rr=[x for x in t if x['direction'].startswith(family)]
        ax[0].semilogy(range(len(rr)),[-x['latent_prior_log_change'] for x in rr],'o-',color=color,label=family)
        ax[1].plot([x['target_RMS'] for x in rr],[x['predicted_DeltaPhi_material_linear'] for x in rr],'o',color=color,label=family)
    ax[0].set(xlabel='Direction × requested distance',ylabel='Loss of latent log-prior density',title='Same actual material distance, very different prior cost');ax[0].legend(fontsize=8)
    ax[1].axhline(2,color='k',ls='--');ax[1].set(xlabel='Actual requested material RMS',ylabel='Cached material-J predicted DeltaPhi',title='Finite decoder changes; NOT new Maxwell predictions');ax[1].legend(fontsize=8);save(fig,'02_decoder_trials')
    with np.load(OUT/'OFFLINE_GT_PAIR_ARRAYS.npz') as z:xa=z['anchor_material'];xg=z['GT_material'];delta=z['white_prediction_difference'];ml=z['material_linear_prediction'];ll=z['latent_linear_prediction']
    with np.load(OUT/'A_ONLINE_GEOMETRY.npz') as z:offline_residual=z['residual'].copy()
    write('OFFLINE_MATERIAL_LINEAR_PHI.json',dict(role='offline postprocessing only',Phi_material_linear=float(.5*np.linalg.norm(offline_residual+ml)**2),Phi_stored_full=float(.5*np.linalg.norm(offline_residual+delta)**2),linear_Phi_error=float(.5*(np.linalg.norm(offline_residual+ml)**2-np.linalg.norm(offline_residual+delta)**2)),relative_finite_Maxwell_remainder=float(np.linalg.norm(delta-ml)/np.linalg.norm(delta)),new_physics=0))
    fig,ax=plt.subplots(1,3,figsize=(13,4),constrained_layout=True)
    for i,(x,title) in enumerate([(xa,'A anchor / measured already'),(xg,'Generating material / OFFLINE only')]):
        im=ax[i].imshow((1+x.real).reshape(24,24).T,origin='lower',extent=[-.5,.5,-.5,.5],vmin=1.1,vmax=2.9,cmap='viridis');ax[i].set(title=title,xlabel='x [m]',ylabel='y [m]')
    fig.colorbar(im,ax=ax[:2],shrink=.8,label='Re epsilon_r')
    ax[2].bar(['Latent J dq','Material J h','Stored full secant'],[np.linalg.norm(ll),np.linalg.norm(ml),np.linalg.norm(delta)],color=['#b96e3f','#246b81','#406e47']);ax[2].set_yscale('log');ax[2].set(ylabel='Whitened response norm',title='Finite decoder distortion dominates this pair');ax[2].tick_params(axis='x',rotation=20)
    fig.suptitle('OFFLINE pair: material RMS .08071; prediction distance 5.637; not a discovered alternate mode');save(fig,'03_offline_pair')
    fig,ax=plt.subplots(1,3,figsize=(14,4),constrained_layout=True);bb=b[:-1]
    ax[0].plot([x['unscaled_GN_Gaussian_step'] for x in bb],'o-',label='Bounded GN before postscale');ax[0].plot([x['actual_Gaussian_step'] for x in bb],'o-',label='Actually executed');ax[0].set(xlabel='B update',ylabel='Gaussian coordinate length',title='Uniform postscale wastes the subproblem solution');ax[0].legend(fontsize=8)
    ax[1].plot([x['trust_agreement_ratio'] for x in bb],'o-');ax[1].axhline(.75,color='k',ls='--');ax[1].set(xlabel='B update',ylabel='Actual / predicted objective decrease',title='First 8 updates support radius adaptation')
    ax[2].plot([x['material_RMS_step'] for x in bb],'o-');ax[2].set(xlabel='B update',ylabel='Actual material RMS step',title='Same latent cap does not fix material distance');save(fig,'04_B_optimizer')
    write('RESULT_SUMMARY.json',dict(target_id='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018',new_Maxwell_states=0,new_complex_RHS=0,new_LU=0,GT_online=False,material_Dg_rank=66,latent_J_condition=a['latent_response_condition_number'],material_normalized_J_condition=a['material_normalized_response_condition_number'],local_linear_noise_scale_directions=38,prior_aware_quadratic_scale_directions=13,new_discovered_worlds=0,GT_pair_prediction_distance=off['whitened_prediction_distance'],fixed_pair_prospective_error=off['fixed_pair_equal_prior_optimal_classification_error'],posterior_event_mass=None,all_new_physical_experiments='NOT_RUN_NOT_AUTHORIZED',next_budget_proposal=dict(A_states=12,A_RHS=288,A_wall=300,B_states=16,B_RHS=640,B_wall=600,total_states=28,total_RHS=928,total_wall=900),NN_GPU_Mie_new3D=0))
    write('COST_AND_SOURCE.json',dict(source_B1_commit='3d01e03a0848e5546e2ff23f9caee22b407fbfa7',online_outputs_frozen_before_offline='48f11e8',physical_target_source='P1 material/physics commits 4f71c0a,117359b,2c97d39',A17_X_source='research/a17x_f3_competition_audit/ at263f3f0',physical_cache_cost=dict(states=26,complex_RHS=976,forward_RHS=208,adjoint_RHS=768,LU=26,wall_seconds=5.476945083006285),new_physics_cost=dict(states=0,complex_RHS=0,LU=0),new_analysis_wall='RESOURCE_LEDGER.jsonl; full human/agent/host time not measured',background_setup='reused paid; not remeasured or declared free',source_cache_release='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p1-physics-qualification-20261018',historical_results_and_ledgers_unchanged=True,no_hash_checks=True))
    write('DEPENDENCIES.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=__import__('scipy').__version__,matplotlib=matplotlib.__version__,physical_backend_imports=0))
    nav='<nav><a href="../../">B1首页</a><a href="index.html">本轮入口</a><a href="report.html">中文报告</a><a href="gallery.html">数值图</a><a href="protocol.html">冻结协议</a><a href="START_HERE.md">GPT入口</a><a href="'+RELEASE+'">证据</a></nav>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P1零物理审计 · 材料弱方向与优化</title><style>'+STYLE+'</style><main>'+nav
    body='''<h1>P1：真实弱方向仍在，但没有重要结构模式见证</h1><p>2026-10-18 · 只用已付费缓存 · 新Maxwell/RHS/LU均0 · 无NN</p><section><h2>四个直接结论</h2><p>A的material-normalized局部弱方向不全是饱和伪影；原prior-aware三条方向在RMS .02有小先验代价，但尚无full-Maxwell验证或显著结构变化。</p><p>A和GT预测白化距离5.637，不能因为都过noise-ball就说不可区分。固定二均值的前瞻误判率约.241%，不是当前数据p值或posterior概率。</p><p>有限latent外推预测1212，完整材料一阶5.52，缓存真实差值5.64；此pair的主要失真来自decoder/坐标，不能把它全部叫高阶Maxwell feedback。</p><p>B原GN步被缩至约6.5%–11.2%；前8步模型下降比>.75。值得测试直接有界adaptive TR，尚不保证能收敛。</p><a href="report.html">完整审计与证据边界 →</a></section><section><h2>新的实验仅冻结，不执行</h2><div class="metrics"><span><b>12 / 288</b>A拟议states / RHS</span><span><b>16 / 640</b>B拟议states / RHS</span><span><b>900 s</b>合计拟议物理wall</span><span><b>0</b>当前已授权新物理</span></div><p>A用material-native和强Full-J比较.02 finite-noise移动；B比较旧LM与adaptive TR。A17-X保留为独立current enrichment竞争机制，不借旧余额，不省略exact T-matrix/goal-ROM，不训练NN。</p><a href="protocol.html">冻结协议、预算与停止 →</a></section><section><h2>坐标几何与原prior必须一起看</h2><img src="figures/01_material_geometry.png" alt="真实材料metric的谱与latent弱方向"><p>数值rank66不等于存在exact fiber；38个线性尺度方向或13个prior-aware二次方向不等于后验模式。</p></section>'''
    (OUT/'index.html').write_text(top+body+'</main></html>');(OUT/'report.html').write_text(top+rendered((OUT/'RESEARCH_REPORT_ZH.md').read_text())+'</main></html>');(OUT/'protocol.html').write_text(top+rendered((OUT/'FROZEN_EXPERIMENT_PROTOCOL.md').read_text())+'</main></html>')
    gallery=top+'<h1>零新增物理：缓存分析图</h1><p>所有有限decoder新材料均未full-Maxwell验证；离线GT不参与在线方向、优化或排名。</p>'
    panels=[('01_material_geometry','A：latent与真实材料尺度','材料归一化消除了大量坐标病态，但仍有非零弱方向；prior费用和boundary必须保留。'),('02_decoder_trials','有限decoder运动与latent prior代价','RMS .02/.04，固定三个方向/家族。预测DeltaPhi来自缓存完整材料J，绝不是新full-forward。'),('03_offline_pair','隔离的A–GT材料与有限响应','双noise-ball通过并不等于统计不可分。GT仅离线，图中未展示新的已发现世界。'),('04_B_optimizer','B：原始GN步与统一缩短','建议adaptive TR仍是未验证算法改进，不以缓存模型预测冒充真实收敛。')]
    for name,title,caption in panels:gallery+='<section><h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+html.escape(title)+'"><p>'+caption+'</p><a href="figures/'+name+'.pdf">PDF</a></section>'
    (OUT/'gallery.html').write_text(gallery+'</main></html>')
    notice='## 最新：P1零新增物理审计与下一最小实验\n\n[中文报告]('+REL+'RESEARCH_REPORT_ZH.md) · [网页与冻结协议](https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL+')\n\n新物理0：材料尺度弱方向仍在；A/GT不是重要模式见证；B支持标准trust-region修复。A 12states/288RHS与B 16states/640RHS仅为建议，未授权、未运行。\n\n'
    for name in ['README.md','START_HERE.md']:
        path=ROOT/name;s=path.read_text()
        if REL not in s:path.write_text(notice+s.replace('## 最新：P1材料先验','## 历史：P1材料先验',1))
    path=ROOT/'index.html';s=path.read_text()
    if REL not in s:path.write_text(s.replace('<div class="grid">','<div class="grid"><section class="card"><h2>最新：P1零物理审计</h2><p>真实弱方向、prior代价、GT pair可区分性与B优化。下一实验只冻结、未授权。</p><a href="'+REL+'">结果与协议 →</a></section>',1).replace('最新：P1材料先验','P1材料先验',1))
    path=ROOT/'research/index.html';s=path.read_text()
    if 'p1_zero_physics_audit_20261018/' not in s:path.write_text(s.replace('<main>','<main><section><h2>最新：P1零物理审计</h2><p>新增physics0；A小型探索与B优化修复仅建议。</p><a href="p1_zero_physics_audit_20261018/">结果、图与协议 →</a></section>',1).replace('最新：P1材料先验','P1材料先验',1))
    archive=ROOT/'publication_local/P1_Zero_Physics_Audit_20261018_Evidence.zip';archive.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for path in sorted(OUT.rglob('*')):
            if path.is_file() and path.name not in {'PUBLICATION_RECEIPT.json','WEB_READBACK.json','EVIDENCE_RELEASE.json'}:z.write(path,path.relative_to(ROOT))
        for name in ['p1_zero_online.py','p1_zero_offline.py','publish_p1_zero_audit.py','publish_lxz.py','publish_research_f.py']:z.write(ROOT/'tools'/name,'tools/'+name)
        for path in (ROOT/'src/p1_material').glob('*.py'):z.write(path,path.relative_to(ROOT))
    write('EVIDENCE_RELEASE.json',dict(asset=archive.name,bytes=archive.stat().st_size,source_physical_cache='download previous P1 physics qualification evidence separately; old cost disclosed',new_physics=0,release=RELEASE))
    print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,new_physics=0)))


if __name__=='__main__':main()

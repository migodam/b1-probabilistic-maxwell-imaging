"""Plot and publish saved P1 material outputs; never execute a physical solve."""
from pathlib import Path
import json,html,platform,sys,zipfile,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from publish_lxz import STYLE
from publish_research_f import rendered

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p1_material_target_20261018'
REL='research/p1_material_target_20261018/'
TAG='p1-material-target-20261018'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/'+TAG


def record(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def save(fig,name):
    fig.savefig(OUT/f'figures/{name}.png',dpi=175,bbox_inches='tight')
    fig.savefig(OUT/f'figures/{name}.pdf',bbox_inches='tight')
    plt.close(fig)


def figures(m,r):
    (OUT/'figures').mkdir(exist_ok=True)
    with np.load(OUT/'inputs/offline_materials/calibration.npz') as z:
        x=z['chi'][:6];f=z['fraction'][:6]
    fig,ax=plt.subplots(3,6,figsize=(16,8.5),constrained_layout=True)
    arrays=[1+x.real,x.imag,f];limits=[(1.1,2.9),(.04,.22),(0,1)]
    for row,(values,lim,cmap,label) in enumerate(zip(arrays,limits,['viridis','magma','Greys'],['Re epsilon_r','Im epsilon_r (loss)','Material fraction'])):
        for col in range(6):
            im=ax[row,col].imshow(values[col].reshape(24,24).T,origin='lower',extent=[-.5,.5,-.5,.5],vmin=lim[0],vmax=lim[1],cmap=cmap)
            if row==0:ax[row,col].set_title(f'Prior draw {col:03}')
            if col==0:ax[row,col].set_ylabel(label+'\ny [m]')
            if row==2:ax[row,col].set_xlabel('x [m]')
        fig.colorbar(im,ax=ax[row,:],shrink=.8)
    fig.suptitle('P1 prior material draws — NOT Maxwell reconstructions / NOT a Top-K posterior',fontsize=16)
    save(fig,'01_material_prior')

    events=m['prior_event_mass_estimates'];keys=list(events)
    names={'multiple_components':'Multiple components','one_component_no_closed_host_pocket':'One component / no host pocket','one_component_one_closed_host_pocket':'One closed host pocket','one_component_multiple_closed_host_pockets':'Multiple closed host pockets','empty_at_fraction_threshold':'No >=4-cell threshold component'}
    y=np.array([events[k]['estimated_prior_mass'] for k in keys]);lo=np.array([events[k]['binomial_Wilson95'][0] for k in keys]);hi=np.array([events[k]['binomial_Wilson95'][1] for k in keys])
    fig,ax=plt.subplots(figsize=(10,4.5),constrained_layout=True)
    ax.barh(np.arange(len(keys)),y,xerr=np.r_[[(y-lo)],[(hi-y)]],color='#25667d',capsize=4)
    ax.set_yticks(np.arange(len(keys)),[names[k] for k in keys]);ax.invert_yaxis();ax.set_xlim(0,.75)
    ax.set_xlabel('Prior event frequency, Wilson 95% interval (n=128)')
    ax.set_title('Events are measured after decoding; no topology labels are inputs')
    save(fig,'02_prior_events')

    fig,ax=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
    for method,label,color in [('frozen_latent','Frozen latent ray','#b57334'),('material_native','Material-native','#25667d')]:
        rr=[a for a in r if a['method']==method]
        ax[0].scatter([a['target_RMS'] for a in rr],[a['actual_RMS'] for a in rr],label=label,alpha=.7,color=color)
        ax[1].plot([a['material_goal_error'] for a in rr],'o-',label=label,color=color)
    ax[0].plot([0,.045],[0,.045],'k--',lw=1);ax[0].set(xlabel='Requested material RMS',ylabel='Actual material RMS');ax[0].legend()
    ax[1].set(xlabel='Fixed case × requested distance',ylabel='Material goal error (no likelihood)');ax[1].legend()
    ax[2].bar(range(6),[a['rank'] for a in m['material_tangent_spectra']],color='#25667d');ax[2].set(ylim=(0,70),xlabel='Fixed calculus point',ylabel='Numerical rank of material Dg')
    fig.suptitle('Decoder-only experiment — no physical compatibility or sampler performance claim')
    save(fig,'03_material_motion')


def main():
    start=time.perf_counter()
    if not (OUT/'MATERIAL_STAGE_CLOSED.json').exists():raise RuntimeError('material stage must already be complete')
    m=json.loads((OUT/'METRICS.json').read_text());r=json.loads((OUT/'results/RETRACTION_RESULTS.json').read_text())
    figures(m,r)
    native=[a for a in r if a['method']=='material_native']
    audit=dict(source='saved RETRACTION_RESULTS.json; no optimization rerun',native_optimizer_success=sum(a['optimizer_success'] for a in native),native_optimizer_iteration_limit=sum(not a['optimizer_success'] for a in native),structure_scale_pass_by_method={k:sum(a['structural_change']['task_structure_pass'] for a in r if a['method']==k) for k in ['frozen_latent','material_native']},Maxwell_validated_endpoints=0,posterior_mass=None)
    record('evidence/RETRACTION_STATUS_AUDIT.json',audit)
    record('DECISION.json',dict(target_id=m['target_id'],material_model='READY',material_motion='RADIUS_MATCHED_BUT_OPTIMIZER_CONVERGENCE_LIMITED',actual_observations='NOT_RUN',independent_physics_anchors='NOT_RUN',ambiguity_witness='NOT_RUN',A17_increment='NOT_RUN',new_physics_authorization=0,next_required='explicit new 32 states / 1024 RHS / 600 s authorization',NN_training=0))
    record('DEPENDENCIES.json',dict(python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=__import__('scipy').__version__,matplotlib=matplotlib.__version__,physical_backend='existing A30; not imported in material stage'))
    record('SOURCE_INVENTORY.json',dict(base_B1_commit='263f3f0d6eaf2f9dde3b3cd5adf75f3af9eeb1f2',model_and_protocol_frozen_commit='4f71c0a826f9db6217e5081bd56300716152db24',paid_interface_implementation_commit='a1b613b91e05ec2d09018364f40951f152936636',A30_backend_commit='717d3d26a290248bd2b9aeb55f3d03694b407b84',background_cache='reused already paid F3 weak background, known host only; historical setup cost not remeasured',historical_F3_commit='a323ef223731870b390f875c55f041aa8cc855b4',old_results_unchanged=True,new_Maxwell_states=0,new_complex_RHS=0,source_snapshots='Git commit history and evidence archive',license_status='existing repository license status retained'))
    nav='<nav><a href="../../">B1首页</a><a href="index.html">P1入口</a><a href="report.html">中文报告</a><a href="gallery.html">先验材料图</a><a href="START_HERE.md">GPT入口</a><a href="'+RELEASE+'">原始数组</a></nav>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>P1 · 相关被动材料先验与Maxwell target</title><style>'+STYLE+'figcaption{color:#516672;font-size:15px}</style><main>'+nav
    body='''<h1>P1已经实现：材料支持就绪，真实Maxwell资格尚未运行</h1><p>2026-10-18 · 独立新prior · 24×24 scalar TM · 无NN训练</p>
<section><h2>这次实际完成了什么？</h2><div class="metrics"><span><b>192</b>已保存的先验材料</span><span><b>66</b>连续latent坐标</span><span><b>24</b>材料侧运动端点</span><span><b>0</b>新Maxwell states / RHS</span></div><p>相关场、被动材料参数、C2混合decoder、正规化latent prior与噪声定义已经写成可运行代码。不同结构由解码后统计，generator没有孔洞数量或topology label。</p><p><strong>这些图是先验样本，不是同一测量的多解或后验模式。</strong>在线观测与独立物理anchor仍为NOT_RUN。</p><a href="report.html">完整中文报告 →</a></section>
<section><h2>实际P1材料，不再只有一个双Gaussian chart</h2><img src="figures/01_material_prior.png" alt="固定前六份P1先验样本，同材料色标"><p>按固定顺序展示前六份calibration draw；没有按形状挑图。宿主区为χ=.1+.04i，不是绝对空气；已知外包络，内部形状未知。</p></section>
<section><h2>材料运动能运行，仍有明确失败</h2><img src="figures/03_material_motion.png" alt="实际材料RMS、目标误差与材料导数rank"><p>两种方法各12次，均达到.02/.04材料RMS。native仅2/12正常收敛，10/12触及迭代上限；两方法各1/12满足预设结构变化尺度。没有Maxwell或likelihood优越性结论。</p></section>
<section><h2>下一步已经具体到可执行入口</h2><p>计费接口、固定两个parent的观测生成、隔离GT的独立Full-J反演与立即持久化均已实现。当前物理预算为0，调用在导入物理backend前终止。</p><p>第一轮资格建议：32个states、1024个complex RHS、600秒CPU。只在新的明确授权后执行；历史预算不借用。合格anchor后再比较F3-native / A17-X，保留exact T-matrix、普通Full-J和goal-oriented ROM。</p><a href="PHYSICS_NEXT_RUN.md">固定任务、预算与停止条件 →</a></section>
<section><h2>证据与限制</h2><p>物理求解0、NN/GPU/Mie/新3D均0。成功材料分析wall .976s，首失败attempt wall未测，不冒充总费用。先验事件频率不是后验概率，所有材料样本现已暴露。</p><p><a href="METRICS.json">实际指标</a> · <a href="RESOURCE_LEDGER.jsonl">资源与失败</a> · <a href="MODEL_AND_PRIOR.md">模型与测度</a> · <a href="REPRODUCE.md">复现入口</a></p></section>'''
    (OUT/'index.html').write_text(top+body+'</main></html>')
    (OUT/'report.html').write_text(top+'<article>'+rendered((OUT/'RESEARCH_REPORT_ZH.md').read_text())+'</article></main></html>')
    gallery=top+'<h1>P1先验图库与材料运动</h1><p>仅材料侧构造。没有散射观测、Maxwell重建、Top-K后验或新物理多解。</p>'
    for name,title,caption in [('01_material_prior','固定前六份先验材料','Re epsilon_r、损耗和混合分数共用色标；x/y按原mesh顺序显示。'),('02_prior_events','先验推前事件频率','128个calibration样本，Wilson95%区间。阈值.5、foreground4连通、host8连通、最小4cells。不是后验质量，也不为展示图赋单点概率。'),('03_material_motion','实际材料距离与有限运动','native只有2/12优化正常收敛，其余迭代上限；不把Dg rank下降当作Maxwell nullspace。')]:
        gallery+='<section><h2>'+title+'</h2><img src="figures/'+name+'.png" alt="'+html.escape(title)+'"><p>'+caption+'</p><a href="figures/'+name+'.pdf">PDF</a></section>'
    (OUT/'gallery.html').write_text(gallery+'</main></html>')
    notice='## 最新：P1材料先验与Maxwell target实现\n\n[中文报告]('+REL+'RESEARCH_REPORT_ZH.md) · [网页与先验材料](https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL+')\n\n192份先验材料、24个材料侧运动端点已保存；真实P1观测和物理anchor尚未运行。新增Maxwell states/RHS均0，不训练NN。\n\n'
    for name in ['README.md','START_HERE.md']:
        path=ROOT/name;s=path.read_text()
        if REL not in s:path.write_text(notice+s.replace('## 最新：A17-X','## 历史：A17-X',1))
    path=ROOT/'index.html';s=path.read_text()
    if REL not in s:path.write_text(s.replace('<div class="grid">','<div class="grid"><section class="card"><h2>最新：P1材料先验实现</h2><p>192份实际先验材料，0新物理；观测与anchor待独立授权。</p><a href="'+REL+'">报告与先验材料 →</a></section>',1).replace('最新：A17-X','A17-X',1))
    path=ROOT/'research/index.html';s=path.read_text()
    if 'p1_material_target_20261018/' not in s:path.write_text(s.replace('<main>','<main><section><h2>最新：P1材料先验与target</h2><p>模型、实际材料数组和付费接口已实现；物理资格尚未运行。</p><a href="p1_material_target_20261018/">报告与材料 →</a></section>',1).replace('最新：A17-X','A17-X',1))
    record('PUBLICATION.json',dict(repository='migodam/b1-probabilistic-maxwell-imaging',new_repository=False,release=RELEASE,old_ledgers_unchanged=True,locked_data=False,private_conversations=False,third_party_fulltext=False,hash_checks=False,scientific_status_independent_of_publication=True))
    archive=ROOT/'publication_local/P1_Material_Target_20261018_Evidence.zip';archive.parent.mkdir(exist_ok=True)
    excluded={'physics.lock','PUBLICATION_RECEIPT.json','WEB_READBACK.json','EVIDENCE_RELEASE.json'}
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file() and p.name not in excluded and '__pycache__' not in p.parts:z.write(p,p.relative_to(ROOT))
        for folder in ['src/p1_material']:
            for p in (ROOT/folder).glob('*.py'):z.write(p,p.relative_to(ROOT))
        for name in ['prepare_p1_material_target.py','run_p1_physics.py','publish_p1_material.py','publish_research_f.py','publish_lxz.py']:
            z.write(ROOT/'tools'/name,'tools/'+name)
    with zipfile.ZipFile(archive) as z:names=z.namelist()
    record('EVIDENCE_RELEASE.json',dict(asset=archive.name,bytes=archive.stat().st_size,files=len(names),material_arrays=sum(n.endswith('.npz') for n in names),physical_endpoints=0,background_cache_reused=True))
    print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,reporting_wall_seconds=time.perf_counter()-start,native_status=audit),ensure_ascii=False))


if __name__=='__main__':main()

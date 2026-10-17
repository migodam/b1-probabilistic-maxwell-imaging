"""Plot saved audit results, prepare existing-site pages and evidence. No physics."""
from pathlib import Path
import sys,json,zipfile,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from publish_lxz import STYLE
from publish_research_f import rendered
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/a17x_f3_competition_audit'
REL='research/a17x_f3_competition_audit/'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/a17x-f3-audit-20261017'

def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
def main():
    rows=json.loads((OUT/'results/CACHED_PORT_AND_SOURCE_AUDIT.json').read_text())['rows'];eligible=[r for r in rows if r['eligible']]
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for r in eligible:
        points=[z for z in r['compression'] if z['ports']<576]
        axes[0].loglog([z['ports'] for z in points],[z['relative_operator_error'] for z in points],marker='o',label=r['method'])
        axes[1].loglog([z['ports'] for z in points],[z['relative_cached_linear_receiver_error'] for z in points],marker='o',label=r['method'])
    for ax,title in zip(axes,['Operator spectral truncation','Cached first-order receiver truncation']):
        ax.axvline(7,color='gray',ls=':');ax.axhline(.01,color='black',ls='--');ax.set(xlabel='Retained physical ports',ylabel='Relative error',title=title);ax.legend(fontsize=8)
    fig.suptitle('Paid F3 endpoints: seven latent coordinates do not imply seven ports. Not a finite-response certificate.')
    for ext in ['png','pdf']:fig.savefig(OUT/f'figures/01_port_compression.{ext}',dpi=175,bbox_inches='tight')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for r in eligible:
        pp=r['current_straight_interpolation'];x=[z['t'] for z in pp]
        axes[0].plot(x,[z['multi_source_row_fit_relative_residual'] for z in pp],marker='o',label=r['method'])
        axes[1].plot(x,[z['fixed_phase_relative_defect'] for z in pp],marker='o',label=r['method'])
    for ax,title in zip(axes,['Shared-material row-fit defect','Fixed-phase material defect']):
        ax.set(xlabel='Straight-current interpolation t',ylabel='Relative defect',title=title);ax.legend(fontsize=8);ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    fig.suptitle('Cached-current interpolation only: no new Maxwell state or candidate validation.')
    for ext in ['png','pdf']:fig.savefig(OUT/f'figures/02_shared_material.{ext}',dpi=175,bbox_inches='tight')
    plt.close(fig)
    nav='<nav><a href="../../">B1</a><a href="index.html">本轮入口</a><a href="report.html">中文报告</a><a href="gallery.html">数值图</a><a href="START_HERE.md">GPT入口</a><a href="'+RELEASE+'">原始证据</a></nav>'
    math='<script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>A17-X × F3 · 端口与共享材料审计</title><style>'+STYLE+'.math{overflow-x:auto;padding:8px}</style>'+math+'<main>'+nav
    body='''<h1>A17-X × F3：先验证计算机制，不预设增量</h1><p>2026-10-17 · 零新增物理审计已完成 · 真实P1竞争尚未启动</p>
<section><h2>七维先验，不是七端口</h2><div class="metrics"><span><b>0</b>新增Maxwell states / RHS</span><span><b>17</b>已付费F3端点审计</span><span><b>575</b>每个存储DeltaX的非零对角项</span><span><b>NOT_RUN</b>真实P1与NN</span></div><p>三个合格weak端点若只保留最大7个端口，算子相对误差仍58%–64%，缓存一阶接收误差103%–115%。没有finite输出证书，不能由低维latent推定低秩材料更新。</p><p><strong>当前裁决：条件代数GO；新物理实验NO-GO，等待真实P1 target及独立预算。A17增量未建立。</strong></p><a href="report.html">完整中文报告 →</a></section>
<section><h2>实际缓存支持了什么？</h2><img src="figures/01_port_compression.png" alt="575秩增量的端口截断误差"><p>全部是历史已付费数组的代数，P0未重跑。图中receiver曲线仅为一阶，不代表截断后完整Maxwell预测。</p></section>
<section><h2>共享材料不会因current插值自动成立</h2><img src="figures/02_shared_material.png" alt="多源共享和固定相位缺陷"><p>三个候选的直线current路径不满足精确共享标量与固定相位条件。这里只作可实现性诊断，没有新材料候选或posterior样本。</p></section>
<section><h2>依据结果，实际修订了一次机制</h2><p>改为固定合法材料，直接按有限状态的source RHS更新current空间，避免预建575列端口传播。12个dense合成增量核验通过；同空间A17与普通Galerkin仍相同，差≤4.62e−16，没有真实速度或发现增量主张。</p><a href="MECHANISM_REVISION.md">修订、对照和数值 →</a></section>
<section><h2>下一轮必须公平竞争</h2><p>proposal轴：F3-native / X-rational / ordinary Full-J。固定材料求值轴：full / exact T-matrix / ordinary goal-ROM / A17。两个独立P1 parents、相同材料距离和全部费用；无目标/预算则不执行。</p><p><a href="COMPETITION_PROTOCOL_FROZEN.md">冻结设计</a> · <a href="EVIDENCE_GAPS_AND_BUDGET.md">缺失证据与条件预算</a> · <a href="GO_NO_GO.md">Go/No-Go</a> · <a href="SOURCE_AND_TARGET_MANIFEST.md">来源与目标</a> · <a href="REPRODUCE.md">复现</a></p></section>'''
    (OUT/'index.html').write_text(top+body+'</main></html>')
    (OUT/'report.html').write_text(top+'<article>'+rendered((OUT/'RESEARCH_REPORT_ZH.md').read_text())+'</article></main></html>')
    (OUT/'gallery.html').write_text(top+'<h1>零新增物理的数值图</h1><p>已付费F3 cache诊断。没有新P1成像，候选概率不估计。</p><section><h2>端口截断</h2><img src="figures/01_port_compression.png"><p>最优算子截断不等于最优goal-oriented current压缩；不排除普通ROM的潜力。</p><a href="figures/01_port_compression.pdf">PDF</a></section><section><h2>直线current路径</h2><img src="figures/02_shared_material.png"><p>不能成为精确共享材料，不代表其他曲线不存在。没有进行新forward。</p><a href="figures/02_shared_material.pdf">PDF</a></section></main></html>')
    notice='## 最新：A17-X × F3 竞争审计\n\n[中文报告]('+REL+'RESEARCH_REPORT_ZH.md) · [网页与数值图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL+')\n\n零新增物理：17份旧F3增量均rank575；小端口假设不能由七维latent推出。合成机制修订完成，真实P1及新预算未就绪，比较NOT_RUN。\n\n'
    for name in ['README.md','START_HERE.md']:
        p=ROOT/name;s=p.read_text()
        if '## 最新：A17-X' not in s:p.write_text(notice+s.replace('## 最新：Research F3','## 历史：Research F3'))
    p=ROOT/'index.html';s=p.read_text()
    if REL not in s:
        card='<section class="card"><h2>最新：A17-X × F3竞争审计</h2><p>零新物理；低维latent不等于低秩端口。真实P1未就绪，A17不预设增量。</p><a href="'+REL+'">报告、缓存诊断与冻结协议 →</a></section>'
        p.write_text(s.replace('<div class="grid">','<div class="grid">'+card,1).replace('最新：Research F3','Research F3'))
    p=ROOT/'research/index.html';s=p.read_text()
    if 'a17x_f3_competition_audit/' not in s:p.write_text(s.replace('<main>','<main><section><h2>最新：A17-X × F3审计</h2><p>0 states / 0 RHS；575秩端口障碍与共享材料诊断。</p><a href="a17x_f3_competition_audit/">报告与条件协议 →</a></section>',1).replace('最新：Research F3','Research F3'))
    ignore='\n# A17-X audit: reused paid arrays and synthetic arrays in the evidence Release.\nresearch/a17x_f3_competition_audit/evidence/reused_f3/*.npz\nresearch/a17x_f3_competition_audit/results/*.npz\n'
    p=ROOT/'.gitignore'
    if 'a17x_f3_competition_audit/evidence/reused_f3' not in p.read_text():p.write_text(p.read_text()+ignore)
    dump(OUT/'DEPENDENCIES.json',dict(python=sys.version,numpy=np.__version__,scipy=__import__('scipy').__version__,matplotlib=matplotlib.__version__,platform=platform.platform(),physical_backend_imported=False,GPU_training=0))
    dump(OUT/'PUBLICATION.json',dict(repository='migodam/b1-probabilistic-maxwell-imaging',new_repository=False,release=RELEASE,license_status_unchanged=True,private_conversations=False,GT_arrays=False,locked_data=False,third_party_fulltexts=False,hash_checks=False,publication_independent_of_scientific_success=True))
    archive=ROOT/'publication_local/A17X_F3_Competition_Audit_20261017_Evidence.zip';archive.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file() and p.name not in ['EVIDENCE_RELEASE.json','PUBLICATION_RECEIPT.json','WEB_READBACK.json'] and '__pycache__' not in p.parts:z.write(p,p.relative_to(ROOT))
        for p in (ROOT/'src/a17x').glob('*.py'):z.write(p,p.relative_to(ROOT))
        for name in ['a17x_zero_audit.py','a17x_source_rhs_revision.py','publish_a17x_audit.py','publish_research_f.py','publish_lxz.py']:z.write(ROOT/'tools'/name,'tools/'+name)
        for name in ['results/RECORDS.json','results/ELIGIBLE_ENDPOINTS.json','results/F15/p250073_weak_a1_snr20_FINAL.json','results/F15/p250073_ring_a1_snr20_FINAL.json']:
            p=ROOT/'research/research_f3_quotient_exploration'/name;z.write(p,p.relative_to(ROOT))
        # Rendering helper imports lxz.proxy; these source files have no physics side effects.
        for p in (ROOT/'src/lxz').glob('*.py'):z.write(p,p.relative_to(ROOT))
    with zipfile.ZipFile(archive) as z:names=z.namelist()
    rec=dict(asset=archive.name,bytes=archive.stat().st_size,files=len(names),reused_F3_array_slices=sum('/reused_f3/' in n and n.endswith('.npz') for n in names),new_Maxwell_endpoints=0,new_physics=0,GT=False,LU_or_pivots=False,release=RELEASE,hash_checks=False)
    dump(OUT/'EVIDENCE_RELEASE.json',rec);print(json.dumps(rec))

if __name__=='__main__':main()

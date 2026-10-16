"""Publish Research F artifacts in the existing B1 site. No physical solves."""
from pathlib import Path
import re,json,html,zipfile,shutil,sys,platform
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from publish_lxz import render,STYLE
OUT=ROOT/'research/research_f_intrinsic_transport'
REL='research/research_f_intrinsic_transport/'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-f-intrinsic-20261017'

def rendered(md):
    blocks=[];imgs=[]
    def mathblock(m):blocks.append(m.group(1));return 'F_MATH_BLOCK_'+str(len(blocks)-1)
    md=re.sub(r'\\\[\s*\n(.*?)\n\\\]',mathblock,md,flags=re.S)
    def picture(m):imgs.append((m.group(1),m.group(2)));return 'F_IMAGE_'+str(len(imgs)-1)
    md=re.sub(r'!\[([^\]]*)\]\(([^)]+)\)',picture,md)
    result=render(md)
    for i,formula in enumerate(blocks):result=result.replace('<p>F_MATH_BLOCK_'+str(i)+'</p>','<div class="math">\\['+html.escape(formula)+'\\]</div>')
    for i,(alt,path) in enumerate(imgs):result=result.replace('F_IMAGE_'+str(i),'<img alt="'+html.escape(alt)+'" src="'+html.escape(path)+'">')
    return result

def main():
    for p in [OUT/'RESEARCH_REPORT_ZH.md']:
        p.write_text(p.read_text().replace(r'\(||Q_{\rm out}C_{\rm prior}||\)', 'prior余核响应范数'))
    nav='<nav><a href="../../">B1</a><a href="index.html">本轮入口</a><a href="report.html">中文报告</a><a href="gallery.html">材料图集</a><a href="START_HERE.md">GPT入口</a><a href="'+RELEASE+'">完整数组</a></nav>'
    math='<script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research F · Maxwell材料几何与物理反演</title><style>'+STYLE+'.math{overflow-x:auto;padding:8px}figcaption{color:#516672;font-size:15px}</style>'+math+'<main>'+nav
    body='''<h1>Research F · Maxwell材料几何与独立物理反演</h1><p>2026-10-17 · 已执行并封账 · 原P0 / 2D scalar-TM</p><section><h2>坐标问题已分离，有限探索仍未建立</h2><p>巨大decoder二阶项中，weak例99.9967%的材料二阶能量来自坐标加速度；ring例85.2285%。真实prior曲率仍非零。完整历史Maxwell二阶向量缺失，保留MISSING。</p><div class="metrics"><span><b>24 / 64</b>新full-material states</span><span><b>960 / 2048</b>complex RHS</span><span><b>5.534 s</b>受控物理wall</span><span><b>NOT_RUN</b>F2有限传输</span></div><p>独立physics-only反演取得两份合法且噪声相容的材料，但未过冻结驻点标准，停止F2，不调整门槛、不训练NN。</p><a href="report.html">完整中文结论与十问 →</a></section><section><h2>原先验内的物理反演，真正运行了</h2><img src="figures/01_material_ring.png" alt="先验初值、独立物理反演、缓存N3及离线真值"><p>同parent的两个已暴露背景条件。物理反演额外已知family/global geometry，不能称与N3同信息量的普遍胜出。GT只用于离线评分。</p></section><section><h2>为什么不继续堆高阶网络？</h2><img src="figures/03_geometry_and_offline_risk.png" alt="坐标与真实prior曲率及离线风险"><p>真实材料曲率、坐标加速度和Maxwell反馈须分开。两例去gauge的响应均rank7，未确认非零exact local fiber；有限噪声区域仍未排除。</p></section><section><h2>继续什么，停止什么</h2><p><strong>CONDITIONAL GO：</strong>合法材料商坐标、真实比例边界和有界物理驻点。<strong>NO-GO：</strong>现在训练learned Explorer或声称发现新posterior modes。</p><p><a href="GO_NO_GO_DECISION.md">逐项裁决</a> · <a href="evidence/INDEPENDENT_MATH_REVIEW.md">独立数学审查</a> · <a href="RESOURCE_SUMMARY.json">完整费用</a> · <a href="REPRODUCE.md">复现边界</a></p></section>'''
    (OUT/'index.html').write_text(top+body+'</main></html>')
    (OUT/'report.html').write_text(top+'<article>'+rendered((OUT/'RESEARCH_REPORT_ZH.md').read_text())+'</article></main></html>')
    gallery=top+'<h1>Research F · 实际材料与数值图</h1><p>F0旧TRAIN数学点 / F1新条件反演分别展示。F2未运行，无传输、posterior样本或新模式图。</p>'
    panels=[('01_material_ring','ring：初值、独立物理解、缓存N3与离线真值','已知家族/全局几何，所有材料使用相同色标。物理残差通过，驻点资格未过。'),('01_material_weak','weak：同一个parent，另一背景','不能将两个背景算作两独立物体；N3 full-forward residual缺失。'),('02_F1_convergence_and_gate','实际12次forward序列与冻结驻点标准','ftol不是KKT；没有虚构中间梯度历史。'),('03_geometry_and_offline_risk','可消除坐标项与真实prior法向项','阴影柱只有历史Maxwell norm，不是缺失向量；区间使用三角不等式，依赖历史norm。')]
    for name,title,caption in panels:gallery+='<section><h2>'+title+'</h2><figure><img src="figures/'+name+'.png" alt="'+title+'"><figcaption>'+caption+'</figcaption></figure><a href="figures/'+name+'.pdf">PDF</a></section>'
    (OUT/'gallery.html').write_text(gallery+'</main></html>')
    notice='## 最新：Research F 材料几何与独立物理反演\n\n[中文报告]('+REL+'RESEARCH_REPORT_ZH.md) · [网页与材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL+')\n\n24 states / 960 RHS封账：分离坐标与真实prior曲率；独立物理反演合法且噪声相容，但驻点资格未过，F2未运行。不训练NN。\n\n'
    for n in ['README.md','START_HERE.md']:
        p=ROOT/n;s=p.read_text()
        if '## 最新：Research F' not in s:p.write_text(notice+s.replace('## 最新诊断：有限Maxwell反馈','## 历史诊断：有限Maxwell反馈'))
    p=ROOT/'index.html';s=p.read_text()
    if REL not in s:
        card='<section class="card"><h2>最新：Research F材料几何</h2><p>24 states：真实先验曲率与坐标项分离；独立反演噪声相容，驻点未达标，F2未运行。</p><a href="'+REL+'">中文报告与材料图 →</a></section>'
        p.write_text(s.replace('<div class="grid">','<div class="grid">'+card,1).replace('最新：有限Maxwell反馈诊断','有限Maxwell反馈诊断'))
    p=ROOT/'research/index.html';s=p.read_text()
    if 'research_f_intrinsic_transport/' not in s:p.write_text(s.replace('<main>','<main><section><h2>最新：Research F材料几何与独立物理反演</h2><p>原P0、24 states / 960 RHS。传输未启动，保留完整负结果。</p><a href="research_f_intrinsic_transport/">报告、原始材料与裁决 →</a></section>',1))
    dep=dict(python=sys.version,platform=platform.platform(),numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,matplotlib=__import__('matplotlib').__version__)
    (OUT/'DEPENDENCIES.json').write_text(json.dumps(dep,indent=2)+'\n')
    old=ROOT.parent/'b1-finite-feedback/research/finite_feedback_fiber_diagnostic_20261016'
    cache=OUT/'evidence/reused_f0'
    for key in ['250033_weak','250034_ring']:
        for source in [old/f'results/{key}_CHART.npz',old/f'results/{key}_JACOBIAN.npz',old/f'inputs/{key}.npz',old/f'inputs/{key}.json']:
            dest=cache/source.relative_to(old);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dest)
    shutil.copy2(old/'METRICS.json',cache/'METRICS.json')
    from prior_constrained import prior
    for bg in ['weak','ring']:
        dest=OUT/f'inputs/background/{bg}.npz';dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(prior.INPUTS/f'background/{bg}.npz',dest)
    pub=dict(repository='migodam/b1-probabilistic-maxwell-imaging',new_repository=False,release=RELEASE,license_status_unchanged=True,third_party_fulltexts=False,private_conversations=False,locked_data=False,hash_checks=False,scientific_success_independent_of_publication=True,original_old_ledgers_unchanged=True)
    (OUT/'PUBLICATION.json').write_text(json.dumps(pub,indent=2)+'\n')
    archive=ROOT/'publication_local/Research_F_Intrinsic_Transport_20261017_Evidence.zip';archive.parent.mkdir(exist_ok=True)
    excluded={'budget.lock','PRIVATE_CONVERSATION.json','EVIDENCE_RELEASE.json','PUBLICATION_RECEIPT.json','WEB_READBACK.json'}
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for source in sorted(OUT.rglob('*')):
            if source.is_file() and source.name not in excluded and '__pycache__' not in source.parts:z.write(source,source.relative_to(ROOT))
        for folder in ['src/research_f','src/prior_constrained','src/finite_feedback','src/lxz']:
            for source in (ROOT/folder).glob('*.py'):z.write(source,source.relative_to(ROOT))
        for name in ['research_f_geometry.py','research_f_figures.py','research_f_stationarity.py','publish_research_f.py']:z.write(ROOT/'tools'/name,'tools/'+name)
    with zipfile.ZipFile(archive) as z:
        names=z.namelist();ns=sum(bool(re.search(r'/F1/.*/state_\d\d\.npz$',n)) for n in names);nj=sum(n.endswith('_JAC.npz') and '/F1/' in n for n in names)
    rec=dict(asset=archive.name,bytes=archive.stat().st_size,files=len(names),full_state_endpoint_arrays=ns,Jacobian_arrays=nj,F0_cache_included=True,paid_background_arrays=True,private_conversation_excluded=True,hash_checks=False,release=RELEASE)
    (OUT/'EVIDENCE_RELEASE.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec))
if __name__=='__main__':main()

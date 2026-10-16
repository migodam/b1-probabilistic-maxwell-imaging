"""Prepare existing B1 site and evidence, without physics or hashes."""
from pathlib import Path
import json,sys,html,re,zipfile,subprocess,platform
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from lxz.proxy import OUT,dump

def inline(s):
    s=html.escape(s)
    s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s)
    return re.sub(r'`([^`]+)`',r'<code>\1</code>',s)

def render(md):
    lines=md.splitlines();blocks=[];i=0
    while i<len(lines):
        s=lines[i]
        if not s.strip():i+=1;continue
        if s.startswith('|'):
            table=[]
            while i<len(lines) and lines[i].startswith('|'):
                a=lines[i].strip('|').split('|')
                if not all(re.match(r'^\s*:?-+:?\s*$',t) for t in a):table.append(a)
                i+=1
            blocks.append('<div class="scroll"><table>'+''.join('<tr>'+''.join(f'<{"th" if n==0 else "td"}>{inline(c.strip())}</{"th" if n==0 else "td"}>' for c in row)+'</tr>' for n,row in enumerate(table))+'</table></div>');continue
        if s.startswith('#'):
            n=len(s)-len(s.lstrip('#'));blocks.append(f'<h{n}>'+inline(s[n:].strip())+f'</h{n}>');i+=1;continue
        blocks.append('<p>'+inline(s)+'</p>');i+=1
    return '\n'.join(blocks)

STYLE='body{font:17px/1.8 system-ui,-apple-system,"PingFang SC",sans-serif;color:#20323e;background:#f6f8fa;margin:0}main{max-width:1100px;margin:auto;padding:42px 24px}h1,h2{line-height:1.4}h2{margin-top:2em}a{color:#126c93}p{overflow-wrap:anywhere}nav{display:flex;gap:18px;flex-wrap:wrap}article,section{background:#fff;border:1px solid #dce4e8;padding:24px;border-radius:12px;margin:22px 0}img{max-width:100%;height:auto}.scroll{overflow-x:auto}table{border-collapse:collapse;font-size:15px;width:100%}td,th{border:1px solid #d3dce1;padding:8px;text-align:left}th{background:#eaf2f5}.metrics{display:flex;gap:24px;flex-wrap:wrap}.metrics b{font-size:28px;display:block}code{background:#edf2f5;padding:2px 4px}'

def main():
    release='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-lxz-20261016'
    nav='<nav><a href="../../">B1首页</a><a href="report.html">中文报告</a><a href="gallery.html">材料图库</a><a href="'+release+'">原始证据下载</a><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_lxz_20261016">源码与审计</a></nav>'
    header='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>LXZ · Maxwell局部–全局多解研究</title><style>'+STYLE+'</style><main>'+nav
    body='''<h1>LXZ · Maxwell 局部–全局多解研究</h1><p>2026-10-16 · 已运行、预算封账。N3强Top-1基线和历史T/G结果保持。</p><section><h2>先看实际结论</h2><div class="metrics"><span><b>0</b>新增显著结构world</span><span><b>0.00436</b>global MH proxy事件MAE</span><span><b>96 / 96</b>二维full states</span><span><b>896 / 3072</b>complex RHS</span></div><p>已有C1提案经正确MH，在本次径向proxy上优于local–global混合。复杂二维只得到公开TRAIN已知材料附近的小幅连续变化；N3附近的合法fiber与分离mode仍未证实。</p><p><a href="report.html">阅读完整中文报告 →</a> · <a href="gallery.html">看材料图和比较图 →</a></p></section><section><h2>概率机制</h2><img src="figures/e1_proxy_event_cdf.png" alt="同调用预算的概率误差"><p>四个暴露观测、四个seed；已排除两次错误geometry结果。它是proxy target机制实验，不是新full-Mie posterior。</p></section><section><h2>真实二维局部材料</h2><img src="figures/e2_materials_250034_ring.png" alt="公开TRAIN校准材料"><p>相同材料色标。小幅材料位移不等于显著结构新解；候选概率为null。首失败任务35个状态原始endpoint数组缺失已披露，费用全部保留。</p></section><section><h2>继续与停止</h2><p>保留normalized proposal+正确MH基准；停止当前复杂混合、A17、Flow优越性主张。不启动E3：没有独立分离mode证据。现在不足以支持强TAP算法论文。</p><p><a href="GO_NO_GO_DECISION.md">机器可读裁决与依据</a> · <a href="TARGET_REGISTRY.json">目标身份</a> · <a href="RESOURCE_SUMMARY.json">完整费用</a> · <a href="START_HERE.md">GPT阅读入口</a></p></section>'''
    (OUT/'index.html').write_text(header+body+'</main></html>')
    (OUT/'report.html').write_text(header+'<article>'+render((OUT/'RESEARCH_REPORT_ZH.md').read_text())+'</article></main></html>')
    gallery=(OUT/'GALLERY.html').read_text();(OUT/'gallery.html').write_text(gallery.replace('<h1>',nav+'<h1>',1))
    notice='## 最新研究：LXZ 局部–全局多解验证\n\n[中文报告](research/research_lxz_20261016/RESEARCH_REPORT_ZH.md) · [网页与材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_lxz_20261016/)\n\n96个二维states封账，未发现显著结构新world；已有C1提案经正确MH为本次proxy基线最佳。N3合法anchor与分离mode仍未证实。\n\n'
    for name in ['README.md','START_HERE.md']:
        p=ROOT/name;s=p.read_text()
        if '## 最新研究：LXZ' not in s:p.write_text(notice+s)
    p=ROOT/'index.html';s=p.read_text()
    if 'research/research_lxz_20261016/' not in s:
        s=s.replace('<div class="grid">','<div class="grid"><section class="card"><h2>最新：LXZ局部–全局多解验证</h2><p>96个二维状态，0显著新结构；C1提案+正确MH胜过当前混合。保留N3，公开失败与完整费用。</p><a href="research/research_lxz_20261016/">看报告、概率误差与材料图 →</a></section>',1)
        s=s.replace('最新：原先验约束竞争解释','历史：原先验约束竞争解释');p.write_text(s)
    p=ROOT/'research/index.html';s=p.read_text()
    if 'research_lxz_20261016/' not in s:p.write_text(s.replace('<main>','<main><section><h2>最新：LXZ 局部–全局Maxwell研究</h2><p>proxy机制与真实二维校准；未建立N3 fiber或分离新模式。</p><a href="research_lxz_20261016/">报告与实际图库 →</a></section>',1))
    ignore='\n# LXZ raw chains and physical endpoints are retained in its evidence Release.\nresearch/research_lxz_20261016/results/**/*.npz\nresearch/research_lxz_20261016/budget.lock\n'
    p=ROOT/'.gitignore'
    if 'research/research_lxz_20261016/results/**/*.npz' not in p.read_text():p.write_text(p.read_text()+ignore)
    dependencies=dict(python=sys.version.split()[0],platform=platform.platform(),
        numpy=__import__('numpy').__version__,scipy=__import__('scipy').__version__,matplotlib=__import__('matplotlib').__version__,
        GPU_training=0,C1_weights='https://github.com/migodam/c1-physics-state-generative-imaging/releases/tag/joint-p2-20261014',
        B1_paid_bank='research/joint_p2/results/COMMON_RESPONSE_BANK.npz',
        A30_backend_commit='717d3d26a290248bd2b9aeb55f3d03694b407b84',A26_prior_commit='2603f3bd16c48e48587fedc4363975779def17f5')
    dump(OUT/'DEPENDENCIES.json',dependencies)
    oldpub=dict(repository='migodam/b1-probabilistic-maxwell-imaging',new_repository=False,license='existing license status unchanged',
        publication_independent_of_scientific_success=True,release=release,hash_checks=False,
        data='authorized synthetic public TRAIN calibration and exposed P2 proxy, no locked data',thirdparty='no copied literature fulltexts',
        raw_first35_failed_endpoint_gap_preserved=True,other_repositories_unmodified=True)
    dump(OUT/'PUBLICATION.json',oldpub)
    archive=ROOT/'publication_local/LXZ_20261016_Evidence.zip';archive.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(OUT.rglob('*')):
            if p.is_file() and p.name not in {'budget.lock','EVIDENCE_RELEASE.json','PUBLICATION_RECEIPT.json'} and '__pycache__' not in p.parts:z.write(p,p.relative_to(ROOT))
        for p in sorted((ROOT/'src/lxz').glob('*.py')):z.write(p,p.relative_to(ROOT))
        for name in ['report_lxz.py','plot_lxz.py','package_lxz_proxy.py','publish_lxz.py']:z.write(ROOT/'tools'/name,'tools/'+name)
    dump(OUT/'EVIDENCE_RELEASE.json',dict(asset=archive.name,bytes=archive.stat().st_size,files=len(zipfile.ZipFile(archive).infolist()),
        release=release,hashes_not_requested_not_computed=True,physics=0,relative_root='existing B1 repo'))
    print(json.dumps(dict(asset=str(archive),bytes=archive.stat().st_size)))

if __name__=='__main__':main()

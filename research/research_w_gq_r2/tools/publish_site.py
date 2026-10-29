"""Present sealed R2 results in the existing B1 website; no physics calls."""
from pathlib import Path
import html,json,re

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
REL='research/research_w_gq_r2/'
GITHUB='https://github.com/migodam/b1-probabilistic-maxwell-imaging/'
BLOB=GITHUB+'blob/main/'+REL
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/'+REL
PAGE='https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL
TAG='research-w-gq-r2-20261029'
RELEASE=GITHUB+'releases/tag/'+TAG
ASSET=GITHUB+'releases/download/'+TAG+'/Research_W_GQ_R2_Results.zip'

STYLE='''*{box-sizing:border-box}body{margin:0;background:#f4f6f8;color:#172633;font:17px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif}main{max-width:1120px;margin:auto;padding:32px 24px 70px}nav{display:flex;gap:10px 22px;flex-wrap:wrap;margin-bottom:32px}a{color:#17668d;text-underline-offset:3px}h1{font-size:clamp(29px,4.5vw,44px);line-height:1.35;letter-spacing:-.03em}h2{font-size:26px;line-height:1.5;margin-top:2.1em}h3{line-height:1.5}article,.panel{background:#fff;border:1px solid #dce3e8;border-radius:14px;padding:28px;margin:22px 0}.eyebrow{font-size:14px;color:#526575;letter-spacing:.04em}.lede{font-size:21px;line-height:1.8}.status{padding:17px 21px;border-left:4px solid #9a7035;background:#fff4df;border-radius:4px;margin:24px 0}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:26px 0}.metric{background:#eaf1f5;border:1px solid #dce5eb;border-radius:12px;padding:18px}.metric strong{display:block;font-size:30px;line-height:1.3;color:#164d68}.metric span{font-size:14px;color:#405766}.muted,figcaption{font-size:14px;color:#506374}figure{margin:30px 0}img{max-width:100%;height:auto;display:block;background:white}figure a{font-size:14px}.table-scroll{overflow-x:auto;margin:20px 0}table{border-collapse:collapse;width:100%;font-size:15px}td,th{border:1px solid #d5dfe6;padding:9px 12px;text-align:left;vertical-align:top}th{background:#edf3f7}code,pre{background:#edf2f5;font-size:.88em}code{padding:2px 4px;border-radius:4px}pre{padding:16px;overflow:auto;border-radius:8px}p,li{overflow-wrap:anywhere}.math{overflow-x:auto;margin:22px 0}footer{border-top:1px solid #d5dfe6;margin-top:38px;padding-top:20px;font-size:14px;color:#526575}blockquote{margin:20px 0;padding:2px 20px;border-left:3px solid #9fb4c1}.two{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:18px}@media(max-width:600px){main{padding:22px 14px 50px}article,.panel{padding:20px}nav{font-size:15px}}'''
MATH=r'''<script>window.MathJax={tex:{inlineMath:[['\\(','\\)']],displayMath:[['\\[','\\]']]},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>'''

def href(s):
    if s.startswith(('http:','https:','#','mailto:')):return s
    if s=='REPORT_ZH.md':return 'report.html'
    if s=='NEXT_PRO_HANDOFF.md':return 'handoff.html'
    if s=='START_HERE.md':return 'index.html'
    f=ROOT/s.split('#')[0]
    if f.exists():
        return s if f.suffix in {'.png','.pdf','.html'} else BLOB+s
    return 'evidence.html'

def inline(s):
    tokens=[]
    def keep(v):
        tokens.append(v);return '\x00'+str(len(tokens)-1)+'\x00'
    s=re.sub(r'`([^`]+)`',lambda m:keep('<code>'+html.escape(m[1])+'</code>'),s)
    s=re.sub(r'!?\[([^\]]+)\]\(([^)]+)\)',lambda m:keep('<a href="'+html.escape(href(m[2]),quote=True)+'">'+html.escape(m[1])+'</a>'),s)
    s=html.escape(s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s)
    return re.sub('\x00(\d+)\x00',lambda m:tokens[int(m[1])],s)

def render(md):
    lines=md.splitlines();out=[];i=0
    while i<len(lines):
        s=lines[i].strip()
        if not s:i+=1;continue
        if s=='---':out.append('<hr>');i+=1;continue
        if s.startswith('```'):
            code=[];i+=1
            while i<len(lines) and not lines[i].strip().startswith('```'):code.append(lines[i]);i+=1
            out.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>');i+=1;continue
        if s==r'\[':
            math=[s];i+=1
            while i<len(lines):
                math.append(lines[i]);i+=1
                if math[-1].strip()==r'\]':break
            out.append('<div class="math">'+html.escape('\n'.join(math))+'</div>');continue
        im=re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',s)
        if im:
            out.append('<figure><img loading="lazy" src="'+html.escape(href(im[2]),quote=True)+'" alt="'+html.escape(im[1],quote=True)+'"><figcaption>'+html.escape(im[1])+'</figcaption></figure>');i+=1;continue
        if s.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=lines[i].strip().strip('|').split('|')
                if not all(re.fullmatch(r'\s*:?-+:?\s*',c) for c in cells):rows.append(cells)
                i+=1
            out.append('<div class="table-scroll"><table>'+''.join('<tr>'+''.join('<'+('th' if n==0 else 'td')+'>'+inline(c.strip())+'</'+('th' if n==0 else 'td')+'>' for c in row)+'</tr>' for n,row in enumerate(rows))+'</table></div>');continue
        h=re.match(r'^(#{1,6})\s+(.*)',s)
        if h:
            n=len(h[1]);out.append(f'<h{n}>'+inline(h[2])+f'</h{n}>');i+=1;continue
        if re.match(r'^(?:[-*]|\d+[.)])\s+',s):
            ordered=bool(re.match(r'^\d+',s));tag='ol' if ordered else 'ul';items=[]
            pattern=r'^\d+[.)]\s+' if ordered else r'^[-*]\s+'
            while i<len(lines) and re.match(pattern,lines[i].strip()):
                items.append(re.sub(pattern,'',lines[i].strip()));i+=1
            out.append('<'+tag+'>'+''.join('<li>'+inline(x)+'</li>' for x in items)+'</'+tag+'>');continue
        para=[s];i+=1
        while i<len(lines) and lines[i].strip() and not re.match(r'^(?:#|\||```|!\[|\\\[|[-*]\s|\d+[.)]\s)',lines[i].strip()):
            para.append(lines[i].strip());i+=1
        out.append('<p>'+inline(' '.join(para))+'</p>')
    return '\n'.join(out)

def page(title,body):
    nav='<nav><a href="../../">B1 首页</a><a href="index.html">本轮概览</a><a href="report.html">完整中文报告</a><a href="gallery.html">实际图表</a><a href="handoff.html">Pro 交接</a><a href="evidence.html">复现与数据</a><a href="'+RAW+'REPORT_ZH.md">匿名 Markdown</a></nav>'
    footer='<footer>2026-10-29 · Research W-GQ-R2 · 二维 scalar-TM，静态、有限开发合同。物理机制证据、在线识别和计算收益分别裁决。历史结果与费用保持原样。<br><a href="'+BLOB+'GO_NO_GO.json">机器裁决</a> · <a href="'+BLOB+'PHYSICS_COST.json">完整费用</a> · <a href="'+BLOB+'FAILURE_LEDGER.jsonl">失败记录</a> · <a href="'+RELEASE+'">完整证据 Release</a></footer>'
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Research W-GQ-R2：真实二维微结构墙体、目标分辨率、有限非线性补偿与墙返回压缩实验。"><title>'+html.escape(title)+'</title><style>'+STYLE+'</style>'+MATH+'</head><body><main>'+nav+body+footer+'</main></body></html>'

def figure(name,caption):
    return '<figure><a href="figures/'+name+'.png"><img loading="lazy" src="figures/'+name+'.png" alt="'+html.escape(caption,quote=True)+'"></a><figcaption>'+html.escape(caption)+'</figcaption><a href="figures/'+name+'.pdf">矢量 PDF</a></figure>'

def main():
    assert (ROOT/'PHYSICS_CLOSED.json').exists()
    body='''<div class="eyebrow">RESEARCH W-GQ-R2 · 实际全波研究 · 2026-10-29</div><h1>墙体材料可以不准，<br>目标仍能被分辨吗？</h1><p class="lede">新的二维 Maxwell 正例与反例同时成立：微观材料很不一样，有限目标族的传播仍可近似等效；但更丰富的墙体补偿，也能掩盖目标位移。</p><div class="status"><strong>裁决：目标分辨率的物理主线值得继续；一般未知墙在线成像和新的 Bayesian 净加速尚未建立。</strong></div><div class="metrics"><div class="metric"><strong>1.13137</strong><span>两纹理墙内材料 RMS 差</span></div><div class="metric"><strong>0.159971</strong><span>held-out 目标族最大白化传播差</span></div><div class="metric"><strong>0.079224 λ₀</strong><span>受控有限参考的位置均值 RMSE</span></div><div class="metric"><strong>722 / 5212</strong><span>实际 FE / complex RHS</span></div></div><p><a href="report.html">阅读完整中文报告 →</a>　<a href="handoff.html">给 GPT Pro 的交接 →</a>　<a href="evidence.html">下载完整复现包 →</a></p>'''
    body+='<section class="panel"><h2>材料不同，目标传播仍接近</h2><p>两个宏观非均匀墙体，其中一个有不规则边界与空腔；两种微结构共享局部平均材料。网格逐条解析微结构界面。材料差保持显著，传播差随尺度下降；最细层覆盖 12 个 held-out 目标、四发三频。独立细化支持有限离散任务的近等效，不是连续目标族的统一证明。</p>'+figure('F02_MICRO_SCALE_AND_ALL_PROBES','真实材料 RMS、共同探针尺度趋势与最细尺度全部目标。')+'</section>'
    body+='<section class="panel"><h2>目标准确，与墙细节仍未知可以共存</h2><p>在明确的 4 墙×12 目标有限 prior 中，位置均值 RMSE 降至 0.079224 λ₀，而墙纹理熵仍为 0.998388 bit。这是受控参考；另一个未知宏观墙在线模型的 RMSE 为 0.68343 λ₀，失败不能被参考成功覆盖。</p>'+figure('F13_CONTROLLED_MICRO_RISK_TEXTURE_AND_TV_TAILS','目标风险下降、纹理仍未知；所有条件共识失败尾部保留。')+'</section>'
    body+='<section class="panel"><h2>反例：新增测量不保证排除重新补偿</h2><p>目标横移 0.05 λ₀，合法墙材料重新优化后，四发三频白化距离在三层离散上为 0.889152 / 0.912855 / 0.913910。它支持有限噪声尺度歧义；不是连续 PDE 的精确相等或全局最优证明。</p>'+figure('F14_FINAL_FINITE_COMPENSATION_AND_THIRD_GRID','真实补偿轨迹、被动材料端点及三层独立网格。')+'</section>'
    body+='<section class="panel"><h2>返回传播可压缩，净加速还没有证据</h2><p>保留精确目标自散射，只压缩墙返回，rank 8 最坏白化误差约 2.05×10⁻⁴。但构建与审核真实消耗 964 RHS / 22 FE；旧 118 墙 bank 的条件共识没有通过。没有把低 rank 或表查询减少写成冷部署加速。</p>'+figure('F08_RETURN_RANK_ERROR_SPECTRUM_AND_PAID_COST','原始失败、复对称配对修复、逐秩误差及全部已付费用。')+'</section>'
    (ROOT/'index.html').write_text(page('Research W-GQ-R2｜未分辨墙体后的目标分辨率',body))
    (ROOT/'report.html').write_text(page('Research W-GQ-R2｜完整中文科研报告','<article>'+render((ROOT/'REPORT_ZH.md').read_text())+'</article>'))
    (ROOT/'handoff.html').write_text(page('Research W-GQ-R2｜GPT Pro 研究交接','<article>'+render((ROOT/'NEXT_PRO_HANDOFF.md').read_text())+'</article>'))
    figs=json.loads((ROOT/'figures/FIGURE_MANIFEST.json').read_text())['figures']
    gallery='<h1>实际材料、复场与目标结果</h1><p>全部图表来自保存的实验数组。点击图片查看原尺寸；原始 CSV、场、后验和作业回执在完整复现包中。各图的范围限制保留如下。</p>'
    for name,meta in sorted(figs.items()):
        gallery+='<section class="panel"><h2>'+name[:3]+'</h2>'+figure(name,meta['explanation_zh'])+'<p class="muted">'+html.escape(meta['scope_zh'])+'</p></section>'
    (ROOT/'gallery.html').write_text(page('Research W-GQ-R2｜全部真实图表',gallery))
    evidence='<article><h1>源码、原始数据与复现</h1><p><strong><a href="'+ASSET+'">下载完整 Research_W_GQ_R2_Results.zip（400,963,935 bytes，约 401 MB）</a></strong></p><p>包内 4,427 项，含原始复响应、必要 P3 场、Jacobian、完整有限后验、材料状态、失败与作业回执、源码及附件理论。没有 locked 数据。解压后入口为 <code>Research_W_GQ_R2/</code>。</p><p><a href="'+RELEASE+'">Release 说明</a> · <a href="'+BLOB+'REPRODUCE.md">完整复现命令</a> · <a href="'+BLOB+'publication/PREPARATION.json">Git 与大数据的分配清单</a></p><h2>结果与费用</h2><ul>'
    for n in ['E1_RESULTS.csv','E2_RESULTS.csv','E3_RESULTS.csv','EVIDENCE_LEDGER.csv','GO_NO_GO.json','PHYSICS_COST.json','PHYSICS_LEDGER.jsonl','QUERY_LEDGER.jsonl','FAILURE_LEDGER.jsonl','CONTRACT.json','THEOREM_STATUS.md']:
        evidence+='<li><a href="'+BLOB+n+'">'+n+'</a></li>'
    evidence+='</ul><h2>读取顺序与状态</h2><p>先读完整报告与 Pro 交接，再沿结果表定位实际数组。仓库直接保留报告、源码、图表、配置、紧凑指标与账本；NPZ/NPY、大型逐点 CSV 只在 Release 中。上轮本地 GQ 的必要源码已作为 <code>reuse/previous_gq</code> 随附，旧 BGGA/W-SLAM 后端来自固定历史版本。</p><p>研究已封账，保存数据回放不做新 PDE。全新物理复现需建立独立输出 scope；不要删除原账本或解除当前封账来覆盖结果。</p><h2>费用的正确版本</h2><p>科研最终费用为 722 FE、5,212 complex RHS、665.087712 wrapper process CPU 秒，峰值 RSS 3,348,283,392 bytes。早期 <code>PACKAGE_RECEIPT.json</code> 是打包过程的历史快照；归档中以 <code>delivery/FINAL_ACCOUNTING.json</code> 和 <code>delivery/PACKAGE_RECEIPT.json</code> 为最终版本。发布的准备、上传与核验另列，不改变科研账本。</p><p>这些是二维、有限任务和开发样本的实验；不是连续 Maxwell 认证、任意未知墙成功或同质量净加速证明。</p></article>'
    (ROOT/'evidence.html').write_text(page('Research W-GQ-R2｜完整证据与复现下载',evidence))
    mark='W_GQ_R2_PUBLICATION'
    summary='微结构材料 RMS 差 1.13137，held-out 目标族最大传播白化距离 0.159971；受控有限参考定位 0.079224 λ₀，宽墙材料仍可补偿 0.05 λ₀ 目标位移。一般未知墙在线与净加速未成立。722 FE / 5212 RHS，全部失败保留。'
    entry='<!-- '+mark+' -->\n## Research W-GQ-R2：墙体未分辨，目标分辨率与合法补偿\n\n'+summary+'\n\n[完整中文网页]('+PAGE+') · [中文 Markdown]('+REL+'REPORT_ZH.md) · [GPT Pro 交接]('+REL+'NEXT_PRO_HANDOFF.md) · [完整原始证据与复现包]('+RELEASE+')\n<!-- /'+mark+' -->\n\n'
    for n in ['README.md','START_HERE.md']:
        p=REPO/n;s=p.read_text()
        if mark not in s:p.write_text(entry+s)
    for n,link in [('index.html',REL),('research/index.html','research_w_gq_r2/')]:
        p=REPO/n;s=p.read_text()
        if mark not in s:
            card='<!-- '+mark+' --><section class="card"><h2>Research W-GQ-R2 · 未分辨墙体后的目标分辨率</h2><p>'+summary+'</p><a href="'+link+'">中文报告、真实图表、反例与完整复现包 →</a></section><!-- /'+mark+' -->'
            p.write_text(s.replace('<main>','<main>'+card,1))
    print(json.dumps(dict(pages=5,figures=len(figs),reports_unchanged=True,new_physics=0),ensure_ascii=False))

if __name__=='__main__':main()

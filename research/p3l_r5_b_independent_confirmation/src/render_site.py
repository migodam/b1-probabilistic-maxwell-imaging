"""Existing B1 presentation, new confirmation pages only. No solver imports."""
from pathlib import Path
import importlib.util,html,re,os,json
P=Path(__file__).resolve().parents[1];ROOT=P.parents[1]
spec=importlib.util.spec_from_file_location('b1_presentation',ROOT/'tools/render_r4_final_site.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/p3l_r5_b_independent_confirmation/RESEARCH_REPORT_ZH.md'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr5-b-confirm-20261027'
def render(path):
    text=path.read_text();body=m.markdown(text)
    body=re.sub(r'<p>!<a href="([^"]+)">([^<]*)</a></p>',r'<figure><img loading="lazy" src="\1" alt="\2"><figcaption>\2</figcaption></figure>',body)
    def link(match):
        url=html.unescape(match[1])
        if re.match(r'^(https?:|#|mailto:)',url):return match[0]
        target=(path.parent/url.split('#')[0]).resolve()
        if url.endswith('.md') and target.is_relative_to(P) and target.is_file():return 'href="'+html.escape(url[:-3]+'.html',quote=True)+'"'
        if url.endswith(('.npz','.pt','.jsonl','.npy')):return 'href="'+RELEASE+'"'
        return match[0]
    body=re.sub(r'href="([^"]+)"',link,body)
    title=next((line.lstrip('# ') for line in text.splitlines() if line.startswith('#')),'R5-B Independent Confirmation')
    nav='<nav>'+''.join('<a href="'+html.escape(os.path.relpath(P/t,path.parent))+'">'+label+'</a>' for t,label in [
        ('report.html','完整报告'),('POOLED_512_STATISTICS.html','独立统计'),('MATERIAL_RISK_COMPARISON.html','材料风险'),
        ('ONLINE_AND_END_TO_END_COST.html','完整成本'),('GO_NO_GO.html','科研裁决'),('EVIDENCE_PACKAGE.html','源码与证据'),('GPT_PRO_HANDOFF/START_HERE.html','Pro交接')])+'<a href="'+RAW+'">匿名Markdown</a></nav>'
    style=m.STYLE+'figure{margin:22px 0}figure img{width:100%;height:auto}figcaption{font-size:14px;color:var(--sub)}.note{border-left:4px solid #a67822;padding:16px;background:#fff9e9}'
    head='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+style+'</style>'
    cards=''
    if path==P/'RESEARCH_REPORT_ZH.md':
        cards='<div class="cards">'+''.join('<div class="card"><b>'+v+'</b><span>'+label+'</span></div>' for v,label in [
           ('512','全新独立parents'),('80.70%','增量保留率点估计'),('HOLD','≥80%统计门尚未通过'),('1,999,504','实际新增complex RHS')])+'</div>'
        cards+='<p class="note">蒸馏相对Strong CNN有风险收益；确定性矩相对K32改善成立。80%保留率区间跨门，全系统推导节省仅0.0117%。模型、prior和物理配置未因结果更改。</p>'
    page=head+'<header>'+nav+'<h1>'+html.escape(title)+'</h1><p>冻结模型 · 新256+256对象 · 同配置TM24物理 · 预测先封存 · 2026-10-27</p></header><main>'+cards+'<article class="report">'+body+'</article><footer>本页只读呈现实际证据；网站发布不改变统计、物理或论文资格。完整原始数组与账本随同仓库Release提供。</footer></main></html>'
    path.with_suffix('.html').write_text(page)
    if path==P/'RESEARCH_REPORT_ZH.md':(P/'report.html').write_text(page);(P/'index.html').write_text(page)
def main():
    files=[f for f in P.rglob('*.md') if not any(x in f.parts for x in ('raw','__pycache__'))]
    for f in files:render(f)
    print(json.dumps(dict(status='RENDERED',pages=len(files)+2,new_RHS=0)))
if __name__=='__main__':main()

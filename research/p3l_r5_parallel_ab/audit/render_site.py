"""Offline R5 report renderer, reusing the existing B1 report presentation."""
from pathlib import Path
import importlib.util,html,re,os,json
P=Path(__file__).resolve().parents[1];ROOT=P.parents[1]
spec=importlib.util.spec_from_file_location('r4_presentation',ROOT/'tools/render_r4_final_site.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/p3l_r5_parallel_ab/RESEARCH_REPORT_ZH.md'
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/p3lr5-parallel-ab-20261026'

def rel(target,fromfile):return os.path.relpath(target,fromfile.parent)
def render(path):
    text=path.read_text();body=m.markdown(text)
    body=re.sub(r'<p>!<a href="([^"]+)">([^<]*)</a></p>',r'<figure><img loading="lazy" src="\1" alt="\2"><figcaption>\2</figcaption></figure>',body)
    def link(match):
        url=html.unescape(match[1])
        if re.match(r'^(https?:|#|mailto:)',url):return match[0]
        target=(path.parent/url.split('#')[0]).resolve()
        if url.endswith('.md') and target.is_relative_to(P) and target.is_file():
            return 'href="'+html.escape(url[:-3]+'.html',quote=True)+'"'
        if url.endswith(('.npz','.pt','.npy')):return 'href="'+RELEASE+'"'
        return match[0]
    body=re.sub(r'href="([^"]+)"',link,body)
    main=path==P/'RESEARCH_REPORT_ZH.md'
    title=next((s.lstrip('# ') for s in text.splitlines() if s.startswith('#')),'Research P3L-R5')
    nav='<nav>'+''.join(f'<a href="{html.escape(rel(P/t,path))}">{label}</a>' for t,label in [
      ('report.html','完整报告'),('stage_a_response_coordinates/RESEARCH_REPORT_ZH.html','A 概率实验'),
      ('stage_b_material_risk/REPORT_ZH.html','B 风险压缩'),('stage_b_material_risk/GALLERY.html','材料图库'),
      ('GO_NO_GO.html','科研裁决'),('EVIDENCE_PACKAGE.html','源码与证据')])+'<a href="'+RAW+'">匿名Markdown</a></nav>'
    style=m.STYLE+'figure{margin:22px 0}figure img{width:100%;height:auto}figcaption{font-size:14px;color:var(--sub)}'
    head='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+style+'</style>'
    cards=''
    if main:
        cards='<div class="cards">'+''.join('<div class="card"><b>'+v+'</b><span>'+l+'</span></div>' for v,l in [
          ('18 / 6','实际Maxwell观测 / 独立parent'),('291,784','A实际新增complex RHS'),
          ('92.96%','增量风险收益保留，开发VAL256'),('2.54%','额外provider成本比，非整个系统')])+'</div>'
        cards+='<p class="note">A全域概率覆盖未通过；B新256+256确认未运行。开发集保留率区间跨80%门槛，不将其写成独立泛化或posterior校准。共同Maxwell反演成本完整保留。</p>'
    page=head+'<header>'+nav+'<h1>'+html.escape(title)+'</h1><p>原P1统计合同 · 实际TM24 · A/B独立实验 · 正负证据与封账</p></header><main>'+cards+'<article class="report">'+body+'</article><footer>本页是研究证据入口；网页发布不改变实验资格。原始数组和模型位于同仓库独立Release。历史研究未覆盖。</footer></main></html>'
    path.with_suffix('.html').write_text(page)
    if main:(P/'report.html').write_text(page);(P/'index.html').write_text(page)

def main():
    files=[x for x in P.rglob('*.md') if '__pycache__' not in x.parts]
    for path in files:render(path)
    print(json.dumps(dict(status='RENDERED',pages=len(files)+2,physics_RHS=0)))
if __name__=='__main__':main()

"""Static report page in existing B1 site; no scientific work or new hosting."""
from pathlib import Path
import html,re
ROOT=Path(__file__).resolve().parents[1]
BASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_w_bgga/'
def inline(s):
    s=html.escape(s)
    def link(m):
        path=html.unescape(m[2]);href=BASE+path if not path.startswith(('http','#')) and path.endswith(('.md','.json','.csv','.jsonl')) else path
        return '<a href="'+html.escape(href,quote=True)+'">'+m[1]+'</a>'
    s=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',link,s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s)
    return re.sub(r'`([^`]+)`',r'<code>\1</code>',s)
def render(md):
    lines=md.splitlines();result=[];i=0
    while i<len(lines):
        s=lines[i]
        if not s.strip():i+=1;continue
        if s.startswith('```'):
            i+=1;code=[]
            while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
            result.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>');i+=1;continue
        if s.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].startswith('|'):
                cells=lines[i].strip('|').split('|')
                if not all(re.match(r'^\s*:?-+:?\s*$',c) for c in cells):rows.append(cells)
                i+=1
            result.append('<div class="scroll"><table>'+''.join('<tr>'+''.join('<'+('th' if n==0 else 'td')+'>'+inline(c.strip())+'</'+('th' if n==0 else 'td')+'>' for c in row)+'</tr>' for n,row in enumerate(rows))+'</table></div>');continue
        image=re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',s.strip())
        if image:
            result.append('<figure><img src="'+html.escape(image[2])+'" alt="'+html.escape(image[1])+'"><figcaption>'+html.escape(image[1])+'</figcaption></figure>');i+=1;continue
        if s.startswith('#'):
            n=min(6,len(s)-len(s.lstrip('#')));result.append(f'<h{n}>'+inline(s[n:].strip())+f'</h{n}>')
        else:result.append('<p>'+inline(s)+'</p>')
        i+=1
    return '\n'.join(result)
def main():
    style='body{background:#f3f6f8;color:#20323f;font:17px/1.85 system-ui,sans-serif;margin:0}main{max-width:1000px;margin:auto;padding:35px 24px}article{background:white;border:1px solid #dce3ea;border-radius:12px;padding:28px}h1{line-height:1.4;font-size:34px}h2{margin-top:2em}a{color:#126d91}nav{display:flex;flex-wrap:wrap;gap:12px 20px;margin-bottom:24px}p{overflow-wrap:anywhere}img{max-width:100%;height:auto}code,pre{background:#edf2f5;font-size:.9em}pre{padding:14px;overflow:auto}.scroll{overflow:auto}table{border-collapse:collapse;width:100%}td,th{border:1px solid #d5dfe7;padding:9px;text-align:left}figcaption{font-size:14px;color:#536371}.status{border-left:4px solid #bc7f25;background:#fff1d8;padding:12px 18px;margin:20px 0}'
    nav='<nav><a href="../">B1研究入口</a><a href="'+BASE+'START_HERE.md">源码与结果</a><a href="'+BASE+'GATE_DECISION.json">分层判决</a><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-bgga-20261028">原始数据Release</a><a href="https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_bgga/RESEARCH_W_BGGA_REPORT_ZH.md">匿名raw</a></nav>'
    page='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research W-BGGA｜精确吸收与计算成本判决</title><style>'+style+'</style></head><body><main>'+nav+'<div class="status">G0真实全波3934 RHS / 482 FE；G1为旧bank重放。精确离散吸收通过，当前加速路线NO-GO；G2未启动。</div><article>'+render((ROOT/'RESEARCH_W_BGGA_REPORT_ZH.md').read_text())+'</article></main></body></html>'
    (ROOT/'index.html').write_text(page)

if __name__=='__main__':main()

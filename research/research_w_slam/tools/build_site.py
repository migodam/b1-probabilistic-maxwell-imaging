"""Static reading pages from actual reports/figures, no physics or new claims."""
from pathlib import Path
import html,re
ROOT=Path(__file__).resolve().parents[1]
BASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_w_slam/'
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
            result.append('<figure><a href="'+html.escape(image[2])+'"><img loading="lazy" src="'+html.escape(image[2])+'" alt="'+html.escape(image[1])+'"></a><figcaption>'+html.escape(image[1])+'</figcaption></figure>');i+=1;continue
        if s.startswith('#'):
            n=min(6,len(s)-len(s.lstrip('#')));result.append(f'<h{n}>'+inline(s[n:].strip())+f'</h{n}>')
        else:result.append('<p>'+inline(s)+'</p>')
        i+=1
    return '\n'.join(result)
STYLE='body{margin:0;background:#f4f7f9;color:#20323f;font:17px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif}main{max-width:1060px;padding:38px 24px 70px;margin:auto}article{background:white;padding:28px;border:1px solid #dbe4eb;border-radius:12px}h1{font-size:clamp(28px,4vw,43px);line-height:1.35}h2{margin:2em 0 .6em;font-size:25px}h3{font-size:20px}nav{display:flex;gap:12px 24px;flex-wrap:wrap;margin:20px 0}a{color:#126d91;text-underline-offset:4px}figure{margin:30px 0}img{width:100%;height:auto}figcaption,footer{color:#5c6b76;font-size:14px}table{border-collapse:collapse;width:100%;font-size:15px}td,th{border:1px solid #d5dfe7;padding:9px;text-align:left}th{background:#eaf2f6}.scroll{overflow:auto}code{font-size:.88em;overflow-wrap:anywhere;background:#edf2f5}pre{padding:15px;background:#edf2f5;overflow:auto}.status{padding:12px 18px;border-left:4px solid #ba7b29;background:#fff2db;margin:22px 0}p{overflow-wrap:anywhere}'
def main():
    title='Research W-SLAM｜Maxwell墙–目标联合反演'
    nav='<nav><a href="../../">B1研究总入口</a><a href="report.html">中文研究报告</a><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-w-slam-20261025">原始数据Release</a><a href="'+BASE+'START_HERE.md">GPT阅读入口</a><a href="'+BASE+'REPRODUCE.md">复现</a><a href="https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_w_slam/RESEARCH_W_SLAM_REPORT_ZH.md">匿名raw</a></nav>'
    prefix='<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+title+'</title><style>'+STYLE+'</style></head><body><main>'+nav
    report=(ROOT/'RESEARCH_W_SLAM_REPORT_ZH.md').read_text()
    page=prefix+'<article>'+render(report)+'</article><footer>二维标量TM；状态与费用以实际报告和账本为准。没有把未执行阶段写成成功。</footer></main></body></html>'
    (ROOT/'index.html').write_text(page);(ROOT/'report.html').write_text(page)
if __name__=='__main__':main()

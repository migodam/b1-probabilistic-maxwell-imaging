"""Render the completed local research documents into the existing Pages tree.

No physics, training, metrics recomputation, Git action, or remote publication.
"""
from pathlib import Path
import html,json,re

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r2_risk_posterior_exploration'


def inline(text):
    text=html.escape(text)
    text=re.sub(r'!\[([^\]]*)\]\(([^)]+)\)',r'<a href="\2"><img loading="lazy" src="\2" alt="\1"></a>',text)
    text=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',text)
    text=re.sub(r'`([^`]+)`',r'<code>\1</code>',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',text)
    return text


def markdown(source):
    lines=source.splitlines();out=[];paragraph=[];code=None;listing=False;i=0
    def flush():
        if paragraph:out.append('<p>'+inline(' '.join(paragraph))+'</p>');paragraph.clear()
    def close_list():
        nonlocal listing
        if listing:out.append('</ul>');listing=False
    while i<len(lines):
        line=lines[i]
        if line.startswith('```'):
            flush();close_list()
            if code is None:code=[]
            else:out.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>');code=None
        elif code is not None:code.append(line)
        elif line.startswith('|') and i+1<len(lines) and re.fullmatch(r'[| :\-]+',lines[i+1]):
            flush();close_list();heads=[x.strip() for x in line.strip('|').split('|')]
            out.append('<div class="table-scroll"><table><thead><tr>'+''.join('<th>'+inline(x)+'</th>' for x in heads)+'</tr></thead><tbody>');i+=2
            while i<len(lines) and lines[i].startswith('|'):
                values=[x.strip() for x in lines[i].strip('|').split('|')]
                out.append('<tr>'+''.join('<td>'+inline(x)+'</td>' for x in values)+'</tr>');i+=1
            out.append('</tbody></table></div>');continue
        elif re.match(r'^#{1,6} ',line):
            flush();close_list();n=len(line)-len(line.lstrip('#'));out.append(f'<h{n}>'+inline(line[n+1:])+f'</h{n}>')
        elif line.startswith('- '):
            flush()
            if not listing:out.append('<ul>');listing=True
            out.append('<li>'+inline(line[2:])+'</li>')
        elif line.startswith('> '):flush();close_list();out.append('<blockquote>'+inline(line[2:])+'</blockquote>')
        elif not line.strip():flush();close_list()
        else:close_list();paragraph.append(line)
        i+=1
    flush();close_list()
    if code is not None:raise ValueError('Unclosed code block in final research document')
    return '\n'.join(out)


STYLE='''*{box-sizing:border-box}body{margin:0;background:#f3f6fa;color:#182b40;font:17px/1.8 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}header,main,footer{max-width:1150px;margin:auto;padding:28px}header{padding-top:56px}h1{font-size:clamp(27px,4vw,43px);line-height:1.35;margin:.4em 0}h2{margin-top:2.3em;border-top:1px solid #dae3eb;padding-top:22px}h3{margin-top:1.8em}a{color:#086ba5;text-decoration:none}a:hover{text-decoration:underline}main{background:white;border-radius:16px;box-shadow:0 8px 32px #1c365510;padding:36px}p{margin:1em 0}nav{display:flex;flex-wrap:wrap;gap:10px 20px;font-size:15px}.badge{display:inline-block;border-radius:20px;padding:4px 14px;background:#dcefe5;color:#1b5b3f;font-size:15px}.scope{background:#edf3f9;border-left:4px solid #739abb;padding:16px 22px}.table-scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:15px;margin:22px 0}th,td{text-align:left;vertical-align:top;border-bottom:1px solid #dce4ed;padding:11px 12px}th{background:#edf3f9}img{display:block;max-width:100%;height:auto;margin:20px auto}code{font:14px/1.7 ui-monospace,monospace;background:#f0f4f8;border-radius:4px;padding:2px 5px}pre{overflow:auto;padding:20px;background:#eef3f8;white-space:pre-wrap}blockquote{background:#fff7e8;border-left:4px solid #c9922a;padding:12px 22px;margin:20px 0}footer{font-size:14px;color:#526579}.hero{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:14px;margin:24px 0}.metric{padding:18px;background:#fff;border-radius:10px;border:1px solid #dce5ed}.metric strong{font-size:30px;display:block;color:#135d71}.gallery{display:grid;gap:25px}.gallery figure{margin:0}.gallery figcaption{font-size:14px;color:#526579}@media(max-width:650px){header,main,footer{padding:20px}main{border-radius:0}th,td{padding:8px}}'''


def document(title,body):
    nav='<nav><a href="index.html">研究入口</a><a href="report.html">完整中文报告</a><a href="gallery.html">成像图库</a><a href="START_HERE.md">匿名 Markdown</a><a href="REPRODUCE.md">复现</a><a href="GPT_PRO_HANDOFF/START_HERE.md">GPT Pro 交接</a></nav>'
    return '<!doctype html><html lang="zh"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+STYLE+'</style></head><body><header>'+nav+'<h1>'+html.escape(title)+'</h1></header><main>'+body+'</main><footer>固定 P1 · 24×24 二维标量 TM/VIE · 真实 Full Maxwell 数据 · 原始负结果与费用保留。网页可访问性与科研结论分别记录。</footer></body></html>'


def main():
    required=['START_HERE.md','RESEARCH_REPORT_ZH.md','GO_NO_GO.md','LOCKED_GENERALIZATION.md']
    for name in required:
        if not (OUT/name).is_file():raise FileNotFoundError('Completed final document required: '+name)
    report=(OUT/'RESEARCH_REPORT_ZH.md').read_text()
    (OUT/'report.html').write_text(document('P3L-R2：材料风险、非 Gaussian 推断与局部几何',markdown(report)))
    intro=markdown((OUT/'START_HERE.md').read_text())
    (OUT/'index.html').write_text(document('P3L-R2 · Maxwell 材料风险研究',intro))
    gallery=['<p class="scope">以下固定展示首批测试 IDs 00000–00003，未按误差挑图。实部与损耗使用统一色标。GT 只用于已完成线上估计后的独立评分；学习图和轨迹平均是风险动作，不是 posterior samples，也不保证是一份 decoder 材料。</p>','<div class="gallery">']
    for i in range(4):
        name=f'R2_locked_{i:05d}'
        gallery.append('<figure><a href="figures/main_v2/'+name+'.png"><img loading="lazy" src="figures/main_v2/'+name+'.png" alt="'+name+'"></a><figcaption>'+name+'：GT、native点、强Probit点、额外16-state普通优化、监督风险动作、曲率轨迹弧长动作。</figcaption></figure>')
    gallery+=['</div>',markdown((OUT/'FAILURE_ANALYSIS.md').read_text())]
    (OUT/'gallery.html').write_text(document('固定对象材料图库与失败分析','\n'.join(gallery)))
    pages={'supervised.html':'SUPERVISED_RISK_RESULTS.md','posterior.html':'NON_GAUSSIAN_POSTERIOR_REPORT.md','geometry.html':'GEOMETRY_TRAJECTORY_RESULTS.md','models.html':'LEARNED_MODEL_COMPARISON.md','locked.html':'LOCKED_GENERALIZATION.md','decision.html':'GO_NO_GO.md'}
    for filename,source in pages.items():
        text=(OUT/source).read_text();title=text.splitlines()[0].lstrip('# ')
        (OUT/filename).write_text(document(title,markdown(text)))
    print('Rendered final documents; no new physics or inference')


if __name__=='__main__':main()

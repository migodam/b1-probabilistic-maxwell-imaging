"""Render the Chinese research report and handoff beside the manuscript."""
from pathlib import Path
import html
import re
import markdown
P=Path(__file__).resolve().parents[1]
RELEASE='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/download/z-alg-closed-loop-20261126/closed_loop_reproduction.zip'
def layout(title,body):
    t=html.escape(title)
    return f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{t}</title><link rel="stylesheet" href="report.css"></head>
<body><nav aria-label="研究导航"><a href="../../../index.html">Maxwell 逆散射研究</a><a href="index.html">本研究</a><a href="report.html">中文报告</a><a href="manuscript_ip.pdf">英文论文</a><a href="handoff.html">研究交接</a></nav>
<main>{body}</main><footer>固定阵列中的介电缺陷定位与照明选择 · 2026 年 11 月</footer></body></html>"""
def document(source):
    s=(P/source).read_text()
    def fmt(m):
        f=m.group(1)
        for a,b in [(r'\preceq','≼'),(r'\le','≤'),(r'\chi','χ'),(r'\ell','ℓ'),(r'\nu','ν'),(r'\times','×'),('^*','*')]:f=f.replace(a,b)
        f=re.sub(r'_([iok])',r'<sub>\1</sub>',f)
        return '<span class="math">'+f+'</span>'
    s=re.sub(r'\\\((.*?)\\\)',fmt,s,flags=re.S)
    for name, page in [('RESEARCH_REPORT_ZH.md', 'report.html'), ('NEXT_PRO_HANDOFF.md', 'handoff.html'), ('REPRODUCE.md', 'reproduce.html'), ('reviews/RESPONSE_ZH.md', 'review.html')]:
        s=s.replace('('+name+')','('+page+')')
    return markdown.markdown(s,extensions=['tables','fenced_code','sane_lists']).replace('<table>','<div class="table"><table>').replace('</table>','</table></div>')
def build():
    for src,out,title in [('RESEARCH_REPORT_ZH.md','report.html','完整中文研究报告'),('NEXT_PRO_HANDOFF.md','handoff.html','固定阵列逆散射研究交接'),('REPRODUCE.md','reproduce.html','复现实验'),('reviews/RESPONSE_ZH.md','review.html','独立评审与修改')]:
        body=document(src)
        if src.startswith('reviews/'):
            body=re.sub(r'href="neon-', 'href="reviews/neon-',body).replace('href="../','href="')
        if out in ('handoff.html','reproduce.html','report.html'):
            body+=f'<p class="actions"><a href="{RELEASE}">下载完整代码与数据</a><a href="IP_Anonymous_LaTeX.zip">下载匿名投稿源码</a></p>'
        (P/out).write_text(layout(title,body),encoding='utf-8')
    home=f"""<div class="rule"></div><h1>能把介电缺陷定位到哪里，<br>下一束波该怎样发？</h1>
<p class="lead">只用一侧的复散射读数，算法逐块排除不可能的缺陷位置，再用完整材料验证留下的位置。对仍能解释旧数据的两类材料，还能在发射前判断需要多少能量、至少要同时控制几个方向。</p>
<div class="actions"><a class="primary" href="manuscript_ip.pdf">阅读英文论文 PDF</a><a href="report.html">完整中文报告</a><a href="handoff.html">ChatGPT 研究交接</a></div>
<div class="twocol"><div><strong>8,243 / 9,216 个中心被排除</strong><p>96² 全视野强散射环，背景和边缘层均可变。余下 973 个中心，一个有完整的相容材料。</p></div><div><strong>最少需要两路相干照明</strong><p>一个原生案例中，15 个单角分别都有共同读数。按选定能量同时发射 67.5° 和 90° 两个方向即可区分全部两类材料。</p></div></div>
<figure><img src="figures/location.png" width="1180" height="794" alt="强散射环和跨边缘双夹杂的真实材料、七束波后的可行位置和补测后的可行位置"><figcaption>左列为生成材料。中列为覆盖全视野的位置判断，右列为追加一束波后的结果。中、右列的星号均已用完整材料通过对应读数的误差检验。橙色表示仍有可能的位置。</figcaption></figure>
<h2>加测的收益取决于要判断什么</h2><p>在 16 件 32² 物体上，七次旧照明后只剩 71 个候选中心。进一步定位时，选得好的单角与相干波束接近。对于两个材料类的最坏情况判别，同能量单角仍可能产生共同读数，而相干发射可以得到严格分离。给普通全波模板拟合核心强度和均匀背景，它也找对 18/18 个中心，其中 1 份材料满足噪声球；允许逐格变化后，18 件都找到了相容材料。这些结论依赖给定轮廓、公差和离散模型。</p>
<p><a href="maps/index.html">逐物体定位与噪声结果</a> · <a href="figures/energy.pdf">能量与硬件图</a> · <a href="review.html">独立审查与修改</a></p>
<h2>研究材料</h2><p>匿名论文包可独立编译。复现包保留模拟原始读数、材料见证、独立数值复核和计时，以及两个暗室的派生统计。实验室原始测量由各发布方提供。</p>
<div class="actions"><a href="IP_Anonymous_LaTeX.zip">匿名 LaTeX 源码</a><a href="{RELEASE}">实验代码与数据 ZIP</a><a href="reproduce.html">复现说明</a><a href="lab/held_ports.csv">暗室逐项数据 CSV</a></div>"""
    (P/'index.html').write_text(layout('介电缺陷定位与下一束照明',home),encoding='utf-8')
if __name__=='__main__':build()

"""Offline static research report and figure gallery; no remote runtime needed."""
from pathlib import Path
import html
import json
import re

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'research/p3l_r4_conditional_transport'

STYLE='''
:root{--ink:#183247;--sub:#506778;--blue:#0a6283;--line:#d5e0e7;--paper:#fff;--bg:#f3f7fa}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:17px/1.85 -apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC",sans-serif}a{color:var(--blue)}
header{background:#102e43;color:#f7fbff;padding:44px max(24px,calc((100vw - 1120px)/2)) 34px}header h1{font-size:clamp(27px,4vw,42px);line-height:1.35;margin:14px 0}header p{max-width:1000px;color:#cee1ed}
nav{display:flex;gap:12px;flex-wrap:wrap}header nav a{color:#edf7ff;text-decoration:none;border:1px solid #668192;border-radius:7px;padding:3px 12px}
main{max-width:1170px;margin:24px auto;padding:0 24px 60px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:12px;margin:22px 0}.card{background:white;border:1px solid var(--line);border-radius:10px;padding:18px}.card b{font-size:27px;display:block}.card span{font-size:14px;color:var(--sub)}.note{border-left:4px solid #b18b38;background:#fff9eb;padding:15px 20px}.report{background:white;border:1px solid var(--line);border-radius:12px;padding:clamp(18px,4vw,46px);overflow-wrap:anywhere}
h1,h2,h3{line-height:1.45}article h1{font-size:29px}article h2{font-size:23px;margin:42px 0 14px;border-top:1px solid var(--line);padding-top:24px}p{margin:15px 0}table{border-collapse:collapse;width:100%;font-size:15px;margin:20px 0;display:block;overflow-x:auto}th,td{padding:11px 13px;border:1px solid var(--line);vertical-align:top;white-space:normal}th{background:#edf3f7;text-align:left}pre{overflow:auto;background:#edf3f7;padding:16px;font-size:13px;line-height:1.6}code{background:#eef3f6;padding:2px 4px;border-radius:3px}pre code{padding:0}li{margin:6px 0}.figure{background:#fff;border:1px solid var(--line);border-radius:12px;padding:18px;margin:24px 0}.figure img{display:block;width:100%;height:auto}.figure p{color:var(--sub);font-size:15px}.figure h2{font-size:23px}.toc{columns:2;gap:25px;font-size:15px}.toc a{display:block;margin:7px 0}footer{border-top:1px solid var(--line);font-size:14px;color:var(--sub);padding-top:24px}button{font:inherit;padding:8px 16px;border:1px solid #bdcdd7;background:white;color:var(--blue);border-radius:7px;cursor:pointer}button[aria-pressed=true]{background:#0a6283;color:#fff}.hidden{display:none}@media(max-width:780px){.cards{grid-template-columns:repeat(2,1fr)}.toc{columns:1}main{padding:0 12px 30px}header{padding:26px 20px}body{font-size:16px}}
'''


def inline(text):
    text=html.escape(text)
    text=re.sub(r'`([^`]+)`',r'<code>\1</code>',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',text)
    text=re.sub(r'\[([^\]]+)\]\(([^\s)]+)\)',r'<a href="\2">\1</a>',text)
    return text


def markdown(text):
    # Minimal renderer for the controlled Markdown produced by the report builder.
    lines=text.splitlines();out=[];paragraph=[];i=0
    def flush():
        if paragraph:out.append('<p>'+inline(' '.join(paragraph))+'</p>');paragraph.clear()
    while i<len(lines):
        line=lines[i]
        if line.startswith('```'):
            flush();i+=1;code=[]
            while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
            out.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>')
        elif line.startswith('#'):
            flush();m=re.match(r'(#{1,6})\s+(.*)',line)
            if m:out.append(f'<h{len(m[1])}>'+inline(m[2])+f'</h{len(m[1])}>')
        elif line.startswith('|'):
            flush();rows=[]
            while i<len(lines) and lines[i].startswith('|'):
                cells=lines[i].strip().strip('|').split('|')
                if not all(re.match(r'^\s*:?-+:?\s*$',c) for c in cells):rows.append(cells)
                i+=1
            out.append('<table>'+''.join('<tr>'+''.join(f'<{"th" if n==0 else "td"}>'+inline(c.strip())+f'</{"th" if n==0 else "td"}>' for c in row)+'</tr>' for n,row in enumerate(rows))+'</table>');i-=1
        elif re.match(r'^\s*(?:[-*]|\d+\.)\s',line):
            flush();items=[]
            while i<len(lines) and re.match(r'^\s*(?:[-*]|\d+\.)\s',lines[i]):
                items.append(re.sub(r'^\s*(?:[-*]|\d+\.)\s+','',lines[i]));i+=1
            out.append('<ul>'+''.join('<li>'+inline(v)+'</li>' for v in items)+'</ul>');i-=1
        elif not line.strip():flush()
        else:paragraph.append(line.strip())
        i+=1
    flush();return '\n'.join(out)


def head(title):
    return '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="P3L-R4: 768个新对象上的Maxwell材料风险、完整条件输运和数值物理计算实测"><title>'+html.escape(title)+'</title><style>'+STYLE+'</style>'


NAV='<nav><a href="report.html">中文报告</a><a href="gallery.html">实际数组图库</a><a href="RESEARCH_REPORT_ZH.md">匿名 Markdown</a><a href="GO_NO_GO.md">五项科研裁决</a><a href="EVIDENCE_PACKAGE.md">证据与源码</a><a href="../">研究目录</a></nav>'


def figure(file,title,caption,group='all'):
    p=S/file
    if not p.exists():raise FileNotFoundError(p)
    return f'<section class="figure {group}" data-group="{group}"><h2>{html.escape(title)}</h2><a href="{file}"><img loading="lazy" src="{file}" alt="{html.escape(title)}"></a><p>{html.escape(caption)}</p></section>'


def main():
    metrics=json.loads((S/'RESULTS_AND_METRICS.json').read_text())
    cards='<div class="cards"><div class="card"><b>768</b><span>两批全新独立测试 parent</span></div><div class="card"><b>13.10% / 11.36%</b><span>相对 R2 CNN 的材料 MSE 改善</span></div><div class="card"><b>66.17%</b><span>固定核数值求值 wall 节省</span></div><div class="card"><b>未合格</b><span>完整后验的 Maxwell 概率校准</span></div></div>'
    intro='<header>'+NAV+'<h1>P3L-R4 · 材料风险与完整条件概率的独立实测</h1><p>原 P1 · 全66维条件输运 · 相阈值矩 · 非线性物理代表 · 决定驱动 Maxwell 计算。真实实验完成，正面与负面结果分别保留。</p></header>'
    note='<p class="note">成像图是监督材料风险动作，不能自动称为 posterior mean。显式 Flow 的 proper score 有相对改善，生成样本仍严重违反原噪声尺度的必要校验。所有新模型在读取这两批 GT 前冻结；现在标签均已 EXPOSED。</p>'
    selected=figure('figures/final_risk/two_cohort_risk_bars_display_v2.png','两批独立材料风险','逐parent平均三个训练seed的loss，不是集成图；固定Gaussian provider几乎追平Spline，因此没有Flow特有风险增益。')
    links=[('MATERIAL_RISK_RESULTS.md','材料风险与强控制'),('CONDITIONAL_FLOW_RESULTS.md','Full-66 条件概率'),('PHASE_MOMENT_RESULTS.md','相阈值条件矩'),('PROBABILITY_CALIBRATION_AUDIT.md','独立概率校准'),('MAXWELL_REUSE_RESULTS.md','有限反馈与决定区间'),('LEARNED_JACOBIAN_RESULTS.md','Learned Jacobian 实测'),('NONHOLONOMIC_GEOMETRY_RESULTS.md','动态非完整几何'),('LOCKED_GENERALIZATION.md','冻结与泛化审计'),('RESOURCE_SUMMARY.md','完整资源语义'),('FAILURE_ANALYSIS.md','失败与 NOT_RUN'),('GPT_PRO_HANDOFF/START_HERE.md','下一轮 Pro 研究接口'),('RESULTS_AND_METRICS.json','机器可读总指标')]
    links.extend([('TRAINING_AND_COST_AUDIT.md','训练、同输入对照与费用'),('THEORY_CALIBRATION_ERROR_BOUND.md','标准校准误差下界'),('PUBLIC_REPLAY.md','可移植离线数值重放')])
    toc='<div class="toc">'+''.join(f'<a href="{p}">{html.escape(n)}</a>' for p,n in links)+'</div>'
    report=head('P3L-R4 完整中文研究报告')+intro+'<main>'+cards+note+selected+toc+'<article class="report">'+markdown((S/'RESEARCH_REPORT_ZH.md').read_text())+'</article><footer>物理/概率资格与网页发布分开。所有图来自保存的真实数组；原始模型、拒绝和失败通过现有 B1 Release 提供。<a href="EVIDENCE_PACKAGE.md">证据包</a></footer></main></html>'
    (S/'report.html').write_text(report);(S/'index.html').write_text(report)
    entries=[]
    for i in range(4):
        entries.append((f'figures/final_risk/gallery_R4L1_{i:05d}.png',f'固定前四对象 R4L1_{i:05d}：真值、独立物理点与风险动作','固定seed1示例，不按GT挑图。全部材料使用相同χ色标；6,418仅为局部head参数，约0.68M全局provider及物理基础图成本另外列出。','risk'))
    entries.extend([
        ('figures/final_risk/two_cohort_risk_bars_display_v2.png','两批新对象的风险','独立parent为统计单位，保留强普通CNN、Gaussian和空间打乱控制。','risk'),
        ('figures/final_risk/anchor_transition_gain_decomposition_display_v2.png','相区风险贡献','过渡相定义来自基础图；约99.9%增益在该区是描述性分解，不是新物理定理。','risk'),
        ('figures/final_probability/proper_NLL_status.png','完整显式概率评分','Spline和三角模型有真实properNLL；CFM数值NLL未合格项不以零填充。','probability'),
        ('figures/final_probability/material_proper_energy_CFM_controls.png','有界物理几何的 CFM 样本评分','几何权重与同谱随机控制相比较；energy score改善不是FULL后验认证。','probability'),
        ('figures/final_probability/fixed0_FullPhi_Gamma256_and_PIT.png','完整生成样本的 Maxwell 校准失败','每对象固定index0，期望Phi均值256；实际大一至两个数量级。','probability'),
        ('figures/final_probability/paired_VI_NLL_Phi_entropy.png','物理训练与概率覆盖的冲突','真实Full梯度能降低response residual，却损害NLL与分布覆盖；成对监督控制保留。','probability'),
        ('figures/final_probability/R4L1_00000_Re_fixed_samples.png','固定生成样本：实部','真实模型样本，不以最接近GT者评分；不宣称这些是已校准posterior samples。','probability'),
        ('figures/final_probability/R4L1_00000_Im_fixed_samples.png','固定生成样本：虚部','与实部同一固定样本；原被动材料和decoder支持不变。','probability'),
        ('figures/final_physics/fig01_locked_chain_cost.png','相同核、相同路径的数值计算费用','安静串行8对象重放；含setup与reduced RHS。没有以较差混合换取求值节省。','physics'),
        ('figures/final_physics/fig02_geometry_gates.png','非完整路径与直接方向的距离门','保留全部计划分母，原teacher1%门没有被事后2%适用门改写。','geometry'),
        ('figures/final_physics/fig03_risk_fit_tradeoff.png','风险图与合法物理代表','固定前64对象，所有原投影预算停止和fallback保留；箭头不是收敛轨迹。','physics'),
        ('figures/final_physics/fig04_fullJ_TR_diagnostics.png','额外 Full-J 控制诊断','GT开放前冻结；data-only更拟合但材料风险更差，prior-aware point接近原baseline。','physics'),
    ])
    buttons='<div id="filters">'+''.join(f'<button data-filter="{g}" aria-pressed="{"true" if g=="all" else "false"}">{name}</button> ' for g,name in [('all','全部'),('risk','材料风险'),('probability','概率'),('physics','物理与成本'),('geometry','几何')])+'</div>'
    gallery=head('P3L-R4 实际数组与实验图库')+intro+'<main>'+note+buttons+''.join(figure(*e) for e in entries)+'''<script>document.querySelectorAll('[data-filter]').forEach(b=>b.onclick=()=>{document.querySelectorAll('[data-filter]').forEach(a=>a.setAttribute('aria-pressed',a===b));document.querySelectorAll('[data-group]').forEach(s=>s.classList.toggle('hidden',b.dataset.filter!=='all'&&s.dataset.group!==b.dataset.filter));});</script></main></html>'''
    (S/'gallery.html').write_text(gallery)
    print(json.dumps(dict(status='STATIC_REPORT_AND_REAL_ARRAY_GALLERY_RENDERED',figure_count=len(entries),RHS=0),ensure_ascii=False))


if __name__=='__main__':main()

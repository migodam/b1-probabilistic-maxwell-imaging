"""Render the post-score public page. Reads scores and parent scientific decision only."""
import html, json, pathlib

P=pathlib.Path(__file__).resolve().parent
S=json.loads((P/'RESULTS_SUMMARY.json').read_text())
D=json.loads((P/'SCIENTIFIC_DECISION.json').read_text())
E=html.escape
def avg(rows,key):
    values=[r[key] for r in rows if r.get(key) is not None]
    return sum(values)/len(values) if values else None
def fmt(value):return '—' if value is None else f'{value:.4g}'
def result_table(backend):
    text='<div class="table"><table><thead><tr><th>A</th><th>噪声</th><th>方法</th><th>有限/计划</th><th>NRMSE 均值</th><th>IoU 均值</th><th>数据门</th><th>材料门</th><th>联合门</th><th>墙钟均值(s)</th></tr></thead><tbody>'
    for a in (.4,2.5,12.):
        for noise in ('clean','noisy'):
            for method in ('native_trf','reciprocal_trf','data_lbfgs','globaldual'):
                rows=[r for r in S['rows'] if r['backend']==backend and r['amplitude']==a and r['noise']==noise and r['method']==method]
                finite=[r for r in rows if r['finite_endpoint']]
                if not rows:continue
                values=[f'{a:g}',noise,method,f'{len(finite)}/{len(rows)}',fmt(avg(finite,'NRMSE')),fmt(avg(finite,'IoU_Otsu'))]
                values += [f"{sum(r[gate]=='PASS' for r in rows)}/{len(rows)}" for gate in ('data_gate','material_gate','joint_gate')]
                values += [fmt(avg(rows,'wall_seconds'))]
                text+='<tr>'+''.join('<td>'+E(v)+'</td>' for v in values)+'</tr>'
    return text+'</tbody></table></div>'

docs=[('RESEARCH_REPORT_ZH.md','中文完整报告'),('THEORETICAL_EXPECTATION.md','低幅度与 Gram 几何'),
      ('ALGORITHM.md','四种算法'),('BLIND_PROTOCOL.md','盲测与归因'),('DIGIT_AND_SHAPE_GENERATION.md','材料生成'),
      ('REPRODUCE.md','复现命令'),('audit/FINAL_SCIENTIFIC_REVIEW.md','独立科学审查'),
      ('RESULTS_SUMMARY.json','逐案例结果'),('COST_LEDGER.csv','成本总账'),('analysis/PAIRED_BY_BACKEND_AMPLITUDE_NOISE.csv','分层配对结果')]
source='https://github.com/migodam/b1-probabilistic-maxwell-imaging/blob/main/research/research_z_alg_r4_diversity/'
doclinks=''.join('<li><a href="'+source+p+'">'+E(label)+'</a> · <a href="'+p+'">原始文件</a></li>' for p,label in docs)
cuda=any(r['backend']=='cuda' and r['finite_endpoint'] for r in S['rows'])
template=r'''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>R4-DIVERSITY · 低中高对比度 Maxwell 成像</title>
<style>*{box-sizing:border-box}body{margin:0;background:#f3f6f5;color:#19342f;font:17px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif}header{background:#164e45;color:white;padding:32px max(24px,calc((100vw - 1260px)/2))}header small{color:#b6d7cd;letter-spacing:.08em}header h1{max-width:1100px;font-size:34px;line-height:1.4}main{max-width:1360px;margin:24px auto;background:white;padding:32px 42px;border:1px solid #dae5df;border-radius:12px}h2{border-top:1px solid #dbe5df;margin-top:38px;padding-top:28px;font-size:25px}h3{font-size:20px}a{color:#176c5e;text-underline-offset:3px}header a{color:white}nav{display:flex;gap:24px;flex-wrap:wrap;font-size:14px}p{max-width:1080px}.status{background:#edf5f1;border-left:4px solid #348170;padding:14px 20px}.controls{display:flex;align-items:center;gap:14px;flex-wrap:wrap;background:#f1f6f3;padding:16px;border-radius:8px}select{font:inherit;padding:6px;background:white;border:1px solid #a8c6b8;border-radius:5px}figure{margin:22px 0}figure img{width:100%;height:auto}figcaption{font-size:14px;color:#566b63}table{border-collapse:collapse;font-size:13px;white-space:nowrap}th,td{border:1px solid #d9e4dd;padding:9px 12px;text-align:left}th{background:#ecf4ef}.table{max-width:100%;overflow:auto;margin:20px 0}details{margin:20px 0}summary{cursor:pointer;font-weight:650}code{background:#eff5f1;padding:2px 5px}li{margin:6px 0}footer{border-top:1px solid #dbe5df;padding-top:24px;color:#5c7168;font-size:14px}.note{font-size:14px;color:#5a7066}@media(max-width:760px){body{font-size:16px}main{margin:10px;padding:22px 16px}header h1{font-size:27px}header{padding:24px}h2{font-size:22px}}</style></head><body>
<header><small>FULL NONLINEAR SCALAR-TM VIE · 24×24 · 12 SHAPES · CONTRAST-PAIRED BLIND STUDY</small><h1>__HEADLINE__</h1><nav><a href="#answer">科学结论</a><a href="#gallery">全部图像</a><a href="#results">分层结果</a><a href="#cost">成本与归因</a><a href="#files">代码与复现</a><a href="../research_z_alg_r4_basin/">原 R4-BASIN</a></nav></header>
<main><section id="answer"><div class="status">__STATUS__</div><p>__SUMMARY__</p><p>本轮以相同形状、相同背景和配对噪声，改变笔画材料增量 A=0.4、2.5、12。CPU 覆盖数字 1–9、环、双岛、弯条，共 72 份观测；CUDA 的预定 18 份子集来自相同的数字 2、7、环。72 个观测属于 12 个形状簇，GPU 不增加独立对象数。</p><p>全部 576 个实材料未知，统一初值 0.1、box [0.02,20]；12 源×48 复接收，k=2π，保留完整多次散射。8→12→24 只改变材料参数化，所有阶段都用 576 阶物理系统。</p></section>
<section id="gallery"><h2>九数字与三种非数字形状：完整保留成功和失败</h2><div class="controls"><label>后端 <select id="backend"><option value="cpu">i9 P 核 · 全 12 形状</option>__CUDA_OPTION__</select></label><label>材料增量 <select id="amplitude"><option value="0.4">低幅度 A=0.4</option><option value="2.5">中幅度 A=2.5</option><option value="12">高幅度 A=12</option></select></label><label>观测 <select id="noise"><option value="clean">无噪声</option><option value="noisy">相对复噪声 10⁻³</option></select></label></div><p class="note">每个幅度层内共用物理材料色阶，所有算法相同。粉色表示超过该显示上限，完整合法材料 box 仍到 20。NRMSE=材料 RMSE/A；IoU 使用重建自身的 Otsu 阈值。图中的 GT 只在预测关闭后用于评分。</p>
<figure><a id="digitsLink" href="figures/cpu_digit_A0.4_clean.png"><img id="digits" src="figures/cpu_digit_A0.4_clean.png" alt="九数字GT与四种算法的在线top1" loading="lazy"></a><figcaption>GT、native TRF、reciprocal TRF、data-L-BFGS、GlobalDual+原目标 correction。点击图查看原图。</figcaption></figure>
<figure><a id="shapesLink" href="figures/cpu_shape_A0.4_clean.png"><img id="shapes" src="figures/cpu_shape_A0.4_clean.png" alt="环、双岛、弯条GT与四种算法" loading="lazy"></a><figcaption>非数字形状采用相同物理采集、初值、先验与预算。</figcaption></figure></section>
<section id="results"><h2>按幅度、噪声、算法分开的结果</h2><p>每个方法上限 45 秒、2000 次完整 LU、60000 个完整系统复 RHS。共同正则与多尺度安排固定。实际消耗和停止原因保留；预算耗尽的有限图像不等于收敛。数据门、材料门和联合门分别计数。</p><details open><summary>i9 P 核：12 个配对形状</summary>__CPU_TABLE__</details><details><summary>4060：预定三形状配对子集（完成状态见表）</summary>__CUDA_TABLE__</details>
<figure><img src="figures/contrast_quality.png" alt="三档幅度的全部材料误差" loading="lazy"><figcaption>每个点保留实际对象结果；材料误差按真实增量归一化。不同幅度/噪声的同形状结果相关。</figcaption></figure></section>
<section id="cost"><h2>GlobalDual 的增量与真实费用</h2><p>直接比较 GlobalDual 与匹配的 data-L-BFGS，才能区分材料相关 Gram 目标与普通伴随优化的贡献。仅比 Full-J 少了接收端 RHS，不足以证明 GlobalDual 有额外收益。</p><figure><img src="figures/cost_vs_image_quality.png" alt="实测墙钟与在线top1材料误差" loading="lazy"><figcaption>实际端点墙钟与材料误差；CPU/CUDA 分开，不以负载不同的总时间宣称硬件加速。</figcaption></figure><figure><img src="figures/paired_wins_and_costs.png" alt="配对优胜平局与计算成本" loading="lazy"><figcaption>预声明 practical tie 是实际差异阈值，不是统计等价检验。分层连续差与形状簇汇总均保存。</figcaption></figure>
<p>低幅度并不自动让 GlobalDual 更优。在稳定的小扰动区域，它趋于固定的源空间加权残差；权重改变优化曲率，不增加测量信息。根附近、无正则且噪声模型正确时，普通最小二乘还有标准局部线性无偏方差最优性。GlobalDual 的价值若存在，必须来自本批实际有限预算的候选路径或条件改善，而不能从 Gram 公式直接推出。</p><p>本批A=0.4时，Born/完整响应差仍约21%–29%；它是低材料增量层，并非认证Born区。中幅度环与弯条的大误差在共用Data粗层阶段已出现，不能单独归罪最后的Gram阶段。</p><figure><img src="figures/coarse_path_attribution.png" alt="封存阶段候选的误差与残差" loading="lazy"><figcaption>仅连接保存的阶段最优端点；GT距离在预测关闭后计算，不是在线触发信号。</figcaption></figure><p>__LIMITATIONS__</p></section>
<section id="files"><h2>全部证据与复现入口</h2><ul>__DOCS__</ul><p><a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-z-alg-r4-diversity-20261030">原始收据与完整复现归档</a> · <a href="https://github.com/migodam/b1-probabilistic-maxwell-imaging/tree/main/research/research_z_alg_r4_diversity/src">算法源码</a></p><p>在线 top1 由原始完整观测目标与共同通用正则选择。best@K 只作离线 oracle 指标。当前保存的是有限优化阶段候选；没有后验概率、概率解集或独立根分支覆盖的主张。</p></section><footer>R4-DIVERSITY · 2026-10-30 · 固定离散 scalar-TM Maxwell 实验，不能推广为三维矢量 Maxwell 或全局唯一/全局收敛定理。保留原 R2/R3/R4 冻结历史。</footer></main>
<script>function gallery(){const b=document.getElementById('backend').value,a=document.getElementById('amplitude').value,n=document.getElementById('noise').value;for(const k of ['digits','shapes']){const p='figures/'+b+'_'+(k==='digits'?'digit':'shape')+'_A'+a+'_'+n+'.png';document.getElementById(k).src=p;document.getElementById(k+'Link').href=p;}}for(const id of ['backend','amplitude','noise'])document.getElementById(id).addEventListener('change',gallery);</script></body></html>'''
replace={'__HEADLINE__':E(D['headline']),'__STATUS__':E(D['status']),'__SUMMARY__':E(D['summary']),
         '__LIMITATIONS__':E(D['limitations']),'__DOCS__':doclinks,'__CPU_TABLE__':result_table('cpu'),
         '__CUDA_TABLE__':result_table('cuda'),'__CUDA_OPTION__':'<option value="cuda">RTX4060 · 三形状子集</option>' if cuda else ''}
for key,value in replace.items():template=template.replace(key,value)
(P/'index.html').write_text(template,encoding='utf-8')

"""Build the static R10 evidence entry from completed score tables."""
from pathlib import Path
import csv,html
ROOT=Path(__file__).resolve().parents[1]
REPO='https://github.com/migodam/b1-probabilistic-maxwell-imaging'
WEB=REPO+'/blob/main/research/research_max_struct_r10/'
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/research_max_struct_r10/CHAT_HANDOFF.md'

def main():
    old=(ROOT.parent/'research_max_struct_r9/index.html').read_text()
    style=old.split('<style>',1)[1].split('</style>',1)[0]
    labels={'cie_single':'CIE','cie_dual':'Dual-CIE','a_plus':'A+ global','a_plus_sourcewise':'A+ sourcewise','matrix_green19':'Green19','matrix_twofold':'Two-fold30','global_dual':'Global-Dual（内部消融）'}
    data=list(csv.DictReader((ROOT/'delegated/twod_score/scored/REGIME_METHOD_SUMMARY.csv').open()))
    table=[]
    for key,label in labels.items():
        r=next(x for x in data if x['regime']=='all' and x['method']==key)
        table.append('<tr><td>'+label+'</td><td>'+r['material_success_count']+'/24</td><td>'+format(float(r['amplitude_NRMSE_mean']),'.5f')+'</td><td>'+format(float(r['external_full_child_CPU_seconds_mean']),'.2f')+'</td></tr>')
    links=''.join('<a href="'+WEB+name+'">'+label+'</a>' for name,label in [('RESEARCH_REPORT_ZH.md','中文完整报告'),('THEORY_DERIVATIONS.md','公式与证明'),('COST_LEDGER.md','费用与失败'),('PAPER_DRAFT.md','英文论文草稿'),('REPRODUCE.md','复现')])
    conditional=(ROOT/'PUBLIC_CONDITIONAL_SUMMARY.html').read_text() if (ROOT/'PUBLIC_CONDITIONAL_SUMMARY.html').exists() else '<p>最后四对象补偿确认的独立结果见完整报告；不与24对象主比较混合。</p>'
    page='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>MAX-STRUCT R10 · Green分母、能量与CIE补偿</title><style>'''+style+'''pre{white-space:pre-wrap;padding:20px;background:#e5eeea;overflow-wrap:anywhere}.small{font-size:14px}</style></head><body><main>
<div class="eyebrow">MAX-STRUCT / R10 · REPRODUCIBLE EVIDENCE</div>
<h1>同一套Maxwell数据，<br>分母改变了哪条搜索路径？</h1>
<p class="lead">新高对比环被救回，开缝壳与桥却出现反向失败。完整比较CIE、Dual-CIE、A+与Green矩阵分母，把表示容量、材料路径、噪声歧义和模型误差分开。</p>
<div class="links">'''+links+'<a class="primary" href="'+RAW+'">Chat读取入口</a><a href="'+REPO+'''/releases/tag/research-max-struct-r10-20261107">完整原始证据</a></div>
<div class="stats"><div class="stat"><b>24 × 7</b><span>独立新二维材料与七种方法</span></div><div class="stat"><b>8 家族</b><span>low / high / mixed 全保留</span></div><div class="stat"><b>4 / 18</b><span>真三维独立材料 / 主反演</span></div><div class="stat"><b>4.41×</b><span>同v08任务CUDA后端实测提速</span></div></div>
<section><h2>先看全部结果，再解释正例</h2><div class="tablewrap"><table><thead><tr><th>方法</th><th>材料描述门</th><th>平均NRMSE</th><th>平均完整CPU秒</th></tr></thead><tbody>'''+''.join(table)+'''</tbody></table></div>
<p class="muted">168条主预测全部封存后评分。不是168个独立材料。NRMSE与IoU联合门只是预声明描述；普通对象的部分成功数变化来自很小的阈值跨线。完整CPU来自Windows子进程，而不是只给材料内循环计时。</p>
<figure><a href="FIGURES/twod_comparison/MATERIALS_high.png"><img loading="lazy" src="FIGURES/twod_comparison/MATERIALS_high.png" alt="八个高对比材料的GT与全部七种算法，包含所有失败"></a><figcaption>环0001：CIE .423 → Green19 .099。桥0022：CIE .097 → Green19 .337。开缝壳0004：CIE .111，而A类约.515。新分母有真实增量，也有实质损害。</figcaption></figure>
<div class="links"><a href="FIGURES/twod_comparison/MATERIALS_low.png">全部8例low</a><a href="FIGURES/twod_comparison/MATERIALS_mixed.png">全部8例mixed</a><a href="FIGURES/twod_comparison/MATERIALS_high.pdf">高对比PDF</a><a href="FIGURES/twod_comparison/COST_PERFORMANCE.png">费用—性能图</a></div></section>
<section><h2>从局部一次传播，改为低秩矩阵分母</h2><pre>R_K = β [I + βD_m − D_mK]⁻¹ [J − D_m(Einc + GJ)]
K = PGP · 全部576个材料未知量保留</pre>
<p>原A+只在标量分母处理自项，再在局部块传播一次。新方法用19或30维小系统处理选定的全域耦合，current、FFT/APCS、prior及原Maxwell corrector保持。对精确盘脉冲模型可证明分母定性可逆；这不保证良好条件数、全局凸性或材料唯一性。</p>
<div class="cards"><article class="card"><h3>不是多给几个Fourier模式</h3><p>新集七臂共享实际FFT/APCS族。末容量的GT current不可表示能量仅占.000225%–.014336%，却有巨大的材料路径分叉。小表示尾不能解释或保证最终材料恢复。</p></article><article class="card"><h3>更好的G逼近也不够</h3><p>去掉左投影得到GP，保留了强输入模式的完整输出。三例实跑仍只保住环，未修复壳，桥更差。没有继续调rank或β来遮掩反例。</p></article></div>
<figure><img loading="lazy" src="FIGURES/twod_comparison/ENERGY_FIXED_TAU_POINT1.png" alt="全部新对象的材料和物理current误差能量分解"><figcaption>绝对误差与百分比、辅助current与物理current、有符号cross分别保留。这些范数不是电磁储能，Green模式也不是材料Jacobian可观测性。</figcaption></figure></section>
<section><h2>只在数据明显未拟合时补偿CIE</h2>'''+conditional+'''<p class="muted">该规则是独立冻结后的四对象小确认，不能回填为原24比较的预注册算法。普通Dual-CIE同触发重试是额外候选/成本控制。未触发的fallback明确NOT_RUN。</p></section>
<section><h2>三维：算法差别小，模型与观测限制更大</h2>
<p>12³个实材料未知量、三分量dyadic Green、60外部偶极源；两对象另用18³生成数据。matched ε9实心块的χ误差约.055–.061，三臂都恢复较好；包覆体与空心体在近噪声拟合下仍有显著材料误差，A+没有独立而明确的恢复优势。</p>
<figure><img loading="lazy" src="FIGURES/vector_comparison/v09_matched12_slices.png" alt="包覆体GT与CIE、Global-Dual及A+三正交材料切片"><figcaption>包覆体GT与CIE估计均在同一实际观测的CN99噪声球内，材料χ误差却约49%。这是有限噪声双材料实例，不是精确多解定理。</figcaption></figure>
<figure><img loading="lazy" src="FIGURES/vector_comparison/v04_matched12_volume.png" alt="高低介电双体的三维GT与三种算法重建"><figcaption>全18主终点保留。双体另有固定80步物理续跑，用来诊断预算；原18不被替换。主三维目标为继承的源归一化Φ3，真实σ白化资格另计。</figcaption></figure>
<p>实心体和包覆体的coarse/fine clean端口差为24.80%和6.41%，大于3%测量噪声。原coarse材料在fine观测的噪声球之外；精确返回它与拟合fine数据不能同时满足。不能将这种模型差异包装成算法失败或算法优势。</p></section>
<section><h2>实际用上i9与4060，也记录慢路与失败</h2><p>二维主168任务使用i9 P核有限队列，完整child CPU合计7.91小时。三维4060改用库提供的retained dense LU，全60源及伴随复用分解；同v08 CIE任务882→200秒，端点残差接近。旧慢运行、工程启动故障、候选费用、评分及打包均入账，GPU时间戳不冒充活跃kernel时间。</p></section>
<footer>论文主张由完整结果决定。本轮保留有价值的条件正例，也公开对应反例和费用；不要求全例获胜，不宣称普遍凸化或连续三维认证。R1–R9冻结证据保持。<br><a href="../research_max_struct_r9/">← R9</a> · <a href="'''+WEB+'''RESEARCH_REPORT_ZH.md">完整中文结论</a> · <a href="'''+RAW+'''">匿名raw交接</a></footer>
</main></body></html>'''
    (ROOT/'index.html').write_text(page)

if __name__=='__main__':main()

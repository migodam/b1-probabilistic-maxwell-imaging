"""Publish presentation of sealed R3 results. No scientific/physics execution."""
from pathlib import Path
import importlib.util,html,json,re
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
spec=importlib.util.spec_from_file_location('r2site',ROOT.parent/'research_w_gq_r2/tools/publish_site.py')
base=importlib.util.module_from_spec(spec);spec.loader.exec_module(base)
REL='research/research_w_gq_r3_joint/'
GITHUB='https://github.com/migodam/b1-probabilistic-maxwell-imaging/'
BLOB=GITHUB+'blob/main/'+REL
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/'+REL
PAGE='https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+REL
TAG='research-w-gq-r3-ji-20261030'
RELEASE=GITHUB+'releases/tag/'+TAG
ASSET=GITHUB+'releases/download/'+TAG+'/RESEARCH_W_GQ_R3_OFFLINE_AUDIT.zip'
base.ROOT=ROOT;base.BLOB=BLOB
MATH=r"""<script>window.MathJax={tex:{inlineMath:[['$','$'],['\\(','\\)']],displayMath:[['$$','$$'],['\\[','\\]']]},options:{skipHtmlTags:['script','noscript','style','textarea','pre','code']}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>"""
def render(md):
    formulas=[]
    def keep(m):
        formulas.append(html.escape(m[0]));return 'MATHPLACEHOLDER'+str(len(formulas)-1)+'END'
    md=re.sub(r'\$\$[\s\S]*?\$\$|(?<!\$)\$(?!\$)[^\n$]*?\$',keep,md)
    s=base.render(md)
    for i,f in enumerate(formulas):
        k='MATHPLACEHOLDER'+str(i)+'END'
        s=s.replace('<p>'+k+'</p>','<div class="math">'+f+'</div>').replace(k,f)
    return s
def page(title,body):
    nav='<nav><a href="../../">B1 首页</a><a href="index.html">R3 概览</a><a href="report.html">中文报告</a><a href="formulas.html">41组公式</a><a href="gallery.html">真实图像</a><a href="digits.html">数字8</a><a href="handoff.html">Pro交接</a><a href="evidence.html">下载与复现</a></nav>'
    note='<p class="muted">发布记录：2026-10-30，按用户“发布”请求公开。科研于2026-10-29 18:25 UTC封账；原报告中的“未发布/未授权”描述当时状态，原始裁决与归档不回写。发布没有新增物理计算。</p>'
    foot='<footer>二维 scalar-TM · 受限未知宏格墙＋连续椭圆目标 · 有限开发样本。<br><a href="'+BLOB+'GO_NO_GO.json">分项裁决</a> · <a href="'+BLOB+'results/COST_AGGREGATION.json">完整费用</a> · <a href="'+BLOB+'FAILURE_LEDGER.jsonl">失败账本</a> · <a href="'+RAW+'REPORT_ZH.md">匿名原始 Markdown</a> · <a href="'+RELEASE+'">完整证据 Release</a></footer>'
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+html.escape(title)+'</title><style>'+base.STYLE+'</style>'+MATH+'</head><body><main>'+nav+note+body+foot+'</main></body></html>'
def fig(path,caption):
    s='<figure><a href="'+path+'.png"><img loading="lazy" src="'+path+'.png" alt="'+html.escape(caption,quote=True)+'"></a><figcaption>'+html.escape(caption)+'</figcaption>'
    if (ROOT/(path+'.pdf')).exists():s+='<a href="'+path+'.pdf">PDF 原图</a>'
    return s+'</figure>'
figs=[
('figures/R0_TARGET_MATERIAL_OVERVIEW','三种方法、三个独立世界的真实目标和重建；不是只展示最好的分支。'),
('figures/R0_WALL_EPS_OVERVIEW','未知墙实介电率图；六个公开宏格材料未知，不是任意未知几何。'),
('figures/R0_WALL_LOSS_OVERVIEW','墙体损耗重建；目标定位成功不意味着复材料恢复成功。'),
('figures/R0_FULLGN_CONTOUR_ZOOM','同尺度连续轮廓与评分网格细化；形状偏差不能归结为展示分辨率。'),
('figures/R0_QUALITY_COST','配对首噪声质量与实际费用；adaptive没有超过普通Full-GN。'),
('figures/R0_ALL_BRANCHES','全部多起点及测量选择；真值误差仅离线展示。'),
('figures/R0_OBJECTIVE_TRAJECTORIES','实际完整Maxwell优化轨迹与停止状态。'),
('figures/R0_MATCHED_CONTROLS','prior-only、先估墙再固定及known-wall诊断对照。'),
('research/delegated/sources/CHART_MECHANISM_RUN/FIGURES/public_2x2_mechanism','公共2×2机制：目标自身变化也能使参考材料forcing chart失配。'),
('exploration/digits/truth_vs_reconstruction','数字8：110个未知材料变量，孔洞未恢复，两个起点均未收敛。')]
def main():
    assert (ROOT/'RUN_CLOSED.json').exists()
    body='<div class="eyebrow">RESEARCH W-GQ-R3-JI · 联合全波成像 · 2026-10-30</div><h1>墙体和目标都未知，<br>这次真正重建出了什么？</h1><p class="lede">从复测量出发，连续估计墙体材料与目标位置、尺寸、方向及复材料，输出真实二维图像。普通联合 Full-GN 得到一个有效盲测样例；自适应 chart 的独立优势尚未成立。</p><div class="status"><strong>可运行联合成像：PASS（受限模型）。盲测准确性：LIMITED PASS。Adaptive增量：NO-GO。数字8：未恢复清楚。</strong></div>'
    metrics=[('0.03234 λ₀','Full-GN 三parent首噪声位置RMSE'),('1 / 3','Full-GN通过位置／形状／实ε三门的parent'),('10 / 18','完成的主比较槽位；控制9/9'),('41','从Green到材料更新的推导公式组')]
    body+='<div class="metrics">'+''.join('<div class="metric"><strong>'+v+'</strong><span>'+t+'</span></div>' for v,t in metrics)+'</div>'
    body+='<p><a href="report.html">完整中文报告 →</a>　<a href="formulas.html">推导公式集 →</a>　<a href="evidence.html">193 MB完整代码和数据 →</a></p>'
    body+='<section class="panel"><h2>真实未知量，真实图像</h2><p>在线只读取测量与公开配置。六个墙宏格的实介电率、损耗和七个连续目标参数共同优化；多起点按测量目标函数与冻结先验选支。三个独立parent全部展示；重复方法和噪声不算新世界。一个连外边界空气缺块世界的位置误差0.01160λ₀、IoU 0.84069、实ε相对误差0.00656，但目标损耗仍不准。</p>'+fig(*figs[0])+'</section>'
    body+='<section class="panel"><h2>不只看IoU：看尺度、边界和材料</h2><p>normal目标偏心且压扁，strong目标放大且偏移。将评价网格从401×301细化至1601×1201，IoU变化不足0.0032，不能用图像像素解释形状失败。物理网格细化另行报告；两例真值的细网格变化超过旧W的0.1白化门，未追认旧资格通过。</p>'+fig(*figs[3])+'</section>'
    body+='<section class="panel"><h2>固定chart为何失效？新算法是否更好？</h2><p>仅改变墙体时，目标Schur相对误差约0.00630；墙体不变而插入强目标后升至4.396。参考载荷子空间遗漏的材料切向不仅来自墙变化，也来自目标局部场与完整反馈。这个例子不否定原生SOM/TSOM。adaptive恢复部分固定chart质量，却用约1.47倍RHS、1.50倍墙钟，没有胜过普通Full-GN。</p>'+fig(*figs[8])+'</section>'
    body+='<section class="panel"><h2>更复杂的数字8：拟合数据仍不能保证看清物体</h2><p>7×7目标复材料格与未知墙共110变量，不使用数字模板或真实笔画支撑。白化测量RMS 0.9801，复contrast相对误差0.7793；两孔未恢复，两个起点均达到迭代上限。1–9中只对8执行了物理反演，其余仅有预定义图样。clean均值距离8.94，不能写成d≤1物理歧义证明。</p>'+fig(*figs[9])+'</section>'
    body+='<section class="panel"><h2>复现与科学边界</h2><p>可移动主入口已实际执行一次独立测量重建；全新依赖安装未运行。完整ZIP含7552项，保留原始复响应、分支、材料数组、失败和账本。数字8未来独立部署封装尚未完成，不冒充与主入口相同的复现状态。A17X已有理论、标准Schur/SOM/TSOM/ROM明确归属，论文级新算法增量未成立。</p><p>总预约12,440 FE / 187,732 complex RHS；两个早期CPU/RSS记录缺失。物理已封账，发布不增加求解。</p></section>'
    (ROOT/'index.html').write_text(page('Research W-GQ-R3-JI｜未知墙与未知目标联合成像',body))
    docs=[('report.html','REPORT_ZH.md','完整中文报告'),('formulas.html','FORMULA_COMPENDIUM_ZH.md','41组推导公式'),('handoff.html','NEXT_PRO_HANDOFF.md','GPT Pro交接'),('run.html','RUN_RECONSTRUCTION.md','单命令复现'),('digits.html','exploration/digits/REPORT_ZH.md','数字8真实实验')]
    for out,source,title in docs:
        md=(ROOT/source).read_text()
        if source.startswith('exploration/digits/'):
            md=re.sub(r'(!?\[[^\]]*\]\()((?!https?:|#)[^)]+)(\))',lambda m:m[1]+'exploration/digits/'+m[2]+m[3],md)
        (ROOT/out).write_text(page('R3｜'+title,'<article>'+render(md)+'</article>'))
    (ROOT/'gallery.html').write_text(page('R3｜全部真实图表','<h1>真实成像与失败画廊</h1><p>全部来自已保存实验数组。点击图片看原尺寸；未运行数字没有重建图。</p>'+''.join('<section class="panel">'+fig(*f)+'</section>' for f in figs)))
    ev='<article><h1>源码、测量和完整证据</h1><p><strong><a href="'+ASSET+'">完整离线审计与复现ZIP：193,037,632 bytes，7552项</a></strong></p><p><a href="'+RELEASE+'">Release说明</a> · <a href="run.html">实际运行入口</a> · <a href="'+RAW+'REPORT_ZH.md">匿名Markdown</a> · <a href="'+RAW+'FORMULA_COMPENDIUM_ZH.md">匿名公式集</a></p><p>解压入口为 <code>research_w_gq_r3_joint/</code>。大型NPZ/NPY数组保留在Release；仓库保留报告、科学源码、图、配置、原始账本、收敛和运行回执。原始离线归档不修改；外部最终费用和发布记录单列。</p><h2>结果与原始记录</h2><ul>'
    for n in ['GO_NO_GO.json','EXPERIMENT_RESULTS.csv','ABLATIONS.csv','CONTRACT.json','METHODS_FROZEN.json','RUN_CLOSED.json','results/COST_AGGREGATION.json','COST_LEDGER.jsonl','PHYSICS_LEDGER.jsonl','FAILURE_LEDGER.jsonl','results/ISOLATION_AUDIT.json','results/PORTABLE_ENTRY_VERIFICATION.json','PACKAGE_COST_SIDECAR.json','PACKAGE_DIRECTORY_READBACK.json','publication/PREPARATION.json']:
        ev+='<li><a href="'+BLOB+n+'">'+n+'</a></li>'
    ev+='</ul><h2>成本与复现边界</h2><p>12,440 FE / 187,732 RHS；已知process CPU 6468.60秒，两个早期CPU/RSS缺失。数字8的276 FE / 2172 RHS独立标记，镜像账本不重复相加。打包CPU 7.056秒在外部sidecar，不回写已封存归档。发布准备、上传和HTTP核验属于交付费用，不是Maxwell实验。</p><p>主portable入口已实际执行，使用内置冻结源码和必要历史后端；全新依赖安装NOT_RUN。数字8可复查保存结果，但其原实验入口沿用已结束的绝对时限，独立新预算封装NOT_IMPLEMENTED。不要解除旧时限或覆盖原账本来制造“重跑”。</p><h2>来源与权利</h2><p>未打包版权教材原文、书页或提取文本。引用与A17X归属保留；本次公开不为既有代码或依赖擅自增加新许可。科研封账与发布授权是不同时间的记录，历史“未发布”字段不更改。</p></article>'
    (ROOT/'evidence.html').write_text(page('R3｜证据与下载',ev))
    pub={'authorization':'User: 发布','date':'2026-10-30','new_physics':0,'scientific_reports_unchanged':True,'historical_publication_fields':'Preserved as scientific close-time state; this record supersedes only publication authorization/status','release':RELEASE,'asset':ASSET,'archive_bytes':193037632,'pages':8}
    (ROOT/'publication/PUBLICATION_NOTE.json').write_text(json.dumps(pub,ensure_ascii=False,indent=2))
    mark='W_GQ_R3_PUBLICATION'
    summary='未知六宏格墙＋连续椭圆目标联合全波成像已运行。三parent首噪声Full-GN位置RMSE 0.03234λ₀；一例通过位置/形状/实ε三门。Adaptive未胜过普通Full-GN；数字8的孔洞未恢复。全部失败、41组公式和完整复现包公开。'
    entry='<!-- '+mark+' -->\n## Research W-GQ-R3-JI：未知墙与未知目标联合成像\n\n'+summary+'\n\n[中文网页]('+PAGE+') · [中文原文]('+REL+'REPORT_ZH.md) · [公式集]('+REL+'FORMULA_COMPENDIUM_ZH.md) · [完整代码和数据]('+RELEASE+')\n<!-- /'+mark+' -->\n\n'
    for n in ['README.md','START_HERE.md']:
        p=REPO/n;s=p.read_text()
        if mark not in s:p.write_text(entry+s)
    for n,link in [('index.html',REL),('research/index.html','research_w_gq_r3_joint/')]:
        p=REPO/n;s=p.read_text()
        if mark not in s:
            card='<!-- '+mark+' --><section class="card"><h2>Research W-GQ-R3-JI · 未知墙–目标联合全波成像</h2><p>'+summary+'</p><a href="'+link+'">真实成像、数字8失败、公式与完整复现包 →</a></section><!-- /'+mark+' -->'
            assert '<main>' in s
            p.write_text(s.replace('<main>','<main>'+card,1))
    print(json.dumps({'pages':8,'figures':len(figs),'new_physics':0}))
if __name__=='__main__':main()


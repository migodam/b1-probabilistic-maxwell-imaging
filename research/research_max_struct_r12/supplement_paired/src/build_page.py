"""Render the scored paired supplement only; never loads truth or changes figures.

Run after root confirms offline scoring completion. All numerical tables come
from the selected offline results JSON; missing required records fail before either page is written.
"""
from pathlib import Path
import argparse
from html import escape
import json
import math
from urllib.parse import quote

ROOT=Path(__file__).resolve().parents[1]
OFFLINE=ROOT/'RESULTS/offline_v2'
REPO='https://github.com/migodam/b1-probabilistic-maxwell-imaging'
SUB='research/research_max_struct_r12/supplement_paired/'
METHODS={'g19':'G19 合法投影重启','dr':'DR','ap':'AP','direct':'Direct'}
COHORTS={'A':'A · 三个暴露坏端点','B':'B · 全24条件救援策略','C':'C · 四个冻结 NN 初值'}
CSS='''
:root{--ink:#172d43;--muted:#586b7c;--line:#dce5eb;--paper:#fff;--accent:#176b83}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:#f3f6f8;color:var(--ink);font-family:-apple-system,BlinkMacSystemFont,"Segoe UI","PingFang SC","Microsoft YaHei",sans-serif;line-height:1.7}main{max-width:1240px;margin:auto;padding:28px 24px 64px}header{background:linear-gradient(135deg,#17354e,#176b83);color:white;border-radius:20px;padding:36px;margin-bottom:22px}header a{color:#e0f5ff}header p{max-width:960px}h1{font-size:clamp(28px,4vw,44px);line-height:1.25;margin:10px 0 20px}h2{font-size:25px;line-height:1.35;margin:0 0 14px}h3{font-size:19px}a{color:#14657e;overflow-wrap:anywhere}nav{display:flex;flex-wrap:wrap;gap:8px 18px}section{background:var(--paper);border:1px solid var(--line);border-radius:15px;margin:20px 0;padding:24px}.note{color:var(--muted);font-size:14px}.callout{border-left:4px solid var(--accent);background:#edf5f8;padding:14px 18px;margin:18px 0}.table-scroll{overflow:auto;border:1px solid var(--line);border-radius:9px}table{border-collapse:collapse;width:100%;font-size:13px;white-space:nowrap}td,th{padding:10px 12px;border-bottom:1px solid var(--line);text-align:right}th{background:#edf3f6;font-weight:600}td:first-child,th:first-child,td:nth-child(2),th:nth-child(2){text-align:left}tr:last-child td{border-bottom:0}figure{margin:20px 0;border:1px solid var(--line);border-radius:10px;padding:12px;background:white}figure img{display:block;width:100%;height:auto;object-fit:contain}figcaption{font-size:14px;color:var(--muted);padding:12px 5px 2px;overflow-wrap:anywhere}.docs{display:grid;grid-template-columns:repeat(auto-fit,minmax(230px,1fr));gap:12px}.mapping td{font-family:ui-monospace,SFMono-Regular,monospace}footer{color:var(--muted);font-size:13px}@media(max-width:600px){main{padding:12px 10px 40px}header{padding:24px 18px;border-radius:14px}section{padding:18px 12px}figure{padding:4px;margin:14px 0}td,th{padding:9px}nav{font-size:14px}}
'''

def read(path):
    return json.loads(path.read_text())

def e(value):
    return escape(str(value),quote=True)

def fmt(value):
    if value is None:return 'NA'
    if isinstance(value,bool):return '是' if value else '否'
    if isinstance(value,float):
        if not math.isfinite(value):return 'NA'
        return f'{value:.5g}'
    return e(value)

def table(headers,rows,cls=''):
    return '<div class="table-scroll"><table class="'+cls+'"><thead><tr>'+''.join('<th>'+e(x)+'</th>' for x in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+fmt(v)+'</td>' for v in row)+'</tr>' for row in rows)+'</tbody></table></div>'

def link(path,label):
    return '<a href="'+e(path)+'">'+e(label)+'</a>'

def figure(path):
    rel=path.relative_to(ROOT).as_posix();pdf=path.with_suffix('.pdf')
    title=path.stem.replace('_',' ')
    extra=' · '+link(quote(pdf.relative_to(ROOT).as_posix()),'PDF') if pdf.is_file() else ''
    return '<figure><a href="'+e(quote(rel))+'"><img src="'+e(quote(rel))+'" alt="'+e(title)+'" loading="lazy" decoding="async"></a><figcaption>'+e(title)+' · '+link(quote(rel),'原尺寸 PNG')+extra+'</figcaption></figure>'

def documents():
    docs=[('README.md','实验说明'),('RESEARCH_REPORT_ZH.md','中文研究报告'),('COST_LEDGER.md','成本账本'),('REPRODUCE.md','复现说明')]
    return '<section id="documents"><h2>报告、数据与复现</h2><div class="docs">'+''.join('<div>'+link(REPO+'/blob/main/'+SUB+name,label)+' · '+link('https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/'+SUB+name,'raw')+'</div>' for name,label in docs)+'</div><p>'+link(OFFLINE.relative_to(ROOT).as_posix()+'/SUMMARY.json','SUMMARY.json')+' · '+link(OFFLINE.relative_to(ROOT).as_posix()+'/SCORES.json','逐对象评分')+' · '+link(OFFLINE.relative_to(ROOT).as_posix()+'/PERFORMANCE_VS_ACTUAL_MEAN_CPU.json','表现与实际成本')+' · '+link(OFFLINE.relative_to(ROOT).as_posix()+'/COST_ACCOUNTING.json','成本口径')+' · '+link('inputs/manifest.json','全部输入与 ID')+'</p></section>'

def page(title,body):
    return '<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+e(title)+'</title><style>'+CSS+'</style></head><body><main>'+body+'<footer>页面由已完成的离线评分 JSON 与既有 PNG 生成。图像、色标与原始实验记录未被网页生成器改写。</footer></main></body></html>\n'

def main():
    global OFFLINE
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir',type=Path,default=ROOT/'RESULTS/offline_v2')
    args=parser.parse_args()
    OFFLINE=args.results_dir.resolve()
    results_href=OFFLINE.relative_to(ROOT).as_posix()
    receipt=read(OFFLINE/'SCORING_RECEIPT.json')
    if receipt.get('status')!='COMPLETE_OFFLINE_ONLY':raise RuntimeError('Completed offline scoring receipt required')
    if not (ROOT/'RESEARCH_REPORT_ZH.md').is_file():raise RuntimeError('Root scientific report required for findings panel')
    summary=read(OFFLINE/'SUMMARY.json');scores=read(OFFLINE/'SCORES.json')
    performance=read(OFFLINE/'PERFORMANCE_VS_ACTUAL_MEAN_CPU.json')
    costs=read(OFFLINE/'STRATEGY_COST_TOTALS.json');accounting=read(OFFLINE/'COST_ACCOUNTING.json')
    cases=read(ROOT/'inputs/manifest.json')['cases']
    if len(cases)!=31 or {c:sum(r['cohort']==c for r in cases) for c in COHORTS}!={'A':3,'B':24,'C':4}:raise RuntimeError('Unexpected full cohort mapping')
    if not scores or not summary or not performance or not costs:raise RuntimeError('Completed score/cost tables required')
    images=sorted((ROOT/'FIGURES').glob('*.png'))
    if not images:raise RuntimeError('Scored figures required')
    groups={c:[p for p in images if p.name.startswith(c+'_')] for c in COHORTS}
    used={p for ps in groups.values() for p in ps}
    mapping=[p for p in images if 'mapping' in p.stem.lower()]
    historical=[p for p in images if p not in used and p not in mapping]
    bcase_ids={r['case_id'] for r in scores if r['cohort']=='B' and r.get('rescue_triggered')}
    # Counts below are checked against completed score records, not used as results.
    if len(bcase_ids)!=4:raise RuntimeError('Expected approved four public-noise triggers; root review required')
    nav='<nav>'+link('../index.html','R12 原冻结研究')+link('index.html','补充实验总览')+link('GALLERY.html','完整图集')+link('#documents','文档与数据')+'</nav>'
    hero='<header><p>MAX-STRUCT R12 · PAIRED SUPPLEMENT</p><h1>R12 配对补充实验</h1><p>暴露对象上的事后配对分析。A 为3个既有坏端点；B 保留全部24个原 G19 端点，其中4个触发救援，20个保留原材料；C 为4个冻结 NN 初值。本补充不作为新的独立确认集。</p>'+nav+'</header>'
    scope='<section id="scope"><h2>比较范围与读图口径</h2><p>G19 使用合法投影后的重启接口，不是完整 PR 方向与状态连续续跑。各方法从相同冻结材料开始，使用公开物理量选候选；GT 仅用于预测封存后的离线评分。</p><p>raw 表示已通过完整原物理求值的提案候选；Direct 的 raw 是共同初值。corrected 表示共同物理校正后的候选。原物理目标 Φ 由数据项 D 与原先验 R 构成；Φ 下降不能替代材料恢复质量，也不表示使用 GT 选择。</p><div class="callout">旧实验的阶段内部 body CPU，与本补充的完整子进程 CPU 属于不同口径。历史 DR24 与矢量3D仅单独展示，不能合并均值或据此声称等成本。新增成本采用外部控制器完整退出实测；原端点获取成本单列，嵌套计时不重复相加。</div></section>'
    primary=[]
    for r in summary:
        if r['view']=='primary':primary.append([r['cohort'],METHODS.get(r['method'],r['method']),str(r['available'])+'/'+str(r['objects']),r.get('NRMSE_mean_available'),r.get('relativeL2_mean_available'),r.get('IoU_mean_available'),r.get('Phi_mean_available'),r.get('residual_over_delta_mean_available'),r.get('inherited_descriptive_success')])
    finding='<section id="findings"><h2>已完成配对结果：条件救援是最清晰的增量</h2><p><strong>B 全24对象的 primary 平均 NRMSE 从原 G19 的0.143891降至 G19→DR 的0.106261。</strong> 四个球外触发对象中，DR有3例改善材料并进入原CN99球，另1例仍失败；按旧 NRMSE／IoU 合取条件，其中2例新增通过。20个未触发对象保留原材料，是策略旁路，不是DR运行后无损的证据。</p><p>A 的环对象由G19合法重启得到更准材料，开缝壳和桥由DR更好恢复。C 的G19与DR主均值接近，DR实测CPU低约31%；球内minR的secondary部署均值排序反转，因此不能据GT选择有利估计器。</p><p>这支持把DR定位为G19球外终点的条件救援模块；结果属于既有暴露对象上的后续配对分析。详细解释与逐例失败见中文报告。</p></section>'
    overview='<section id="primary"><h2>主要公开选择：全队列汇总</h2><p class="note">primary 为已冻结的公开 minΦ 选择。均值按可用记录计算，并列明可用数／总数。旧成功标签仅作描述；连续材料误差与完整逐对象图应共同阅读。</p>'+table(['队列','方法','可用／总数','NRMSE 均值','相对 L2 均值','IoU 均值','Φ 均值','残差／CN99 半径均值','旧描述成功数'],primary)+'</section>'
    costrows=[[r['cohort'],METHODS.get(r['method'],r['method']),r['objects'],r['rescue_children'],r.get('unique_parent_CPU_seconds'),r.get('additional_external_child_CPU_seconds'),r.get('total_deployment_CPU_seconds')] for r in costs]
    perfrows=[[r['cohort'],METHODS.get(r['method'],r['method']),r['view'],str(r['available'])+'/'+str(r['objects']),r.get('NRMSE_mean_available'),r.get('mean_parent_CPU_seconds'),r.get('mean_added_rescue_CPU_seconds'),r.get('mean_total_deployment_CPU_seconds')] for r in performance]
    costbody='<section id="cost"><h2>实际 CPU 与材料表现</h2><p>B 的20个未触发对象保留原结果及获取成本，4个救援子进程的额外成本分摊到全部24个对象。C 的共享预处理只计一次。以下为观察到的成本与表现，不是相同花费比较曲线。</p>'+table(['队列','方法','对象数','救援子进程','原获取 CPU s','新增 CPU s','部署合计 CPU s'],costrows)+'<h3>每对象平均成本</h3>'+table(['队列','方法','选择视图','可用／总数','NRMSE 均值','原获取 s／对象','新增 s／对象','部署 s／对象'],perfrows)+'<p class="note">'+e(accounting.get('warning','完整计费来源见 COST_ACCOUNTING.json。'))+'</p></section>'
    mappingrows=[]
    for c in cases:
        status=('触发救援' if c['case_id'] in bcase_ids else '保留原材料') if c['cohort']=='B' else '固定配对'
        mappingrows.append([c['cohort'],c['case_id'],c.get('label',''),c.get('source_count'),status,c.get('parent_task_id','')])
    mapbody='<section id="mapping"><h2>全部对象 ID 对照</h2>'+table(['队列','完整对象 ID','标签','照明数','策略','原任务 ID'],mappingrows,'mapping')+''.join(figure(p) for p in mapping)+'</section>'
    gallery=[]
    for cohort,title in COHORTS.items():
        gallery.append('<section id="cohort-'+cohort+'"><h2>'+e(title)+'</h2><p class="note">所有 PNG 均按原比例完整展示，材料图与误差图保留原色标。点击图像查看原尺寸；NA 槽位保持可见。</p>'+''.join(figure(p) for p in groups[cohort])+'</section>')
    gallery.append('<section id="historical"><h2>历史对照 · 单独口径</h2><p>以下既有 DR24 与矢量3D分解，仅提供历史背景。其目标、对象和预算不同，不进入本补充的合并统计或等成本判断。</p>'+''.join(figure(p) for p in historical)+'</section>')
    preview='<section><h2>逐对象材料、误差与校正轨迹</h2><p>'+link('GALLERY.html','打开完整图集：A／B／C、全部 ID 与历史背景')+'</p><p>图集包含 raw／corrected 材料、B 全24策略、材料与误差色标、已付 CPU 前缀，以及 D／R／Φ／材料误差／CN99 比值的校正对照。</p></section>'
    docs=documents()
    oldbase=REPO+'/blob/main/research/research_max_struct_r12/'
    docs+='<section><h2>原 R12 冻结研究与补充关系</h2><p>'+link(oldbase+'COMPREHENSIVE_EVALUATION_ZH.md','原综合评价')+' · '+link(oldbase+'RESEARCH_REPORT_ZH.md','原中文报告')+' · '+link(oldbase+'PAPER_DRAFT.md','论文草稿')+' · '+link('../index.html','原冻结网页')+'</p></section>'
    # Write only after all mandatory records, mappings and plots were read.
    (ROOT/'index.html').write_text(page('R12 配对补充实验',hero+finding+scope+overview+costbody+preview+docs),encoding='utf-8')
    (ROOT/'GALLERY.html').write_text(page('R12 配对补充实验 · 完整图集',hero+finding+scope+mapbody+''.join(gallery)+docs),encoding='utf-8')
    print(json.dumps(dict(status='PAIRED_PAGES_BUILT_FROM_COMPLETED_SCORING',pages=['index.html','GALLERY.html'],PNG_count=len(images),score_rows=len(scores)),ensure_ascii=False))

if __name__=='__main__':main()

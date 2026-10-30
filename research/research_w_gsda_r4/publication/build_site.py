"""Publish the closed R4 artifact set without changing any scientific output."""
from pathlib import Path
import argparse,csv,datetime,html,json,re,shutil,subprocess,time
from zoneinfo import ZoneInfo

SLUG='research_w_gsda_r4'
REPO='https://github.com/migodam/b1-probabilistic-maxwell-imaging'
TAG='research-w-gsda-r4-20261030'
RELEASE=REPO+'/releases/tag/'+TAG
FULL=REPO+'/releases/download/'+TAG+'/RESEARCH_W_GSDA_R4_20261030.zip'
WEB=REPO+'/releases/download/'+TAG+'/RESEARCH_W_GSDA_R4_WEB_CHAT.zip'
RAW='https://raw.githubusercontent.com/migodam/b1-probabilistic-maxwell-imaging/main/research/'+SLUG+'/'
CSS='''*{box-sizing:border-box}body{margin:0;background:#f4f6f8;color:#192b37;font:17px/1.85 system-ui,-apple-system,"PingFang SC",sans-serif}main{max-width:1160px;margin:auto;padding:30px 24px 70px}nav{display:flex;gap:9px 21px;flex-wrap:wrap;margin-bottom:25px}a{color:#17668d;text-underline-offset:3px}h1{font-size:clamp(29px,4.5vw,43px);line-height:1.4;letter-spacing:-.025em}h2{font-size:26px;line-height:1.5;margin-top:1.8em}h3{line-height:1.5}.panel,article{background:white;border:1px solid #dce3e8;border-radius:14px;padding:28px;margin:23px 0}.lede{font-size:21px}.status{background:#fff1eb;border-left:4px solid #b55c43;padding:17px 21px;margin:24px 0;border-radius:4px}.note{background:#edf3f7;border-left:3px solid #7c9bad;padding:13px 19px}.muted,figcaption{font-size:14px;color:#526575}.metrics{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:12px;margin:27px 0}.metric{padding:18px;background:#eaf1f5;border:1px solid #dce5eb;border-radius:12px}.metric strong{display:block;font-size:29px;line-height:1.4;color:#164d68}.metric span{font-size:14px}.table-scroll{overflow-x:auto;margin:20px 0}table{border-collapse:collapse;width:100%;font-size:14px;display:block;overflow-x:auto}th,td{border:1px solid #d7e0e6;padding:9px 11px;vertical-align:top;text-align:left}th{background:#edf3f7}img{width:100%;height:auto;display:block;background:white}figure{margin:30px 0}figure>a{display:block}figure img{border:1px solid #dce3e8}figcaption{padding-top:10px}code,pre{background:#edf2f5;font-size:.88em}code{padding:2px 4px;border-radius:4px}pre{padding:16px;overflow-x:auto;border-radius:8px}pre code{padding:0}p,li{overflow-wrap:anywhere}blockquote{margin:20px 0;padding:2px 20px;border-left:3px solid #9fb4c1}footer{margin-top:38px;padding-top:20px;border-top:1px solid #d5dfe6;font-size:14px;color:#526575}.downloads{display:flex;gap:12px 20px;flex-wrap:wrap}.downloads a{display:inline-block;padding:9px 14px;background:#e6eff4;border:1px solid #cbdce6;border-radius:8px}@media(max-width:620px){main{padding:20px 13px 45px}.panel,article{padding:18px}nav{font-size:15px}h2{font-size:23px}.lede{font-size:19px}}'''

def jread(p):return json.loads(p.read_text())
def dump(p,obj):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def fig(path,caption):return '<figure><a href="'+html.escape(path)+'"><img loading="lazy" src="'+html.escape(path)+'" alt="'+html.escape(caption)+'"></a><figcaption>'+html.escape(caption)+' 点击查看原始分辨率。</figcaption></figure>'
def card(title,content):return '<section class="panel"><h2>'+title+'</h2>'+content+'</section>'

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--repository',type=Path,required=True);p.add_argument('--node',required=True);p.add_argument('--marked-cli',required=True);a=p.parse_args()
    source=a.source.resolve();repo=a.repository.resolve();dest=repo/'research'/SLUG;dest.mkdir(parents=True,exist_ok=True)
    copied=[];omitted=[]
    for f in sorted(source.rglob('*')):
        if not f.is_file():continue
        rel=f.relative_to(source)
        if '__pycache__' in rel.parts or any(x.startswith('.') for x in rel.parts):continue
        if f.suffix in {'.npz','.npy','.zip','.pyc','.pyo'} or rel.as_posix()=='PHYSICS_LEDGER.jsonl':
            omitted.append(dict(path=str(rel),bytes=f.stat().st_size,location='existing repository Release, full ZIP'))
            continue
        target=dest/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,target)
        copied.append(dict(path=str(rel),bytes=f.stat().st_size))
    dump(dest/'publication/GIT_FILE_LAYOUT.json',dict(git_files=copied,arrays_and_large_journal_in_release=omitted,full_asset=FULL,web_asset=WEB,source_science_bytes_unchanged=True,sha256_not_used=True))
    ledger=source/'PHYSICS_LEDGER.jsonl'
    dump(dest/'PHYSICS_LEDGER_LOCATION.json',dict(path_in_full_archive='research_w_gsda_r4/PHYSICS_LEDGER.jsonl',bytes=ledger.stat().st_size,full_asset=FULL,web_asset=WEB,warning='Restore archive before running summarize_cost.py; the compact Git tree does not contain this authoritative full journal'))
    receipt=jread(source/'delivery/DELIVERY_RECEIPT.json');cost=jread(source/'PHYSICS_COST.json');summary=jread(source/'evaluation/FINAL_METHOD_SUMMARY.json')
    now=datetime.datetime.now(ZoneInfo('Asia/Singapore')).isoformat()
    note=dict(authorization='User: r4发布一下',publication_prepared_at_singapore=now,timezone='Asia/Singapore',new_physics_FE=0,new_complex_RHS=0,scientific_gates_changed=False,scientific_reports_unchanged=True,closed_time_publication_fields='Original NOT_AUTHORIZED/NOT_DONE fields describe scientific close-time state; this new publication record supersedes publication authorization only',release=RELEASE,full_asset=FULL,full_bytes=receipt['archive_bytes'],web_asset=WEB,web_bytes=receipt['web_chat_archive_bytes'],source_checkout='R4 local closed source; source commits in ENVIRONMENT.json')
    dump(dest/'publication/PUBLICATION_NOTE.json',note)
    write=(dest/'publication/PUBLICATION_NOTE.md').write_text
    write('# R4 发布记录\n\n用户已明确授权“r4发布一下”。本次只公开已封账结果，新增物理FE/RHS为0，不改变任何科学裁决。\n\n原报告、GO_NO_GO和ZIP中的“未发布/未授权”是科研封账时的历史状态，原件保留。本记录及网页发布栏说明之后的授权与发布；不将失败改写为成功。\n\n完整ZIP保留全部数组与大传播缓存；网页复盘ZIP省去117份大缓存，但保留全部49个非开发终点、观测、预测、Jacobian快照、物理核验数组、代码与账本。Git保留报告、图、源码、配置与紧凑记录，大数组和完整PHYSICS_LEDGER在Release。\n\n解压两种ZIP时，其顶层是research_w_gsda_r4；在仓库research/目录解压以恢复相同路径。未恢复完整账本前，不要在瘦Git副本运行成本重汇总；已有PHYSICS_COST.json是原研究记录。\n\n[完整数据Release]('+RELEASE+') · [匿名原始报告]('+RAW+'REPORT_ZH.md)\n')
    banner='<p class="muted">发布于2026-10-30（新加坡时间），经用户明确授权。科研已于2026-10-30 06:02封账；本次新增0 FE / 0 RHS。原报告中“未发布”保留为封账历史状态。<a href="publication/PUBLICATION_NOTE.md">发布记录</a></p>'
    nav='<nav><a href="../../">B1首页</a><a href="index.html">R4概览</a><a href="report.html">完整报告</a><a href="gallery.html">实际图像</a><a href="methods.html">方法实现</a><a href="physics.html">物理与误差</a><a href="handoff.html">网页端交接</a><a href="evidence.html">下载与复现</a></nav>'
    footer='<footer>静态二维 scalar-TM；未知墙限3×2材料格，目标为7×7自由复材料格。有限pilot，无连续PDE或概率安全证书。<br><a href="GO_NO_GO.json">分项裁决</a> · <a href="PHYSICS_COST.json">费用汇总</a> · <a href="FAILURE_LEDGER.jsonl">失败</a> · <a href="'+RAW+'REPORT_ZH.md">匿名原始Markdown</a> · <a href="'+RELEASE+'">完整证据Release</a></footer>'
    def page(name,title,body):
        (dest/name).write_text('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="description" content="Research W-GSDA-R4 已完成的全波实验、九数字实际失败、精确吸收核验与成本边界。"><title>'+html.escape(title)+'</title><style>'+CSS+'</style></head><body><main>'+nav+banner+body+footer+'</main></body></html>')
    def render(md):
        result=subprocess.run([a.node,a.marked_cli,'--gfm'],input=md,text=True,capture_output=True,check=True)
        return result.stdout
    downloads='<div class="downloads"><a href="'+WEB+'">网页复盘包 · 232 MB</a><a href="'+FULL+'">完整数据包 · 632 MB</a><a href="'+RAW+'REPORT_ZH.md">原始报告Markdown</a></div>'
    rows=[]
    labels={'B0_FULL_FWI':'普通Full-FWI + TV','FULL_GN_TV_SUPPLEMENT':'B0 + Full-GN补充','B4_GSDA':'原生CSI + 风险筛选B4'}
    for s in summary:
        rows.append('<tr><td>'+labels[s['campaign']]+'</td><td>9 / 9</td><td>'+f"{s['mean_IoU']:.3f}"+'</td><td>'+f"{s['mean_complex_relative_L2']:.3f}"+'</td><td>'+f"{s['mean_train_white_RMS']:.3f}"+'</td><td>'+f"{s['mean_heldout_white_RMS']:.3f}"+'</td><td>0 / 9</td></tr>')
    metrics_table='<table><thead><tr><th>方法</th><th>实际数字运行</th><th>IoU均值</th><th>复材料相对L2</th><th>训练白化RMS</th><th>留出白化RMS</th><th>完整成像通过</th></tr></thead><tbody>'+''.join(rows)+'</tbody></table>'
    hero='<h1>数字1–9都真正重建了，<br>但还没有达到准确穿墙成像</h1><p class="lede">固定观测Green、完整材料吸收、内部谱预条件和可逆REOPEN均有实际代码与数据。九数字运行完成；准确材料/形状成像与同质量加速没有成立。</p><div class="status"><strong>可运行：受限模型PASS。九数字完整成像：0/9。当前GSDA算法增量：NO-GO。</strong><br>不能把同离散吸收恒等成立或电流残差较小，写成目标图像已恢复。</div><div class="metrics"><div class="metric"><strong>27 + 21 + 1</strong><span>九数字三方法、成对控制、TV0补充的实际终点</span></div><div class="metric"><strong>0 / 9</strong><span>每个主方法的完整材料／形状门通过数</span></div><div class="metric"><strong>12</strong><span>B2/B3真实定时REOPEN；自动触发为0</span></div><div class="metric"><strong>+432 RHS</strong><span>B4与普通B5三个输出相同，额外诊断费用</span></div></div>'+downloads
    body=hero+card('所有方法与失败一起展示',metrics_table+fig('DIGITS_1_TO_9_CONTRAST.png','真实7×7材料图；真值、普通FWI、B0+Full-GN与B4固定列和统一色轴。图像来自实际数组，无逐数字真值选支。'))
    body+=card('不只看IoU：孔洞、笔画、材料与损耗', '<p>普通FWI和Full-GN留下一些笔画，但6/8/9的闭环与孔洞仍错，损耗大多接近零或弥散。B4多为弱材料云团。每个数字都实际执行，但B4数字3的最终物理目标选择了第0次公开均匀起点，不能说每例都优于初始化。</p>'+fig('DIGITS_1_TO_9_LOSS.png','全部九数字的损耗系数ell；共享色轴，任意频率虚介电率为ell/f。'))
    body+=card('同数据TV消融有部分收益，不能归给Green吸收','<p>数字8在相同数据、起点和其余先验下：TV0→TV.2，IoU 0.343→0.500，复材料相对误差0.943→0.756；质心误差0.075→0.106 λ0变差。两种都只恢复1个孔洞，没有恢复完整8。不同TV总目标不可直接比较大小。</p>'+fig('figures/tv_ablation/DIGIT8_TV0_VS_TV02.png','真实同数据数字8消融；孔洞、损耗与顶端伪热点的失败保留。'))
    body+=card('机制可算对，不等于优化和成像有增量','<p>精确吸收包含背景数据、有效照明及返回接收路径，同离散完整Maxwell核验通过。B2/B3实际吸收后重开，但目标图并未准确。三个独立控制中B4与普通同谱B5的最终材料和复预测逐元素完全相同，反而多用432 RHS；风险筛选全部KEEP_ACTIVE，没有付费目标J筛查。</p><p>48模态投影的全场误差仍约0.89–0.93，投影内Schur精确不能修复全空间信息丢失。主B3/B4只作保留全部电流方向的预条件，未实现安全删除完整弱空间。</p>'+fig('figures/paired_controls_with_B0/NATIVE_CONTROLS_CONTRAST.png','同三个非数字世界的全部方法；B0含强伪热点，CSI各版本为弱云团，均不包装为准确物体。'))
    body+=card('数字背后的范围与账本','<p>13个生成世界中1个开发，12个测试parent，每个1份噪声。49个终点是同批parent的方法与补充运行，不是49个独立物理样本。模型为4发16收、三频、实虚每坐标σ=.01，墙12与目标98个实变量共同优化；无字形模板或真实墙Green。</p><p>本轮预约13,158 FE / 191,940 complex RHS；物理在158.80分钟封账。已记录Book CPU约6,928秒，但早期Python启动未覆盖；不能把它当完整端到端CPU。仅2/5粗细网格真值差满足旧0.1白化门，无连续PDE认证。多数基线仍未达驻点，未宣称成熟FWI的上限。</p><p><a href="report.html">完整中文实验报告 →</a>　<a href="handoff.html">把真实障碍交给网页端分析 →</a></p>')
    page('index.html','Research W-GSDA-R4｜九数字真实全波实验与失败',body)
    page('report.html','R4完整中文实验报告','<article>'+render((dest/'REPORT_ZH.md').read_text())+'</article>')
    page('methods.html','R4算法实现与范围','<article>'+render((dest/'ALGORITHM_IMPLEMENTATION.md').read_text())+'</article>'+card('从测量到材料图',render((dest/'RUN_ALGORITHM.md').read_text())))
    page('physics.html','R4物理验证与数值误差','<article>'+render((dest/'PHYSICS_VALIDATION.md').read_text())+'</article>'+fig('figures/physics/GD_MODES_AND_FEEDBACK.png','实际GD模态、交叉反馈消去误差和完整场误差；没有把小矩阵恒等误差当成物理成功。')+fig('figures/physics/INDEPENDENT_MESH_CHECK.png','全部预声明粗细网格结果；0.1门以上失败保留。')+fig('figures/physics/ABSORPTION_REOPEN_BORN.png','精确吸收/重开与完整目标Born诊断；物理一致性不替代成像质量。'))
    gallery='<h1>实际图像与失败边界</h1><p>所有图来自封存数组。7×7是物理材料分辨率，放大图片不是新的物理成像分辨率。</p>'
    for path,caption in [('DIGITS_1_TO_9_CONTRAST.png','全部数字对比，固定列，不挑最好方法'),('DIGITS_1_TO_9_LOSS.png','全部损耗材料图'),('figures/tv_ablation/DIGIT8_TV0_VS_TV02.png','同数据TV消融'),('figures/paired_controls_with_B0/NATIVE_CONTROLS_CONTRAST.png','普通FWI与六个原生电流控制'),('figures/quality_cost/CSI_PSEUDODATA_VS_PHYSICAL_CHECKPOINTS.png','CSI电流数据拟合和对应材料完整FEM目标的差距'),('figures/quality_cost/ENDPOINT_QUALITY_VS_CPU.png','实际终点质量与已记录CPU；GN使用B0+GN累计费用'),('figures/quality_cost/HELDOUT_ERROR_VS_FITTING_CPU.png','独立留出频率误差与拟合CPU')]:gallery+=fig(path,caption)
    page('gallery.html','R4实际九数字图像',gallery)
    page('handoff.html','R4网页端研究交接','<article>'+render((dest/'WEB_CHAT_HANDOFF.md').read_text())+'</article><article>'+render((dest/'QUESTIONS_FOR_WEB_CHAT.md').read_text())+'</article>')
    evidence='<h1>代码、数据与复现</h1>'+downloads+'<div class="note">完整ZIP保留所有缓存与原始数组；232 MB网页包省去117个大型传播缓存，保留全部49个非开发终点和物理核验数组。Git中的瘦副本不含NPZ/NPY与46 MB原始物理日志。先解压数据包再复查，切勿在缺失账本时重汇总费用。</div>'
    evidence+=card('实际可读的原始入口','<ul><li><a href="'+RAW+'REPORT_ZH.md">匿名中文报告Markdown</a></li><li><a href="'+RAW+'WEB_CHAT_HANDOFF.md">匿名网页端交接Markdown</a></li><li><a href="'+RAW+'QUESTIONS_FOR_WEB_CHAT.md">实验暴露的具体问题</a></li><li><a href="DIGITS_1_TO_9_METRICS.csv">27个主方法数字行</a>；<a href="EXPERIMENT_RESULTS.csv">49个非开发终点</a>；<a href="METHOD_ABLATIONS.csv">21个成对控制</a></li><li><a href="GO_NO_GO.json">机器可读分项裁决</a>；<a href="COST_ACCOUNTING_SCOPE.md">费用测量缺口</a>；<a href="FAILURE_LEDGER.jsonl">全部失败</a></li><li><a href="publication/GIT_FILE_LAYOUT.json">Git与Release文件分配</a>；<a href="delivery/ARCHIVE_MANIFEST.json">完整包清单</a>；<a href="delivery/WEB_CHAT_SCOPE.json">精简包省略清单</a></li></ul>')
    evidence+=card('如何恢复目录与重放已有结果','<p>把ZIP解压到仓库的 <code>research/</code>：ZIP内顶层是 <code>research_w_gsda_r4/</code>。例如：</p><pre><code>cd research\nunzip /path/to/RESEARCH_W_GSDA_R4_WEB_CHAT.zip\ncd research_w_gsda_r4\npython tools/score_all.py --results results\npython tools/summarize_cost.py</code></pre><p>重放评分不调用PDE。新物理部署必须显式给独立预算；原研究已封账。本次发布没有新PDE。全新依赖环境安装和新部署封装的再次PDE运行仍为NOT_RUN，不冒称做过。</p>')
    evidence+='<article>'+render((dest/'RUN_ALGORITHM.md').read_text())+'</article><article>'+render((dest/'COST_ACCOUNTING_SCOPE.md').read_text())+'</article>'
    page('evidence.html','R4下载与复现',evidence)
    page('failures.html','R4全部数字失败与范围','<article>'+render((dest/'DIGITS_1_TO_9_FAILURES.md').read_text())+'</article><article>'+render((dest/'SOURCE_CHANGE_RECORD.md').read_text())+'</article>')
    marker='W_GSDA_R4_PUBLICATION'
    desc='九数字全部实际全波反演；普通FWI、Full-GN与GSDA均0/9完整成像门。精确吸收核验通过，但风险版与普通同谱对照输出相同、额外432RHS，当前算法增量NO-GO。原始数组、失败、成本和两个复现包公开。'
    for filename,prefix in [('index.html','research/'),('research/index.html','')]:
        f=repo/filename;text=f.read_text()
        if '<!-- '+marker+' -->' in text:raise RuntimeError('publication marker already exists; do not duplicate')
        section='<!-- '+marker+' --><section class="card"><h2>Research W-GSDA-R4 · 九数字真实全波反演与成像失败</h2><p>'+desc+'</p><p><a href="'+prefix+SLUG+'/">九数字全部图像、实验报告与下载 →</a></p></section><!-- /'+marker+' -->'
        f.write_text(text.replace('<main>','<main>'+section,1))
    f=repo/'START_HERE.md';f.write_text('<!-- '+marker+' -->\n## Research W-GSDA-R4：九数字全部运行，准确成像与算法增量未成立\n\n'+desc+'\n\n[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/'+SLUG+'/) · [原始报告](research/'+SLUG+'/REPORT_ZH.md) · [网页端交接](research/'+SLUG+'/WEB_CHAT_HANDOFF.md) · [完整数据Release]('+RELEASE+')\n<!-- /'+marker+' -->\n\n'+f.read_text())
    (dest/'.gitignore').write_text('__pycache__/\n*.pyc\n*.npz\n*.npy\n*.zip\n*.lock\nPHYSICS_LEDGER.jsonl\n')
    notes='''Research W-GSDA-R4 已封账实验公开，新增物理FE/RHS为0。\n\n- 数字1–9均真正生成全波数据、在未知墙下反演；普通FWI、B0+Full-GN、B4各九例，完整材料/形状门均0/9。\n- 固定GS contrast-source和精确材料吸收在同离散核验通过；未建立连续Maxwell认证。\n- 三个成对控制中B4与普通B5材料及复预测完全相同，额外432RHS；当前算法增量与同质量加速NO-GO。\n- 同数据数字8的TV.2有部分支持/材料改善，但两孔仅恢复一个；不得归因于Green吸收。\n- 完整包包含原始复场/观测、材料、全部缓存、源码、配置、失败和账本。网页包仅省去大缓存，范围清单明确。\n- 科研封账时“未发布”字段保留；publication/PUBLICATION_NOTE说明本次用户授权，不回写科学历史。\n\n只读网页和raw报告、下载不需要登录。\n'''
    (dest/'publication/RELEASE_NOTES.md').write_text(notes)
    dump(dest/'publication/PREPARATION.json',dict(source=str(source),repository=str(repo),prepared_singapore=now,source_scientific_files_changed=False,new_FE=0,new_RHS=0,html_pages=['index.html','report.html','methods.html','physics.html','gallery.html','handoff.html','evidence.html','failures.html'],copied_file_count=len(copied),omitted_to_release_count=len(omitted),raw_report=RAW+'REPORT_ZH.md',release=RELEASE,full_asset_bytes=receipt['archive_bytes'],web_asset_bytes=receipt['web_chat_archive_bytes']))
    print(json.dumps(dict(copied=len(copied),omitted=len(omitted),pages=8,destination=str(dest)),ensure_ascii=False))

if __name__=='__main__':main()

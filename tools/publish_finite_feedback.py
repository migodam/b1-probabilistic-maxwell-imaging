"""Prepare the existing B1 site and a no-hash evidence archive. Zero physics."""
from pathlib import Path
import json,sys,zipfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from finite_feedback.common import OUT,dump
from publish_lxz import render,STYLE

def main():
    rel='research/finite_feedback_fiber_diagnostic_20261016/'
    release='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/finite-feedback-fiber-20261016'
    nav='<nav><a href="../../">B1首页</a><a href="index.html">本轮入口</a><a href="report.html">中文报告</a><a href="gallery.html">数值与材料图</a><a href="'+release+'">完整证据</a><a href="START_HERE.md">GPT入口</a></nav>'
    top='<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Maxwell有限反馈与局部解集诊断</title><style>'+STYLE+'</style><main>'+nav
    body='''<h1>Maxwell 有限反馈与GT局部连接诊断</h1><p>2026-10-16 · 独立小型CPU实验，已完成并封账。</p><section><h2>弱方向很多，是否就能走到真实材料？</h2><p><strong>没有得到这种连接证据。</strong> N3–GT差量的98.09%–99.06%在所测阈值下属于局部弱方向，但N3不在原P0合法先验及原noise球内。Full-feedback能准确预测有限变化，却不会创造可行材料世界。</p><div class="metrics"><span><b>45 / 48</b>新二维states</span><span><b>1,880 / 2,048</b>complex RHS</span><span><b>11.146 s</b>受控wall</span><span><b>0</b>显著新结构</span></div><p><a href="report.html">阅读独立中文报告 →</a> · <a href="gallery.html">看原始响应与材料变化 →</a></p></section><section><h2>本轮最清楚的失败归因</h2><p>一例合法TRAIN chart的finite响应达到冻结切向预测的580.75倍；同一Maxwell Jacobian对实际材料增量的预测误差却仅0.446%。材料decoder曲率与Maxwell反馈曲率必须分开。</p><img src="figures/06_chart_vs_maxwell.png" alt="材料参数化曲率与Maxwell反馈的拆分"></section><section><h2>高阶反馈预测确实改善</h2><img src="figures/02_gt_feedback.png" alt="一至二十阶与完整有理反馈比较"><p>这是离线GT直线路径的finite-response核验，不能当作通向GT的在线算法或等数据fiber。</p></section><section><h2>合法TRAIN材料路径</h2><img src="figures/04_material_250034_ring.png" alt="已知TRAIN材料的合法小幅变化"><p>保持原P0和noise阈值的小幅连续变化，未出现显著新结构；已知材料不称未知成像成功。</p></section><section><h2>继续哪条方向？</h2><p>先得到合法且敏感分量合格的physics-only重建anchor，并改善材料decoder/有限retraction，再讨论tangent field或finite transport。不启动新NN、Flow、Mie或三维campaign。</p><p><a href="THEORY_QUESTIONS_FOR_PRO.md">给下一轮Pro的理论问题</a> · <a href="PRIOR_SUPPORT_AND_ANCHOR_AUDIT.md">先验与anchor审计</a> · <a href="RESOURCE_SUMMARY.json">资源</a> · <a href="REPRODUCE.md">复现</a></p></section>'''
    (OUT/'index.html').write_text(top+body+'</main></html>')
    (OUT/'report.html').write_text(top+'<article>'+render((OUT/'DIAGNOSTIC_REPORT_ZH.md').read_text())+'</article></main></html>')
    gallery=(OUT/'gallery.html').read_text()
    if nav not in gallery:gallery=gallery.replace('<h1>',nav+'<h1>',1)
    (OUT/'gallery.html').write_text(gallery)
    notice='## 最新诊断：有限Maxwell反馈与GT差量\n\n[中文报告]('+rel+'DIAGNOSTIC_REPORT_ZH.md) · [网页与数值图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/'+rel+')\n\n45个二维states、1880个complex RHS；高阶响应已核验，GT合法局部路径未建立。材料decoder曲率是当前关键限制。\n\n'
    for name in ['README.md','START_HERE.md']:
        pp=ROOT/name;s=pp.read_text()
        if '## 最新诊断：有限Maxwell反馈' not in s:pp.write_text(notice+s.replace('## 最新研究：LXZ','## 历史研究：LXZ'))
    card='<section class="card"><h2>最新：有限Maxwell反馈诊断</h2><p>45个二维状态：GT差量多数在弱方向，但没有合法连接路径。分离decoder曲率与Maxwell反馈。</p><a href="'+rel+'">中文报告与原始数值图 →</a></section>'
    pp=ROOT/'index.html';s=pp.read_text()
    if rel not in s:pp.write_text(s.replace('<div class="grid">','<div class="grid">'+card,1).replace('最新：LXZ局部–全局多解验证','LXZ局部–全局多解验证'))
    pp=ROOT/'research/index.html';s=pp.read_text()
    if 'finite_feedback_fiber_diagnostic_20261016/' not in s:pp.write_text(s.replace('<main>','<main><section><h2>最新：有限反馈与GT差量</h2><p>45新states；没有合法GT路径，关键限制是decoder曲率与anchor。</p><a href="finite_feedback_fiber_diagnostic_20261016/">报告与原始图 →</a></section>',1))
    dump(OUT/'PUBLICATION.json',dict(existing_repository='migodam/b1-probabilistic-maxwell-imaging',new_repository=False,release=release,
        historical_results_unchanged=True,synthetic_exposed_development_and_TRAIN_only=True,private_Dyson_conversation_excluded=True,
        third_party_fulltexts_included=False,hashes_computed=False,license_status_unchanged=True,
        scientific_status_independent_of_publication=True))
    archive=ROOT/'publication_local/Finite_Feedback_Fiber_Diagnostic_20261016_Evidence.zip';archive.parent.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=5) as z:
        for pp in sorted(OUT.rglob('*')):
            if pp.is_file() and pp.name not in {'budget.lock','DYSON_CONVERSATION.json','EVIDENCE_RELEASE.json','PUBLICATION_RECEIPT.json'} and '__pycache__' not in pp.parts:z.write(pp,pp.relative_to(ROOT))
        for pp in (ROOT/'src/finite_feedback').glob('*.py'):z.write(pp,pp.relative_to(ROOT))
        for name in ['analyze_finite_feedback.py','report_finite_feedback.py','plot_finite_feedback.py','publish_finite_feedback.py']:z.write(ROOT/'tools'/name,'tools/'+name)
    dump(OUT/'EVIDENCE_RELEASE.json',dict(asset=archive.name,bytes=archive.stat().st_size,files=len(zipfile.ZipFile(archive).infolist()),
        release=release,all_45_paid_endpoint_arrays=True,private_conversation_excluded=True,hash_checks=False))
    print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size)))

if __name__=='__main__':main()

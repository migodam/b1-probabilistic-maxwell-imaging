# coding: utf-8
"""Build static publication around frozen scientific results; no physics."""
from pathlib import Path
import json,platform
P=Path(__file__).resolve().parents[1]
REPO=P.parents[1]
def main():
    s=(P.parent/'research_gn_green_error_20261106/index.html').read_text()
    substitutions={
      'GN–GT · Green强弱与Two-fold误差能量':'GN–Green R2 · 有效材料搜索空间',
      'GN error geometry · 2026-11-06':'GN–Green R2 · 2026-11-07',
      'GN–GT误差落在哪里？':'往 GT 搜索的材料空间有多大？',
      '144个冻结GN/TRF结果；材料与真实电流误差；Gin/Gob强、弱、补空间；交集、串联与TSOM投影分别计量。':'共享实材料、多照明约束、有限步与校正；52端点空间诊断、96搜索任务，保留失败和普通优化的增量。',
      'MATHEMATICAL_DEFINITIONS.md':'MATHEMATICAL_THEORY.md',
      'results/REQUESTED_TABLES_ZH.md':'CODEX_TO_PRO_HANDOFF.md','完整请求表':'下一轮理论交接',
      'results/SUMMARY.json':'results/RESULTS_SUMMARY.json','results/requested_per_case_statistics.csv':'results/RESULTS_BY_CASE.csv',
      'research_gn_green_error_20261106':'research_gn_search_space_r2',
      '主强门为σ/σmax≥.1；能量为平方范数。真交集不一定覆盖全域，TSOM投影空间可以重叠；静态G的材料坐标分解不等于非线性材料Jacobian。':'已披露历史开发对象；能量为材料误差平方范数。当前空间不等于全局解集；离线 best@K 不等于在线 top1；NN与TSOM因果消融未运行。'
    }
    for a,b in substitutions.items():s=s.replace(a,b)
    s=s.replace('displayMath:[["$$","$$"]]','displayMath:[["$$","$$"]],inlineMath:[["$","$"]]')
    (P/'index.html').write_text(s)
    manifest=dict(source_commit_at_start='32f7bfb308fd583da51bf3414a8aaad8d9da5ba5',publication_base_commit='fd970123dfbbba38118301b89d60d26c5c9903dd',physics='research/research_z_alg_r4_diversity/src/physics.py unchanged',historical_spectra='gn-green-error-geometry-20261106 Release full-reproduction.zip',data_status='HISTORICAL_DISCLOSED_DEVELOPMENT',hash_checks='NOT_PERFORMED_PER_USER_PREFERENCE',python='3.13.13',numpy='2.5.3',scipy='1.18.1',platform=platform.platform(),threads_per_worker=1,max_numerical_workers=2,NN_training=False)
    (P/'audit/SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    report=P/'RESEARCH_REPORT_ZH.md';text=report.read_text()
    if '中位覆盖77.63%'not in text:text=text.replace('完整8/16/32/64/128/256/576表','1%泄漏宽空间在低对比端点有466–469维，中位覆盖77.63%的GT误差。完整8/16/32/64/128/256/576表')
    report.write_text(text)
    for name in ['results/e2/SCORES.json','results/e3/SCORES.json']:
        f=P/name;d=json.loads(f.read_text());d.pop('truth_disclosed_only_after_seal',None);d['historical_truth_already_disclosed']=True;d['formal_scoring_after_prediction_seal']=True;f.write_text(json.dumps(d,indent=2)+'\n')
    f=P/'src/score_search.py';t=f.read_text().replace('truth_disclosed_only_after_seal=True','historical_truth_already_disclosed=True,formal_scoring_after_prediction_seal=True');f.write_text(t)
    for name in ['README.md','START_HERE.md']:
        f=REPO/name;t=f.read_text()
        if '<!-- GN_GREEN_R2 -->'in t:continue
        entry='<!-- GN_GREEN_R2 -->\n## GN–Green R2：共享材料搜索空间与真实小型搜索\n\n509维复电流暗空间不等于材料latent；最弱64维仅覆盖约3%GT误差。52端点空间、96搜索任务全部完成；Green过滤未改善在线top1，严格停止零正则TRF四例有实质增量。\n\n[中文网页](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/research_gn_search_space_r2/) · [复现与下载](research/research_gn_search_space_r2/REPRODUCE.md) · [下一轮理论交接](research/research_gn_search_space_r2/CODEX_TO_PRO_HANDOFF.md)\n<!-- /GN_GREEN_R2 -->\n\n'
        pos=t.find('\n')+1;f.write_text(t[:pos]+'\n'+entry+t[pos:])
    f=REPO/'index.html';t=f.read_text()
    if 'research/research_gn_search_space_r2/'not in t:
        card='<section class="card"><h2>GN–Green R2：往GT搜索的材料空间</h2><p>52端点空间诊断、96搜索任务。最弱64维仅覆盖约3%误差；过滤搜索未改善在线top1，普通优化停止与正则偏差已分离。完整失败、成本和NN下一轮建议保留。</p><a href="research/research_gn_search_space_r2/">中文报告、真实图像与完整复现 →</a></section>'
        f.write_text(t.replace('<main>','<main>'+card,1))
if __name__=='__main__':main()

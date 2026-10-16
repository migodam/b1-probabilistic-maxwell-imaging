"""Publish the sealed F3 results in the existing B1 site; no physics calls."""
from pathlib import Path
import html, json, platform, shutil, sys, zipfile
from publish_lxz import STYLE
from publish_research_f import rendered

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'research/research_f3_quotient_exploration'
REL = 'research/research_f3_quotient_exploration/'
RELEASE = 'https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/research-f3-quotient-20261017'
OLD = ROOT.parent / 'b1-research-f/research/research_f_intrinsic_transport'


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def main():
    if not (OUT / 'STAGE_CLOSED.json').exists():
        raise RuntimeError('Only publish sealed research')
    nav = '<nav><a href="../../">B1首页</a><a href="index.html">本轮入口</a><a href="report.html">中文报告</a><a href="gallery.html">材料与路径图</a><a href="START_HERE.md">GPT入口</a><a href="' + RELEASE + '">原始证据</a></nav>'
    math = '<script>window.MathJax={tex:{inlineMath:[["\\\\(","\\\\)"]],displayMath:[["\\\\[","\\\\]"]]}};</script><script defer src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>'
    top = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Research F3 · Maxwell商坐标与有限噪声探索</title><style>' + STYLE + '.math{overflow-x:auto;padding:8px}figcaption{color:#516672;font-size:15px}</style>' + math + '<main>' + nav
    body = '''<h1>Research F3 · 商坐标解决驻点障碍，有限探索给出停止边界</h1>
<p>2026-10-17 · 已执行并封账 · 原P0 / 24×24 scalar TM/VIE · 无神经网络训练</p>
<section><h2>这轮真正得到了什么？</h2><div class="metrics"><span><b>2</b>新资格下的物理驻点</span><span><b>3</b>很接近的竞争端点</span><span><b>0</b>新增拓扑类型</span><span><b>55 / 64</b>新增材料状态</span></div>
<p>消除共同权重尺度并保留正确非均匀先验后，只需4次新状态求解，两份旧F1材料就达到新注册的边界驻点标准。旧Research F的Gate2失败记录保留。</p>
<p>weak条件发现RMS≈0.005276、ΔΦ≤2的合法候选；ring没有。三个端点相互RMS仅0.000639–0.001100，没有拓扑变化，不能称三个结构世界或后验模式。</p>
<p><strong>裁决：结束原P0局部Explorer研究；停止当前高阶反馈增量和NN训练主张。</strong>普通材料retraction已改善有限步漂移，高阶方法没有同成本优势。下一阶段另立有物理依据的P1，本轮未启动。</p>
<p><a href="report.html">完整中文报告与五项决定 →</a> · <a href="GO_NO_GO_DECISION.md">逐项裁决</a></p></section>
<section><h2>材料真的变了，但没有变成另一种结构</h2><img src="figures/01_material_weak.png" alt="旧F1、新物理anchor、native和高阶候选及离线真值，共同色标"><p>同一已暴露parent的两个背景条件，额外已知family/global geometry。GT只用于离线评分。候选是candidate-only，概率质量为null；不与N3作同信息量的胜负比较。</p><a href="gallery.html">查看全部材料、差分与路径图 →</a></section>
<section><h2>真正的增量来自材料坐标与retraction</h2><img src="figures/06_matched_distance.png" alt="同材料RMS的四方法likelihood漂移"><p>weak在相同材料RMS下，冻结切向ΔΦ=10.644，普通native=1.141，高阶反馈=1.633。D每条件56 RHS，对照native前7端点冷启动等价也是56 RHS；weak合格数分别1和2。没有高阶增量。</p></section>
<section><h2>合格端点不等于合格传输路径</h2><img src="figures/04_path_obstruction.png" alt="9个固定参数中点的likelihood障碍"><p>全部9个中点仍在原噪声球，但全部ΔΦ&gt;2。只否定所测参数直线路径的竞争likelihood保持性，不证明其他路径不存在，也不证明分离模式。</p></section>
<section><h2>完整费用与概率边界</h2><div class="metrics"><span><b>712 / 2048</b>complex RHS</span><span><b>17.641 s</b>受控作业wall</span><span><b>55 / 55</b>立即保存的端点</span><span><b>0</b>NN / GPU / Mie / 3D</span></div><p>包含数值失败、恢复、拒绝和中点诊断，账本未重置；方法独立cold/warm速度排名未测。MH每条件仅8步、接受3步，不能提供后验质量、混合或credible interval结论。</p><p><a href="MATCHED_COST_RESULTS.md">同成本比较</a> · <a href="RESOURCE_LEDGER.jsonl">逐项账本</a> · <a href="POSTERIOR_TARGET_AUDIT.md">概率目标</a> · <a href="FAILURE_ANALYSIS.md">失败记录</a> · <a href="REPRODUCE.md">复现与封账保护</a></p></section>'''
    (OUT / 'index.html').write_text(top + body + '</main></html>')
    (OUT / 'report.html').write_text(top + '<article>' + rendered((OUT / 'RESEARCH_REPORT_ZH.md').read_text()) + '</article></main></html>')
    gallery = top + '<h1>Research F3 · 实际材料、竞争likelihood与路径</h1><p>所有候选来自真实full-Maxwell求值；非后验模式图库。GT仅离线评分，材料图共用实部/损耗色标。</p>'
    panels = [
        ('01_material_weak', 'weak：基础材料与候选', '有任务尺度竞争端点，但只是同一小材料区域。native代表按Φ选择，不按GT误差选择。'),
        ('01_material_ring', 'ring：反例保留', '原基础材料已经很好；当前候选没有达到ΔΦ≤2。不能用noise-ball通过代替概率重要性。'),
        ('05_candidate_differences', '三个端点的真实材料差分', '数值不同，彼此距离小于预注册任务尺度；多尺度拓扑计数没有改变，不是三个worlds。'),
        ('06_matched_distance', '相同材料RMS的四臂比较', '普通材料retraction优于冻结切向；没有稳定高阶增量。第二目标的距离失败在原始表单独记录。'),
        ('02_RMS_likelihood', '所有端点与任务阈值', 'MH散点是proposal，不能把被拒绝的点算成保留链或teacher；没有residual-softmax概率。'),
        ('03_MH_chain', '真实保留链，包含重复状态', '冻结对角反射proposal与完整likelihood/prior接受；8步不足以评估混合或后验质量。'),
        ('04_path_obstruction', '端点之间的固定参数线段', '后续冻结的9个中点诊断，无GT、无优化。全部noise-feasible、全部不competitive；不排除其他路径。'),
    ]
    for name, title, caption in panels:
        gallery += '<section><h2>' + title + '</h2><figure><img src="figures/' + name + '.png" alt="' + html.escape(title) + '"><figcaption>' + caption + '</figcaption></figure><a href="figures/' + name + '.pdf">PDF</a></section>'
    gallery += '<section><h2>原权重诱导的商先验</h2><img src="quotient_prior.png" alt="原权重抽样与解析商密度"><p>原8维Uniform不是7维Uniform；比例边界密度为0，材料闭支持仍合法。</p><a href="QUOTIENT_PRIOR_DERIVATION.md">推导与独立数值核验</a></section>'
    (OUT / 'gallery.html').write_text(gallery + '</main></html>')

    notice = '## 最新：Research F3 商坐标与有限噪声探索\n\n[中文报告](' + REL + 'RESEARCH_REPORT_ZH.md) · [网页与实际材料图](https://migodam.github.io/b1-probabilistic-maxwell-imaging/' + REL + ')\n\n55 states / 712 RHS封账：两个新资格物理驻点，weak小幅竞争端点，无新拓扑、高阶方法无同成本增量。结束原P0局部Explorer路线，不训练NN。\n\n'
    for name in ['README.md', 'START_HERE.md']:
        path = ROOT / name
        s = path.read_text()
        if '## 最新：Research F3' not in s:
            path.write_text(notice + s.replace('## 最新：Research F 材料几何', '## 历史：Research F 材料几何'))
    path = ROOT / 'index.html'
    s = path.read_text()
    if REL not in s:
        card = '<section class="card"><h2>最新：Research F3商坐标与有限噪声探索</h2><p>55 states / 712 RHS。两个新驻点；小幅竞争端点，无新拓扑。高阶反馈无同成本增量，结束原P0 Explorer。</p><a href="' + REL + '">看报告与真实材料图 →</a></section>'
        path.write_text(s.replace('<div class="grid">', '<div class="grid">' + card, 1).replace('最新：Research F材料几何', 'Research F材料几何'))
    path = ROOT / 'research/index.html'
    s = path.read_text()
    if 'research_f3_quotient_exploration/' not in s:
        path.write_text(s.replace('<main>', '<main><section><h2>最新：Research F3商坐标与有限噪声探索</h2><p>已完成并封账；停止原P0局部Explorer，不将小变化包装成多模式。</p><a href="research_f3_quotient_exploration/">报告、材料与路径 →</a></section>', 1).replace('最新：Research F材料几何', 'Research F材料几何'))

    source = dict(
        baseline_B1='bb23534dd04983217668d301af6857f0db61885f',
        historical_F_baseline='8b620b315b2b4b212d0f1e4c143ced49f22f02e0',
        executed_first_attempt='b87413bbe06100600988ed87f69b48e21dfe3230',
        executed_recovery='4d45c3b7db2d17ca8f4c27280d75c719ed8491a4',
        executed_fixed_path_audit='bd60ee61afd5e17f2e4246582b6d7196abba1f1c',
        A30_backend='717d3d26a290248bd2b9aeb55f3d03694b407b84',
        A26_generator='2603f3bd16c48e48587fedc4363975779def17f5',
        final_changes='reporting, evidence packaging and legacy/background path portability only; no new physical evaluations',
        reused_F1_states=['state_00 and state_11 plus paid Jacobians for two exposed background conditions'],
        historical_costs_unchanged=True, source_snapshots='evidence/source_executed_attempt1; executed git commits',
        probability_semantics='candidate_only, no posterior mass estimates',
    )
    write_json(OUT / 'SOURCE_INVENTORY.json', source)
    write_json(OUT / 'DEPENDENCIES.json', dict(python=sys.version, platform=platform.platform(), numpy=__import__('numpy').__version__, scipy=__import__('scipy').__version__, matplotlib=__import__('matplotlib').__version__, GPU_training=0, source_pins='SOURCE_INVENTORY.json'))
    for bg in ['weak', 'ring']:
        case = f'p250073_{bg}_a1_snr20'
        for j in [0, 11]:
            for suffix in ['.npz', '_JAC.npz']:
                path = OLD / f'results/F1/{case}/state_{j:02}{suffix}'
                dest = OUT / 'inputs/legacy_cache' / path.relative_to(OLD)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
        for folder in ['online', 'offline_only']:
            for suffix in ['.json', '.npz']:
                path = OLD / f'inputs/{folder}/{case}{suffix}'
                dest = OUT / 'inputs/legacy_cache' / path.relative_to(OLD)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, dest)
        sys.path.insert(0, str(ROOT / 'src'))
        from prior_constrained import prior
        dest = OUT / f'inputs/background/{bg}.npz'
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(prior.INPUTS / f'background/{bg}.npz', dest)
    shutil.copy2(OLD / 'results/F1/SUMMARY.json', OUT / 'inputs/legacy_cache/results/F1/SUMMARY.json')
    write_json(OUT / 'PUBLICATION.json', dict(repository='migodam/b1-probabilistic-maxwell-imaging', new_repository=False, release=RELEASE, license_status_unchanged=True, private_conversations=False, third_party_fulltexts=False, locked_data=False, hash_checks=False, old_ledgers_unchanged=True, scientific_decision_independent_of_publication=True))
    archive = ROOT / 'publication_local/Research_F3_Quotient_Exploration_20261017_Evidence.zip'
    archive.parent.mkdir(exist_ok=True)
    excluded = {'budget.lock', 'EVIDENCE_RELEASE.json', 'PUBLICATION_RECEIPT.json', 'WEB_READBACK.json'}
    with zipfile.ZipFile(archive, 'w', zipfile.ZIP_DEFLATED, compresslevel=5) as z:
        for path in sorted(OUT.rglob('*')):
            if path.is_file() and path.name not in excluded and '__pycache__' not in path.parts:
                z.write(path, path.relative_to(ROOT))
        for path in sorted((ROOT / 'src').rglob('*.py')):
            z.write(path, path.relative_to(ROOT))
        for path in sorted((ROOT / 'research/prior_constrained_20261015/source_snapshots').glob('*')):
            if path.is_file():
                z.write(path, path.relative_to(ROOT))
        for name in ['research_f3_prior_audit.py', 'research_f3_results.py', 'publish_research_f3.py', 'publish_research_f.py', 'publish_lxz.py']:
            z.write(ROOT / 'tools' / name, 'tools/' + name)
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
    write_json(OUT / 'EVIDENCE_RELEASE.json', dict(asset=archive.name, bytes=archive.stat().st_size, files=len(names), new_endpoint_arrays=sum('/results/endpoints/' in n and n.endswith('.npz') for n in names), new_Jacobian_arrays=sum('/results/jacobians/' in n and n.endswith('.npz') for n in names), legacy_state_arrays=4, legacy_Jacobians=4, background_arrays=2, release=RELEASE, zero_new_physics=True, hash_checks=False))
    print((OUT / 'EVIDENCE_RELEASE.json').read_text())


if __name__ == '__main__':
    main()

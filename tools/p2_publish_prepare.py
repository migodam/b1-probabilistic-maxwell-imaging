"""Prepare existing-repository publication and one self-contained evidence ZIP.

This never calls a solver or training. Historical assets stay in place. Raw
arrays and dense evaluation rows go to Release, summaries to the existing Git.
"""
from pathlib import Path
import os,sys,json,zipfile,platform,html
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
C1=Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'))
from src.b1_joint.accounting import dump

def prepend(path,text):
    old=path.read_text();path.write_text(text+'\n\n'+old)

def main():
    out=ROOT/'research/joint_p2';results=out/'results';tag='joint-p2-20261014';release='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/'+tag;crelease='https://github.com/migodam/c1-physics-state-generative-imaging/releases/tag/'+tag
    # Portable C1 location for readers cloning the two existing repositories.
    for name in ['src/b1_joint/run.py','tools/p2_correct_proposals.py','tools/p2_final_analysis.py','tools/p2_accept_reference.py']:
        p=ROOT/name;s=p.read_text()
        if 'import os\n' not in s:s='import os\n'+s if not s.startswith('"""') else s.replace('from __future__ import annotations','from __future__ import annotations\nimport os',1) if 'from __future__ import annotations' in s else s.replace('from pathlib import Path','import os\nfrom pathlib import Path',1)
        s=s.replace("ROOT.parent.parent/'C1/c1-joint-p2'","Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'))")
        p.write_text(s)
    source=dict(stage='B1-C1-JOINT-P2-20261013',publication_date='2026-10-14',B1_readonly_baseline='1cd474cbf6e72c6bf9a070a78293664149a472a5',C1_readonly_baseline='9f2c51882429817db95e5801f5937257eb373b2a',new_repository=False,artifact_release=release,C1_release=crelease,raw_data_scope='only new authorized P2 synthetic observations; old locked data not read',license='existing repository license status unchanged',external_source='existing B1 vector-Mie; third-party papers linked/paraphrased, no copied full text',numerical_status='GATE_DECISION.json independent from public availability')
    dump(out/'PUBLICATION.json',source)
    import numpy,scipy,matplotlib,torch
    dump(out/'DEPENDENCIES.json',dict(python=sys.version.split()[0],platform=platform.platform(),numpy=numpy.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,torch=torch.__version__,GPU_training=0))
    ignore='\n# Joint P2 raw data and checkpoints are retained in the evidence Release.\nresearch/joint_p2/results/forward_cache/\nresearch/joint_p2/results/C1_*DATA*.npz\nresearch/joint_p2/results/FRESH_CONTINUOUS_REFERENCE*.npz\npublication_local/\n'
    with (ROOT/'.gitignore').open('a') as f:f.write(ignore)
    with (C1/'.gitignore').open('a') as f:f.write('\n# Full P2 evaluation rows and weights are published as Release evidence.\np2/results/evaluation/\np2/results/models/*.pt\np2/results/continuous/CONTINUOUS_COMPARISON.json\np2/results/continuous_partition/CONTINUOUS_COMPARISON.json\np2/results/continuous_accepted/CONTINUOUS_COMPARISON.json\npublication_local/\n')
    notice='# 最新研究：B1 × C1 Joint P2\n\n[中文报告](research/joint_p2/RESEARCH_REPORT_ZH.md) · [交互Top-K图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/joint_p2/) · [科学裁决](research/joint_p2/GATE_DECISION.md) · [完整证据下载]('+release+')\n\n同一连续材料prior保留实心/中心/环形空腔解释，提供事件质量、条件不确定性、full-Mie补偿与独立传感器更新。新增3459/4000次物理attempt；C1未达到完整连续posterior目标。下面原有内容保留为历史研究。'
    prepend(ROOT/'START_HERE.md',notice);prepend(ROOT/'README.md',notice)
    for name in ['RESEARCH_REPORT_ZH.md','RESULTS.md','DECISION.md','THEORY.md']:
        prepend(ROOT/name,'> 最新 Joint P2 的'+('理论与范围' if name=='THEORY.md' else '研究结果')+'见 [独立中文报告](research/joint_p2/RESEARCH_REPORT_ZH.md) 和 [新阶段入口](research/joint_p2/START_HERE.md)。以下历史正文保留。')
    cnotice='# 最新：B1 × C1 Joint P2\n\n[共同图库](https://migodam.github.io/b1-probabilistic-maxwell-imaging/research/joint_p2/) · [C1独立报告](p2/C1_P2_RESEARCH_REPORT_ZH.md) · [裁决](p2/C1_P2_DECISION.md) · [权重与原始记录]('+crelease+')\n\n共同四维vector-Mie posterior上，小型蒸馏MDN的事件MAE约.0132；完整协方差并未稳定胜出事件质量，连续条件CDF约.142，完整posterior目标未达到。下面P1等原结果保留。'
    prepend(C1/'START_HERE.md',cnotice);prepend(C1/'README.md',cnotice)
    for folder in ['evaluation','models','continuous','continuous_partition','continuous_accepted']:
        (C1/'p2/results'/folder/'DATA_ACCESS.md').write_text('# 原始数据获取\n\n原文件保留在研究工作目录；大数组、逐观测原始JSON和模型权重在 [P2证据Release]('+crelease+')。下载ZIP，按原相对目录解压至对应repo根目录。最终连续结果为 `continuous_accepted`；旧节点/direct/partition近似原样保留。\n')
    for p in (C1/'p2').glob('*.md'):
        s=p.read_text();s=s.replace('(results/evaluation/EVALUATION_SUMMARY.json)','('+crelease+')');p.write_text(s)
    rootpage=ROOT/'index.html';s=rootpage.read_text();s=s.replace('<div class="grid">','<div class="grid"><section class="card"><h2>最新：Joint P2材料歧义</h2><p>同一观测保留多种材料解释与真实事件质量；独立几何传感器可进一步更新判断。完整连续posterior精度尚未达到。</p><a href="research/joint_p2/">浏览候选图、概率与补偿 →</a></section>',1);s=s.replace('最新：E0 连续摘要修复','历史：E0 连续摘要修复');rootpage.write_text(s)
    overview=ROOT/'research/index.html';s=overview.read_text();s=s.replace('<main>','<main><section><h1>最新：B1 × C1 Joint P2</h1><p>材料歧义的概率保留、补偿解释与下游交付。</p><a href="joint_p2/">打开实际图库与研究报告 →</a></section>',1);overview.write_text(s)
    # Shared relative-directory archives make all replays independent of this
    # workstation path. No hashes, historical copies or new physical bank.
    local=ROOT/'publication_local/joint_p2';local.mkdir(parents=True,exist_ok=True)
    bzip=local/'B1_Joint_P2_Evidence.zip';czip=local/'C1_Joint_P2_Models_and_Evaluation.zip'
    with zipfile.ZipFile(bzip,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file() and p.name!='physics.lock':z.write(p,p.relative_to(ROOT))
        for p in sorted((ROOT/'src/b1_joint').glob('*.py')):z.write(p,p.relative_to(ROOT))
        for p in sorted((ROOT/'tools').glob('p2_*.py')):z.write(p,p.relative_to(ROOT))
    with zipfile.ZipFile(czip,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted((C1/'p2').rglob('*')):
            if p.is_file() and '__pycache__' not in p.parts:z.write(p,p.relative_to(C1))
    dump(local/'ASSET_LIST.json',dict(B1_asset=dict(name=bzip.name,bytes=bzip.stat().st_size,source_relative_root='B1 existing repo'),C1_asset=dict(name=czip.name,bytes=czip.stat().st_size,source_relative_root='C1 existing repo'),physics=0,old_results_untouched=True))
    print(json.dumps(dict(B1_evidence=str(bzip),C1_evidence=str(czip))))
if __name__=='__main__':main()

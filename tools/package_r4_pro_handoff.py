"""Small GPT Pro reading/code package; full raw evidence is a separate release."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1];S=ROOT/'research/p3l_r4_conditional_transport'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=S/'packaging/pro_handoff')
    p.add_argument('--name',default='Research_P3LR4_Results_and_Pro_Handoff.zip')
    a=p.parse_args()
    start=time.perf_counter();out=a.out
    out.mkdir(parents=True,exist_ok=False)
    files=set(S.glob('*.md'))|set((S/'GPT_PRO_HANDOFF').glob('*'))
    for pattern in ('*.json','*.csv'):files.update((S/'evidence_summary').glob(pattern))
    for folder in ('final_risk','final_probability','final_physics'):
        for pattern in ('*.png','*.pdf','*.md','*.csv'):
            files.update((S/'figures'/folder).glob(pattern))
    for name in ('SOURCE_INVENTORY.json','ENVIRONMENT.json','RUNTIME_VERSIONS.json','GATE_SUMMARY.json',
                 'RESEARCH_CLOSURE.json','RESULTS_AND_METRICS.json','RESOURCE_LEDGER.json','FAILURE_INDEX.json'):
        files.add(S/name)
    for name in ('JOINT_SIMULATION_PROPER_SCORE_INTERPRETATION.md','RISK_CONTEXT_AND_METRIC_THEORY.md',
                 'PHYSICS_TRANSPORT_THEORY_REVISION.md','GEOMETRY_MATCHED_DISTANCE_GATE_REVIEW.md',
                 'GEOMETRIC_CANDIDATE_ARCHIVE_REVIEW.md','ADAPTIVE_CFM_DENSITY_PROTOCOL.md'):
        files.add(S/'AUDIT'/name)
    # Small post-pack additions are distinct from the full source/archive snapshot.
    for folder in ('post_score_calibration_bound','delivery_layout','public_replay'):
        base=S/'AUDIT'/folder
        for pattern in ('*.md','*.json','*.csv','*.tsv'):
            files.update(base.glob(pattern))
    files.update((ROOT/'src/p3r4').glob('*.py'))
    files.update((ROOT/'src/p1_material').glob('*.py'))
    files.update((ROOT/'tools').glob('*r4*.py'))
    files={p for p in files if p.is_file()};manifest=[]
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bundle=out/a.name
    with zipfile.ZipFile(bundle,'x',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        z.writestr('READ_THIS_FIRST.md',
            '# P3L-R4 实验结果与下一轮 Pro 阅读包\n\n'+
            '先读workspace/research/p3l_r4_conditional_transport/START_HERE.md和GPT_PRO_HANDOFF/START_HERE.md。\n\n'+
            '本包包含真实研究结论、机器指标、全部关键图表和R4源码；不是全量原始物理缓存或模型权重包。'+
            '完整原始证据保存在同一B1仓库Release分卷，未用本阅读包替代。\n\n'+
            '源码版本：'+revision+'\n')
        for p in sorted(files):
            v=p.stat();name='workspace/'+str(p.relative_to(ROOT));z.write(p,name)
            manifest.append(dict(path=name,bytes=v.st_size,mtime_ns=v.st_mtime_ns))
        z.writestr('READING_MANIFEST.json',json.dumps(dict(revision=revision,files=manifest,hashes_used=False),ensure_ascii=False,indent=2))
    receipt=dict(status='PRO_READING_PACKAGE_COMPLETE',bundle=str(bundle),bytes=bundle.stat().st_size,
                 files=len(files),revision=revision,wall_seconds=time.perf_counter()-start,RHS=0,
                 full_raw_evidence=False,hashes_used=False)
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))


if __name__=='__main__':main()

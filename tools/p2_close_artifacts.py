"""Freeze actual receipts and prepare public evidence without further physics."""
from pathlib import Path
import os,sys,json,csv,re,zipfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
C1=Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'))
from tools.p2_resources import main as resources
from src.b1_joint.accounting import dump
def main():
    resources();out=ROOT/'research/joint_p2';results=out/'results';receipt=json.loads((results/'RESOURCE_SUMMARY.json').read_text());wall=receipt['controlled_job_wall_seconds']
    jobs=[json.loads(a) for a in (results/'JOB_LEDGER.jsonl').read_text().splitlines()];correction=json.loads((results/'JOB_ATTRIBUTION_CORRECTION.json').read_text());override={a['job_id']:a['actual_own_new_physics'] for a in correction['overrides']};cpu_children={a['path'].split('/')[1]:a['child_cpu_seconds'] for a in receipt['child_cpu_receipts']}
    with (results/'COST_TABLE.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['job','status','new_physics','wall_seconds','CPU_parent_seconds','CPU_child_seconds','error']);writer.writeheader()
        for a in jobs:writer.writerow(dict(job=a['id'],status=a['status'],new_physics=override.get(a['id'],a['new_physics']),wall_seconds=a['wall_seconds'],CPU_parent_seconds=a['cpu_seconds'],CPU_child_seconds=cpu_children.get(a['id'],0),error=a.get('error')))
    receipt['job_failure_count']=sum(a['status']=='FAIL' for a in jobs);receipt['CPU_parent_plus_child_seconds']=sum(a['cpu_seconds']+cpu_children.get(a['id'],0) for a in jobs);receipt['scientific_run_closed']=True;receipt['further_physics_requires_explicit_scope_and_remaining_budget']=True;dump(results/'RESOURCE_SUMMARY.json',receipt)
    p=out/'RESOURCE_SUMMARY.md';s=p.read_text();s=s.replace('受控jobwall约 **1723s**',f'最终受控jobwall **{wall:.3f}s（{wall/60:.2f}分钟）**');s=s.replace('截至交付数值作业快照','截至最终归档');s=s.replace('最终精确数值随gallery和必要artifact作业记入','已包含gallery、全部失败与必要artifact作业，精确数值记入');s+='\n\n研究数值作业已关闭。完整控制wall不含对话/模型服务等待、人工写作与公开网络传输，不能当作整个chat的历时。最终阶段总账与逐job表见 [COST_TABLE.csv](results/COST_TABLE.csv)。五次cold/warm测量覆盖proposal acquisition；重复整套full-Mie inference的五轮cold/warm没有在剩余配额内运行，记NOT_RUN，不拿毫秒proposal时间冒充完整posterior速度。\n';p.write_text(s)
    p=out/'GATE_DECISION.json';d=json.loads(p.read_text());d['resource'].update(controlled_job_wall_seconds=wall,job_failures=receipt['job_failure_count']);d['publication_date']='2026-10-14';dump(p,d)
    p=out/'THEORY.md';p.write_text(p.read_text()+'\n\n材料事件context更新的精确期望、TV误差界和latent/material KL链式分解见 [THEORY_HANDOFF.md](THEORY_HANDOFF.md)。这些标准恒等式解释交付信息需求，不被申报为新原创定理。\n')
    p=ROOT/'START_HERE.md';p.write_text(p.read_text().replace('# 最新研究入口：B1 E0','# 历史入口：B1 E0'))
    for name in ['START_HERE.md','README.md']:
        p=C1/name;s=p.read_text().replace('最新：[C1-P1','历史：[C1-P1').replace('**最新研究：[C1-P1','**历史P1研究：[C1-P1');p.write_text(s)
    # Git retains the accepted aggregate; large rows remain Release evidence.
    dump(C1/'p2/results/REFINED_FRESH_EVENT_COMPARISON.json',json.loads((results/'REFINED_FRESH_EVENT_COMPARISON.json').read_text()))
    p=out/'FAILURES.md';p.write_text(p.read_text()+'\n\n发布整理还保留`p2-positive-partition-CDF-reference`语法括号失败、accepted summary中parent字符串导致NaN及C1 metadata字段名错误的0物理cost receipts。后续分别修复；accepted参考数组已保存，汇总恢复没有重复积分或物理。服务器丢失已完成会话句柄时，以持久化job receipt为准，未再次启动已完成作业。\n')
    # Rebuild the two existing-repository evidence files from finalized records.
    local=ROOT/'publication_local/joint_p2';local.mkdir(parents=True,exist_ok=True)
    for dest,root,scope in [('B1_Joint_P2_Evidence.zip',ROOT,[out,ROOT/'src/b1_joint']),('C1_Joint_P2_Models_and_Evaluation.zip',C1,[C1/'p2'])]:
        with zipfile.ZipFile(local/dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for base in scope:
                for p in sorted(base.rglob('*')):
                    if p.is_file() and '__pycache__' not in p.parts and p.name!='physics.lock':z.write(p,p.relative_to(root))
            if root==ROOT:
                for p in sorted((ROOT/'tools').glob('p2_*.py')):z.write(p,p.relative_to(ROOT))
    print(json.dumps(dict(physics=receipt['used'],wall=wall,release_assets=[a.name for a in local.glob('*.zip')])))
if __name__=='__main__':main()

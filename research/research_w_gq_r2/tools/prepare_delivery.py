"""Local reproducible source/input bundle. No physics or hash checks."""
from pathlib import Path
import json,shutil,sys,platform,importlib.metadata
from r2.accounting import BudgetBook
ROOT=Path(__file__).resolve().parents[1]
def main():
    old=ROOT.parent/'research_w_gq';reuse=ROOT/'reuse/previous_gq'
    for area,names in [('src/gq',['__init__.py','physics.py','contrast.py']),('tools',['p3_bank.py'])]:
        for name in names:
            dest=reuse/area/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(old/area/name,dest)
    shutil.copytree(old/'vendor',ROOT/'reuse/fixed_vendor',dirs_exist_ok=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    shutil.copytree(old/'inputs/old_finite',ROOT/'inputs/old_finite',dirs_exist_ok=True)
    for name in ['RESEARCH_W_GQ_REPORT_ZH.md','GPT_PRO_HANDOFF.md','GATE_DECISION.json','COST_SUMMARY.json']:
        dest=ROOT/'source_review/previous_local_gq'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(old/name,dest)
    versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','scikit-fem','triangle','matplotlib']}
    env=dict(python=sys.version,executable=sys.executable,platform=platform.platform(),machine=platform.machine(),dependencies=versions,BLAS_threads=1,used_GPU=False,old_fixed_commit='73ba6eec3746548bf21259f627b635678709df00',new_commit='NOT_CREATED',inputs_only_exposed_development=True,no_hash_checks=True)
    (ROOT/'ENVIRONMENT.json').write_text(json.dumps(env,indent=2));(ROOT/'requirements.txt').write_text('\n'.join(f'{n}=={v}' for n,v in versions.items())+'\n')
    (ROOT/'SOURCE_REUSE.json').write_text(json.dumps(dict(reused_source=str(reuse.relative_to(ROOT)),fixed_vendor='reuse/fixed_vendor',old_local_GQ_not_committed=True,old_prior='inputs/old_finite',purpose='required standalone reproduction inputs, not a backup workflow',source_files=[str(p.relative_to(ROOT)) for p in (ROOT/'reuse').rglob('*.py')]),indent=2))
    b=BudgetBook();marks=ROOT/'results/SCIENTIFIC_FAILURES_RECORDED.json'
    if not marks.exists():
        b.failure('E1','unknown macro finite imaging','MODEL_MISMATCH: A2 position RMSE .683435; not accurate online imaging')
        b.failure('E1','continuous FWI trajectories','All90 actual calls retained, nine capped starts NONCONVERGED; objective-selected ordinary joint branch remains wrong')
        b.failure('E3','original Hermitian low-rank projection','COMPRESSION_FAILURE: rank32 maxwhite .23094, preserved before transpose correction')
        b.failure('E3','conditional consensus','NOT_CERTIFIED: old bank sampled/full conditional diameter near1; physical micro groups have five TV>.02 tails')
        marks.write_text(json.dumps(dict(recorded=True,science_failures_not_runtime_exceptions=True)))
    (ROOT/'PHYSICS_PROTOCOL.md').write_text('本轮完整物理与统计合同见 CONTRACT.json；源/载荷、网格和边界实现见 docs/PHYSICS_API.md，结果范围见 REPORT_ZH.md。\n')
    print(json.dumps(dict(source_files=len(list((ROOT/'reuse').rglob('*.py'))),dependencies=versions,physical_new_calls=0)))
if __name__=='__main__':main()

"""Export the fixed completed observation000 as a self-contained path mirror.

Copies paid evidence bytes only, no physics/Book, inference or GT array access.
The native schema is unchanged: referenced paths resolve under the export root.
"""
from pathlib import Path
import datetime as dt
import json
import shutil

ROOT=Path(__file__).resolve().parents[1]

def main():
    gate=json.loads((ROOT/'results/P1B/offline_scores/COMPLETE.json').read_text())
    if gate.get('status')!='OFFLINE_SCORES_READY_FOR_PARENT_REVIEW':
        raise RuntimeError('Offline completion required before public export')
    source=ROOT/'results/P1B/handoff/observation_000/EXAMPLE.json'
    example=json.loads(source.read_text())
    if example.get('observation')!=0 or example.get('truth_access') is not False:
        raise RuntimeError('Only the frozen non-GT observation000 package is exported')
    refs=set()
    def visit(x):
        if isinstance(x,str) and x.startswith(('results/','inputs/','configs/')):
            p=ROOT/x
            if p.is_file():refs.add(x)
        elif isinstance(x,dict):
            for v in x.values():visit(v)
        elif isinstance(x,list):
            for v in x:visit(v)
    visit(example)
    refs.add('results/P1B/handoff/observation_000/EXAMPLE.json')
    refs.update(['configs/W_CONTROLLED_PRIOR.json','configs/P1_PROBABILITY_PROTOCOL.json',
                 'configs/ORIGINAL_P0.json','inputs/W_FAMILY_WALLS.npz'])
    out=ROOT/'MaterialHypothesisPackage'
    if out.exists():raise RuntimeError('Public export exists; inspect instead of overwriting')
    total=0
    for rel in sorted(refs):
        p=ROOT/rel
        if p.is_symlink() or not p.resolve().is_relative_to(ROOT.resolve()):raise RuntimeError('Unsafe reference')
        if 'TRUTH' in rel.upper() or 'GENERATION_PLAN' in rel.upper():raise RuntimeError('No GT/plans in hypothesis export')
        size=p.stat().st_size;total+=size
        if size>32_000_000 or total>128_000_000:raise RuntimeError('Unexpected size; inspect reference before publication')
    out.mkdir()
    copied=[]
    for rel in sorted(refs):
        d=out/rel;d.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/rel,d)
        copied.append(dict(source=rel,public_path=rel,bytes=d.stat().st_size,copy='whole-file byte stream'))
    # Main entry is also at export root for a readable, stable machine URL.
    shutil.copyfile(source,out/'EXAMPLE.json')
    manifest=dict(schema='W_PUBLIC_PATH_MIRROR_V1',utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        native_schema_preserved=True,entry='EXAMPLE.json',path_resolution='relative to this export root',
        fixed_example=0,outside_family_mass=None,original_minimum_feature='UNRESOLVED',
        new_physics=0,GT_read=False,files=copied,total_referenced_bytes=total,
        source_revision='Per-job frozen SOURCE_SNAPSHOT/CONFIG_SNAPSHOT; final publication revision separate',
        scope='Two fixed public wall fields; continuous branch-event masses with finite joint particles; not general unknown-wall posterior')
    (out/'MANIFEST.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    (out/'README.md').write_text('''# 固定 observation000 的概率材料交接包

入口 `EXAMPLE.json` 保留原生 W_PROBABILITY_HANDOFF_V1。所有文件引用相对本目录，按原路径镜像；MANIFEST.json列出逐文件字节数与来源。包含两个事件定义、连续参考事件质量、prior/后验odds、有限条件不确定性、完整实际Full响应、材料数组、联合有限权重和8份SMC轨迹。只有最终粒子权重属于最终有限posterior，不能把不同source阶段的轨迹当同一分布样本。

这是两个固定公开墙体及相关目标位置带的proper prior；家族外质量未量化。材料均值和位置均值是Bayes action，不是第三个合法世界。参数化/表示和shared-gain规则见配置；精度资格是经验求积交叉验证，不是连续PDE或未知支持全覆盖证书。GT未进入本交接包；GT只在单独离线图库与评分证据中出现。

完整campaign所有观测与粒子、其余原始物理响应及第一joint witness的全频复场见同一仓库的 W Release。此网页小包不声称包含全部本地原始场。新增物理0；所有响应原始费用保留在根账本。
''')
    print(json.dumps(dict(status='PUBLIC_HYPOTHESIS_COPY_COMPLETE',files=len(copied),bytes=total,new_physics=0)))

if __name__=='__main__':main()

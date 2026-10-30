"""One-time documentation append for R12; no physics, scoring or tests."""
from pathlib import Path
import csv
import json
import re

ROOT = Path(__file__).resolve().parents[3]
AT = ROOT / 'research/maxwell_formula_atlas'
CHAPTER = AT / '12_BILATERAL_MATERIAL_MODES_AND_FIXED_POINTS.md'

text = CHAPTER.read_text()
entries = [dict(id=m[1], title=m[2], chapter=CHAPTER.name,
                anchor=m[1].lower(),
                source='R12 explicit finite-model derivations; per-item source and validity in chapter12')
           for m in re.finditer(r'^## (F-BF-\d+) · (.+)$', text, re.M)]
registry = json.loads((AT / 'FORMULA_REGISTRY.json').read_text())
registry['version'] = '2026-10-30 R12'
registry['formulas'].extend(entries)
registry['count'] = len(registry['formulas'])
(AT / 'FORMULA_REGISTRY.json').write_text(json.dumps(registry, ensure_ascii=False, indent=2) + '\n')

combined = (AT / 'FORMULA_ATLAS.md').read_text()
combined = combined.replace('P1—R11主链', 'P1—R12主链', 1)
(AT / 'FORMULA_ATLAS.md').write_text(combined.rstrip() + '\n\n---\n\n' + text)

changelog = AT / 'CHANGELOG.md'
with changelog.open('a') as handle:
    handle.write('\n2026-10-30 R12：追加 F-BF-001—009，当前双边材料导数与端口对、原 Gaussianized66 实白化拉回、材料局部商与原先验增广、同J嵌套SVD冗余、完整非根固定点导数、数据/MAP/RML区别、固定弱坐标及有限secant经验边界。已有158条稳定编号与历史源文保留，总计167条。新增Maxwell RHS=0、LU=0；不新增数值恒等式重放、评分、测试或文献调用。P1现有八臂覆盖均0/4、在线选择仍为原锚点，N1边界失败及文献Tier C身份范围原样保留。\n')

sources = [
    ('shared/physics.py', 'Physics.evaluate and BILATERAL_MODE_DERIVATION section1', '', ''),
    ('agent_a/modes.py', 'bilateral_gram / pair_modes', 61, 92),
    ('shared/physics.py', 'white_jac / Problem.decode; BILATERAL_MODE_DERIVATION section4', 35, 74),
    ('agent_a/modes.py', 'material_canonical_basis / proposal_basis', 37, 196),
    ('NESTED_SVD_REDUNDANCY_AUDIT.md', 'same-J singular-index restriction', '', ''),
    ('FIXED_POINT_SPECTRUM_AUDIT.md', 'sections1,3,4 full nonroot derivative and coordinate caveats', '', ''),
    ('FIXED_POINT_SPECTRUM_AUDIT.md', 'section2 data damping versus MAP/RML', '', ''),
    ('agent_b/engine.py', 'weak_frame / transverse_correct', 89, 380),
    ('agent_c/p1_diagnostics/FINITE_SECANTS.json', 'existing finite-node diagnostic; no recomputation', '', ''),
]
with (AT / 'SOURCE_FORMULA_INDEX.csv').open('a', newline='') as handle:
    writer = csv.writer(handle)
    for entry, (path, section, start, end) in zip(entries, sources):
        writer.writerow([entry['id'], entry['title'],
                         'research/p3l_r12_bilateral_fixedpoint/' + path,
                         section, start, end,
                         'R12 added derivation; existing 158 source rows preserved',
                         'DERIVED_FINITE_MODEL_WITH_DECLARED_IMPLEMENTATION_OR_EXISTING_EVIDENCE',
                         'Original Gaussianized66; full-grid RMS576; standard identity not novelty; finite DEV evidence only'])

print(json.dumps(dict(version=registry['version'], formulas=registry['count'],
                      chapter=CHAPTER.name, new_physics_RHS=0,
                      new_physics_LU=0, new_numerical_runs=0)))

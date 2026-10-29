"""Public dependency metadata inventory. Never opens NPZ/truth or scientific runners."""
from pathlib import Path
import csv,json,time,sys,platform,importlib.metadata
P=Path(__file__).resolve().parents[1];ROOT=P.parents[4]
B=ROOT/'Gaussian/B1';R4=B/'b1-p3l-r4-conditional-transport'
R6=B/'b1-p3l-r6-task-optimal';R7=B/'b1-p3l-r7-material-fibers'
A30=ROOT/'Gaussian/A30/a30-maxwell-material-manifold'
OUT=P/'audit/reproducibility'
DENY={'private','offline_truth','offline_only','.git','__pycache__','.venv','.venv_nn'}
rows={}
def add(path,branches,role,required='YES',mapping='PRESERVE_GAUSSIAN_SIBLING_LAYOUT'):
    path=Path(path)
    if any(x in DENY for x in path.parts):raise RuntimeError('Private path forbidden in inventory')
    rel=str(path.relative_to(ROOT))
    if rel in rows:
        old=rows[rel];old['branches']=';'.join(sorted(set(old['branches'].split(';'))|set(branches.split(';'))));return
    stat=path.stat() if path.exists() else None
    rows[rel]=dict(original_absolute_path=str(path),archive_path=rel,branches=branches,
        role=role,required=required,exists=path.exists(),bytes=None if stat is None else stat.st_size,
        content_opened=False if path.suffix=='.npz' else 'SOURCE_OR_METADATA_ONLY',
        relocation=mapping,release_tier='FULL_PUBLIC_EVIDENCE_AND_REPLAY')
def tree(path,branches,role):
    path=Path(path)
    if not path.exists():add(path,branches,role);return
    for q in sorted(path.rglob('*')):
        if q.is_file() and not q.is_symlink() and not any(x in DENY for x in q.parts):add(q,branches,role)

# Original P1, target/acquisition, and code imported by the historical Target.
for name in ['src/p1_material','src/p3r4']:
    tree(R4/name,'E0;E1;E2;E3;E3_REV1;E4;E5','Original P1/core source closure')
for name in ['research/p3l_r4_conditional_transport/configs/TARGET.json',
             'research/p3l_r4_conditional_transport/inputs/background/weak.npz',
             'research/p1_material_target_20261018/inputs/geometry.json']:
    add(R4/name,'E0;E1;E2;E3;E3_REV1;E4;E5','Frozen target, geometry or public Green/acquisition arrays')
# Historical physics module imports its own sibling file paths; preserve names.
for name in ['physics.py','material_jets.py','decoder_geometry.py','run_e2.py']:
    add(R7/'research/p3l_r7_material_solution_fibers/src'/name,'E0;E1;E2;E3;E3_REV1','Imported historical R7 source')
for i in range(8):
    hp=R7/'research/p3l_r7_material_solution_fibers/E2/parents'/f'R7_LOCKED_{i:03d}'
    for name in ['ANCHOR.npz','BASE_GEOMETRY.npz','OBSERVATION.npz']:
        add(hp/name,'E0;E1;E2;E3','Exposed DEV raw anchor/geometry/observation; no truth')
    for radius in [.01,.02]:
        add(hp/'B_FULLJ'/f'RMS_{radius:.4f}'/'RESULT.json','E3','Historical native reproducibility reference')
    for t in [0,.25,.5,.75,1]:
        add(hp/'paths/B_FULLJ_CHORD'/f'NODE_{t:.2f}.npz','E1','Exposed historical path/J arrays')
    add(hp/'B_FULLJ/RMS_0.0200/ENDPOINT.npz','STABLE_ROI','Matched fixed-state exact-ROI benchmark input')
for i in range(4):
    add(R7/'research/p3l_r7_material_solution_fibers/E2/parents'/f'R6_CDEV_{i:04d}'/'ANCHOR.npz',
        'E0;E3','Original DEV4 global weak-frame source')
for name in ['OBSERVATIONS.npz','RESULTS.json','raw/REGULAR_material_GH21.npz','raw/PURE_OVERLAP_material_GH21.npz']:
    add(R6/'research/p3l_r6_task_optimal_maxwell/track_a'/name,'E4','Controlled exposed R6 observation/reference imports')
# Direct576 historical adapter transitive closure includes A24 -> vendored A2.
for name in ['src/a30','vendor/frozen_a29/src/a24','vendor/frozen_a29/src/a20','vendor/frozen_a29/vendor/a2']:
    tree(A30/name,'E2;E3;E3_REV1','Unchanged A30/A24/A2 CPU solver closure')
# Exact source citations for formulas, including old F/F3 and R6 anchor rules.
index=P/'formula_collection/SOURCE_FORMULA_INDEX.csv'
if index.exists():
    with index.open() as f:
        for r in csv.DictReader(f):
            if r.get('source_path'):add(ROOT/r['source_path'],'FORMULA_ATLAS','Exact cited source',required='FORMULA_REPLAY')
# E5 runtime directly needs original prior/background and this already-paid cache.
add(P/'results/stable_roi/SETUP.npz','E5','Exact nonuniform ROI setup; already-paid264 full RHS, no recharge',
    mapping='PROTOCOL_PATH_REBASE_AND_CONTENT_IDENTITY')
# E4 has another independent scalar-family Schur cache.
add(P/'results/probability/transfer_setup.npz','E4','Scalar-amplitude Schur cache; scope differs from arbitrary nonuniform ROI')

OUT.mkdir(parents=True,exist_ok=True)
ordered=[rows[k] for k in sorted(rows)]
with (OUT/'EXTERNAL_DEPENDENCY_INVENTORY.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(ordered[0]));w.writeheader();w.writerows(ordered)
branch_expected={
'E0':['audit/exploration/E0_SOURCE_AUDIT.json'],
'E1':['audit/exploration/E1_HISTORICAL_SUBSPACE_STABILITY.json'],
'E2_E3':['results/exploration/SUMMARY.json','results/exploration/E3_ANALYSIS.json','results/exploration/RESOURCE_LEDGER.jsonl'],
'E3_REV1':['results/exploration_rev1/SUMMARY.json','results/exploration_rev1/RESOURCE_LEDGER.jsonl'],
'STABLE_ROI':['results/stable_roi/RESULTS.json','results/stable_roi/SETUP.npz','results/stable_roi/COST_LEDGER.jsonl'],
'E4':['results/probability/STATUS.json','results/probability/final_numerical_comparison_v3.json','results/probability/original1d_results_v2.json','results/probability/aggregate_costs.json','audit/probability/ledger.jsonl'],
'E5':['results/possibility_volume/GENERATION_SEALED.json','results/possibility_volume/PREDICTIONS_SEALED.json','results/possibility_volume/offline_statistics/SCORING_COMPLETE.json','results/possibility_volume/offline_statistics/SUMMARY.json']}
status={k:[dict(path=v,exists=(P/v).exists()) for v in vv] for k,vv in branch_expected.items()}
runtime={}
for name in ['numpy','scipy','matplotlib','threadpoolctl']:
    try:runtime[name]=importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:runtime[name]=None
payload=dict(created_utc=time.time(),scope='PUBLIC metadata and source dependency inventory; E5 may still be active',
    no_private_or_truth_content_read=True,npz_contents_opened=False,new_Maxwell_RHS=0,
    dependency_count=len(ordered),missing_dependencies=[r for r in ordered if not r['exists']],
    runtime=dict(python_version=sys.version,invoked_executable=sys.executable,
                 required_project_venv_path=str(ROOT/'Gaussian/.venv_nn/bin/python'),
                 venv_symlink_resolved=False,platform=platform.platform(),package_versions=runtime),
    expected_outputs=status,
    archive_layout='Archive root contains Gaussian/B1 and Gaussian/A30 sibling trees; keep existing repository directory names',
    public_study_tree=str(P.relative_to(ROOT)),
    public_exclusions=sorted({'.git','__pycache__','.venv','.venv_nn'}),
    default_online_only_relocation_exclusions=sorted(DENY),
    full_synthetic_archive_offline_evidence='audit/reproducibility/OFFLINE_SYNTHETIC_EVIDENCE_INVENTORY.csv; include with role OFFLINE_ONLY after all32 seal/score; not personal data',
    transport_caveat='Absolute runtime paths are handled only in a NEW transport derivative by relocate_r9_replay.py; original identities retained.',
    private_regeneration='Public seed protocol permits future explicitly authorized replication. This inventory neither reads nor regenerates private truth.')
(P/'audit/reproducibility/DEPENDENCY_STATUS.json').write_text(json.dumps(payload,indent=2)+'\n')
(P/'audit/reproducibility/PYTHON_REQUIREMENTS_OBSERVED.txt').write_text(
    '# Observed installed packages, not a new installation or universal lockfile\n'+
    '\n'.join(f'{k}=={v}' for k,v in runtime.items() if v)+'\n')
print(json.dumps(dict(dependency_count=len(ordered),missing=len(payload['missing_dependencies']),
    E5_prediction_seal=status['E5'][1]['exists'],E5_scoring_complete=status['E5'][2]['exists'],new_Maxwell_RHS=0)))

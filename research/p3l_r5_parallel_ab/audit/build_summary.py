"""Build compact read-only closeout from actual A/B receipts; no physics."""
from pathlib import Path
from datetime import datetime,timezone
import json,csv,hashlib,zipfile

P=Path(__file__).resolve().parents[1]
def read(f):return json.loads((P/f).read_text())
def put(f,x):(P/f).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

a=read('stage_a_response_coordinates/COST_LEDGER_FINAL.json')
b=read('stage_b_material_risk/RESULTS.json')
bc=read('stage_b_material_risk/COST_LEDGER.json')
cost=read('stage_b_material_risk/OPTIMIZED_COST_BENCHMARK.json')
fam={x['family']:x for x in b['summary']}
base=fam['strong_CNN']['seed_mean_parent_MSE'];old=fam['original_K32_frozen']['seed_mean_parent_MSE']
student=fam['distilled_w128_n32768']['seed_mean_parent_MSE']
rows=list(csv.DictReader((P/'stage_a_response_coordinates/OBSERVATION_SUMMARY.csv').open()))
metric=dict(campaign='P3L-R5 parallel A/B',status='CLOSED',date_local='2026-10-26',
 target='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018',historical_source_commit='85b7bcbf39f0eeab7062ecce17531538151082bf',
 A=dict(physical_observations=len(rows),independent_parents=6,noise_levels=[1,.5,.25],
        statistical_target='normalized truncated 2D submodels; NOT full66 marginal',
        numerically_converged_reference_moment_checks=sum(r['numerical_reference_moment_gate']=='True' for r in rows),
        coverage='NOT_CERTIFIED',global_chart_TV='NOT_CERTIFIED',full66='NOT_RUN',
        nonlinear_RTO='LOCAL_NUMERICAL_ROOTS_ONLY_GLOBAL_NORMALIZATION_UNQUALIFIED',
        ordinary_spline='NOT_RUN',manifold_lifting='NOT_RUN',cost=a,
        seam=read('audit/seam_maxwell/RESULTS.json'),
        curvature_probe=read('audit/curvature_bounds/RESULTS.json')),
 B=dict(statistical_target='supervised risk action and functional-context distillation; not calibrated posterior',
        training_joint_parents=32768,validation_joint_parents=2048,expensive_context_train=1024,
        evaluation_parents=256,evaluation_role='EXPOSED_VALIDATION_USED_FOR_SELECTION',
        new_locked_parents=0,new_locked_status='RESOURCE_LIMITED_NOT_RUN',new_physics_RHS=0,
        strong_CNN_MSE=base,Gaussian_K32_MSE=old,deterministic_matched_MSE=fam['deterministic_matched']['seed_mean_parent_MSE'],
        compact_MSE=student,retained_incremental_gain=(base-student)/(base-old),
        retention_exploratory_98_333_percent_interval=[.7160,1.3563],
        original_provider_parameters=682943,compact_provider_parameters=214272,
        risk_head_parameters=6418,training_seed_count=3,
        CPU_MPS_replay=read('stage_b_material_risk/INFERENCE_REPLAY.json'),
        probability_calibration='NOT_CLAIMED',Maxwell_specific_gain='NOT_ESTABLISHED',
        cost_source='stage_b_material_risk/OPTIMIZED_COST_BENCHMARK.json',
        cost=bc,summary=b['summary']),
 decisions=dict(A='REDIRECT',B1='CONDITIONAL_GO_DEVELOPMENT_ONLY',
                B2='STOP_CURRENT_LINEAR_CONTEXT',AB='NOT_FUSED',strong_IEEE_TAP='NOT_YET_SUPPORTED'))
put('RESULTS_AND_METRICS.json',metric)
put('GATE_SUMMARY.json',metric['decisions'])
put('RESOURCE_LEDGER.json',dict(status='CLOSED',physical_ledger='shared/cost_ledger/LEDGER.jsonl',
   A=a,B_training_cost=bc,B_physics=read('shared/cost_ledger/CLOSURE.json')['stage_B_physics'],
   historical_cost='Reused, not billed as new; cold-source generation/anchor costs remain separately disclosed.',
   delivery_cost='Separate DELIVERY_LEDGER.jsonl; never added to physics savings.',
   no_budget_transfer=True,unused_A_states=32768-a['states'],unused_B_confirmation_RHS=32768))
# Root journal indexes authoritative immutable sources rather than double billing.
journal=P/'RESOURCE_LEDGER.jsonl'
if journal.stat().st_size==0:
    events=[dict(event='AUTHORITATIVE_SOURCE',role='physical',path='shared/cost_ledger/LEDGER.jsonl'),
            dict(event='AUTHORITATIVE_SOURCE',role='training',path='stage_b_material_risk/COST_LEDGER.json'),
            dict(event='CLOSEOUT_SUMMARY_NOT_ADDITIONAL_CHARGE',A=a,B_new_physics_RHS=0,
                 GPU_fit_seconds=bc['GPU_training_fit_seconds'],
                 GPU_conservative_job_seconds=bc['GPU_reserved_whole_job_seconds_conservative'])]
    with journal.open('a') as f:
        for x in events:f.write(json.dumps(x,ensure_ascii=False)+'\n')

configs=[]
for f in sorted(P.rglob('*.json')):
    if any(part in ('__pycache__','v1_failed_preserved') for part in f.parts):continue
    if (f.name.startswith(('FROZEN','PROTOCOL','SELECTION')) or f.name in ('DISTILL_REVISION_PROTOCOL.json','DISTILL_SCALE_RELEASE.json','DISTILL_FINAL_CAPACITY_PROTOCOL.json','NESTED_RRR_PROTOCOL.json')) and f.stat().st_size<100000:
        try:v=json.loads(f.read_text())
        except (ValueError,UnicodeError):continue
        raw=json.dumps(v,sort_keys=True,separators=(',',':'),ensure_ascii=True).encode()
        configs.append(dict(path=str(f.relative_to(P)),config_id_blake2b=hashlib.blake2b(raw,digest_size=16).hexdigest()))
put('CONFIG_IDENTITIES.json',dict(method='canonical JSON BLAKE2b-128; no data SHA-256 checks',configs=configs))
zp=Path("/Volumes/migodam's-external-brain/archive/user-home/downloads/Research_P3LR5_Complete_Research_and_Codex_Handoff.zip")
with zipfile.ZipFile(zp) as z:
    entries=[dict(path=x.filename,bytes=x.file_size) for x in z.infolist() if not x.is_dir()]
put('SOURCE_INVENTORY.json',dict(theory_input=dict(filename=zp.name,bytes=zp.stat().st_size,entries=entries,
    reading='Formal report text/source read; selected equation pages rendered; mathematical notes and references read; no input SHA manifest verification'),
    R4_source='85b7bcbf39f0eeab7062ecce17531538151082bf',
    historical_assets='audit/HISTORICAL_ASSET_AUDIT.json',
    dependencies='audit/EVIDENCE_INPUT_REFERENCES.json',
    code_snapshot='Final Git science revision plus archived as-run versions; publication receipt records commit',
    no_historical_edits=True))
print(json.dumps(dict(status='CLOSED',observations=len(rows),A_RHS=a['all_rhs'],B_MSE=student,
                     risk_retention=(base-student)/(base-old),configs=len(configs))))

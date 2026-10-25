"""Stat-only R5 inventory and exact external input references. No digest/payload scan."""
from pathlib import Path
from collections import defaultdict
import os,json,time,subprocess
ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parents[1];AUDIT=OUT/'audit'
HIST=ROOT.parent/'b1-p3l-r4-conditional-transport';R4=HIST/'research/p3l_r4_conditional_transport'
A30=ROOT.parents[1]/'A30/a30-maxwell-material-manifold'
SELF={'NEW_R5_FILE_INVENTORY.tsv','NEW_R5_INVENTORY.json','EVIDENCE_INPUT_REFERENCES.json','REPRO_PATH_LAYOUT.md'}
def write(p,data):p.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
def git(p,*args):
 try:return subprocess.check_output(['git','-C',str(p),*args],text=True,stderr=subprocess.STDOUT).strip()
 except Exception as e:return repr(e)
def ref(path,role,needed_for,kind='file'):
 path=Path(path);r=dict(path=str(path),role=role,needed_for=needed_for,kind=kind,exists=path.exists())
 if path.is_file():s=path.stat();r.update(bytes=s.st_size,mtime=s.st_mtime)
 elif path.is_dir():r.update(read_scope='directory reference only; no recursive payload scan')
 return r

def main():
 start=time.time();ext=defaultdict(lambda:dict(files=0,bytes=0));groups=defaultdict(lambda:dict(files=0,bytes=0));files=[];links=[]
 for folder,dirs,names in os.walk(OUT,followlinks=False):
  dirs[:]=sorted(d for d in dirs if d not in ('__pycache__','.git','.DS_Store'))
  for name in sorted(names):
   path=Path(folder)/name;rel=path.relative_to(OUT)
   if name=='.DS_Store' or (path.parent==AUDIT and name in SELF):continue
   if path.is_symlink():links.append(dict(path=str(rel),target=os.readlink(path)));continue
   s=path.stat();suffix=path.suffix.lower() or '<none>';group=rel.parts[0]
   for target in (ext[suffix],groups[group]):target['files']+=1;target['bytes']+=s.st_size
   files.append((str(rel),s.st_size,s.st_mtime,suffix))
 files.sort()
 with (AUDIT/'NEW_R5_FILE_INVENTORY.tsv').open('w') as f:
  f.write('relative_path\tbytes\tmtime_epoch\textension\n')
  for path,size,mtime,suffix in files:f.write(f'{path}\t{size}\t{mtime}\t{suffix}\n')
 raw={k:v for k,v in ext.items() if k in ('.npz','.npy','.pt','.jsonl','.csv','.tsv','.log')}
 write(AUDIT/'NEW_R5_INVENTORY.json',dict(status='STAT_ONLY_SNAPSHOT',scan_start_utc_epoch=start,scan_end_utc_epoch=time.time(),root=str(OUT),regular_files=len(files),regular_bytes=sum(x[1] for x in files),by_extension=dict(sorted(ext.items())),by_top_level=dict(sorted(groups.items())),raw_data_extensions=raw,raw_data_files=sum(x['files'] for x in raw.values()),raw_data_bytes=sum(x['bytes'] for x in raw.values()),symlinks_not_followed=links,excluded=['__pycache__','.git','.DS_Store',*sorted(SELF)],manifest='audit/NEW_R5_FILE_INVENTORY.tsv',no_SHA256=True,no_payload_reads=True,concurrency_note='Parent may add final reports/package metadata after this snapshot; packaging should enumerate actual final payload again. No historical tree was copied.'))
 references=[
  ref(R4/'inputs/background/weak.npz','same paid known-background TM24 arrays; no geometry rebuild',['A_physics']),
  ref(ROOT/'research/p1_material_target_20261018/configs/ONLINE_TARGET.json','authoritative original P1 target/prior/noise',['A_physics','B_preparation']),
  ref(ROOT/'research/p1_material_target_20261018/inputs/geometry.json','original576 points for P1 MaterialPrior',['A_physics','B_preparation']),
  ref(ROOT/'src/p1_material/prior.py','unchanged P1 decoder source',['A_physics','B_preparation']),
  ref(ROOT/'src/p1_material/target.py','authoritative original packing/noise reference',['provenance']),
  ref(A30/'src','A30 adapter source tree',['A_physics'],'directory'),
  ref(A30/'vendor/frozen_a29/src','frozen A29/A24/A2 backend source tree imported by A30',['A_physics'],'directory'),
  ref(R4/'models/density/extended_diag_gaussian_s2026102501/BEST.pt','frozenGaussian682943parameter provider',['B_original_K32_inference','B_deterministic_inference','B_teacher_regeneration']),
  ref(R4/'inputs/joint/train32768.npz','historical32768 joint training observations/material labels; no newTRAIN generation',['B_retraining_or_teacher_regeneration']),
  ref(R4/'inputs/joint/validation.npz','historical2048 joint validation parents',['B_retraining_or_teacher_regeneration']),
  ref(R4/'inputs/strong_context_1024_256/anchor_train.npz','original96+64strong optimized1024TRAIN context',['B_retraining','B_nested_OOF']),
  ref(R4/'inputs/strong_context_1024_256/anchor_validation.npz','original96+64strong optimized256VAL context',['B_inference_replay']),
  ref(HIST/'src','historical R4/R2 risk, density, feature and compatibility source imports',['B_all_inference'],'directory'),
 ]
 for seed in (2026102501,2026102502,2026102503):
  references.extend([ref(R4/f'models/risk/strongrisk_spectral_s{seed}/BEST.pt','frozen5842parameter strong reference head',['B_reference_comparison','B_RRR_inference']),ref(R4/f'models/risk/riskfusion_diag_gaussian_s{seed}/BEST.pt','frozen6418parameter originalK32 fusion head',['B_original_K32_inference'])])
 local=[]
 for name in ('prepared_train.npz','prepared_validation.npz','geometry.pt','ASSET_INVENTORY.json','inference.py','replay_inference.py','INFERENCE_REPLAY.json','SELECTION_FROZEN.json','OPTIMIZED_COST_BENCHMARK.json'):
  local.append(ref(OUT/'stage_b_material_risk'/name,'newR5 prepared context/API/replay or audit artifact',['B_replay']))
 for seed in (2026102601,2026102602,2026102603):
  for name in ('PROVIDER.pt','HEAD.pt'):
   local.append(ref(OUT/f'stage_b_material_risk/learned_context/distilled_w128_n32768_s{seed}/{name}','newR5 width128 finalized student/provider head',['B_student_inference']))
  local.append(ref(OUT/f'stage_b_material_risk/experiments/deterministic_matched_s{seed}/BEST.pt','newR5 matched deterministic fusion head',['B_deterministic_inference']))
 replay=json.loads((OUT/'stage_b_material_risk/INFERENCE_REPLAY.json').read_text())
 write(AUDIT/'EVIDENCE_INPUT_REFERENCES.json',dict(status='EXACT_PATH_SIZE_REFERENCE_MANIFEST',no_SHA256=True,no_historical_payload_copy=True,indispensable_historical_references=references,new_R5_required_artifacts=local,source_identity=dict(R5_base_commit=git(ROOT,'rev-parse','HEAD'),R4_historical_source_commit='85b7bcbf39f0eeab7062ecce17531538151082bf',R4_checkout_HEAD=git(HIST,'rev-parse','HEAD'),R4_source_vs_delivery_diff=git(ROOT,'diff','--name-only','85b7bcbf','092da74e'),A30_checkout_HEAD=git(A30,'rev-parse','HEAD'),A30_embedded_frozen_SOURCE_COMMIT='5d5e31fcfc56ac8ef9d64c46b4b7bbdedef96bb4',A30_relevant_source_worktree_status=git(A30,'status','--short','--','src/a30/physics.py','vendor/frozen_a29/src')),path_layout='audit/REPRO_PATH_LAYOUT.md',B_entrypoint='stage_b_material_risk/inference.py:Action(method,seed,device=cpu).predict(y,online,parents,batch_size=64)',B_replay_script='stage_b_material_risk/replay_inference.py',B_replay_receipt_status=replay['status'],B_replay_declared_tolerance=replay['tolerance'],B_replay_max_abs=max(r['max_abs_prediction_difference'] for r in replay['records']),B_replay_note='Preserve actual PASS/FAIL receipt and CPU/MPS drift; no automatic tolerance widening.',scope='Reference only; a filesystem path is not a public download URL. Parent release must link actual external input assets or document retrieval separately.'))
 (AUDIT/'REPRO_PATH_LAYOUT.md').write_text('''# Reproduction path layout\n\nKeep the existing sibling layout under Gaussian/: B1/b1-p3l-r5-parallel-ab (newR5), B1/b1-p3l-r4-conditional-transport (historical R4), and A30/a30-maxwell-material-manifold (frozen paid solver). The R5 shared adapter resolves these sibling locations directly. Exact files, byte sizes and source Git identities are in EVIDENCE_INPUT_REFERENCES.json. Historical R4 data are referenced, not copied wholesale.\n\nThe original P1 decoder/config/points remain inside the R5 repository at src/p1_material/ and research/p1_material_target_20261018/. A30 requires its src/ and vendor/frozen_a29/src/ code trees plus the named paid weak.npz background. Importing background arrays is not a new solve. Any fresh forward/inverse action must use the shared cumulative R5 ledger; published reproduction is not authorization to reset or silently replay paid physics.\n\nUse the existing Gaussian/.venv_nn/bin/python with NumPy/SciPy/Torch and VECLIB_MAXIMUM_THREADS=1, OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1. B's final observation-only entry is stage_b_material_risk/inference.py:Action. y is finite real source-major512; online keys are observation,residual,q0,chi0,J from the unchanged historical96+64 strong anchor. Original K32 additionally requires parent IDs for its frozen sample-seed law. Student inference needs only its R5 PROVIDER.pt/HEAD.pt plus historical feature-source code and online anchor input.\n\nThe executable replay script uses retained VAL data and therefore tests serialization/device replay, not fresh independent generalization. Preserve INFERENCE_REPLAY.json's declared tolerance and actual status. It may contain small CPU/MPS float32 drift while material risk agrees.\n\nA5 labels are directional: global decoder constants are along one fixed unit strong W segment at a fixed weak center. The legacy RESULTS key L_on_rho_ball does not establish a full2D/66-dimensional ball. Norm values are Float64 evaluations of exact-arithmetic inequalities, without directed-rounding enclosure. Parent captions must state this scope.\n''')
 print(json.dumps(dict(files=len(files),bytes=sum(x[1] for x in files),raw_files=sum(x['files'] for x in raw.values()),raw_bytes=sum(x['bytes'] for x in raw.values()),historical_references=len(references),newR5_references=len(local),B_replay=replay['status'])))
if __name__=='__main__':main()

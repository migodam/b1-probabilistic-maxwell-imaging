#!/usr/bin/env python3
"""Inert prospective R4L1 same-cache validation; explicit parent release only.

--prepare is zero-physics. --preflight checks frozen online files without Book
or Target. No command waits for data then launches itself. Pair starts once.
"""
from pathlib import Path
import argparse,json,sys,subprocess,os,time
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
import numpy as np
import run_r4_interval_chain_global_gmres as common
from p3r4.evaluation import checkpoint_identity
from p3r4.core import OUT,write,from_probit

engine=common.engine
EXPERIMENT='interval_chain_locked_v1';PARENTS=tuple(f'R4L1_{i:05d}' for i in range(8));SEED=20261025171
ARMS=common.ARMS;CAPS=dict(full=65536,richardson=78154,global_gmres=78155,global_gmres_accuracy=78155)
PREPARED=OUT/'AUDIT/locked_samecache/FROZEN_PROTOCOL.json';RESULT=OUT/'results/IV'/EXPERIMENT
MARKER=OUT/'results/IV/interval_chain_locked_pair';LOCKED=OUT/'configs/FROZEN_LOCKED.json'
engine.EXPERIMENT=EXPERIMENT;engine.PARENTS=PARENTS;engine.SEED=SEED;engine.RESULT=RESULT;engine.CAP=CAPS
engine.LEDGER_STAGE='IV';engine.CONTEXT_DIR=OUT/'results/contexts/locked';engine.CACHE_DIR=OUT/'results/contexts/cache'
engine.FROZEN_DESIGN_EXTENSION=dict(common.engine.FROZEN_DESIGN_EXTENSION,experiment=EXPERIMENT,ledger_stage='IV',
 cohort='fresh first8R4L1, after parent FROZEN_LOCKED and closed online context acquisition',seed=SEED,
 no_automatic_locked_unseal=True,source_protocol='AUDIT/global_gmres/PROTOCOL_ZH.md',source_target_unchanged=True)

def protocol():
 return dict(status='PREPARED_NOT_RELEASED',experiment=EXPERIMENT,campaign_id='P3LR4-20261025',ledger_stage='IV',parents=PARENTS,
  steps=1024,beta=.002,seed=SEED,arms=ARMS,combined_complex_RHS_cap=300000,arm_complex_RHS_caps=CAPS,arm_wall_seconds=1200,arm_RSS_cap_bytes=3*1024**3,
  target_config=json.loads((OUT/'configs/TARGET.json').read_text()),target_identity=checkpoint_identity(OUT/'configs/TARGET.json'),
  common_randomness='each arm independently resets SeedSequence([20261025171,parent_index]);66 normal thenuniform every step',
  rights='each approximate arm independently cold pays576inversecolumns perparent; same eight-source R0/G/S matrix rights,stored canonical endpoint-once arrays',
  reduced_LS='actual complex small LS RHS charged1 before each new solve;cached same-endpoint solution replay0;no RHS debit for cached R0 matrix multiply',
  accuracy_arm='both endpoints white_error_upper<=.05 AND Phiwidth<=.05 before sameU MH, max4 refinement then paidFull fallback',
  initial_context_cost='new locked96+64 context owner separately disclosed,never treated as free cold acquisition',
  required_inputs=dict(locked_marker=str(LOCKED),contexts=str(engine.CONTEXT_DIR),cache=str(engine.CACHE_DIR),observations=str(OUT/'inputs/online')),
  release_required_fields=['released:true','experiment:interval_chain_locked_v1','campaign_id:P3LR4-20261025','prepared_protocol_identity','context_receipt'],
  source_modules=['tools/run_r4_interval_chain_locked.py','tools/run_r4_interval_chain_global_gmres.py','tools/run_r4_interval_chain_fast.py','src/p3r4/decision_refresh_global_gmres.py','src/p3r4/decision_refresh_fast.py'],
  no_GT=True,new_physics_during_preparation=0,NN_training=0,posterior_reference='NONE',floating_certification=False,
  failure_policy='missing marker/context/cache or release fails before Book; no substitute parent,rerun,adaptive beta,iteration orbudget change',
  comparison='allarms close before independent Full audit; flags/xpaths/ALL potential proofs/final intervals/fixedaccuracy;ordinary Krylov benefit not excluded',
  novelty_or_GO_claim=False)

def preflight(release=None):
 missing=[];identities=[];anchors=[]
 for p in (PREPARED,LOCKED,OUT/'configs/TARGET.json'):
  if not p.is_file():missing.append(str(p))
 if missing:return dict(status='NOT_READY_MISSING_PREREQUISITES',missing=missing,new_RHS=0)
 prepared=json.loads(PREPARED.read_text());freeze=json.loads(LOCKED.read_text())
 if not isinstance(freeze,dict) or not freeze:raise ValueError('nonempty FROZEN_LOCKED object required')
 if prepared['experiment']!=EXPERIMENT or prepared['parents']!=list(PARENTS) or prepared['target_config']!=json.loads((OUT/'configs/TARGET.json').read_text()):raise ValueError('prepared locked roster/target changed')
 identities.extend(checkpoint_identity(p) for p in (PREPARED,LOCKED,OUT/'configs/TARGET.json'))
 if release is not None:
  r=json.loads(release.read_text())
  if r.get('released') is not True or r.get('experiment')!=EXPERIMENT or r.get('campaign_id')!='P3LR4-20261025':raise PermissionError('explicit matching parent release required')
  if r.get('prepared_protocol_identity')!=checkpoint_identity(PREPARED):raise ValueError('release must name exact prepared protocol identity')
  receipt=Path(r.get('context_receipt','')).resolve()
  if OUT/'IV/jobs' not in receipt.parents or receipt.name!='CLOSED.json' or not receipt.is_file():raise ValueError('context acquisition must have closed StageIV receipt')
  if json.loads(receipt.read_text())['status']!='COMPLETE':raise ValueError('context acquisition not COMPLETE')
  source=json.loads((receipt.parent/'SOURCE_REVISION.json').read_text())
  if 'locked' not in source.get('argv',[]):raise ValueError('context receipt must be locked context acquisition')
  identities.extend([checkpoint_identity(release),checkpoint_identity(receipt)])
 for parent in PARENTS:
  paths=[engine.CONTEXT_DIR/f'{parent}.npz',engine.CACHE_DIR/f'{parent}.npz',OUT/'inputs/online'/f'{parent}.npz']
  if any(not p.is_file() for p in paths):missing.extend(str(p) for p in paths if not p.is_file());continue
  with np.load(paths[0],allow_pickle=False) as z:x=z['x'].copy();q0=z['q0'].copy();y=z['observation'].copy();phi=float(z['phi'])
  with np.load(paths[1],allow_pickle=False) as z:
   if not all(k in z.files for k in ('q','prediction','field','current','LU','pivots','receiver_T')):raise ValueError('incomplete paid anchor cache; no rebuilding')
   q=z['q'].copy()
  with np.load(paths[2],allow_pickle=False) as z:online=z['observation'].copy()
  if x.shape!=(66,) or not np.isfinite(x).all() or not np.array_equal(q,q0) or np.max(abs(from_probit(x)-q))>1e-12 or not np.array_equal(y,online):raise ValueError('locked online context/cache/observation mismatch')
  identities.extend(checkpoint_identity(p) for p in paths);anchors.append(dict(parent=parent,anchor_Phi=phi,role='online context, no stationarity/posterior claim'))
 pinned=MARKER/'INPUT_PREFLIGHT.json'
 if release is not None and pinned.is_file() and identities!=json.loads(pinned.read_text())['identities']:
  raise ValueError('locked inputs/release differ from coordinator preflight; stop before Book')
 return dict(status='READY_EXPLICIT_RELEASE_REQUIRED' if not missing else 'NOT_READY_MISSING_PREREQUISITES',missing=missing,identities=identities,anchors=anchors,new_RHS=0)

def pair(release):
 checked=preflight(release)
 if checked['missing']:raise FileNotFoundError(json.dumps(checked))
 if MARKER.exists() or RESULT.exists():raise FileExistsError('preserve first attempt; never paid replay')
 MARKER.mkdir(parents=True);write(MARKER/'INPUT_PREFLIGHT.json',checked);timings={};records=[];env=dict(os.environ,R4_ARRAY_CODEC='stored')
 for arm in (*ARMS,'compare'):
  start=time.perf_counter();p=subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),arm,'--release',str(release.resolve())],cwd=ROOT,env=env)
  timings[arm]=dict(seconds=time.perf_counter()-start,returncode=p.returncode);records.append(dict(arm=arm,returncode=p.returncode))
  if arm==ARMS[-1] and p.returncode==0:write(RESULT/'ARM_PROCESS_WALL.json',timings)
  if p.returncode:write(MARKER/'FAILED.json',dict(status='FAILED',records=records));raise SystemExit(p.returncode)
 write(MARKER/'COMPLETE.json',dict(status='COMPLETE',records=records,comparison=str(RESULT/'COMPARISON.json')))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','preflight','pair',*ARMS,'compare'));p.add_argument('--release',type=Path);a=p.parse_args()
 if a.action=='prepare':
  if PREPARED.exists():raise FileExistsError('preserve frozen preparation')
  write(PREPARED,protocol());print(json.dumps(dict(status='PREPARED_NOT_RELEASED',protocol=str(PREPARED),new_RHS=0)));return
 if a.action=='preflight':print(json.dumps(preflight(a.release)));return
 if a.action=='compare':common.compare();return
 if a.release is None:p.error('parent release required; no automatic launch')
 if a.action=='pair':pair(a.release)
 else:
  checked=preflight(a.release)
  if checked['missing']:raise FileNotFoundError(json.dumps(checked))
  engine.run(a.action,a.release)

if __name__=='__main__':main()

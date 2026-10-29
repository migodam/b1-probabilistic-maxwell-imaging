"""Freeze public GENERAL6 protocol, then privately generate observations after solver freeze.
This is a generator process and must never be imported by solver code.
"""
import argparse,json,os,secrets,sys,time
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from scipy.stats import qmc
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from physics import Model

def freeze():
 cp=P/'experiments/GENERAL6_CONTRACT.json'
 if cp.exists():raise RuntimeError('GENERAL6 frozen contract already exists; do not overwrite')
 original=json.loads((P/'experiments/contract.json').read_text());e1=original['presets']['E1_LOW6'];g={'centers':e1['centers'],'sources':e1['sources']+[[-.67,.42]],'receivers':e1['receivers'],'radius':.04,'k':2*np.pi};model=Model(None,**g);incident=np.linalg.svd(model.E0,compute_uv=False);starts=[np.full(6,15.25)]
 for x in qmc.Sobol(6,scramble=False).random_base2(7):
  m=.5+29.5*x
  if all(np.linalg.norm(m-y)>1e-12 for y in starts):starts.append(m)
  if len(starts)==64:break
 c={'schema':'Z-ALG-R2-GENERAL6-v1','frozen_utc':datetime.now(timezone.utc).isoformat(),'scope':'Independent NEW six-cell/three-source balanced quotient study; all6 independent real material cells unknown; not a retuning/replacement of prior studies, not Full66','geometry':g,'material_box':[.5,30.],'observations':{'sources':3,'complex_receivers':1,'real_channels':6,'packing':'source major,receiver major,real then imaginary','noise':'none'},'geometry_only_audit':{'incident_singular_values':incident.tolist(),'incident_rank_rtol1e12':int(sum(incident>incident[0]*1e-12)),'setup_counts':dict(model.counts),'setup_timings':dict(model.timings)},'hidden_targets':{'count':3,'distribution':'independent product uniform on[.5,30] for all6cells','seed_policy':'fresh secrets.randbits(128); private only,never logged publicly or passed to solver','generation_precondition':'core solver implementation+parameters frozen in an external receipt before generation; no GT exposure to engineer'},'initialization':{'box_midpoint':15.25,'ReducedFullJ64_starts':[x.tolist() for x in starts],'prefix8':'first8 of the identical64starts; optional reporting only','start_policy':'midpoint then first63 distinct unscrambled Sobol6 points scaled to box; midpoint duplicate skipped','forbidden':'GT-derived support,coordinate,initialization,hyperparameter,root selection'},'methods':{'main':'general balanced bilinear quotient for6unknowns/3sources; prescribed20eigenvalue endpoint candidates; numerical method, no rigorous completeness claim','baseline':'ReducedFullJ64 exact same fixed-active SMWkernel and analytic full material Jacobian; all6free, no oracle'},'optimizer':{'max_nfev_per_start':120,'ftol':1e-12,'xtol':1e-12,'gtol':1e-12,'loss':'linear','x_scale':1.,'residual_scaling':'all observations divided by fixed observed target norm'},'validation':{'independent':'fresh numpy full Maxwell solve at every unique legal candidate; fitted current never reused','response_relative':1e-8,'state_relative':1e-10,'resolvent_sigma_min_floor':1e-6,'merge_material_RMS':1e-7,'contains_GT_material_RMS':1e-5,'box_tolerance':1e-8},'accounting':{'generation_actual_complex_RHS':9,'generation':'one fresh3source full forward per target, setup separately timed','RHS':'full source/adjoint RHS separate from reduced algebraic solves','all_paths_roots_starts_failures':'preserve without truth-dependent filtering','single_BLAS_thread':True,'method_case_wall_seconds':120,'method_case_actual_complex_RHS_cap':15000},'blindness':'computational observation-only process contract; shared filesystem is not OS sandbox; separate private generator/scorer','score_policy':'only after parent or solver engineer confirms all predictions frozen; retain every root error; nearest truth describes root-set coverage only, not GT-selected reconstruction'}
 cp.write_text(json.dumps(c,indent=2)+'\n');print('GENERAL6 public protocol frozen; no hidden targets generated; incident rank',c['geometry_only_audit']['incident_rank_rtol1e12'])

def generate(receipt):
 cp=P/'experiments/GENERAL6_CONTRACT.json';c=json.loads(cp.read_text());receipt=Path(receipt).resolve()
 if not receipt.exists():raise RuntimeError('Need engineer/parent core solver freeze receipt before generating observations')
 private=P/'private/general6';public=P/'experiments/general6_public';private.mkdir(exist_ok=True);public.mkdir(exist_ok=True)
 if (private/'sealed_targets.npz').exists() or (public/'manifest.json').exists():raise RuntimeError('GENERAL6 observations already generated; no overwrite')
 seed=secrets.randbits(128);rng=np.random.default_rng(seed);t=time.perf_counter();cpu=time.process_time();model=Model(None,**c['geometry']);setup_wall=time.perf_counter()-t;setup_cpu=time.process_time()-cpu;targets={};rows=[];forward_cost=[]
 for i in range(3):
  case=f'GENERAL6_{i:02d}';m=rng.uniform(.5,30,6);targets[case]=m;tick=time.perf_counter();ct=time.process_time();state=model.forward(m,independent=True);forward_cost.append({'case':case,'complex_RHS':3,'wall_seconds':time.perf_counter()-tick,'CPU_seconds':time.process_time()-ct});np.savez(public/f'{case}.npz',y=state['y']);rows.append({'case':case,'observation_file':f'{case}.npz'})
 np.savez(private/'sealed_targets.npz',**targets);secret=private/'generator_secret.json';secret.write_text(json.dumps({'seed':str(seed),'created_utc':datetime.now(timezone.utc).isoformat()},indent=2)+'\n');os.chmod(secret,0o600);os.chmod(private/'sealed_targets.npz',0o600)
 (public/'manifest.json').write_text(json.dumps({'schema':c['schema'],'cases':rows,'allowed':'public geometry,box,observations only','forbidden':'private material truth/seed'},indent=2)+'\n')
 r={'status':'OBSERVATIONS_GENERATED_AFTER_SOLVER_FREEZE','UTC':datetime.now(timezone.utc).isoformat(),'solver_freeze_receipt':str(receipt),'contract':str(cp.resolve()),'targets':3,'actual_complex_RHS':9,'geometry_setup_wall_seconds':setup_wall,'geometry_setup_CPU_seconds':setup_cpu,'forward_cost':forward_cost,'total_wall_seconds':time.perf_counter()-t,'total_CPU_seconds':time.process_time()-cpu,'model_counts':dict(model.counts),'model_timings':dict(model.timings),'GT_values_disclosed':False,'seed_disclosed':False};(P/'results/general6_generation_receipt.json').write_text(json.dumps(r,indent=2)+'\n');print('GENERAL6 three observation-only cases generated after solver freeze;9actual complex RHS;truth/seed suppressed.')

def main():
 a=argparse.ArgumentParser();group=a.add_mutually_exclusive_group(required=True);group.add_argument('--freeze-contract',action='store_true');group.add_argument('--generate',action='store_true');a.add_argument('--solver-freeze-receipt');args=a.parse_args()
 if args.freeze_contract:freeze()
 else:
  if not args.solver_freeze_receipt:raise RuntimeError('--solver-freeze-receipt required for generation')
  generate(args.solver_freeze_receipt)
if __name__=='__main__':main()

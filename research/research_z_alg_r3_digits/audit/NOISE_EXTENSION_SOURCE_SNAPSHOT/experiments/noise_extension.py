"""Paired 1e-3 noise extension, prefrozen methods, same nine materials.

Never call this a fresh material holdout. Main material truth is already
released before this optional stage; no method tuning is permitted.
"""
import argparse,json,os,subprocess,sys,time,uuid
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import score_predictions as frozen_scorer
ROOT=Path(__file__).resolve().parents[1]
METHODS=['fullj_native','fullj_reciprocal','globaldual']
def dump(p,x):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n')
def generate(parent_go):
 if not parent_go or not (ROOT/'audit/NOISE_EXTENSION_SOURCE_FREEZE.json').exists():raise RuntimeError('Parent GO and full noise source freeze required')
 if not (ROOT/'audit/EXPERIMENT_ALL_PREDICTIONS_FROZEN_GT_RELEASE.json').exists():raise RuntimeError('Core63 must finish first')
 out=ROOT/'public_data/noise_extension';truth=ROOT/'private_truth/noise_extension'
 if out.exists() and list(out.glob('*.npz')):raise RuntimeError('Refusing overwrite existing noise observations')
 contract=json.loads((ROOT/'configs/NOISE_EXTENSION_CONTRACT.json').read_text());level=float(contract['relative_complex_noise_norm']);mapping=json.loads((ROOT/'truth_after_freeze/holdout/mapping_seedless.json').read_text());records=mapping['records']
 if len(records)!=9:raise RuntimeError('All same nine target materials required')
 start=time.perf_counter();entropy=int.from_bytes(os.urandom(16),'big');rng=np.random.default_rng(entropy);rng.shuffle(records);rows=[];private=[];out.mkdir(parents=True,exist_ok=True);truth.mkdir(parents=True,exist_ok=True)
 for old in records:
  original=old['observation_id'];ident='noise_'+uuid.uuid4().hex[:12];clean=np.load(ROOT/'public_data/holdout'/f'{original}.npz')['Y'];direction=rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape);noise=direction*(level*np.linalg.norm(clean)/np.linalg.norm(direction));Y=clean+noise;sigma=float(np.sqrt(np.mean(abs(noise)**2)))
  np.savez_compressed(out/f'{ident}.npz',Y=Y,sigma=sigma,observation_id=np.array(ident))
  material=np.load(ROOT/'truth_after_freeze/holdout'/f'{original}.npz');np.savez_compressed(truth/f'{ident}.npz',m=material['m'],mask=material['mask'],Y_clean=clean)
  params=dict(old['parameters']);params.update(paired_source_observation_id=original,noise_relative_norm=level)
  private.append(dict(observation_id=ident,digit=old['digit'],parameters=params))
  rows.append(dict(observation_id=ident,actual_relative_complex_noise_norm=float(np.linalg.norm(Y-clean)/np.linalg.norm(clean)),sigma=sigma))
 dump(truth/'mapping_PRIVATE.json',dict(entropy=str(entropy),records=private))
 receipt=dict(split='noise_extension',scope='paired same nine material targets; not a fresh material holdout',count=9,relative_complex_noise_norm=level,rows=rows,generation_seconds=time.perf_counter()-start,additional_PDE_LU=0,additional_explicit_PDE_RHS=0,reason_no_new_PDE='reuses exact noiseless public main Y; original physics generation paid separately',noise_distribution='independent circular complex Gaussian direction per case, radially normalized to fixed global norm; not iid Gaussian after normalization',main_material_GT_already_disclosed=True,method_tuning_after_main_GT=False)
 dump(ROOT/'audit/EXPERIMENT_GENERATION_NOISE_EXTENSION.json',receipt);print(json.dumps(receipt))
def run():
 if not (ROOT/'audit/NOISE_EXTENSION_SOURCE_FREEZE.json').exists():raise RuntimeError('Noise source freeze required')
 cases=sorted((ROOT/'public_data/noise_extension').glob('*.npz'))
 if len(cases)!=9:raise RuntimeError('All nine paired noise observations required')
 env=os.environ.copy()
 for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:env[k]='1'
 start=time.perf_counter();rows=[]
 def work(job):
  case,method=job;out=ROOT/'results/noise_extension'/method/f'{case.stem}.json';out.parent.mkdir(parents=True,exist_ok=True)
  if out.exists():raise FileExistsError(out)
  command=[sys.executable,str(ROOT/'run_solver.py'),'--public',str(case),'--operators',str(ROOT/'public_data/operators.npz'),'--method',method,'--config',str(ROOT/'configs/noisy_methods.json'),'--output',str(out)]
  begun=time.perf_counter();p=subprocess.run(command,cwd=ROOT,env=env,capture_output=True,text=True);out.with_suffix('.log').write_text(p.stdout+'\nSTDERR\n'+p.stderr);row=dict(observation_id=case.stem,method=method,process_exit=p.returncode,process_wall_seconds=time.perf_counter()-begun,output=str(out.relative_to(ROOT)))
  if out.exists():
   r=json.loads(out.read_text());row.update(status=r.get('status'),wall_seconds=r.get('wall_seconds'),cpu_seconds=r.get('cpu_seconds'),all_unknowns_released=r.get('all_unknowns_released'),relative_residual=(r.get('top1') or {}).get('relative_residual'),accounting=r.get('accounting'))
  else:row['status']='PROCESS_FAILED_NO_OUTPUT'
  print(json.dumps({k:v for k,v in row.items() if k!='accounting'}),flush=True);return row
 with ThreadPoolExecutor(max_workers=2) as executor:
  for future in as_completed([executor.submit(work,(case,method)) for case in cases for method in METHODS]):rows.append(future.result())
 receipt=dict(split='noise_extension',scope='paired same nine materials, not18independent materials',workers=2,BLAS_threads_per_worker=1,rows=rows,total_batch_wall_seconds=time.perf_counter()-start,configuration=json.loads((ROOT/'configs/noisy_methods.json').read_text()),methods_frozen_before_main_GT=True,no_GT_input_to_solver=True,all_processes_succeeded=all(r['process_exit']==0 for r in rows));dump(ROOT/'audit/EXPERIMENT_BATCH_NOISE_EXTENSION.json',receipt)
 print(json.dumps(dict(outputs=len(rows),wall_seconds=receipt['total_batch_wall_seconds'],success=receipt['all_processes_succeeded'])))
def main():
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['generate','run','freeze','score']);ap.add_argument('--parent-go',action='store_true');a=ap.parse_args()
 if a.action=='generate':generate(a.parent_go)
 elif a.action=='run':run()
 elif a.action=='freeze':
  frozen_scorer.freeze('noise_extension');p=ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_NOISE_EXTENSION.json';r=json.loads(p.read_text());r.update(truth_read_before_freeze=True,material_GT_already_disclosed=True,scope='paired noise extension on same9 known materials; methods/source frozen before mainGT, predictions frozen before extension scoring');dump(p,r)
 else:
  if not (ROOT/'audit/EXPERIMENT_PREDICTIONS_FREEZE_NOISE_EXTENSION.json').exists():raise RuntimeError('All27 noisy predictions must freeze first')
  frozen_scorer.score('noise_extension')
if __name__=='__main__':main()

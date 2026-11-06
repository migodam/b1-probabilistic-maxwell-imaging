"""Generate the root-frozen bothblocks16 observations once; no inverse or scoring."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import argparse,json,time,sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,required=True);ap.add_argument('--execute',action='store_true');ap.add_argument('--cpu',type=int);ap.add_argument('--vendor-path');a=ap.parse_args()
 if a.cpu is not None:
  import psutil
  psutil.Process().cpu_affinity([a.cpu])
  if os.name=='nt':psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
 if a.vendor_path:sys.path.insert(0,a.vendor_path)
 cfg=json.loads(a.config.read_text())
 if cfg.get('seed')!=202611071116:raise ValueError('Frozen root seed differs')
 if not a.execute or cfg['status'] not in ('FROZEN_FOR_GENERATION','FROZEN_BEFORE_GENERATION_AND_CONFIRMATION_GT_ACCESS'):raise ValueError('Root freeze then explicit --execute required')
 expected=[(f,r,12. if (i+j)%2==0 else 18.) for i,f in enumerate(('annulus','slotted_shell','multi_island','curved_bars','comb','coated_inclusion','perforated_star','bridge')) for j,r in enumerate(('high','mixed'))]
 if [(r['family'],r['regime'],r['high_chi']) for r in cfg['cases']]!=expected:raise ValueError('Frozen16 balance differs')
 if (cfg['generation_grid'],cfg['inverse_grid'],cfg['sources'],cfg['receivers'],cfg['material_background'],cfg['material_box'],cfg['relative_complex_RMS_noise'])!=(24,24,12,48,.1,[.02,20.],.01):raise ValueError('Physics/noise contract differs')
 pub=(REPO/cfg['public_output']).resolve();private=(REPO/cfg['private_output']).resolve();receiptpath=ROOT/'delegated/diagnostics/bothblocks_confirmation/GENERATION_RECEIPT.json'
 if not pub.is_relative_to(ROOT/'data_public') or not private.is_relative_to(ROOT/'sealed_truth'):raise ValueError('Owned paths required')
 for p in (pub,private,receiptpath):
  if p.exists():raise FileExistsError(p)
 pub.mkdir(parents=True);private.mkdir(parents=True);receiptpath.parent.mkdir(parents=True,exist_ok=True)
 receipt=dict(status='RUNNING',requested=16,completed=0,truth_access='GENERATOR_ONLY',inversions=0,operator_setup=[],objects=[],prior='Same inherited physical Regularizer, no inverse here',CPU_scope='Absolute whole generator child since creation, includes imports/setup/forward/noise/serialization')
 wall=time.perf_counter();model=None;records=[];metadata=[];score=[]
 try:
  import numpy as np
  from scipy.stats import chi2
  from threadpoolctl import threadpool_limits
  from research.research_max_struct_r10.src.twod_materials import draw_material
  from research.research_fbg_r3_experiments.src.physics import MaxwellModel
  dump(private/'FROZEN_GENERATOR_CONFIG.json',cfg)
  with threadpool_limits(limits=1):
   c,t=time.process_time(),time.perf_counter();source=(REPO/cfg['operator_source']).resolve();destination=pub/'operators.npz';shutil.copy2(source,destination)
   with np.load(destination,allow_pickle=False) as z:public={k:z[k].copy() for k in z.files}
   if (public['G'].shape,public['E0'].shape,public['S'].shape)!=((576,576),(576,12),(48,576)):raise ValueError('Operator shape')
   if {'truth','m_true','Y_clean','family','seed'}&set(public):raise ValueError('Private operator payload')
   model=MaxwellModel(public);receipt['operator_setup'].append(dict(CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,operator_source=cfg['operator_source'],original_generation='REUSED_NOT_RERUN',physical_cost=dict(model.cost)))
   children=np.random.SeedSequence(cfg['seed']).spawn(16)
   for index,row in enumerate(cfg['cases']):
    c,t=time.process_time(),time.perf_counter();before=dict(model.cost);geometry_child,noise_child=children[index].spawn(2);geometry_rng=np.random.default_rng(geometry_child);noise_rng=np.random.default_rng(noise_child);m,parameters=draw_material(geometry_rng,row['family'],row['regime'],row['high_chi'])
    if not np.isfinite(m).all() or np.min(m)<.02 or np.max(m)>20:raise ValueError('Illegal material; no reroll')
    state=model.state(m);clean=state['Y'].copy();sigma=float(.01*np.linalg.norm(clean)/np.sqrt(clean.size));noise=sigma/np.sqrt(2)*(noise_rng.normal(size=clean.shape)+1j*noise_rng.normal(size=clean.shape));obs=pub/(row['case_id']+'.npz')
    np.savez_compressed(obs,Y=clean+noise,sigma=sigma,noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2)),noise_model=np.array('iid_circular_complex'),grid_shape=np.array([24,24]),material_bounds=np.array([.02,20.]))
    np.savez_compressed(private/(row['case_id']+'.npz'),m=m,Y_clean=clean)
    metadata.append(dict(**row,parameters=parameters,independent_child_spawn_key=list(children[index].spawn_key),geometry_spawn_key=list(geometry_child.spawn_key),noise_spawn_key=list(noise_child.spawn_key),realized_relative_noise=float(np.linalg.norm(noise)/np.linalg.norm(clean))))
    score.append(dict(case_id=row['case_id'],amplitude=float(np.max(m)-.1)))
    records.append(dict(case_id=row['case_id'],public_files=[str(destination.relative_to(REPO)),str(obs.relative_to(REPO))]))
    after=dict(model.cost);receipt['objects'].append(dict(case_id=row['case_id'],CPU_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,physical_cost={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)}));receipt['completed']=len(records);dump(receiptpath,receipt)
   dump(private/'generator_metadata.json',dict(seed=cfg['seed'],records=metadata));dump(private/'score_metadata.json',dict(records=score,scope='OFFLINE_ONLY'))
   dump(pub/'manifest.json',dict(status='GENERATED_PUBLIC_OBSERVATIONS',records=records,inverse_input_policy='Only public_files; no geometry/regime/seed/amplitude/GT'))
   candidates=dict(fullG19=('full_material_readout','prepare_full'),fullGreenBothBlocks=('full_material_readout','prepare_full'))
   groups={'baseline':[],'candidate':[]}
   for record in records:
    for method in ('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks'):
     ident=record['case_id']+'__'+method+'__480CPU';task=dict(task_id=ident,case_id=record['case_id'],method=method,arm='a_plus' if method in candidates else method,dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,public_files=record['public_files'],config='research/research_max_struct_r11/configs/twod_bothblocks16_solver.json',output_dir='research/research_max_struct_r11/RESULTS/bothblocks16/'+ident,experiment_scope='FRESH_INDEPENDENT_CONFIRMATION')
     if method in candidates:
      call,prepare=candidates[method];task['material_readout']=dict(module_file='research/research_max_struct_r11/src/material_readouts.py',callable=call,prepare=prepare,label=method)
     task['runner']='diagnose_bothblocks_run.py' if method=='fullGreenBothBlocks' else 'twod_run.py'
     groups['candidate' if method in candidates else 'baseline'].append(task)
   methods=('cie_single','cie_dual','global_dual','fullG19','fullGreenBothBlocks')
   groups['all']=sorted(groups['baseline']+groups['candidate'],key=lambda t:(t['case_id'],(methods.index(t['method'])-int(t['case_id'][-4:]))%5))
   for group,tasks in groups.items():dump(ROOT/('configs/twod_bothblocks16_'+group+'_tasks.json'),dict(status='PUBLIC_TASKS_NO_GT',input_base='.',tasks=tasks))
   receipt['status']='GENERATED_SEPARATE_TRUTH'
 except BaseException as exc:receipt.update(status='FAILED_PARTIAL_GENERATION_PRESERVED',error=repr(exc));raise
 finally:
  receipt.update(whole_child_CPU_seconds=time.process_time(),wall_seconds=time.perf_counter()-wall,physical_cost=dict(model.cost) if model else {});dump(receiptpath,receipt)
 print(json.dumps(dict(status=receipt['status'],completed=receipt['completed'],CPU_seconds=receipt['whole_child_CPU_seconds'])))
if __name__=='__main__':main()

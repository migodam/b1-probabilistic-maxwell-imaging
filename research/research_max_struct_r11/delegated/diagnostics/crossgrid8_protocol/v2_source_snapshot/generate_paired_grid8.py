"""Frozen8 geometries, paired matched24/cross48to24 data; generator only."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import sys,json,argparse,time,shutil,traceback
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--execute',action='store_true');a=ap.parse_args()
 cfg=json.loads((HERE/'GENERATOR_CONFIG_V2.json').read_text());solver=json.loads((HERE/'SOLVER_CONFIG_V2.json').read_text())
 if not a.execute or cfg['status']!='FROZEN_FOR_GENERATION' or solver['status']!='FROZEN_FOR_INVERSION':raise PermissionError('Root freeze and explicit execution required')
 if cfg['seed']!=202611071118 or cfg['independent_materials']!=8 or cfg['conditions']!=['matched24','cross48to24'] or cfg['methods']!=['cie_single','cie_dual','global_dual','fullG19']:raise ValueError('Frozen paired contract differs')
 pub=REPO/cfg['public_output'];private=REPO/cfg['private_output'];receipt=HERE/'GENERATION_V2_RECEIPT.json';manifest=HERE/'TASKS_V2.json'
 for p in (pub,private,receipt,manifest):
  if p.exists():raise FileExistsError('Preserve prior generation attempt: '+str(p))
 pub.mkdir(parents=True);private.mkdir(parents=True)
 for condition in cfg['conditions']:(pub/condition).mkdir()
 wall=time.perf_counter();fine=coarse_model=None;records=[];metadata=[]
 rec=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),independent_materials=8,grid_conditions=cfg['conditions'],inverse_tasks=64,completed_materials=0,truth_access='GENERATOR_ONLY',new_inverse_runs=0,per_object=[],CPU_scope='Absolute generator child process_time; external wait4 authoritative and nested timers not added')
 dump(receipt,rec)
 try:
  import numpy as np
  from scipy.stats import chi2
  from threadpoolctl import threadpool_limits
  from research.research_max_struct_r10.src.twod_materials import draw_material
  from research.research_z_alg_r4_diversity.src.physics import Model
  with threadpool_limits(limits=1):
   source=REPO/cfg['coarse_operator_source'];shutil.copy2(source,pub/'operators.npz')
   with np.load(source,allow_pickle=False) as z:coarse={k:z[k].copy() for k in z.files}
   if (coarse['G'].shape,coarse['E0'].shape,coarse['S'].shape)!=((576,576),(576,12),(48,576)):raise ValueError('Original acquisition shape differs')
   if {'truth','gt','m_true','m_truth','Y_clean','family','seed','mask'} & set(coarse):raise ValueError('Private payload in public operators')
   if abs(float(coarse['k'])-2*np.pi)>1e-14 or abs(float(coarse['radius'])-.018)>1e-14:raise ValueError('Inherited wavelength/radius differs')
   axis=(np.arange(48)+.5)*.96/48-.48;yy,xx=np.meshgrid(axis,axis,indexing='ij');centers=np.column_stack((xx.ravel(),yy.ravel()))
   geom=dict(centers=centers,sources=coarse['sources'],receivers=coarse['receivers'],radius=.009,k=float(coarse['k']))
   fine=Model(geom);coarse_model=Model(coarse);np.savez_compressed(private/'FINE_OPERATOR_GEOMETRY.npz',**geom)
   dump(private/'FROZEN_GENERATOR_CONFIG.json',cfg)
   children=np.random.SeedSequence(cfg['seed']).spawn(8)
   for row,child in zip(cfg['cases'],children):
    beg=time.process_time();bf=dict(fine.counts);bc=dict(coarse_model.counts);geometry_child,noise_child=child.spawn(2)
    m,params=draw_material(np.random.default_rng(geometry_child),row['family'],row['regime'],12.,side=48)
    if not np.isfinite(m).all() or m.min()<.02 or m.max()>20:raise ValueError('Invalid fixed draw; no reroll')
    blocks=m.reshape(24,2,24,2);target=blocks.mean(axis=(1,3)).ravel()
    pure_low=np.all(blocks==2.5,axis=(1,3)).ravel();pure_high=np.all(blocks==12.,axis=(1,3)).ravel();pure_bg=np.all(blocks==.1,axis=(1,3)).ravel();boundary=~(pure_low|pure_high|pure_bg)
    clean_cross=fine.forward(m,label='PAIRED8_FINE_GENERATOR')['y'].copy()
    clean_matched=coarse_model.forward(target,label='PAIRED8_MATCHED_GENERATOR')['y'].copy()
    rng=np.random.default_rng(noise_child);standard=(rng.normal(size=clean_cross.shape)+1j*rng.normal(size=clean_cross.shape))/np.sqrt(2)
    private_arrays=dict(m=target,m_fine=m,standard_complex_noise=standard,pure_low=pure_low,pure_high=pure_high,pure_background=pure_bg,partial_volume_boundary=boundary)
    cid=row['case_id']
    for condition,clean in [('matched24',clean_matched),('cross48to24',clean_cross)]:
     sigma=float(.01*np.linalg.norm(clean)/np.sqrt(clean.size));obs=pub/condition/(cid+'.npz')
     np.savez_compressed(obs,Y=clean+sigma*standard,sigma=sigma,noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2)),noise_model=np.array('iid_circular_complex'),grid_shape=np.array([24,24]),material_bounds=np.array([.02,20.]))
     private_arrays['Y_clean_'+condition]=clean;private_arrays['sigma_'+condition]=np.array(sigma)
     records.append(dict(case_id=cid,grid_condition=condition,public_files=[(pub/'operators.npz').relative_to(REPO).as_posix(),obs.relative_to(REPO).as_posix()]))
    np.savez_compressed(private/(cid+'.npz'),**private_arrays)
    metadata.append(dict(**row,parameters=params,case_spawn_key=list(child.spawn_key),geometry_spawn_key=list(geometry_child.spawn_key),noise_spawn_key=list(noise_child.spawn_key),shared_standard_noise=True,target='Same2x2 projected m for both conditions; not exact fine/coarse scattering equivalence'))
    af=dict(fine.counts);ac=dict(coarse_model.counts);rec['per_object'].append(dict(case_id=cid,CPU_seconds=time.process_time()-beg,fine_dimension=2304,fine_counts={k:af.get(k,0)-bf.get(k,0) for k in set(af)|set(bf)},coarse_dimension=576,coarse_counts={k:ac.get(k,0)-bc.get(k,0) for k in set(ac)|set(bc)}));rec['completed_materials']=len(metadata);dump(receipt,rec)
   dump(private/'generator_metadata.json',dict(seed=cfg['seed'],records=metadata,independent_materials=8,noise_pairing='Conditions share same standard draw; not independent noise replicates'))
   dump(pub/'manifest.json',dict(status='GENERATED_PAIRED_PUBLIC_OBSERVATIONS',independent_materials=8,records=records,inverse_input_policy='Only coarse operators+noisyY/rawsigma/legal inverse metadata; no fine material/phase masks/family/seed'))
   tasks=[]
   for i,row in enumerate(records):
    methods=cfg['methods'][i%4:]+cfg['methods'][:i%4]
    for method in methods:
     tid=row['case_id']+'__'+row['grid_condition']+'__'+method+'__480CPU';t=dict(task_id=tid,case_id=row['case_id'],grid_condition=row['grid_condition'],arm='a_plus' if method=='fullG19' else method,dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,public_files=row['public_files'],config=(HERE/'SOLVER_CONFIG_V2.json').relative_to(REPO).as_posix(),output_dir=(Path(cfg['inverse_output'])/tid).as_posix(),experiment_scope='SUPPLEMENTARY_PAIRED_GRID_CONTROL_8_INDEPENDENT_MATERIALS')
     if method=='fullG19':t['material_readout']=dict(module_file='research/research_max_struct_r11/src/material_readouts.py',callable='full_material_readout',prepare='prepare_full',label='fullG19')
     tasks.append(t)
   dump(manifest,dict(status='PUBLIC_TASKS_NO_GT',input_base='.',independent_materials=8,tasks=tasks));rec['status']='GENERATED_PAIRED_SEPARATE_TRUTH'
 except BaseException:
  rec.update(status='FAILED_PARTIAL_GENERATION_PRESERVED',traceback=traceback.format_exc());raise
 finally:
  rec.update(finished_utc=datetime.now(timezone.utc).isoformat(),whole_child_CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,fine_dimension=2304,fine_counts=dict(fine.counts) if fine else {},coarse_dimension=576,coarse_counts=dict(coarse_model.counts) if coarse_model else {},future_model_error_FW=0);dump(receipt,rec)
 print(json.dumps(dict(status=rec['status'],completed_materials=rec['completed_materials'],independent_materials=8,inverse_tasks=64,new_inverse_runs=0)))
if __name__=='__main__':main()

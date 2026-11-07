"""Root-gated8 fine disk-pulse data cases; no inverse or coarse-truth forward."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import sys,json,argparse,time,shutil,traceback
HERE=Path(__file__).absolute().parent;ROOT=HERE.parents[2];REPO=ROOT.parents[1];sys.path.insert(0,str(REPO))
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def portable(p):return Path(str(p).replace('\\','/'))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--config',type=Path,default=HERE/'GENERATOR_CONFIG.json');ap.add_argument('--solver-config',type=Path,default=HERE/'SOLVER_CONFIG.json');ap.add_argument('--execute',action='store_true');ap.add_argument('--cpu',type=int);ap.add_argument('--vendor-path');a=ap.parse_args()
 cfg=json.loads(a.config.read_text());solver=json.loads(a.solver_config.read_text())
 if not a.execute or cfg['status']!='FROZEN_FOR_GENERATION' or solver['status']!='FROZEN_FOR_INVERSION':raise PermissionError('Root must review/freeze both configs before explicit --execute; no arrays generated')
 expected=[(f,r) for f in ('annulus','slotted_shell','bridge','multi_island') for r in ('high','mixed')]
 if [(r['family'],r['regime']) for r in cfg['cases']]!=expected or len({r['case_id'] for r in cfg['cases']})!=8:raise ValueError('Exactly8 frozen independent identities required')
 if (cfg['seed'],cfg['generation_grid'],cfg['inverse_grid'],cfg['sources'],cfg['receivers'],cfg['ROI_width'],cfg['wavelength'],cfg['radius_over_spacing'],cfg['relative_complex_RMS_noise'])!=(202611071118,48,24,12,48,.96,1.,.45,.01):raise ValueError('Fixed cross-grid acquisition differs')
 if any(r['high_chi']!=12. or r['low_chi']!=2.5 for r in cfg['cases']):raise ValueError('No contrast sweep')
 if cfg['material_background']!=.1 or cfg['material_box']!=[.02,20.] or cfg['methods']!=['cie_single','cie_dual','global_dual','fullG19']:raise ValueError('Frozen material/method contract differs')
 if solver['whole_child_CPU_budget_seconds']!=480 or solver['corrector_maxiter']!=3000 or solver['common']['host']!=.1 or solver['common']['material_box']!=[.02,20.]:raise ValueError('Inherited inverse budget/prior initialization contract differs')
 if a.cpu is not None:
  import psutil
  psutil.Process().cpu_affinity([a.cpu])
  if os.name=='nt':psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
 if a.vendor_path:sys.path.insert(0,a.vendor_path)
 pub=REPO/portable(cfg['public_output']);private=REPO/portable(cfg['private_output']);receipt=HERE/'GENERATION_RECEIPT.json';manifest=HERE/'TASKS.json'
 if not pub.absolute().is_relative_to(ROOT/'data_public') or not private.absolute().is_relative_to(ROOT/'sealed_truth'):raise ValueError('Owned8-only output paths required')
 for p in (pub,private,receipt,manifest):
  if p.exists():raise FileExistsError('Preserve prior attempts: '+str(p))
 pub.mkdir(parents=True);private.mkdir(parents=True);wall=time.perf_counter();fine=None;records=[];metadata=[]
 rec=dict(status='RUNNING',started_utc=datetime.now(timezone.utc).isoformat(),requested=8,completed=0,GT_access='GENERATOR_ONLY_NO_INVERSE',new_inverse_runs=0,coarse_projected_truth_forward=0,physical_matrix_dimension=2304,per_object=[],CPU_scope='Absolute generator child process_time includes imports/setup/forward/export; final external OS wrapper authoritative')
 dump(receipt,rec)
 try:
  import numpy as np
  from scipy.stats import chi2
  from threadpoolctl import threadpool_limits
  from research.research_max_struct_r10.src.twod_materials import draw_material
  from research.research_z_alg_r4_diversity.src.physics import Model
  with threadpool_limits(limits=1):
   source=REPO/portable(cfg['coarse_operator_source']);shutil.copy2(source,pub/'operators.npz')
   with np.load(source,allow_pickle=False) as z:coarse={k:z[k].copy() for k in z.files}
   if {'truth','gt','m_true','m_truth','Y_clean','family','seed','mask'} & set(coarse):raise ValueError('Unexpected private operator payload')
   if (coarse['G'].shape,coarse['E0'].shape,coarse['S'].shape)!=( (576,576),(576,12),(48,576)):raise ValueError('Original coarse acquisition shape')
   if not np.isclose(float(coarse['k']),2*np.pi,rtol=0,atol=1e-14) or not np.isclose(float(coarse['radius']),.018,rtol=0,atol=1e-14):raise ValueError('Original wavelength/disk area differs')
   axis=(np.arange(48)+.5)*.96/48-.48;yy,xx=np.meshgrid(axis,axis,indexing='ij');centers=np.column_stack((xx.ravel(),yy.ravel()))
   geom=dict(centers=centers,sources=coarse['sources'].copy(),receivers=coarse['receivers'].copy(),radius=.45*(.96/48),k=float(coarse['k']))
   fine=Model(geom) # Existing true finite disk-pulse operator, one shared setup.
   dump(private/'FROZEN_GENERATOR_CONFIG.json',cfg)
   np.savez_compressed(private/'FINE_OPERATOR_GEOMETRY.npz',**geom)
   children=np.random.SeedSequence(cfg['seed']).spawn(8)
   for row,child in zip(cfg['cases'],children):
    beg=time.process_time();before=dict(fine.counts);geometry_child,noise_child=child.spawn(2)
    m,params=draw_material(np.random.default_rng(geometry_child),row['family'],row['regime'],12.,side=48)
    if not np.isfinite(m).all() or m.min()<.02 or m.max()>20:raise ValueError('Invalid fixed draw; no reroll')
    projected=m.reshape(24,2,24,2).mean(axis=(1,3)).ravel()
    clean=fine.forward(m,label='CROSSGRID8_GENERATOR_ONLY')['y'].copy();sigma=float(.01*np.linalg.norm(clean)/np.sqrt(clean.size));rng=np.random.default_rng(noise_child)
    noise=sigma/np.sqrt(2)*(rng.normal(size=clean.shape)+1j*rng.normal(size=clean.shape));cid=row['case_id'];observation=pub/(cid+'.npz')
    np.savez_compressed(observation,Y=clean+noise,sigma=sigma,noise_radius=float(np.sqrt(chi2.ppf(.99,2*clean.size)/2)),noise_model=np.array('iid_circular_complex'),grid_shape=np.array([24,24]),material_bounds=np.array([.02,20.]))
    np.savez_compressed(private/(cid+'.npz'),m=projected,m_fine=m,Y_clean=clean,fine_shape=np.array([48,48]),coarse_shape=np.array([24,24]))
    metadata.append(dict(**row,parameters=params,case_spawn_key=list(child.spawn_key),geometry_spawn_key=list(geometry_child.spawn_key),noise_spawn_key=list(noise_child.spawn_key),projection='Exact2x2 arithmetic fine-cell average, not exact coarse forward equivalence'))
    records.append(dict(case_id=cid,public_files=[(pub/'operators.npz').relative_to(REPO).as_posix(),observation.relative_to(REPO).as_posix()]))
    after=dict(fine.counts);rec['per_object'].append(dict(case_id=cid,CPU_seconds=time.process_time()-beg,physical_matrix_dimension=2304,physical_cost={k:after.get(k,0)-before.get(k,0) for k in set(before)|set(after)}));rec['completed']=len(records);dump(receipt,rec)
   dump(private/'generator_metadata.json',dict(seed=cfg['seed'],records=metadata,representation='Projected material target only; fine disk-data/coarse disk-model mismatch'))
   dump(pub/'manifest.json',dict(status='GENERATED_PUBLIC_OBSERVATIONS',records=records,inverse_input_policy='Only coarse operators, fine-generated noisyY, rawsigma/CN radius and legal24-grid metadata; no fine material or family/seed'))
   tasks=[]
   for i,row in enumerate(records):
    methods=cfg['methods'][i%4:]+cfg['methods'][:i%4]
    for method in methods:
     tid=row['case_id']+'__'+method+'__480CPU';t=dict(task_id=tid,case_id=row['case_id'],arm='a_plus' if method=='fullG19' else method,dpcs_rule='noise_ball_99',initialization='uniform',initial_material=.1,public_files=row['public_files'],config=a.solver_config.absolute().relative_to(REPO).as_posix(),output_dir=(portable(cfg['inverse_output'])/tid).as_posix(),experiment_scope='SUPPLEMENTARY_INDEPENDENT_FINE_DISK_DATA_COARSE_DISK_MODEL')
     if method=='fullG19':t['material_readout']=dict(module_file='research/research_max_struct_r11/src/material_readouts.py',callable='full_material_readout',prepare='prepare_full',label='fullG19')
     tasks.append(t)
   dump(manifest,dict(status='PUBLIC_TASKS_NO_GT',input_base='.',tasks=tasks));rec['status']='GENERATED_SEPARATE_FINE_TRUTH_AND_PROJECTED_TARGET'
 except BaseException:
  rec.update(status='FAILED_PARTIAL_GENERATION_PRESERVED',traceback=traceback.format_exc());raise
 finally:
  rec.update(finished_utc=datetime.now(timezone.utc).isoformat(),whole_child_CPU_seconds_near_write=time.process_time(),wall_seconds=time.perf_counter()-wall,physical_cost=dict(fine.counts) if fine else {},no_coarse_truth_forward_before_prediction_seal=True);dump(receipt,rec)
 print(json.dumps(dict(status=rec['status'],completed=rec['completed'],new_inverse_runs=0)))
if __name__=='__main__':main()

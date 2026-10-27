"""Freeze -> zero-RHS predictions -> seal -> CDEV-only offline score."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
import sys,json,time,csv,hashlib,datetime,resource
from pathlib import Path
import numpy as np
HERE=Path(__file__).resolve().parent;OUT=HERE.parents[1];sys.path.insert(0,str(OUT));sys.path.insert(0,str(HERE))
from closure_action import Action
from track_c.anchors import CONFIGS,REVISION1_CONFIGS
from track_c.frozen_heads import packed_measurement,SEEDS
from p3r2.residual_risk import online_features
from p3r2.risk_models import observation_features
from p3r4.risk import feature_image
NAMES=list(CONFIGS)+list(REVISION1_CONFIGS);MODES=['q','zero','negativeq'];METHOD='deterministic_matched';BASE='native96_probit64';CPU0=time.process_time();WALL0=time.perf_counter()
def write(name,v):(HERE/name).write_text(json.dumps(v,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else x.tolist())+'\n')
def guard():
 if time.process_time()-CPU0>112:raise TimeoutError('112CPU seconds conservative guard for120 requested')
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>1_900_000_000:raise MemoryError('under2GB memory guard')
def contexts(parents,name):
 values=[];reports=[]
 for parent in parents:
  folder=OUT/'track_c/contexts/CDEV'/parent/name
  with np.load(folder/'CONTEXT.npz') as z:values.append({k:z[k] for k in ['observation','residual','q0','chi0','J']})
  rr=json.loads((folder/'REPORT.json').read_text());assert rr['GT_read'] is False and rr['path_complete'];reports.append(rr)
 return values,reports

def diagnostic(action,values,name):
 rows=[];sparse=[];g=action.geometry;cfg=action.cfg;A=action.A;scale=action.normalization['image_scale']
 for i,c in enumerate(values):
  aa={k:v[None] for k,v in c.items()};rw=observation_features(aa['residual'])*(np.sqrt(2)/cfg['noise']['sigma']);j=np.asarray(aa['J'],np.float32)
  if not g['j_is_whitened']:j=j*(np.sqrt(2)/cfg['noise']['sigma'])
  grad=np.einsum('nri,nr->ni',j,rw,optimize=True)[:,:64];rawtrue=action.raw_image(aa,grad);rawq=action.raw_image(aa);diff=rawq[:,4:6]-rawtrue[:,4:6];z=np.asarray(aa['q0'],np.float32)[0,:64];err=float(np.linalg.norm(grad[0]-z));f=cfg['field_mean']+A@z
  from p1_material.prior import saturation_prime
  df=saturation_prime(f/cfg['phase_transition_halfwidth'])/cfg['phase_transition_halfwidth']*action.mask;q=aa['q0'][0];delta=complex(cfg['inclusion_increment_real_interval'][0]+np.diff(cfg['inclusion_increment_real_interval'])[0]*q[64],cfg['inclusion_increment_imag_interval'][0]+np.diff(cfg['inclusion_increment_imag_interval'])[0]*q[65]);B=abs(delta)*df[:,None]*A;op=float(np.sqrt(max(np.linalg.eigvalsh(B.T@B)[-1],0)));bound=op*err
  rows.append(dict(config=name,parent_index=i,KKT_g_minus_z_norm=err,g_norm=float(np.linalg.norm(grad[0])),z_norm=float(np.linalg.norm(z)),raw_back_difference_norm=float(np.linalg.norm(diff)),normalized_back_difference_norm=float(np.linalg.norm(diff/scale[:,4:6])),decoder_operator_norm=op,Lipschitz_bound=bound,bound_slack=bound-float(np.linalg.norm(diff)),cosine_g_z=float(grad[0]@z/max(np.linalg.norm(grad[0])*np.linalg.norm(z),1e-30))))
  if i in [0,17,63]:
   yy=np.array([packed_measurement(c['observation'])]);original=feature_image(yy,g,aa,action.normalization);reconstructed=action.image(aa,grad);error=float(np.max(abs(original-reconstructed)));assert np.allclose(original,reconstructed,rtol=2e-5,atol=2e-6);sparse.append(dict(config=name,parent_index=i,max_abs_normalized_image_error=error,rtol=2e-5,atol=2e-6,status='PASS',scope='exact originalgradient supplied to analyticalimage; noJproxy'))
  if i%16==0:guard()
 return rows,sparse

def bootstrap(v):
 result=np.empty((20000,*v.shape[1:]));rng=np.random.default_rng(202610270608)
 for begin in range(0,20000,256):
  n=min(256,20000-begin);ids=rng.integers(0,64,(n,64));result[begin:begin+n]=v[ids].mean(1)
 return result

def main():
 parents=[r['parent'] for r in json.loads((OUT/'configs/parents_CDEV.json').read_text())['parents']];assert len(parents)==64 and all(p.startswith('R6_CDEV_') for p in parents)
 if (HERE/'ALL_PREDICTIONS_SEALED.json').exists():raise RuntimeError('already sealed; noautomatic rerun or overwrite')
 write('CONFIG_FROZEN.json',dict(freeze_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),method=METHOD,supported_inference_methods=[METHOD,'distilled_w128_n32768'],configs=NAMES,seeds=list(SEEDS),closure_modes=MODES,replacement='only normalized residual-adjoint fielddirection; g≈+z',unchanged=['anchor','decoderweak/strongsensitivity','u2','normalization','frozenhead','frozenprovider'],noise_ratio=1.,future_noise_formula='g_feature≈(sigma_new/sigma_featurecfg)*z',known_approximation='exactonlyatfieldstationarity; noJconstructedormutated; inferenceusesq0/chi0/y',data='CDEV64 alreadyexposed development cohort; all newpredictions sealedbeforecurrentofflineGTread',baseline='previouslysealed native96_probit64 truegradient fullDet, reused',bootstrap_replicates=20000,bootstrap_seed=202610270608,risk_upper_one_sided_confidence=1-.05/3,NI_ratio_upper_max=1.01,wall_saving_lower_min=.20,physics_RHS=0,new_training_jobs=0))
 oldseal=json.loads((OUT/'track_c/evaluation_CDEV/ALL_PREDICTIONS_SEALED.json').read_text());baseline=[]
 for seed in SEEDS:
  row=next(r for r in oldseal['artifacts'] if r['config']==BASE and r['method']==METHOD and r['seed']==seed);path=OUT/'track_c/evaluation_CDEV'/row['path'];payload=path.read_bytes();assert len(payload)==row['bytes'] and hashlib.blake2b(payload).hexdigest()==row['blake2b']
  with np.load(path) as z:assert z['parent'].tolist()==parents;baseline.append({k:z[k] for k in ['prediction','inference_wall_seconds']})
 model_loads=[];artifacts=[];diag=[];audit=[];registry=[];anchorwall=[]
 # Load8contexts oneconfigurationatatime to keep memorybounded.
 for name in NAMES:
  vals,reports=contexts(parents,name);anchorwall.append([r['end_to_end_anchor_wall_seconds'] for r in reports]);registry.append(dict(config=name,source_context_template=str(OUT/'track_c/contexts/CDEV/<parent>'/name/'CONTEXT.npz'),n=64,reports_truth_read=False,registry='original5' if name in CONFIGS else 'revision1_3'))
  for seed in SEEDS:
   load=time.perf_counter();positive=Action(METHOD,seed,'q');model_loads.append(dict(config=name,seed=seed,wall_seconds=time.perf_counter()-load))
   if seed==SEEDS[0]:dd,aa=diagnostic(positive,vals,name);diag+=dd;audit+=aa
   for mode in MODES:
    action=positive if mode=='q' else Action(METHOD,seed,mode,_base=positive.base);pred=[];wall=[];cpu=[]
    for parent,c in zip(parents,vals):
     online={k:v[None] for k,v in c.items()};yy=np.array([packed_measurement(c['observation'])]);t=time.perf_counter();u=time.process_time();pred.append(action.predict(yy,online,[parent],batch_size=1)[0]);cpu.append(time.process_time()-u);wall.append(time.perf_counter()-t)
    filename=f'{name}__{mode}_s{seed}.npz';np.savez_compressed(HERE/filename,parent=np.array(parents),prediction=np.array(pred),inference_wall_seconds=wall,inference_CPU_seconds=cpu);payload=(HERE/filename).read_bytes();artifacts.append(dict(path=filename,config=name,mode=mode,seed=seed,bytes=len(payload),blake2b=hashlib.blake2b(payload).hexdigest()));guard()
   del positive
  print(json.dumps(dict(event='CONFIG_PREDICTED',config=name,CPU_seconds=time.process_time()-CPU0)),flush=True)
 write('SOURCE_CONTEXT_REGISTRY.json',registry);write('FEATURE_DIAGNOSTICS.json',dict(records=diag,sparse_original_formula_checks=audit,physics_RHS=0,GT_read=False));write('ALL_PREDICTIONS_SEALED.json',dict(status='SEALED_BEFORE_CURRENT_SCORE',CDEV_already_exposed=True,current_GT_read=False,artifacts=artifacts,baseline_original_seal=str(OUT/'track_c/evaluation_CDEV/ALL_PREDICTIONS_SEALED.json'),parents=parents,model_loads=model_loads,physics_RHS=0,created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat()))
 # Revalidate every sealed byte before the FIRST current truth open.
 for row in artifacts:payload=(HERE/row['path']).read_bytes();assert len(payload)==row['bytes'] and hashlib.blake2b(payload).hexdigest()==row['blake2b']
 truth=[]
 for parent in parents:
  with np.load(OUT/'raw/offline_truth'/f'{parent}.npz') as z:assert str(z['parent'])==parent;truth.append(z['chi'])
 truth=np.array(truth);baseerr=np.mean([np.mean(abs(r['prediction']-truth)**2,1) for r in baseline],0);baseanchor=anchorwall[0];basewall=np.array(baseanchor)+np.mean([r['inference_wall_seconds'] for r in baseline],0);risk=[];whole=[];labels=[];seedwise=[]
 for i,name in enumerate(NAMES):
  for mode in MODES:
   rs=[];tw=[]
   for seed in SEEDS:
    with np.load(HERE/f'{name}__{mode}_s{seed}.npz') as z:rs.append(np.mean(abs(z['prediction']-truth)**2,1));tw.append(z['inference_wall_seconds'])
   risk.append(np.mean(rs,0));whole.append(np.array(anchorwall[i])+np.mean(tw,0));labels.append((name,mode));seedwise.append(rs)
 risk=np.array(risk).T;whole=np.array(whole).T;bb=bootstrap(np.column_stack((baseerr,basewall,risk,whole)));summaries=[];m=len(labels)
 for j,(name,mode) in enumerate(labels):
  ratio=bb[:,2+j]/bb[:,0];saving=1-bb[:,2+m+j]/bb[:,1];ru=float(np.quantile(ratio,1-.05/3));sl=float(np.quantile(saving,.05/3));row=dict(config=name,closure_mode=mode,mean_seed_MSE=float(risk[:,j].mean()),baseline_truegradient_fullDet_MSE=float(baseerr.mean()),MSE_ratio=float(risk[:,j].mean()/baseerr.mean()),MSE_ratio_upper98p333=ru,whole_wall_seconds=float(whole[:,j].mean()),wall_saving=float(1-whole[:,j].mean()/basewall.mean()),wall_saving_lower98p333=sl,NI1pct=ru<=1.01,wall20pct=sl>=.2,development_joint_gate=ru<=1.01 and sl>=.2);summaries.append(row)
 np.savez_compressed(HERE/'PARENT_LOSSES.npz',parent=parents,configs=NAMES,modes=MODES,seeds=SEEDS,seedwise_MSE=np.array(seedwise),meanseed_MSE=risk,wholewall=whole,baseline_MSE=baseerr,baseline_wholewall=basewall)
 write('RESULTS.json',dict(status='DEVELOPMENT_KKT_CLOSURE_COMPLETE',all_arms=summaries,baseline_truegradient_fullDet_MSE=float(baseerr.mean()),baseline_whole_wall_seconds=float(basewall.mean()),bootstrap_replicates=20000,unit='independentparent; seedlossesaveragedfirst; notensembleimage',bootstrap_seed=202610270608,thresholds='unchanged NIupper<=1.01; wallsavinglower>=.20, one-sided98.333%',scope='alreadyexposed CDEVdevelopment, multiplefixedcontrols; noindependentconfirmation, noJacobianidentityawayfromstationarity',wall_scope='actualcachedanchorreceipts+singleparentinference; nofictionalrefundforJ; currentloadtiming',physics_RHS=0,new_training_jobs=0,oldtest_or_Dtruth_read=False))
 write('JOB.json',dict(status='COMPLETE',CPU_seconds=time.process_time()-CPU0,wall_seconds=time.perf_counter()-WALL0,CPU_threads=1,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physics_RHS=0,new_training_jobs=0));print(json.dumps(summaries,indent=2),flush=True)
if __name__=='__main__':
 try:main()
 except BaseException as e:write('JOB.json',dict(status='STOPPED',error=repr(e),CPU_seconds=time.process_time()-CPU0,wall_seconds=time.perf_counter()-WALL0,physics_RHS=0));raise

"""Freeze independent R4 authority and reusable P1 data without a physics query."""
from pathlib import Path
import sys, json, shutil, subprocess, os, platform, time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,CAPS,STAGE_WALL,STAGE_GPU,write,save,Book,pack
import numpy as np
from scipy.special import ndtri,ndtr

def main():
    old=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
    a=OUT/'configs/AUTHORIZATION.json'
    if a.exists():raise FileExistsError('authorization already frozen')
    write(a,dict(campaign_id='P3LR4-20261025',explicit_user_authorization=True,
       authority='Direct user Research P3L-R4 task section0.1 and13; overrides package static-only text',
       RHS_caps=CAPS,stage_controlled_wall_seconds=STAGE_WALL,stage_GPU_allocation_seconds=STAGE_GPU,
       total_RHS=58000000,total_controlled_wall_seconds=338400,total_GPU_allocation_seconds=345600,
       cpu_wall_definition='sum controlled job wall; actual CPU/core time separate; parallel jobs each charged',
       gpu_definition='conservative whole MPS training job wall; no CUDA or paid cloud',
       old_balances_used=False,new_3D=False,new_Mie=False,paid_cloud=False))
    write(OUT/'configs/PROTOCOL.json',dict(
      campaign_id='P3LR4-20261025',target_id='P1-CORRELATED-HOST-MIXTURE-TM24-v1-20261018',
      statistical_target='FULL latent posterior dx66; material law is decoder pushforward including plateaus',
      coordinate='x=(z64,ndtri(u2)), exact N(0,I66) prior; no material-volume factor',
      likelihood='Phi=complex squared residual/sigma^2=real whitened squared residual/2',
      observed_packing='8 source rows, each Re32 then Im32, 512 real entries',
      historical_R2_commit='6eadb4690c962aaa95821fc4056b677f0f924a6b',
      old_panels='all historical DEV and both R2 256 panels EXPOSED',
      joint_training_sizes=[4096,16384,32768],validation_parents=2048,
      initial_anchor_context=[512,128],expanded_anchor_context_max=[2048,512],
      new_development=dict(n=128,seed=20261025010,prefix='R4DEV'),
      new_locked=dict(n=512,seed=20261025020,prefix='R4L1',generate_after_method_freeze=True),
      independent_confirmation=dict(n=256,seed=20261025021,prefix='R4L2',generate_after_separate_freeze=True),
      training_seeds=[2026102501,2026102502,2026102503],
      density=dict(families=['full_cov_gaussian','spline_full66','triangular_full66','CFM_full66'],
        input='all raw512 measurements, train-only scaling; paid anchor conditioning separate',
        initial_updates=6000,batch=256,early_stop='validation only including update0',
        capacity_max=2000000,initial_capacity_target=[300000,800000],CFM_NFE=[16,32,64]),
      risk=dict(families=['frozen_R2_CNN','same_capacity_local_CNN','phase8','phase16','fullfield_moment','affine','spatial','monotone_calibration'],
        initial_updates=3000,batch=64,capacity_max=30000,phase_target='joint E[delta(u)*1{s>t}|context], never factorized'),
      generation=dict(decoded_samples=128,fixed_sample_index=0,full_likelihood_samples=64,
        methods_max=3,likelihood_SBC_parent_subset=128),
      correction=dict(development_parents=8,candidates=1024,expand_only_if_ESS_at_least=50,
        max_weight_below=.2,no_GT_initialization=True),
      local_geometry=dict(ranks=[2,4,8],material_RMS=[.02,.04],stress=.08,
        finite_difference_scales=[.001,.002,.004],compare_direct_bracket=True,
        candidate_only_unless_full_rank_MH_corrected=True),
      evaluation=dict(primary_metric='full576 complex material MSE; parent is statistical unit',
        primary_comparisons=['validation-frozen risk winner vs frozen R2 CNN','full66 density vs same-information conditional Gaussian','Maxwell context vs same-capacity decoder/random control'],
        bootstrap_draws=4000,bootstrap_seed=20261025099,primary_CI=.9833333333333333,
        positive_risk_relative_target=.02,physics_speed_target=.2,
        probability=['full66 NLL','material energy score','fixed material projections and phase-event proper scores','data-dependent SBC','one fixed draw per parent Gamma256 residual law'],
        posterior_certification='diagnostics only; no global certificate from tests'),
      GT_isolation='TRAIN supervision permitted; online loads allowlisted y/context only; locked truth scored in separate process after inference closed',
      memory=dict(physical_workers=1,training_workers=1,transfer_cache=2,numerical_threads=1,
        physical_process_RSS_limit_bytes=3221225472,training_RSS_limit_bytes=6442450944),
      freeze_changes='development amendments written before locked observations; no threshold/noise/prior relaxation'))
    shutil.copyfile(old/'configs/TARGET.json',OUT/'configs/TARGET.json')
    (OUT/'inputs/background').mkdir(parents=True,exist_ok=True)
    shutil.copyfile(old/'inputs/background/weak.npz',OUT/'inputs/background/weak.npz')
    write(OUT/'SOURCE_INVENTORY.json',dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
       R2_fixed_commit='6eadb4690c962aaa95821fc4056b677f0f924a6b',R2_cache_root=str(old),
       theory_zip='/Users/migodam/Downloads/Research_P3LR3_Complete.zip',
       package_files=[dict(path=str(p.relative_to(OUT/'source_package')),bytes=p.stat().st_size) for p in (OUT/'source_package').rglob('*') if p.is_file()],
       checksum_checks=False,history_modified=False,old_physics_cost='R2 5,956,328 RHS sealed; paid background provenance retained',
       A30_root=str(ROOT.parents[1]/'A30/a30-maxwell-material-manifold')))
    book=Book('0','prepare_and_restore_data',max_wall=1800,analysis_only=True)
    try:
      from p1_material.prior import MaterialPrior
      cfg=json.loads((OUT/'configs/TARGET.json').read_text())
      pts=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points']
      prior=MaterialPrior(cfg,pts)
      rows=[]
      for name,source in [('train32768',old/'A/data/train32768.npz'),('validation',old/'A/data/validation.npz')]:
        with np.load(source,allow_pickle=False) as d:
          q=d['q'];x=np.c_[q[:,:64],ndtri(q[:,64:])];obs=d['observation'];parent=d['parent']
          y=np.concatenate((obs.real.transpose(0,2,1),obs.imag.transpose(0,2,1)),axis=2).reshape(len(q),512)
          s=(cfg['field_mean']+q[:,:64]@prior.A[prior.mask].T)/prior.temperature
          error=float(np.max(abs(ndtr(x[:,64:])-q[:,64:])))
          if not np.isfinite(x).all() or error>1e-14:raise ValueError('exact probit input contract')
          path=OUT/'inputs/joint'/f'{name}.npz'
          save(path,x=x.astype('float32'),y=y.astype('float32'),material=d['chi'].astype('complex64'),s=s.astype('float32'),parent=parent,
               source=np.array(str(source)),target_id=np.array(cfg['target_id']))
          if name=='train32768':
            for n in (4096,16384):
              save(OUT/'inputs/joint'/f'train{n}.npz',x=x[:n].astype('float32'),y=y[:n].astype('float32'),material=d['chi'][:n].astype('complex64'),s=s[:n].astype('float32'),parent=parent[:n],target_id=np.array(cfg['target_id']))
          rows.append(dict(name=name,n=len(q),source=str(source),bytes=source.stat().st_size,probit_CDF_roundtrip_max=error,restored_path=str(path)))
      write(OUT/'AUDIT/DATA_RESTORE.json',dict(status='RESTORED_FROM_PAID_LOCAL_CACHE',datasets=rows,
           full_context_subset_only=[512,128],new_Maxwell_RHS=0,network_features_float32=True,physical_observation_source='original complex128 online files, not rounded NN features'))
      write(OUT/'ENVIRONMENT.json',dict(platform=platform.platform(),logical_CPU=os.cpu_count(),memory_bytes=25769803776,
            free_disk_bytes=shutil.disk_usage(ROOT).free,numerical_threads=1,MPS_available=True,CUDA=False))
      book.close()
    except BaseException as e:
      book.close('FAILED',repr(e));raise
    (OUT/'START_HERE.md').write_text('# Research P3L-R4：完整条件输运、材料风险与物理计算\n\n状态：RUNNING。独立授权，原 P1 统计合同保持不变。完整概率、材料风险、Maxwell 代表和几何候选分别判定；历史结果保持封账。真实训练及新物理实验在独立追加式账本下执行。\n')
    print('R4 independent protocol frozen; paid data restored; 0 new Maxwell RHS')

if __name__=='__main__':main()

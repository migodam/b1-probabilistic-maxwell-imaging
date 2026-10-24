#!/usr/bin/env python3
"""Read-only paid-response audit. No Target, backend, solve, fit or truth read."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import sys,json,csv,time,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import ndtri,ndtr
from p1_material.prior import MaterialPrior
from p3r4.core import OUT,Book,write,save,pack

GROUPS=('gaussian_s1','spline_s1','triangular_s1','cfm_s1__nfe32')
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(16))
SOURCE=OUT/'results/probability_physics_first/physics'
DEST=OUT/'AUDIT/probability_response_mechanism'

def csv_write(path,rows):
 if not rows:return
 with path.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def white_rows(values,sigma):
 return np.array([pack(z)*np.sqrt(2)/sigma for z in values])

def main():
 if DEST.exists():raise FileExistsError('fresh audit only; do not overwrite closed evidence')
 gate=json.loads((SOURCE/'PHYSICS_COMPLETE.json').read_text())
 if gate['status']!='COMPLETE' or gate['truth_read'] is not False:raise ValueError('closed observation-only physical source required')
 if json.loads(Path(gate['receipt']).read_text())['status']!='COMPLETE':raise ValueError('source Book not closed')
 b=Book('I','probability_response_mechanism_v1',max_wall=600,analysis_only=True);status='COMPLETE';error=None;rows=[];sample_rows=[];missing=[]
 try:
  DEST.mkdir();cfg=json.loads((OUT/'configs/TARGET.json').read_text());sigma=cfg['noise']['sigma']
  points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];prior=MaterialPrior(cfg,points)
  write(DEST/'FROZEN_PROTOCOL.json',dict(groups=GROUPS,parents=PARENTS,sample_indices=list(range(64)),uniform_original_index_weights=True,
   response_center='mean original64 physical complex predictions',decomposition='DF=C0 Dg_x0 dx + C0(dchi-Dg_x0 dx) + (DF-C0 dchi)',
   C0='receiver_T0[q,p]*field0[p,s], receiver/source packing explicitly converted by core.pack',coordinates='full66 raw probit, analytic native derivative times property normal PDF',
   identity='mean Phi=Phi(mean F)+mean ||F-mean F||²/sigma²',source_gate=str(SOURCE/'PHYSICS_COMPLETE.json'),source_receipt=gate['receipt'],
   source_physics_RHS=gate['costs']['RHS'],anchor_role='online qualified context only, not posterior truth',GT_read=False,new_RHS=0,new_LU=0,new_states=0,NN_fits=0,
   caveats=['finite decoder secant includes coordinate/decoder nonlinearity, not intrinsic second derivative','signed cross terms retained; component norm ratios are not causal percentages','sample mean response is not Full response at a material or latent mean','all development data exposed; no posterior certification']))
  for parent in PARENTS:
   anchor_path=OUT/'results/contexts/cache'/f'{parent}.npz';anchor=None
   if anchor_path.exists():
    with np.load(anchor_path,allow_pickle=False) as z:
     required=('q','chi','prediction','field','receiver_T')
     absent=[key for key in required if key not in z.files]
     if not absent:anchor={key:z[key].copy() for key in required}
     else:missing.append(dict(parent=parent,status='MISSING',reason='anchor fields '+','.join(absent),path=str(anchor_path)))
   else:missing.append(dict(parent=parent,status='MISSING',reason='anchor cache',path=str(anchor_path)))
   if anchor is not None:
    q0=anchor['q'];x0=q0.copy()
    if np.any((q0[64:]<=0)|(q0[64:]>=1)):raise ValueError('nonfinite anchor probit; no clipping')
    x0[64:]=ndtri(q0[64:]);decoded,_=prior.decode(q0)
    if not np.array_equal(decoded,anchor['chi']):raise ValueError('cached anchor decoder mismatch')
    D=prior.derivative(q0)*np.r_[np.ones(64),np.exp(-.5*x0[64:]**2)/np.sqrt(2*np.pi)][None,:]
    C=np.einsum('qp,ps->qsp',anchor['receiver_T'],anchor['field']).reshape(256,576)
   for group in GROUPS:
    folder=SOURCE/group/parent;record=json.loads((folder/'SUMMARY.json').read_text())
    if record['sample_indices']!=list(range(64)) or record['successful']!=64 or record['failed']:raise ValueError('all predeclared64 samples required; no survivor filtering')
    predictions=[];xs=[];chis=[];phis=[];paths=[];y=None
    for index in range(64):
     path=folder/f'sample_{index:03d}.npz';paths.append(str(path))
     with np.load(path,allow_pickle=False) as z:
      if int(z['index'])!=index:raise ValueError('sample index mismatch')
      f=z['prediction'].copy();x=z['x'].copy();q=z['q'].copy();chi=z['chi'].copy();obs=z['observation'].copy();phi=float(z['phi'])
     if not np.allclose(q[64:],ndtr(x[64:]),rtol=0,atol=2e-7) or not np.array_equal(x[:64],q[:64]):raise ValueError('sample coordinates mismatch')
     if y is None:y=obs
     elif not np.array_equal(obs,y):raise ValueError('different online observations')
     if f.shape!=(32,8) or chi.shape!=(576,) or not np.isfinite(f).all():raise ValueError('sample arrays invalid')
     predictions.append(f);xs.append(x);chis.append(chi);phis.append(phi)
    F=np.array(predictions);xs=np.array(xs);chis=np.array(chis);phis=np.array(phis);fbar=F.mean(axis=0)
    phi_computed=np.sum(abs(F-y)**2,axis=(1,2))/sigma**2;phi_bar=float(np.sum(abs(fbar-y)**2)/sigma**2)
    variance=float(np.mean(np.sum(abs(F-fbar)**2,axis=(1,2)))/sigma**2);mean_phi=float(phis.mean());identity_error=mean_phi-phi_bar-variance
    if np.max(abs(phi_computed-phis))>1e-8*max(1.,float(phis.max())):raise ValueError('saved Phi packing/units mismatch')
    if abs(identity_error)>1e-10*max(1.,mean_phi):raise ValueError('bias/dispersion identity failed')
    row=dict(model=group,parent=parent,n=64,mean_Phi=mean_phi,Phi_response_center=phi_bar,response_variance=variance,identity_error=identity_error,
     anchor_Phi=float(np.sum(abs(anchor['prediction']-y)**2)/sigma**2) if anchor is not None else None,telescope_status='COMPLETE' if anchor is not None else 'MISSING',
     latent_white_norm_mean=None,decoder_secant_white_norm_mean=None,Maxwell_remainder_white_norm_mean=None,exact_DF_white_norm_mean=None,
     latent_decoder_cross2_mean=None,latent_Maxwell_cross2_mean=None,decoder_Maxwell_cross2_mean=None,triangle_cancellation_gap_mean=None,
     response_energy_closure_max=None,telescope_vector_closure_max=None)
    arrays=dict(sample_indices=np.arange(64),source_sample_paths=np.array(paths),source_anchor_path=np.array(str(anchor_path)),prediction_mean=fbar,observation=y,
     Phi=phis,Phi_response_center=np.array(phi_bar),response_variance=np.array(variance),mean_Phi=np.array(mean_phi),identity_error=np.array(identity_error))
    if anchor is not None:
     dx=xs-x0;dchi=chis-anchor['chi'];material_linear=dx@D.T
     linear=(material_linear@C.T).reshape(64,32,8);decoder=((dchi-material_linear)@C.T).reshape(64,32,8)
     exact=F-anchor['prediction'];maxwell=exact-(dchi@C.T).reshape(64,32,8)
     components=np.stack([white_rows(linear,sigma),white_rows(decoder,sigma),white_rows(maxwell,sigma)],axis=1)
     exact_white=white_rows(exact,sigma);norms=np.linalg.norm(components,axis=2);exact_norm=np.linalg.norm(exact_white,axis=1)
     gram=np.einsum('nip,njp->nij',components,components);cross=2*gram[:,[0,0,1],[1,2,2]]
     closure=np.linalg.norm(exact_white-components.sum(axis=1),axis=1);energy=exact_norm**2-(np.sum(norms**2,axis=1)+cross.sum(axis=1))
     cancellation=norms.sum(axis=1)-exact_norm
     residual0=pack(anchor['prediction']-y)*np.sqrt(2)/sigma
     allcomponents=np.concatenate([np.broadcast_to(residual0,(64,1,512)),components],axis=1)
     phi_gram=np.einsum('nip,njp->nij',allcomponents,allcomponents);phi_closure=phis-.5*phi_gram.sum(axis=(1,2))
     if closure.max()>1e-9*max(1.,float(exact_norm.max())) or abs(phi_closure).max()>1e-8*max(1.,float(phis.max())):raise ValueError('telescope/phase energy audit failed')
     row.update(latent_white_norm_mean=float(norms[:,0].mean()),decoder_secant_white_norm_mean=float(norms[:,1].mean()),Maxwell_remainder_white_norm_mean=float(norms[:,2].mean()),exact_DF_white_norm_mean=float(exact_norm.mean()),
      latent_decoder_cross2_mean=float(cross[:,0].mean()),latent_Maxwell_cross2_mean=float(cross[:,1].mean()),decoder_Maxwell_cross2_mean=float(cross[:,2].mean()),triangle_cancellation_gap_mean=float(cancellation.mean()),
      response_energy_closure_max=float(abs(energy).max()),telescope_vector_closure_max=float(closure.max()))
     arrays.update(x0=x0,dx=dx,anchor_prediction=anchor['prediction'],anchor_chi=anchor['chi'],Dg_x0=D,exact_DF_white=exact_white,
      component_names=np.array(['latent_linear','finite_decoder_secant','finite_Maxwell_remainder']),component_vectors_white=components,component_norms_white=norms,
      component_gram_white=gram,cross2_order=np.array(['latent_decoder','latent_Maxwell','decoder_Maxwell']),signed_cross2_white=cross,
      exact_DF_norm_white=exact_norm,triangle_cancellation_gap=cancellation,telescope_vector_closure=closure,response_energy_closure=energy,
      residual0_white=residual0,Phi_component_gram_white=phi_gram,Phi_gram_component_names=np.array(['anchor_residual','latent_linear','finite_decoder_secant','finite_Maxwell_remainder']),Phi_energy_closure=phi_closure)
     for i in range(64):sample_rows.append(dict(model=group,parent=parent,sample_index=i,Phi=float(phis[i]),exact_DF_norm_white=float(exact_norm[i]),latent_norm_white=float(norms[i,0]),decoder_secant_norm_white=float(norms[i,1]),Maxwell_remainder_norm_white=float(norms[i,2]),
      latent_decoder_cross2=float(cross[i,0]),latent_Maxwell_cross2=float(cross[i,1]),decoder_Maxwell_cross2=float(cross[i,2]),triangle_cancellation_gap=float(cancellation[i]),vector_closure=float(closure[i]),energy_closure=float(energy[i]),Phi_energy_closure=float(phi_closure[i]),source=str(paths[i])))
    save(DEST/'vectors'/group/f'{parent}.npz',**arrays);rows.append(row)
   print(parent,flush=True)
  csv_write(DEST/'PER_PARENT.csv',rows);csv_write(DEST/'PER_SAMPLE_TELESCOPE.csv',sample_rows)
  means=[]
  for group in GROUPS:
   selected=[r for r in rows if r['model']==group];r=dict(model=group,parents=len(selected),samples=sum(z['n'] for z in selected))
   for key in rows[0]:
    if key not in ('model','parent','n','telescope_status'):r[key]=float(np.mean([z[key] for z in selected])) if all(z[key] is not None for z in selected) else None
   means.append(r)
  csv_write(DEST/'MODEL_MEANS.csv',means)
  result=dict(status='COMPLETE',groups=means,planned_parent_model_cells=64,completed_parent_model_cells=len(rows),samples=4096,telescope_samples=len(sample_rows),missing=missing,
   max_abs_center_identity_error=max(abs(r['identity_error']) for r in rows),max_abs_telescope_energy_closure=max(abs(r['energy_closure']) for r in sample_rows) if sample_rows else None,
   max_telescope_vector_closure=max(r['vector_closure'] for r in sample_rows) if sample_rows else None,source_physics_RHS=gate['costs']['RHS'],new_complex_RHS=0,new_Maxwell_states=0,new_physical_LU=0,neural_optimizer_steps=0,GT_read=False)
  write(DEST/'RESULTS.json',result)
  lines=['# 已付费 Full 响应的中心与离散度审计','', '这里审计四个原 full66 模型的前 16 个固定 DEV 父对象，每对象保留原索引 0–63 的全部 Full 样本，包括重复样本。没有读取真值、启动求解器或训练。','',
   '响应空间的恒等式是 mean Phi = Phi(Fbar) + mean ||Fk−Fbar||²/sigma²；Fbar 是物理响应均值，不是潜变量均值或材料均值的 Full 响应。中心项表示相对在线观测的响应失配，离散项表示这 64 个响应的经验离散度；两者都不能单独证明分布“过宽”或后验正确性。锚点 Phi 仅列作在线上下文。','',
   '|模型|mean Phi|Phi(Fbar)|响应离散项|锚点 Phi|','|---|---:|---:|---:|---:|']
  for r in means:lines.append(f"|{r['model']}|{r['mean_Phi']:.6f}|{r['Phi_response_center']:.6f}|{r['response_variance']:.6f}|{r['anchor_Phi'] if r['anchor_Phi'] is not None else 'MISSING'}|")
  lines += ['',f"64 个父对象×模型单元的恒等式最大绝对误差为 {result['max_abs_center_identity_error']:.3g}。",'',
   '有限响应伸缩分解逐样本保存 DF = C0 Dg_x0 dx + C0[(chi−chi0)−Dg_x0 dx] + [DF−C0(chi−chi0)]。C0 直接由已付费 T0 和总场 E0 数组相乘，Dg 使用原始解析解码器，并对两维 probit 属性乘正常密度链式因子。没有构造后端、读 LU 或求新逆列。三个向量依次是潜坐标线性响应、有限解码器割线余项、有限 Maxwell 响应余项。有限解码器割线包含坐标与解码非线性，不能称固有二阶曲率。','',
   '保存了白化向量、范数、带符号交叉项、三角抵消差、三分量 Gram 以及包含锚点残差的四分量 Phi Gram。平方范数的精确重建包含交叉项，因此范数大小或范数占比不是因果贡献比例；强抵消时尤其不能这样解释。', '',
   f"完成 {len(sample_rows)} 个有限分解；缺失锚点/字段 {len(missing)} 个，缺失项不会触发补求解。最大向量闭合误差 {result['max_telescope_vector_closure']}，最大响应能量绝对闭合误差 {result['max_abs_telescope_energy_closure']}。",'',
   '机器结果：PER_PARENT.csv、MODEL_MEANS.csv、PER_SAMPLE_TELESCOPE.csv；vectors/<model>/<parent>.npz 保存完整向量和逐样本来源。输入为 results/probability_physics_first/physics/ 的已关闭样本和 results/contexts/cache/ 的原在线锚点；来源总物理成本 36,352 RHS，本审计不再次计费。新 RHS、状态、物理 LU、神经优化均为 0；本次受控/CPU/RSS 成本见 I/jobs/probability_response_mechanism_v1/CLOSED.json。解释与最终科学判断由主线程保留。']
  (DEST/'REPORT_ZH.md').write_text('\n'.join(lines)+'\n')
 except BaseException as exc:status='ERROR';error=repr(exc);write(DEST/'ERROR.json',dict(error=error,traceback=traceback.format_exc()));raise
 finally:b.close(status,error)

if __name__=='__main__':main()

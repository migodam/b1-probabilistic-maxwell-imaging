#!/usr/bin/env python3
"""Final no-fit first64TRAIN softmax conditioning diagnostic, no phase rescue."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import argparse,json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np,torch
from p3r4.risk_anchorbin import AnchorBinRisk,EPSILON
from p3r4.risk_bin import bin_phase_labels
from p3r4.core import Book,OUT,write,save


def summarize(values,mask):
    count=mask.sum(1);means=np.divide((values*mask).sum(1),count,out=np.full(len(count),np.nan),where=count>0)
    finite=np.isfinite(means)
    return dict(cells_total=int(count.sum()),parents_nonempty=int(finite.sum()),parent_mean=float(means[finite].mean()) if finite.any() else None,
        pixel_mean_descriptive=float(values[mask].mean()) if mask.any() else None,parent_means=means.tolist())


def clean(v):
    if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
    if isinstance(v,list):return [clean(x) for x in v]
    if isinstance(v,float) and not np.isfinite(v):return None
    return v


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--job',default='phase_head_conditioning_first64_train_v1');a=p.parse_args()
    result=OUT/'AUDIT/risk/phase_conditioning_first64_v1'
    if result.exists():raise FileExistsError('one final diagnostic, no repeated tuning')
    book=Book('I',a.job,max_wall=120,analysis_only=True);status='COMPLETE';error=None;start=time.monotonic()
    try:
        torch.set_num_threads(1)
        old=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
        with np.load(old/'A/anchor_data/train.npz',allow_pickle=False) as z:parents=z['parent'][:64].astype(str);q0=z['q0'][:64];chi0=z['chi0'][:64]
        if parents.tolist()!=[f'R2_train_{i:05d}' for i in range(64)]:raise ValueError('fixed first64TRAIN only')
        saved=torch.load(old/'A/residual_training/n512/local_geometry_cnn_s33/BEST.pt',map_location='cpu',weights_only=False);g=saved['geometry'];cfg=g['config'];roi=g['roi']
        with np.load(OUT/'inputs/joint/train32768.npz',allow_pickle=False) as z:
            lookup={v:i for i,v in enumerate(z['parent'].astype(str))};ids=np.array([lookup[v] for v in parents]);x_GT=z['x'][ids];s_GT=z['s'][ids]
        from p1_material.prior import MaterialPrior
        points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];prior=MaterialPrior(cfg,points)
        s0=(cfg['field_mean']+q0[:,:64]@prior.A[roi].T)/prior.temperature
        masks=dict(host=s0<=-1,transition=(s0>-1)&(s0<1),pure=s0>=1,all=np.ones_like(s0,dtype=bool))
        result.mkdir(parents=True);allrows=[];arrays={}
        for k in (8,16):
            model=AnchorBinRisk('phase'+str(k),g)
            with torch.no_grad():anchor_p,_,_=model.baseline(torch.tensor(q0,dtype=torch.float32),torch.tensor(chi0,dtype=torch.complex64))
            labels=bin_phase_labels(x_GT,g,model.thresholds.numpy(),s_GT).astype('float64');B=model.bounds.numpy().astype('float64');local={}
            for name,p0 in [('uniform',np.full((64,2,k+1,256),1/(k+1))),('anchor',anchor_p.numpy().astype('float64'))]:
                p0/=p0.sum(2,keepdims=True);tail=np.flip(np.cumsum(np.flip(p0[:,:,1:],axis=2),axis=2),axis=2);M=B[None,:,None,None]*tail
                indicator=(np.arange(k+1)[None,:]>np.arange(k)[:,None]).astype('float64')
                J=B[None,:,None,None,None]*p0[:,:,None,:,:]*(indicator[None,None,:,:,None]-tail[:,:,:,None,:])
                jac=np.sqrt(np.sum(J*J,axis=(2,3))) # [parent,component,cell]
                integrated=(J/k).sum(2);integrated_norm=np.sqrt(np.sum(integrated*integrated,axis=2))
                # Identical actual weighted squared moment objective on this64
                # TRAIN batch; no optimizer/parameter update is performed.
                lossgrad=2*np.sum((M-labels)[:,:,:,None,:]*J/k,axis=2)/(64*2*256)
                grad=np.sqrt(np.sum(lossgrad*lossgrad,axis=2));local[name]=dict(jacobian=jac,integrated=integrated_norm,loss_gradient=grad)
                for metric,values in local[name].items():arrays[f'k{k}_{name}_{metric}']=values
            row=dict(K=k,formula='B*p_k*(1{k>j}-tail_j)',moment_rule='exact_equal_h_mass_bins',epsilon=EPSILON,strata={})
            for phase,mask in masks.items():
                detail={}
                for metric in ('jacobian','integrated','loss_gradient'):
                    detail[metric]={}
                    for component,c in [('real',0),('imag',1),('combined',None)]:
                        u=local['uniform'][metric][:,c] if c is not None else np.sqrt(np.sum(local['uniform'][metric]**2,axis=1))
                        b=local['anchor'][metric][:,c] if c is not None else np.sqrt(np.sum(local['anchor'][metric]**2,axis=1))
                        us,bs=summarize(u,mask),summarize(b,mask);ratio=None if us['parent_mean'] in (None,0) else bs['parent_mean']/us['parent_mean']
                        detail[metric][component]=dict(uniform=us,anchor=bs,ratio_of_parent_mean_anchor_to_uniform=ratio)
                row['strata'][phase]=detail
            allrows.append(row)
        save(result/'PER_PARENT_SCALES.npz',parent=parents,anchor_standardized_phase=s0,**arrays)
        write(result/'RESULTS.json',clean(dict(status='COMPLETE_FINAL_DIAGNOSTIC_NO_RESCUE',train_parents=parents.tolist(),rows=allrows,
            labels='same fixed first64TRAIN jointGT only, neveronline/test features',physics_RHS=0,training_updates=0,
            interpretation='finite gradients do not imply well-conditioned learning; nearonehot logits can suppress moment JacobianO(eps); finaldiagnostic not another initialization search',
            phase_family='STOP_EXPANSION_AFTER_SIX_REVISION2_FITS',wall_seconds=time.monotonic()-start)))
        report=['# 最终相矩诊断：有限梯度不等于易优化','',
            '固定最初64个TRAIN父样本，不拟合、不读取DEV/locked标签，不新增物理。原512/128训练context为48+48，与实际revision2六fits一致。比较uniform exact-bin与eps=1e−4 anchor-mass初始化，保持同一精确箱标签。','',
            '对每个阈值箱j及质量logit k，∂M_j/∂l_k=B p_k(1{k>j}−tail_j)。当anchor平台给出近one-hot质量时，多数质量p_k为O(eps)，tail_j接近0或1，故有限但很小的导数可能限制更新。该机制是softmax条件数现象，不是梯度NaN；也不能仅由它证明所有实际风险差来自初始化。','',
            '|K|anchor相区|联合moment-Jacobian幅度比|积分材料头导数幅度比|实际箱loss logit梯度幅度比|','|---|---|---|---|---|']
        for row in allrows:
            for phase in ('host','transition','pure'):
                values=[row['strata'][phase][metric]['combined']['ratio_of_parent_mean_anchor_to_uniform'] for metric in ('jacobian','integrated','loss_gradient')]
                report.append('|'+str(row['K'])+'|'+phase+'|'+'|'.join(f'{v:.6g}' if v is not None else 'NOT_EVALUATED' for v in values)+'|')
        report+=['','比值先按每parent非空相区取均值，再在parent上平均，未把像素当独立试验。第一项是完整K×(K+1)Jacobian Frobenius，第二项是1/K积分后材料头导数L2，第三项为同一64TRAIN batch真实加权箱MSE对logits的梯度L2。保留逐parent尺度和空相区，不用GT调eps。','',
            '结合父调度器已观察的revision2验证失败，停止phase家族扩展。此诊断不授权新初始化、学习率或额外训练救结果；finite-gradient PASS只证明导数数值有限，不能证明学习良态。']
        (result/'REPORT_ZH.md').write_text('\n'.join(report)+'\n');print(json.dumps(dict(status='COMPLETE_FINAL_DIAGNOSTIC_NO_RESCUE',out=str(result),physics_RHS=0,training_updates=0)))
    except BaseException as exc:status='ERROR';error=repr(exc);raise
    finally:book.close(status,error)

if __name__=='__main__':main()

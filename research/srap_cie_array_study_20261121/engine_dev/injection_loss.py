"""Resolve the saved first-cycle reinjection error into orthogonal current spaces."""
import argparse
import csv
import json
from pathlib import Path
import time
from check_interventions import (CASE_NAMES,HERE,HISTORICAL,engine,model_for,np,read_public)
from scipy.linalg import qr

ARMS=('srap','graph','full_ap')

def energy(x):
    return float(np.vdot(x,x).real)

def orthogonal_bases(groups):
    return [(ids,qr(V,mode='economic',check_finite=False)[0]) for ids,V in groups]

def project(Z,bases):
    out=np.zeros_like(Z)
    for ids,V in bases:
        out[:,ids]=V@(V.conj().T@Z[:,ids])
    return out

def read_one_case(name,public):
    p,c=read_public(name,public)
    model=model_for(p,c)
    modes,_,rank=engine.schedule(c)
    jd,groups=engine.dpcs(model,p['Y'],rank,False)
    _,cap_groups=engine.dpcs(model,p['Y'],rank,True)
    D=orthogonal_bases(groups)
    Dc=orthogonal_bases(cap_groups)
    U=engine.FourierCurrents(model.shape,modes[0],groups,model.clock.counts)
    with np.load(HERE/'checkpoints'/f'{name}__shared.npz') as s:
        j0=s['before_J']
    rows=[]
    for arm in ARMS:
        with np.load(HERE/'checkpoints'/f'{name}__{arm}.npz') as f:
            first=6*f['AP_component_x']
            jt=6*f['AP_sx']
            ji=f['AP_J']
        target=jt-j0
        loss=jt-ji
        first_D=project(first,D)
        target_D=project(target,D)
        loss_D=project(loss,D)
        loss_other=loss-loss_D
        loss_Dc=project(loss,Dc)
        total=energy(loss)
        change=energy(target)
        inn=np.vdot(loss_other,loss_D)
        # Check against the actual historical numerical Fourier operator, not an ideal FFT mask.
        normal=U.adjoint(loss)
        # Test membership of each kept current in the affine DPCS lock.
        lock_error=energy(project(ji-jd,D))/max(energy(ji),1e-300)
        rows.append(dict(case=name,method=arm,stage=0,step=16,
             source_DPCS_ranks='/'.join(str(b.shape[1]) for _,b in D),
             source_cap_ranks='/'.join(str(b.shape[1]) for _,b in Dc),
             correction_energy=change,loss_energy=total,
             locked_energy=energy(loss_D),other_unrepresentable_energy=energy(loss_other),
             loss_per_target_energy=total/change,
             locked_per_target_energy=energy(loss_D)/change,
             other_per_target_energy=energy(loss_other)/change,
             locked_share_of_loss=energy(loss_D)/total,
             other_share_of_loss=energy(loss_other)/total,
             droppable_locked_share_of_loss=energy(loss_D-loss_Dc)/total,
             raw_routed_direction_locked_share=energy(first_D)/max(energy(first),1e-300),
             target_direction_locked_share=energy(target_D)/change,
             normalized_normal_equation=float(np.sqrt(energy(normal)/total)),
             orthogonality_relative=float(abs(inn)/total),
             pythagoras_relative=abs(total-energy(loss_D)-energy(loss_other))/total,
             returned_lock_error_relative=float(np.sqrt(lock_error))))
    return rows

def table(head,rows):
    return '\n'.join('| '+' | '.join(map(str,r))+' |' for r in [head,['---']*len(head),*rows])

def main(public=HISTORICAL):
    start=time.perf_counter();cpu=time.process_time()
    rows=[r for n in CASE_NAMES for r in read_one_case(n,public)]
    wall=time.perf_counter()-start;cpused=time.process_time()-cpu
    with (HERE/'injection_loss.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    selected=[r for r in rows if r['method']=='srap']
    names=['180° Austria 旋转 90°','120° Austria','180° 三体','Fresnel Twin 8 GHz 半发射']
    def pc(x):
        return f'{100*x:.2f}%'
    summary=[]
    raw=[]
    for name,r in zip(names,selected):
        summary.append([name,pc(r['loss_per_target_energy']),pc(r['locked_per_target_energy']),pc(r['other_per_target_energy']),pc(r['locked_share_of_loss']),pc(r['other_share_of_loss'])])
        raw.append([name,pc(r['raw_routed_direction_locked_share']),pc(r['target_direction_locked_share']),pc(r['droppable_locked_share_of_loss'])])
    residual=max(r['normalized_normal_equation'] for r in rows)
    report=f"""# 傅里叶回注丢失了哪些电流

在四个开发问题的首个 SRAP 周期，DPCS 锁定空间占总丢失能量的 {min(r['locked_share_of_loss'] for r in selected)*100:.1f}%–{max(r['locked_share_of_loss'] for r in selected)*100:.1f}%。四例分别损失了目标电流改变量平方能量的 {', '.join(pc(r['loss_per_target_energy']) for r in selected)}。下表中的读数由已保存的双场和电流直接计算。能量指对所有空间单元和发射源求复数模平方和。

{table(['公开状态','总损失 / 目标变化','DPCS 锁定 / 目标变化','其余不可表示 / 目标变化','DPCS 占总损失','其余占总损失'],summary)}

## 两个不同的电流位置

每次分叉前的旧电流为 J₀，已在合法仿射空间 J_D + range(U) 内。先对交替投影给出的电流变化做路由，得到 **6ΠΔX**，相应场记忆为 ΠΔQ。第二次跨源材料投影要求每个单元共享实标量材料，把 (J₀/6+ΠΔX,Q₀+ΠΔQ) 投到该约束，给出新的 **J_target=6X_C**。这一非线性操作会改变电流方向。实际要回注的变化是 J_target−J₀。傅里叶最小二乘回注给出 J_injected=J_D+U a，损失是 L=J_target−J_injected。

对每组接收掩码单独正交化 DPCS 右奇异向量，构造 P_D。损失分解为 L=P_D L+(I−P_D)L。第一项是锁定的接收子空间，第二项位于其正交补；固定当前 Fourier 模式后，后者也无法由合法增量表示。12 个保存的分支中，实际算子的正规方程误差 ‖U*L‖/‖L‖ 最大为 {residual:.2e}，平方和分解相对误差最大为 {max(r['pythagoras_relative'] for r in rows):.2e}。这项检查直接使用源特定的 (I−V_D V_D*)F，未将所有高频分量当作丢失量。

{table(['公开状态','路由电流 6ΠΔX 中的 DPCS 能量占比','二次材料投影后目标改变量中的 DPCS 能量占比','0.1σ₁ 截断可解锁的空间占总损失'],raw)}

末列给出的空间占比采用旧秩上限与 0.1σ₁ 阈值保留集合的交集。计算时将原有锁定投影 P_D 再分为保留锁定的 P_cap 与被阈值释放的 P_D−P_cap。分母为原 SRAP 已丢失的能量，值为 ‖(P_D−P_cap)L‖²/‖L‖²。120° Austria 中，这部分占总损失 39.71%；同一案例另有 45.36% 落在 DPCS 正交补内。Fresnel 的 5 个既有 DPCS 方向全部满足阈值，所以这次截断可解锁的份额为零。

**因果范围。** 解除 DPCS 锁定的电流还须落入新的 Fourier 增量空间才可能回注。改变截断也会改变初始确定性电流 J_D、全程 PR 轨迹、第二次材料投影的目标和截止时间内可选材料。由此，表中的份额不能直接换算成最终材料误差改善。只比较同一目标方向的空间含量，也不能将全程截断实验的增益单独归因于回注。

共读取 4 份公开训练包、4 个共同检查点以及每例 SRAP、传播相容路由、全 AP 的三份分支，均在第 0 阶段第 16 步。这里没有执行反演迭代、重选材料或读取真值和留出场。数值计算与读取文件共耗墙钟 {wall:.3f} s、进程 CPU {cpused:.3f} s，不含 Python 与依赖库启动。[逐例逐分支绝对能量、归一化能量与残差](injection_loss.csv) · [复现代码](injection_loss.py)
"""
    (HERE/'INJECTION_LOSS.md').write_text(report)
    print(table(['case','total','lock/T','other/T','lock/L','other/L'],summary))
    print(table(['case','first lock','target lock','droppable/L'],raw))
    print(wall,cpused,residual)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--public',type=Path,default=HISTORICAL)
    main(p.parse_args().public)

"""One native96 exact-mean rather than flat-core postmeasurement refit."""
from fractions import Fraction
import json

from attack import HERE, BANDS, BareOpen, elapsed, np, read, real, strict_replay, timer, write_json
from expand import AllColumnJacobian, additional_candidates
from scipy.optimize import lsq_linear


SIDE=96
SAMPLE="specimen_013"
MEAN=2
MAX_ITER=16
TARGET=.975
DYADIC=2**20


def shift_to_sum(x,lo,hi,target):
    """Euclidean projection onto an interval box and one sum constraint."""
    left=float(np.min(lo-x))
    right=float(np.max(hi-x))
    for _ in range(52):
        shift=(left+right)/2
        current=np.clip(x+shift,lo,hi)
        if current.sum()<target:
            left=shift
        else:
            right=shift
    return np.clip(x+(left+right)/2,lo,hi)


def round_core(x):
    """Exact dyadic material sum with legal represented values."""
    ticks=np.clip(np.rint(x*DYADIC).astype(np.int64),0,int(4.5*DYADIC))
    remaining=int(len(x)*MEAN*DYADIC-int(ticks.sum()))
    if remaining>0:
        ix=int(np.argmin(ticks))
    elif remaining<0:
        ix=int(np.argmax(ticks))
    else:
        ix=0
    ticks[ix]+=remaining
    if ticks.min()<0 or ticks.max()>int(4.5*DYADIC) or ticks.sum()!=len(x)*MEAN*DYADIC:
        raise ArithmeticError("Could not represent the constrained mean exactly")
    return ticks.astype(float)/DYADIC,ticks


def main():
    total=timer()
    _,roi,_,readings,labels,r,energy,sources=read(SIDE)
    observed=readings[SAMPLE]
    initial=HERE/"expanded_96_specimen_013_2.npz"
    low=HERE/"candidate_96_specimen_013_0.npz"
    with np.load(initial,allow_pickle=False) as f:
        x=f["material"]
        assert np.all(x[roi]==MEAN)
    with np.load(low,allow_pickle=False) as f:
        zero=f["material"]
        assert np.array_equal(observed[0],f["observed_400"]) and np.array_equal(observed[1],f["observed_460"])
        assert np.all(zero[roi]==0.)
    before=json.loads(initial.with_suffix(".json").read_text())
    original_trace=json.loads((HERE/"candidate_96_specimen_013_2.json").read_text())
    with np.load(HERE/"candidate_96_specimen_013_2.npz",allow_pickle=False) as z:
        source_support=z["original_support"]
    added=np.concatenate([np.array(d["added_sites"],dtype=int) for d in before["steps"]])
    public_sites=np.unique(np.r_[source_support,added])
    models=[AllColumnJacobian(SIDE,b,e) for b,e in zip(BANDS,energy)]
    first=[w.evaluate(x) for w in models]
    initial_defect=[real(p-y)/r for (p,j),y in zip(first,observed)]
    j_start=[j/r for p,j in first]
    newsites=additional_candidates(SIDE,x,j_start,initial_defect,public_sites,80)
    support=np.unique(np.r_[roi,public_sites,newsites])
    local=np.searchsorted(support,roi)
    # Numerical predictions share a single material at both bands.
    forward=[BareOpen(SIDE,b,support,e) for b,e in zip(BANDS,energy)]
    actions=[]
    threshold=0.90
    radius=1.0
    minimum=float("inf")
    best=x.copy()
    for it in range(MAX_ITER):
        now=timer()
        outcomes=first if it==0 else [w.evaluate(x) for w in models]
        d=[real(p-y)/r for (p,j),y in zip(outcomes,observed)]
        js=[j[:,support]/r for p,j in outcomes]
        base=[float(np.linalg.norm(e)) for e in d]
        if max(base)<minimum:
            minimum=max(base)
            best=x.copy()
        print(it,"base",base,"core range",x[roi].min(),x[roi].max(),flush=True)
        if max(base)<TARGET:
            break
        # Put more local emphasis on the farther ball. The mean row is
        # enforced exactly by projection before any nonlinear trial.
        bias=float(np.clip((base[0]/base[1])**2,.4,2.5))
        row=np.zeros(len(support))
        row[local]=18./np.sqrt(len(local))
        robust=np.vstack((js[0]*bias,js[1],row))
        goal=-np.r_[d[0]*bias,d[1],0.]
        delta=lsq_linear(robust,goal,bounds=(np.maximum(-radius,-x[support]),
                                           np.minimum(radius,4.5-x[support])),
                         method="trf",lsq_solver="lsmr",lsmr_tol=8e-5,lsmr_maxiter=160,
                         max_iter=18,tol=2e-4)
        trials=[]
        for ts in (1.,.5,.25):
            candidate=np.clip(x[support]+ts*delta.x,0.,4.5)
            candidate[local]=shift_to_sum(candidate[local],np.zeros(len(local)),np.full(len(local),4.5),MEAN*len(local))
            fields=[m.fields(candidate,jac=False)[0] for m in forward]
            resid=[float(np.linalg.norm(p-y)/r) for p,y in zip(fields,observed)]
            trials.append((max(resid),candidate,resid))
            if max(resid)<threshold:
                break
        use=min(trials,key=lambda z:z[0])
        actions.append(dict(iteration=it,baseline=base,accepted=use[0]<max(base),
                           beta_trust=radius,old_ball_weight=bias,
                           trial_ratios=[p[2] for p in trials],linear_iterations=int(delta.nit),
                           best_new=use[2],**elapsed(now)))
        print(it,"trial",use[2],flush=True)
        if use[0]>=max(base)-1e-5:
            radius*=.52
            if radius<.08:
                break
        else:
            x=np.zeros(SIDE**2)
            x[support]=use[1]
            if use[0]<minimum:
                minimum=use[0]
                best=x.copy()
            if use[0]<TARGET:
                break
            if use[0]<max(base)*.85:
                radius=min(radius*1.3,2.)
    best[roi],ticks=round_core(best[roi])
    mean=sum(Fraction(v) for v in best[roi])/len(roi)
    if mean!=2 or not (best.min()>=0 and best.max()<=4.5):
        raise RuntimeError("ROI-mean equality or physical box failed")
    # The old point verifier checks flatness separately. Use a single marker
    # for that guard: the actual mean is tested exactly above. Its physical
    # calculation sees the complete material on the original 96² grid.
    certificate,fields=strict_replay(SIDE,best,observed,r,energy,roi[:1])
    low_certificate,low_fields=strict_replay(SIDE,zero,observed,r,energy,roi)
    lim=[certificate[str(b)]["strict_residual_over_r_upper"] for b in BANDS]
    lim0=[low_certificate[str(b)]["strict_residual_over_r_upper"] for b in BANDS]
    success=max(lim)<1
    low_ok=max(lim0)<1
    prefix="mean_task"
    with (HERE/(prefix+".npz")).open("wb") as f:
        np.savez_compressed(f,high_mean_material=best,low_mean_material=zero,roi_indices=roi,
                            exact_dyadic_core_ticks=ticks,dyadic_denominator=DYADIC,
                            allowed_trial_sites=support,additional_trial_sites=newsites,
                            observed_400=observed[0],observed_460=observed[1],**fields,
                            **{"low_"+k:v for k,v in low_fields.items()})
    receipt=dict(side=SIDE,sample=SAMPLE,earlier_calibrated_class=labels[SAMPLE],actual_observations=sources,
                 search="Mean-constrained, bounded Gauss–Newton along one continuation from the saved constant-core trial. ROI shape is free. Shared exterior supported on public and data-selected trial cells.",
                 start_file=initial.name,prior_paid_start_cpu_s=original_trace["cpu_s"]+before["cpu_s"],
                 exact_roi_mean=str(mean),dyadic_denominator=DYADIC,dyadic_tick_sum=int(ticks.sum()),
                 core_cells=len(roi),physical_core_area_m2=len(roi)*(2/SIDE)**2,
                 roi_range=[float(min(best[roi])),float(max(best[roi]))],
                 energies=energy,noise_radius=r,trial_steps=actions,strict_high=certificate,
                 strict_low_zero=low_certificate,high_status="verified in both balls" if success else "unresolved",
                 low_status="verified in both balls" if low_ok else "unresolved",
                 task_diameter_lower=2 if low_ok and success else None,
                 cost_full_columns={str(b):m.cost for b,m in zip(BANDS,models)},
                 cost_numerical_trials={str(b):m.cost for b,m in zip(BANDS,forward)},**elapsed(total))
    write_json(prefix+".json",receipt)
    hi_title="找到核心均值恰为 2 的共享双频材料，两个残差球均严格通过。" if success else "这次搜索得到一个均值恰为 2 的最佳点，高均值连续材料类仍未决。"
    conclusion="同一观测的可行核心均值直径至少为 2。" if success and low_ok else ""
    nums=np.ceil(np.array([*lim0,*lim])*1e7)/1e7
    text=rf"""# 96² 已采 specimen_013 的核心均值再攻击

{hi_title}{conclusion}

只复用两份实际观测。400 MHz 旧总波读数来自已封存的共同中点，460 MHz 使用原带噪读数与 11/16 发射功率。未访问离线材料。继续使用原名义 25 端开路、8 发 8 收以及两个独立绝对硬 Frobenius 球。两个频率均使用同一份材料。

在上一轮常数核心试算上，允许 164 个核心像素在 \([0,4.5]\) 内重新分配。一次均值约束 Gauss–Newton 轨迹同时调整外区。核心面积 \(0.0711805556\,\mathrm{{m^2}}\)。最后对核心用步长 \(2^{{-20}}\) 的二进制分数表示，使 164 格的总和精确等于 328。高候选核心格值范围为 [{receipt['roi_range'][0]:.6f}, {receipt['roi_range'][1]:.6f}]。按独立裸 VIE、25 端加载和开路闭合复核后，两端严格上界如下。

| 核心均值 | 400 MHz 严格残差上界 / \(r\) | 460 MHz 严格残差上界 / \(r\) |
|---:|---:|---:|
| 0 | {nums[0]:.7f} | {nums[1]:.7f} |
| 2 | {nums[2]:.7f} | {nums[3]:.7f} |

新增试算与两份严格复核合计 {receipt["cpu_s"]:.3f} 秒 CPU。原常数核心结果和已封存材料保持原样。搜索支持集及完整费用见收据。

[收据](mean_task.json) · [材料与原波幅](mean_task.npz) · [复现入口](mean_task.py)
"""
    (HERE/"mean_task_REPORT_ZH.md").write_text(text)
    print("STRICT",lim0,lim,"exact core mean",mean,"CPU",receipt["cpu_s"],flush=True)


if __name__=="__main__":
    main()

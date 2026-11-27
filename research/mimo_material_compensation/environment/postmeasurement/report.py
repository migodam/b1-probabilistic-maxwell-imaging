"""Verify the readout contract and summarize only actually attained bounds."""
from attack import HERE, ENV, R, RX, TX, np, fro, read, up, rd, at_load, strict_replay, elapsed, timer, write_json
from fractions import Fraction
import json


def summarize():
    start=timer()
    selected={}
    original={f.stem.removeprefix("candidate_"):json.loads(f.read_text())
              for f in HERE.glob("candidate_*.json")}
    extensions={f.stem.removeprefix("expanded_"):json.loads(f.read_text())
                for f in HERE.glob("expanded_*.json")}
    costs=[]
    for key,old in sorted(original.items()):
        side=old["side"]
        sample=old["specimen"]
        core=old["core_fixed"]
        verify=old["independent_audit"]
        case_old=HERE/f"candidate_{key}.npz"
        case_new=HERE/f"expanded_{key}.npz"
        ex=extensions.get(key)
        if ex and (verify is None or max(v["strict_residual_over_r_upper"] for v in ex["independent_audit"].values())
                   <max(v["strict_residual_over_r_upper"] for v in verify.values())):
            use=case_new
            verify=ex["independent_audit"]
        else:
            use=case_old
        with np.load(use,allow_pickle=False) as q:
            x=q["material"]
            roi=q["roi_indices"]
            ys=[q["observed_400"],q["observed_460"]]
        base=read(side)
        if not (np.array_equal(ys[0],base[3][sample][0]) and np.array_equal(ys[1],base[3][sample][1])
                and np.all(x[roi]==core) and np.min(x)>=0 and np.max(x)<=4.5):
            raise ArithmeticError("Candidate does not satisfy the originally loaded input contract")
        if not verify:
            verify,raw=strict_replay(side,x,ys,base[5],base[6],roi)
            costs.append(dict(sample_key=key,paid_extra_audit=verify))
            with (HERE/f"extra_audit_{key}.npz").open("wb") as f:
                np.savez_compressed(f,**raw)
        exact_distances={str(int(k)):max(abs(Fraction(a)-Fraction(b)) for a,b in zip(v,x))
                         for k,v in base[2].items()}
        beta_reference={k:[float(np.nextafter(float(v),-np.inf)),up(float(v))]
                        for k,v in exact_distances.items()}
        outside_both=all(v>Fraction(1,100) for v in exact_distances.values())
        strict=[verify[str(b)]["strict_residual_over_r_upper"] for b in (400,460)]
        selected[key]=dict(side=side,sample=sample,core=core,opposite=old["is_opposite_to_earlier_class"],
                           sealed_class=old["declared_earlier_class"],independent_audit=verify,
                           both_balls=bool(max(strict)<1.),file=use.name,reference_max_distance_interval=beta_reference,
                           outside_both_calibration_boxes=outside_both,
                           core_area_m2=old["area_m2"],core_cells=old["core_cells"],
                           effort_sources=[f"candidate_{key}.json"]+([f"expanded_{key}.json"] if ex else []))
    check=[]
    for side in (32,96):
        ix,roi,centers,reads,labels,r,energy,path=read(side)
        if side==32:
            witness=ENV/"verified_responses_32_2.0_phase_32_two_shells_400.npz"
        else:
            witness=ENV/"verified_responses_96_2.0_wave_96_2.0_1_400.npz"
        y_old=next(iter(reads.values()))[0]
        with np.load(witness,allow_pickle=False) as f:
            p=at_load(f["s1"],1.)[np.ix_(RX,TX)]
            m=at_load(f["s0"],1.)[np.ix_(RX,TX)]
        midpoint=(p+m)/(2*np.sqrt(8))
        diff=float(np.linalg.norm(midpoint-y_old))
        if diff>1e-13:
            raise ArithmeticError("Actual old midpoint has a different absolute calibration")
        check.append(dict(side=side,old_total_midpoint_error_fro=diff,old_ref_source=witness.name,
                          candidate_pair_y_read_from=path,old_total_power=energy[0],new_total_power=energy[1],
                          old_per_shot_amplitude=np.sqrt(energy[0]/8),
                          new_per_shot_amplitude=np.sqrt(energy[1]/8),radius=r))
    write_json("REPORT_RECEIPT.json",dict(arms=selected,contract_normalization=check,
               paid_extra_independent_audits=costs,**elapsed(start)))
    return selected,original,extensions


def read_consumed(d):
    base=d["cpu_s"]
    if "optimizer_cost" in d:
        lu=sum(x["material_lu"] for x in d["optimizer_cost"].values())
        rhs=sum(x["material_forward_rhs"]+x["material_adjoint_rhs"] for x in d["optimizer_cost"].values())
    else:
        lu=sum(x["material_lu"] for x in d["whole_grid_jacobian_cost"].values())
        rhs=sum(x["forward_right_hand_sides"]+x["transpose_right_hand_sides"] for x in d["whole_grid_jacobian_cost"].values())
        for i in d["steps"]:
            lu+=sum(x["material_lu"] for x in i["physics_cost_accepted"].values())
            rhs+=sum(x["material_forward_rhs"]+x["material_adjoint_rhs"] for x in i["physics_cost_accepted"].values())
    audit=d.get("independent_audit")
    if audit:
        lu+=len(audit)
        rhs+=sum(e["strict_inverse_rhs"]+e["strict_physical_rhs"] for e in audit.values())
    return base,lu,rhs


def report(data,original,extension):
    rates=[]
    for f in sorted(data):
        d=data[f]
        lo=d["independent_audit"]["400"]["strict_residual_over_r_upper"]
        hi=d["independent_audit"]["460"]["strict_residual_over_r_upper"]
        lo,hi=np.ceil(np.array([lo,hi])*1e7)/1e7
        prior="高" if d["sealed_class"]=="high" else "低"
        status="两球内" if d["both_balls"] else "该候选未通过，材料类未决"
        rates.append(f"| {d['side']}² | {d['sample']} | {prior} | {int(d['core'])} | {lo:.7f} | {hi:.7f} | {status} |")
    ledger={32:[0.,0,0],96:[0.,0,0]}
    for d in [*original.values(),*extension.values()]:
        side=d["side"]
        vals=read_consumed(d)
        ledger[side]=[a+b for a,b in zip(ledger[side],vals)]
    table="\n".join(rates)
    example=next(v for v in data.values() if v["both_balls"] and v["opposite"])
    coef=example["reference_max_distance_interval"]
    prev=example["independent_audit"]
    a=np.ceil(prev["400"]["strict_residual_over_r_upper"]*1e7)/1e7
    b=np.ceil(prev["460"]["strict_residual_over_r_upper"]*1e7)/1e7
    outgoing=[d for d in data.values() if d["opposite"] and d["both_balls"]]
    text=rf"""# 已采双频读数的宽材料盒再攻击

同一份 32² 低规格读数允许核心均值 **0 和 2** 的两份全波材料。此前闭环把 batch_leaf 判为低规格。用这份已经记录的旧 400 MHz 与新 460 MHz 读数，在原 \([0,4.5]\) 全材料盒内搜索共享双频材料，找到核心固定为 2 的反例。独立裸 VIE 与 25 端电路复核的残差上界分别为 \( {a:.7f}r \) 和 \( {b:.7f}r \)。因此宽材料盒下，可行核心均值跨度至少为 2。该反例超出此前按两份参考定义的 \(\pm0.01\) 校准邻域。

32² 固定核心为 16 格、\(0.0625\,\mathrm{{m^2}}\)，96² 固定核心为 164 格、\(0.0711805556\,\mathrm{{m^2}}\)。两频率的材料完全共享，外区允许落入 \([0,4.5]\)。搜索先从公开的两规格启动，对环带划分 4 与 8 扇区，再逐格放开公共参考支撑。未找到同球候选的分支再根据全图导数增加外部像素。搜索没有读取离线样品或目标标签。单独保留旧、新读数的绝对 Frobenius 球。源幅使用 \(\sqrt{{E/8}}\)，旧频率双方 \(E=1\)，新频率 32² 为 \(E=1/8\)，96² 为 \(E=11/16\)。两份旧读数均对记录的总波幅中点做精确归一化回查。

| 网格 | 记录样品 | 旧闭环规格判定 | 强制核心均值 | 400 MHz 严格残差上界 / \(r\) | 460 MHz 严格残差上界 / \(r\) | 结论 |
|---|---|---|---:|---:|---:|---|
{table}

表中未通过只表示该点没有通过外包复算。对应连续材料类未决。表内对向同时通过数为 {len(outgoing)}。32² 反例距公开高参考的最大像素差下界为 {np.floor(coef["2"][0]*1e5)/1e5:.5f}，远大于闭环的校准预算。该反例限定了闭环可解释的材料类，原来两份校准规格之间的一侧区分证书仍按其已声明的材料邻域解释。

搜索和独立点态核验的 32² / 96² 单线程 CPU 分别约为 {ledger[32][0]:.2f} / {ledger[96][0]:.2f} 秒。过程合计材料 LU 分别为 {ledger[32][1]} / {ledger[96][1]} 次，前向、转置及逆核验列总数分别为 {ledger[32][2]} / {ledger[96][2]}。所有实际计算费用在各分支收据中列出。复算见 attack.py 和 expand.py，汇总可通过 report.py 再生成。对于表内压力测试中仅为汇总新核查的严格点，额外费用在汇总收据。

[逐例严格残差与归一化回查](REPORT_RECEIPT.json) · [32² 决定性错规格材料](candidate_32_batch_leaf_2.npz) · [32² 同读数低核心材料](candidate_32_batch_leaf_0.npz) · [同色标材料配对图](paired_32.pdf)
"""
    (HERE/"REPORT_ZH.md").write_text(text)


def figure(data):
    import os
    os.environ["MPLBACKEND"]="Agg"
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(7.0,3.5),layout="constrained")
    for core,axis in zip((0,2),axes):
        key=f"32_batch_leaf_{core}"
        q=data[key]
        with np.load(HERE/q["file"],allow_pickle=False) as datafile:
            x=datafile["material"].reshape(32,32)
            mask=np.zeros(x.size)
            mask[datafile["roi_indices"]]=1.
        m=axis.imshow(x,extent=(-1,1,-1,1),origin="lower",vmin=0,vmax=4.5,cmap="viridis")
        a=(np.arange(32)+.5)/16-1
        axis.contour(a,a,mask.reshape(32,32),[.5],colors=["white"],linewidths=.9)
        r400=q["independent_audit"]["400"]["strict_residual_over_r_upper"]
        r460=q["independent_audit"]["460"]["strict_residual_over_r_upper"]
        axis.set_title(f"Core = {core} | 400: {r400:.3f}r, 460: {r460:.3f}r",fontsize=10)
        axis.set_xlabel("x (m)")
        axis.set_ylabel("y (m)")
    fig.colorbar(m,ax=axes,label="Real material contrast",shrink=.78)
    fig.savefig(HERE/"paired_32.pdf")
    fig.savefig(HERE/"paired_32.png",dpi=170)
    plt.close(fig)


if __name__=="__main__":
    data,original,extension=summarize()
    report(data,original,extension)
    figure(data)

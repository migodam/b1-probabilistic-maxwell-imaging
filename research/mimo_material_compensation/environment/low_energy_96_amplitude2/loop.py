"""Freeze, acquire, decide and only then score the same-amplitude 96² run."""
from pathlib import Path
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import json
from io import BytesIO
import os
import secrets
import sys
import zipfile

for name in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","MKL_NUM_THREADS"):
    os.environ[name]="1"
sys.dont_write_bytecode=True

ROOT=Path(__file__).resolve().parent
ENV=ROOT.parent
sys.path.insert(0,str(ENV))
from frequency_classes import STUDY, VerifiedWave, down, elapsed, fro, now, np, rd, up, write_json
from matched_amplitude import read_pair
from passive_replay import load_cover, stored_response
from probe import Circuit, at_load
from flint import arb


TARGET_RADIUS=0.01
CENTER_CORE=2.0
GRID=96
NOISE_FRACTION=.85
SAFETY_THRESHOLD=1.04
CAPTURES=("specimen_013","specimen_074")


def stamp():
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def read(fname):
    return json.loads((ROOT/fname).read_text())


def require_absent(*names):
    for name in names:
        if (ROOT/name).exists():
            raise FileExistsError(f"Sealed packet {name} already exists")


def upper(x):
    return float(np.nextafter(float(x)+rd(3)*abs(x),np.inf))


def real(a):
    return np.r_[a.real.ravel(),a.imag.ravel()]


def response_unit(a):
    return a[8:16,:8]


def safe_sqrt(energy):
    center=energy**.5
    error=up(float(abs(arb(energy).sqrt()-arb(center)).upper()))
    return center,error


def save_array(filename,**arrays):
    with (ROOT/filename).open("xb") as out:
        np.savez_compressed(out,**arrays)


def freeze():
    start=now()
    require_absent("proposal.json","selector.npz","readings.npz","decision.json")
    packet=json.loads((ENV/"matched_amplitude.json").read_text())
    opposite=json.loads((ENV/"kyp_verified_responses_96_2.0_wave_96_2.0_1_400.json").read_text())
    if not opposite["feasible"] or packet["nominal_core_mean_gap"]!=CENTER_CORE or packet["side"]!=GRID:
        raise ValueError("Different design packet")
    old=opposite["upper_over_2r"]
    r=packet["radius_sqrtW"]
    band=packet["outputs"]["460"]
    bound=next(b for b in band["budgets"] if b["beta"]==TARGET_RADIUS)
    g1=bound["open_lower_over_two_r"]
    with np.load(ENV/"matched_amplitude.npz",allow_pickle=False) as raw:
        direction=raw[bound["lambda_key"]]
        empty=raw["open_empty_offset_460"]
        A_old=raw["midpoint_400_total"]
        stock=[raw[f"stock_460_{tag}"] for tag in ("plus","minus")]
    if fro(direction)>1.:
        raise ArithmeticError("The frozen realified direction exceeds unit norm")
    endpoints=[]
    for n in (0,1):
        base=real(stock[n]+empty)
        subtract=bound["endpoint_stock_terms"][n]
        transfer=bound["endpoint_open_uniform_terms"][n]["response_error_fro_upper"]
        v=float(direction@base)
        arith=up(rd(len(base)+5)*fro(direction)*fro(base)+rd(4)*(fro(stock[n])+fro(empty)))
        err=up(band["anchor_audits"][n]["response_error_upper"]+
               band["open_coefficient_transfer"]["empty_total_Fro_error_upper"]+arith+transfer)
        delta=up(subtract["safe_support"]+subtract["derivative_term_error"]+subtract["nonlinear"]+err)
        endpoints.append(dict(center_score=v,support=subtract["safe_support"],
                              j_derivative_error=subtract["derivative_term_error"],
                              nonlinear=subtract["nonlinear"],physical_open_transfer=transfer,
                              other_round_and_field_errors=err,total_paid_deduction=delta))
    fw_hi=down(endpoints[0]["center_score"]-endpoints[0]["total_paid_deduction"])
    fw_lo=upper(endpoints[1]["center_score"]+endpoints[1]["total_paid_deduction"])
    if down((fw_hi-fw_lo)/(2*r)) <= 1.:
        raise ArithmeticError("The two endpoints do not yield a common separating direction")
    actions=[]
    for n in range(1,17):
        e=n/16
        amp,err_amp=safe_sqrt(e)
        lock=up(err_amp*(abs(fw_hi)+abs(fw_lo))+rd(6)*(amp*(abs(fw_hi)+abs(fw_lo))+r))
        high=down(amp*fw_hi-r-lock)
        low=upper(amp*fw_lo+r+lock)
        pair_lower=down((amp-err_amp)*g1)
        actions.append(dict(level=n,energy=e,amplitude=amp,amplitude_error_upper=err_amp,
                        proven_pair_ratio_over_two_r=pair_lower,high_noisy_score_lower=high,
                        low_noisy_score_upper=low,positive_gap=high>low,
                        meets_frozen_margin=pair_lower>SAFETY_THRESHOLD and high>low))
    sel=next(x for x in actions if x["meets_frozen_margin"])
    assert all(not e["meets_frozen_margin"] for e in actions[:sel["level"]-1])
    noisylow,noisyhigh=sel["low_noisy_score_upper"],sel["high_noisy_score_lower"]
    threshold=(noisylow+noisyhigh)*.5
    if not noisylow<threshold<noisyhigh:
        raise ArithmeticError("No rounded threshold inside the verified interval")
    exact_beta=Fraction(1,100)
    interval_low=[0.,float(np.nextafter(float(exact_beta),np.inf))]
    interval_high=[float(np.nextafter(float(Fraction(2)-exact_beta),-np.inf)),
                   float(np.nextafter(float(Fraction(2)+exact_beta),np.inf))]
    assert Fraction(interval_low[1])>=exact_beta
    assert Fraction(interval_high[0])<=2-exact_beta and Fraction(interval_high[1])>=2+exact_beta
    save_array("selector.npz",direction=direction,threshold=threshold,old_reading_total=A_old)
    record=dict(stage="freeze",sealed_utc=stamp(),side=GRID,core_cells=packet["core_cells"],core=2.,
                core_area_m2=packet["core_area_m2"],beta=TARGET_RADIUS,radius=r,
                observation_id=packet["radius_source_id"],frequency_MHz=460,old_frequency_MHz=400,
                receive_layout="A with 25th probe open; original eight hold ports terminated",
                convention="Total (including empty field) 8x8 realified response. Real frequency-independent contrast.",
                scope="Uniform discrimination only for the union of the two declared reference neighborhoods.",
                input_contract=["matched_amplitude.json","matched_amplitude.npz",
                                "kyp_verified_responses_96_2.0_wave_96_2.0_1_400.json"],
                upper400=old,lower460=g1,necessary_400=down(1./old**2),
                theoretical_sufficient_460=up(1./g1**2),energy_ratio_lower=down((g1/old)**2),
                operating_energy=sel["energy"],operating_level_num=sel["level"],
                actual_energy_ratio_lower=down(1./(old*old*sel["energy"])),
                engineering_margin=SAFETY_THRESHOLD,actions=actions,
                pure_gap_minimum_power=next(x["energy"] for x in actions if x["positive_gap"]),
                high_clean_score_lower=fw_hi,low_clean_score_upper=fw_lo,
                high_noisy_score_lower=noisyhigh,low_noisy_score_upper=noisylow,
                decision_threshold=threshold,endpoint_certificates=endpoints,
                interval_high=interval_high,interval_low=interval_low,display_intervals={
                "high":["1.99","2.01"],"low":["0","0.01"]},
                shots=8,complex_reads=64,noise_fraction_of_radius=NOISE_FRACTION,
                acquisition_design={"amplitude_cap":0.008,"outside_sites":96,
                    "family":"Independent seeded smooth material perturbations with sparse additional exterior changes. One specimen per class. Independent complex noise orientation at prescribed fraction of fixed absolute radius."},
                **elapsed(start))
    write_json(ROOT/"proposal.json",record)
    print("frozen",sel["level"],"/16",threshold,noisylow,noisyhigh,flush=True)


def extra_sites(ix):
    outside=np.setdiff1d(np.arange(GRID**2),ix)
    return outside[np.linspace(13,len(outside)-14,96,dtype=int)]


def synthesize(material,ix,sel,rng):
    x=material.copy()
    jj,ii=np.divmod(ix,GRID)
    films=[]
    for k in range(3):
        dy,dx=rng.uniform(.08,.40,size=2)
        phase=rng.uniform(0,2*np.pi)
        films.append(np.cos(dy*jj+dx*ii+phase))
    field=.0065*np.mean(films,axis=0)+.0015*rng.uniform(-1.,1.,size=len(ix))
    x[ix]=np.clip(x[ix]+field,0,4.5)
    x[sel]=np.clip(x[sel]+rng.uniform(0.001,0.0075,size=len(sel)),0,4.5)
    if up(float(np.max(np.abs(x-material))))>TARGET_RADIUS or np.array_equal(x,material):
        raise ArithmeticError("Specimen does not meet the frozen generator contract")
    return x


def power_stats(s,energy):
    ports=tuple(range(24))
    amp,_=safe_sqrt(energy/8)
    loaded=at_load(s,1.,ports=ports)[:,:8]*amp
    diagonals=loaded[np.arange(8),np.arange(8)]
    reflect=float(np.sum(abs(diagonals)**2))
    all_out=float(np.sum(abs(loaded)**2))
    passive=all_out-reflect
    return dict(total_available_W=energy,total_accepted_W=energy-reflect,driven_reflected_W=reflect,
                other_terminations_W=passive,received_A_W=float(np.sum(abs(loaded[8:16])**2)),
                radiated_or_dissipated_W=energy-all_out)


def acquisition():
    start=now()
    p=read("proposal.json")
    require_absent("readings.npz","acquisition.json","offline_specimens.npz","offline_keys.json")
    materials,ix,roi,r,obs=read_pair()
    sel=extra_sites(ix)
    ids=np.unique(np.r_[ix,sel])
    wave=Circuit(GRID,460,ids)
    rng=secrets.SystemRandom()
    assign=[0,1]
    rng.shuffle(assign)
    meta={}
    measured={}
    actual={}
    scores={}
    amp,_=safe_sqrt(p["operating_energy"])
    for tag,k in zip(CAPTURES,assign):
        specimen_seed=rng.getrandbits(128)
        noise_seed=rng.getrandbits(128)
        genome=synthesize(materials[k],ix,sel,np.random.default_rng(specimen_seed))
        s=wave.solve(genome[ids])
        clean=response_unit(at_load(s,1.))*amp/np.sqrt(8)
        noise_rng=np.random.default_rng(noise_seed)
        direction=noise_rng.normal(size=clean.shape)+1j*noise_rng.normal(size=clean.shape)
        recorded=clean+((NOISE_FRACTION*r)/np.linalg.norm(direction))*direction
        measured[tag]=recorded
        actual[tag+"_material"]=genome
        actual[tag+"_clean_numerical"]=clean
        actual[tag+"_scattering_460"]=s
        actual[tag+"_noise_numerical"]=recorded-clean
        scores[tag]=power_stats(s,p["operating_energy"])
        meta[tag]=dict(specimen_seed_hex=f"{specimen_seed:032x}",noise_seed_hex=f"{noise_seed:032x}",
                       reference_index=k,label="high" if k==0 else "low")
    save_array("offline_specimens.npz",active_union=ids,roi=roi,**actual)
    write_json(ROOT/"offline_keys.json",dict(stage="sealed generator",sealed_utc=stamp(),samples=meta))
    save_array("readings.npz",**measured)
    write_json(ROOT/"acquisition.json",dict(stage="acquire",sealed_utc=stamp(),frequency_MHz=460,open_probe=True,
                read_count_per_sample=64,shots_per_sample=8,readings=list(CAPTURES),
                selected_available_W=p["operating_energy"],noise_fraction_of_fixed_ball=NOISE_FRACTION,
                power_numerical=scores,counts=dict(wave.count),**elapsed(start)))
    print("acquired",len(measured),"blind observations",flush=True)


def decide():
    start=now()
    require_absent("decision.json","offline_scores.json")
    p=read("proposal.json")
    require=[ROOT/"selector.npz",ROOT/"readings.npz"]
    if not all(f.is_file() for f in require):
        raise FileNotFoundError("Selection and acquired data must already be sealed")
    threshold=p["decision_threshold"]
    with np.load(ROOT/"selector.npz",allow_pickle=False) as target:
        direction=target["direction"]
    if fro(direction)>1.:
        raise ArithmeticError("Unverified direction")
    answers={}
    with np.load(ROOT/"readings.npz",allow_pickle=False) as observations:
        for tag in sorted(observations.files):
            y=real(observations[tag])
            score=float(direction@y)
            error=up(rd(len(y))*fro(y)*fro(direction))
            lo,hi=down(score-error),upper(score+error)
            if lo>threshold:
                label="high"
            elif hi<threshold:
                label="low"
            else:
                label="unresolved"
            answers[tag]=dict(label=label,score=score,score_lower=lo,score_upper=hi,
                threshold=threshold,ROI_interval=p.get("interval_"+label,[0.,4.5]),
                certified_gap_from_threshold=max(0.,down(lo-threshold),down(threshold-hi)))
    write_json(ROOT/"decision.json",dict(stage="decide",sealed_utc=stamp(),
                         accessed=["proposal.json","selector.npz","readings.npz"],results=answers,
                         stopping="One verified scalar separates both known material neighborhoods.",
                         **elapsed(start)))
    print("predictions sealed")


def certified_total(wave,ver,contrast,energy):
    s,error,receipt=stored_response(wave,ver,contrast)
    p,open_err,_=load_cover(s,error,1.,0.,0.)
    act=response_unit(p)
    amp,amp_err=safe_sqrt(energy/8.)
    estimate=act*amp
    guard=up((amp+amp_err)*open_err+amp_err*fro(act)+rd(4)*amp*fro(act))
    return estimate,guard,s,receipt


def arb_mean(roi):
    x=sum(Fraction(y) for y in roi)/len(roi)
    a=float(x)
    return x,a


def score():
    start=now()
    p=read("proposal.json")
    acquired=read("acquisition.json")
    decision=read("decision.json")
    require_absent("offline_scores.json","offline_forward_audit.npz")
    if not (p["sealed_utc"]<acquired["sealed_utc"]<decision["sealed_utc"]):
        raise RuntimeError("Sealed order is inconsistent")
    released=read("offline_keys.json")
    references,ix,roi,r,_=read_pair()
    start_replay=now()
    with np.load(ROOT/"offline_specimens.npz",allow_pickle=False) as f:
        ids=f["active_union"]
        specimens={tag:f[tag+"_material"] for tag in CAPTURES}
        nominal={tag:f[tag+"_clean_numerical"] for tag in CAPTURES}
        noise={tag:f[tag+"_noise_numerical"] for tag in CAPTURES}
    additional=extra_sites(ix)
    if not np.array_equal(ids,np.unique(np.r_[ix,additional])):
        raise ArithmeticError("Generator support receipt does not agree")
    with np.load(ROOT/"selector.npz",allow_pickle=False) as f:
        old=f["old_reading_total"]
    with np.load(ROOT/"readings.npz",allow_pickle=False) as f:
        measurements={k:f[k] for k in f.files}
    rows={}
    arrays={}
    for tag in CAPTURES:
        truth=released["samples"][tag]
        ref=references[truth["reference_index"]]
        x=specimens[tag]
        exact_max_change=max(abs(Fraction(a)-Fraction(b)) for a,b in zip(x,ref))
        if exact_max_change>Fraction(1,100):
            raise ArithmeticError("Outward contract check failed")
        regenerated=synthesize(ref,ix,additional,np.random.default_rng(int(truth["specimen_seed_hex"],16)))
        if not np.array_equal(regenerated,x):
            raise ArithmeticError("Sealed material is not reproducible")
        draw=np.random.default_rng(int(truth["noise_seed_hex"],16))
        z=draw.normal(size=nominal[tag].shape)+1j*draw.normal(size=nominal[tag].shape)
        y=nominal[tag]+(NOISE_FRACTION*r/np.linalg.norm(z))*z
        if not np.array_equal(y,measurements[tag]):
            raise ArithmeticError("Sealed reading is not reproducible")
        roi_exact,roi_value=arb_mean(x[roi])
        predicted=decision["results"][tag]["ROI_interval"]
        covers=Fraction(predicted[0])<=roi_exact<=Fraction(predicted[1])
        rows[tag]=dict(label_revealed=truth["label"],label_prediction=decision["results"][tag]["label"],
                       correct=truth["label"]==decision["results"][tag]["label"],
                       region_mean=roi_value,mean_in_closed_preferred_interval=covers,
                       roi_interval=predicted,max_change_to_own_reference_upper=up(float(exact_max_change)),
                       new_outside_nonzero_cells=int(np.count_nonzero(x[additional])),
                       measured_noise_numerical_over_r=up(fro(noise[tag])/r))
        if not rows[tag]["correct"] or not covers:
            raise ArithmeticError("Blind decision or mean interval failed")
    costs={"sealed_specimen_and_noise_replay":elapsed(start_replay)}
    for band,energy,observations in ((400,1.,dict.fromkeys(CAPTURES,old)),
                                     (460,p["operating_energy"],measurements)):
        work=now()
        w=Circuit(GRID,band,ids)
        basis=dict(spectrum=w.symbol,grid_shape=(GRID,GRID),S=w.c,E0=w.b,
                   low_rank_left=np.zeros((GRID**2,1)),low_rank_right=np.zeros((1,GRID**2)))
        verifier=VerifiedWave(basis)
        for tag in CAPTURES:
            pred,err,s,proof=certified_total(w,verifier,specimens[tag][ids],energy)
            d=pred-observations[tag]
            exact_r=up((fro(d)+err+rd(3)*(fro(pred)+fro(observations[tag])))/r)
            rows[tag][f"strict_{band}_residual_over_r"]=exact_r
            rows[tag][f"strict_{band}_error_fro_upper"]=err
            rows[tag][f"strict_{band}_audit"]=proof
            if exact_r>=1:
                raise ArithmeticError(f"Reading inconsistent with sealed material at {band}")
            arrays[f"{tag}_response_{band}"]=pred
            arrays[f"{tag}_scattering_{band}"]=s
        costs[str(band)]=dict(**elapsed(work),strict_material_lu=2,
                              strict_inverse_columns=2*len(ids),physical_rhs=50,
                              explicit_open_circuit_inversions=2)
    save_array("offline_forward_audit.npz",old_common_reading=old,**arrays)
    commit=stamp()
    write_json(ROOT/"offline_scores.json",dict(stage="score",revealed_utc=commit,read_after_prediction=True,
               results=rows,score_cost_detail=costs,**elapsed(start)))
    write_report()
    print([(name,s["label_prediction"],s["region_mean"],s["strict_400_residual_over_r"],
            s["strict_460_residual_over_r"]) for name,s in rows.items()],flush=True)


def write_report():
    p=read("proposal.json")
    a=read("acquisition.json")
    prediction=read("decision.json")
    s=read("offline_scores.json")
    data=s["results"]
    table="\n".join(
        f"| {name} | {'高' if d['label_prediction']=='high' else '低'} | {d['region_mean']:.9f} | "
        f"{np.ceil(d['strict_400_residual_over_r']*1e7)/1e7:.7f} | "
        f"{np.ceil(d['strict_460_residual_over_r']*1e7)/1e7:.7f} | "
        f"{np.ceil(d['max_change_to_own_reference_upper']*1e7)/1e7:.7f} |"
        for name,d in data.items())
    freeze_pct=p["operating_level_num"]
    old_display=np.ceil(p["upper400"]*1e8)/1e8
    new_display=np.floor(p["lower460"]*1e8)/1e8
    ratio_display=np.floor(p["actual_energy_ratio_lower"]*1e7)/1e7
    energy_required=np.floor(p['necessary_400']*1e7)/1e7
    sufficient=np.ceil(p['theoretical_sufficient_460']*1e7)/1e7
    text=rf"""# 96² 同核心增幅的低功率闭环

两件新仿真样品各用一次 460 MHz 开路 A 读数，在解盲前分别判为高、低两种规格，解盲后均正确。两件样品也都与同一个 400 MHz 旧读数球相容。

核心增幅 2，区域为 164 个原生像素、{p['core_area_m2']:.10f} 平方米。两类允许在各自参考的全图所有像素独立改变最多 0.01，且材料保持在 [0, 4.5]。采集前冻结两个均值区间的向外取整版本，低规格为 [0, 0.01]，高规格为 [1.99, 2.01]。两新材料的生成器使用独立随机数，并在两份样品中均加入 96 个原参考为零的外区像素。新测量使用相同的绝对硬 Frobenius 半径，噪声取半径的 85%。未读取旧观测的目标或测量标签。

两类在 460 MHz 的严格全类间距至少为 {new_display:.8f} 倍 \(2r\)。400 MHz 的 KYP 上界为 {old_display:.8f} 倍 \(2r\)，适用于该两个参考在本合同内的相干发射和被动负载。460 MHz 冻结阈值为 {p['decision_threshold']:.12g}。沿用预先声明的 1.04 归一化分离余量，最低可选档为 **{freeze_pct}/16**。只要求判别区间严格分开的数学最低档为 {p["pure_gap_minimum_power"]:.4f}。归一化的理论 460 MHz 充分量上界为 {sufficient:.7f}。400 MHz 保证区分所需能量至少 {energy_required:.7f}，相当于实际采用 460 MHz 档位的至少 {ratio_display:.7f} 倍。该比值来自全类证书与 400 MHz 基准反例，两次样品分类用于独立检查这次动作。

| 样品 | 解盲前判定 | 解盲后的真实区域均值 | 400 MHz 旧球严格残差 / \(r\) | 460 MHz 新球严格残差 / \(r\) | 相对自身规格的最大像素改动上界 |
|---|---|---:|---:|---:|---:|
{table}

每件样品发射 8 次、获得 64 个复数读数。总可用功率为 {p["operating_energy"]:.4f} W，按两频率等曝光换算归一化能量。总计 16 次发射和 128 个复数读数。发射持续时间未记录，能量以同曝光基准为单位。采集模拟使用 {a["counts"]["lu"]} 次裸材料 LU、{a["counts"]["rhs"]} 个右端，决策为每样品一次 128 维实投影。400 MHz 复算计入离线评分，逐频核验两次 LU。冻结、采集、决策和评分的实际单线程 CPU 用时分别为 {p["cpu_s"]:.4f}、{a["cpu_s"]:.4f}、{prediction["cpu_s"]:.4f}、{s["cpu_s"]:.4f} 秒。

[冻结动作与阈值](proposal.json) · [解盲前结果](decision.json) · [离线评分与成本](offline_scores.json) · [400 MHz KYP 证书](../kyp_verified_responses_96_2.0_wave_96_2.0_1_400.json) · [460 MHz 全类证书](../matched_amplitude.json)

本目录保留已封存样品与种子。要新建同合同轨迹，在同一环境下的空目录依次执行 loop.py freeze、acquisition、decide、score、archive。评分阶段从封存种子重生材料和测量，逐项核对原读数。
"""
    (ROOT/"REPORT_ZH.md").write_text(text)


def archive():
    start=now()
    if not (ROOT/"offline_scores.json").is_file():
        raise FileNotFoundError("The archive requires the sealed score")
    own=("loop.py","proposal.json","selector.npz","readings.npz","acquisition.json","decision.json",
         "offline_keys.json","offline_specimens.npz","offline_forward_audit.npz","offline_scores.json","REPORT_ZH.md")
    src=("matched_amplitude.json","matched_amplitude.npz","wave_96_2.0_1.npz",
         "kyp_verified_responses_96_2.0_wave_96_2.0_1_400.json",
         "verified_responses_96_2.0_wave_96_2.0_1_400.npz",
         "verified_responses_96_2.0_wave_96_2.0_1_460.npz","matched_amplitude.py",
         "frequency_classes.py","frequency_classes_open.py","frequency_classes_independent.py",
         "probe.py","passive_replay.py","bounded_real.py")
    workspace=ENV.parents[2]
    external=[STUDY/"diagnostics/verified_wave.py",STUDY/"diagnostics/aperture_control.py",
              STUDY/"diagnostics/aperture_control_freeze.json",STUDY/"imaging/shared.py",
              STUDY/"imaging/FREEZE.json",STUDY/"ports/model.py",
              STUDY/"ports/public/operators/n96_f400_a08.npz",
              STUDY/"ports/public/operators/n96_f460_a08_c048.npz"]
    external.extend(sorted((STUDY/"vendor_r8").rglob("*.py")))
    external.extend(ENV.parent/"budgets"/f for f in ("bound.py","budget_proof.py","quick.py","freeze.json"))
    path=ROOT/"matched_amplitude_energy_archive.zip"
    pending=path.with_suffix(".pending")
    entries={str(f.relative_to(workspace)):f for f in [*(ROOT/n for n in own),*(ENV/n for n in src),*external]}
    old=STUDY/"ports/public/cases/n96_c01.npz"
    header=BytesIO()
    with np.load(old,allow_pickle=False) as public:
        np.savez_compressed(header,observation_id=public["observation_id"],radius=public["radius"])
    metadata=str(old.relative_to(workspace))
    readme="""# 96² 同核心低功率闭环复算包

解压后保留 research/ 下相对位置。Python 3.13 依赖 NumPy、SciPy、python-flint、threadpoolctl。
子目录 low_energy_96_amplitude2 的 REPORT_ZH.md 链接报告、冻结动作、盲判和评分。
原 n96_c01 公开观测的档只保留了代码实际读取的 observation_id 与 radius，不包括任何目标数据。
要新建采集轨迹，把 loop.py 复制到 environment 下一个空子目录，再按报告顺序执行。
"""
    with zipfile.ZipFile(pending,"x",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as out:
        for name,f in entries.items():
            out.write(f,name)
        out.writestr(metadata,header.getvalue())
        out.writestr("README.md",readme)
    with zipfile.ZipFile(pending) as out:
        assert set(out.namelist())==set(entries)|{metadata,"README.md"}
        if out.testzip() is not None:
            raise ArithmeticError("Archive failed CRC")
    pending.replace(path)
    write_json(ROOT/"archive_receipt.json",dict(sealed_utc=stamp(),archive=path.name,
                entries=len(entries)+2,byte_count=path.stat().st_size,**elapsed(start)))
    print(path,flush=True)



if __name__=="__main__":
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("stage",choices=("freeze","acquisition","decide","score","archive"))
    args=arg.parse_args()
    globals()[args.stage]()

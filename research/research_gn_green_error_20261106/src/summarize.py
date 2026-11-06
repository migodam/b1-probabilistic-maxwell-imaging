"""Readable requested energy tables from saved diagnostic coefficients; no physics."""
from pathlib import Path
import csv,json,math,time
from collections import defaultdict
import numpy as np
HERE=Path(__file__).resolve().parents[1]
R=HERE/"results"

def loadcsv(name):
    with (R/name).open() as f:return list(csv.DictReader(f))

def writecsv(name,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with (R/name).open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)

def fnum(x):
    x=float(x)
    return f"{x:.3g}" if 0<abs(x)<.001 else f"{x:.3f}"

def main():
    t=time.perf_counter()
    summary=json.loads((R/"SUMMARY.json").read_text())
    cases=loadcsv("case_metrics.csv")
    byid={r["task_id"]:r for r in cases}
    rows=loadcsv("single_operator.csv")
    conditional=defaultdict(lambda:dict(component=0.,total=0.,n=0))
    for r in rows:
        if float(r["tau"])!=.1:continue
        key=tuple(r[k] for k in ("method","amplitude","noise","quantity","operator","band","joint_gate"))
        d=conditional[key];d["component"]+=float(r["component_energy"]);d["total"]+=float(r["total_energy"]);d["n"]+=1
    cr=[]
    for key,d in conditional.items():
        r=dict(zip(("method","amplitude","noise","quantity","operator","band","joint_gate"),key))
        r.update(cases=d["n"],energy=d["component"],total_energy=d["total"],pooled_energy_pct=100*d["component"]/d["total"])
        cr.append(r)
    writecsv("conditional_success_energy.csv",cr)
    rawstats=[]
    for r in rows:
        if float(r["tau"])!=.1:continue
        frac=float(r["fraction"])
        rawstats.append({k:r[k] for k in ("task_id","method","amplitude","noise","quantity","operator","band","joint_gate")}|
                        dict(energy_pct=100*frac,norm_pct=100*math.sqrt(frac),
                             component_RMS=float(r["component_RMS"]),component_energy=float(r["component_energy"]),
                             total_energy=float(r["total_energy"]),
                             output_energy_pct=100*float(r["output_fraction"])))
    writecsv("requested_per_case_statistics.csv",rawstats)
    spec=np.load(R/"fixed_spectra.npz")
    modes=np.load(R/"cumulative_mode_energies.npz")
    threshold=defaultdict(lambda:dict(strong=0.,total=0.,count=0))
    # This is a response curve of stored coefficients, not another fit or a tuned threshold.
    cuts=(.001,.01,.03,.1,.3)
    for tid,q,op,e in zip(modes["task_id"],modes["quantity"],modes["operator"],modes["mode_energy"]):
        c=byid[str(tid)]
        s=spec[("real_" if q=="material_real" else "")+("in_s" if op=="G_in" else "ob_s")]
        rel=s/s[0]
        for tau in cuts:
            key=(c["method"],c["amplitude"],c["noise"],str(q),str(op),tau)
            d=threshold[key];d["strong"]+=float(e[rel>=tau].sum());d["total"]+=float(e.sum());d["count"]+=1
    tr=[]
    for k,d in threshold.items():
        r=dict(zip(("method","amplitude","noise","quantity","operator","tau"),k))
        r.update(cases=d["count"],strong_energy_pct=100*d["strong"]/d["total"])
        tr.append(r)
    writecsv("threshold_response_from_saved_coefficients.csv",tr)
    cohort=[]
    for method in ("native_trf","reciprocal_trf"):
        for noise in ("clean","noisy"):
            for amplitude in (.4,2.5,12.):
                rs=[r for r in cases if r["method"]==method and r["noise"]==noise and float(r["amplitude"])==amplitude]
                me=sum(float(r["material_error_energy"]) for r in rs)
                je=sum(float(r["current_error_energy"]) for r in rs)
                cohort.append(dict(method=method,noise=noise,amplitude=amplitude,cases=len(rs),
                    joint_PASS=sum(r["joint_gate"]=="PASS" for r in rs),
                    pooled_material_RMSE=math.sqrt(me/(len(rs)*576)),
                    pooled_current_RMS=math.sqrt(je/(len(rs)*576*12)),
                    material_error_energy=me,current_error_energy=je,
                    internal_field_error_energy=sum(float(r["internal_field_error_energy"]) for r in rs),
                    true_data_error_energy=sum(float(r["true_data_error_energy"]) for r in rs),
                    median_relative_true_data_error=float(np.median([float(r["relative_true_data_error"])for r in rs]))))
    writecsv("cohort_absolute_energy.csv",cohort)
    lines=["# 请求的误差比例表","",
           "主阈值 sigma/sigma_max≥0.1 为强，1e-12<sigma/sigma_max<0.1 为弱，余下为右域补空间。全部百分数是平方范数能量；每组按12个对象的误差能量合并。",""]
    for quantity in ("material","current"):
        lines += [f"## Native clean：{quantity}，单算子","",
                  "| A | G | 强 % | 弱 % | 补 % |","|---|---|---:|---:|---:|"]
        for a in (.4,2.5,12.):
            for op in ("G_in","G_ob"):
                rs=[r for r in summary["primary_single"] if r["method"]=="native_trf" and r["noise"]=="clean" and r["amplitude"]==a and r["quantity"]==quantity and r["operator"]==op]
                es={r["band"]:100*r["pooled_fraction"] for r in rs}
                lines.append(f"| {a:g} | {op} | {fnum(es['strong'])} | {fnum(es['weak'])} | {fnum(es['complement'])} |")
        lines += ["",f"## Native clean：{quantity}，先内部后接收串联投影","",
                  "| A | 强强 % | 强弱 % | 弱强 % | 弱弱 % | 含补空间 % |",
                  "|---|---:|---:|---:|---:|---:|"]
        for a in (.4,2.5,12.):
            rs=[r for r in summary["primary_twofold"] if r["method"]=="native_trf" and r["noise"]=="clean" and r["amplitude"]==a and r["quantity"]==quantity and r["definition"]=="in_then_ob"]
            es={(r["in_band"],r["ob_band"]):100*r["pooled_fraction"]for r in rs}
            comp=sum(v for (i,j),v in es.items() if "complement" in (i,j))
            values=[es[i,j] for i,j in (("strong","strong"),("strong","weak"),("weak","strong"),("weak","weak"))]+[comp]
            lines.append("| "+f"{a:g} | "+" | ".join(fnum(v)for v in values)+" |")
        lines+=[""]
    lines+=["上表强弱组合是P_ob,b P_in,a的串联能量，严格交集和正交化TSOM投影空间另存aggregate_twofold.csv。它们不相同，不能互相替换。","",
            "三个严格交集/投影定义、顺序反转、噪声、Reciprocal、实材料限制控制、绝对能量及PASS/FAIL条件统计都在对应CSV中。"]
    (R/"REQUESTED_TABLES_ZH.md").write_text("\n".join(lines)+"\n")
    (R/"POSTPROCESS_RECEIPT.json").write_text(json.dumps(dict(wall_seconds=time.perf_counter()-t,new_physics_LU=0,
        new_GN_optimizations=0,posthoc_threshold_curve="read-only stored modal coefficients; primary0.1andprotocolsensitivitiesunchanged"),indent=2))
    print(json.dumps(dict(status="REQUESTED_TABLES_WRITTEN",seconds=time.perf_counter()-t)))
if __name__=="__main__":main()


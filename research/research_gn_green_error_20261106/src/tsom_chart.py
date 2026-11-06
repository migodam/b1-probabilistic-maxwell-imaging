"""Energy in the literal legacy [V_ob,strong, (I-P_ob,strong)V_in,strong] chart."""
from pathlib import Path
import csv,json,time
import numpy as np
from scipy.linalg import svd
from threadpoolctl import threadpool_limits
P=Path(__file__).resolve().parents[1];R=P/"results"
def main():
    z=np.load(R/"fixed_spectra.npz");e=np.load(R/"error_arrays.npz")
    with (R/"case_metrics.csv").open() as f:meta={r["task_id"]:r for r in csv.DictReader(f)}
    vo=z["ob_V"][:,z["ob_s"]>=.1*z["ob_s"][0]]
    vi=z["in_V"][:,z["in_s"]>=.1*z["in_s"][0]]
    start=time.perf_counter()
    zz=vi-vo@(vo.conj().T@vi)
    q,d,_=svd(zz,full_matrices=False,check_finite=False)
    q=q[:,d>1e-10]
    rows=[]
    for tid,dm,dj in zip(e["task_id"],e["delta_m"],e["delta_j"]):
        m=meta[str(tid)]
        for quantity,x in (("material",dm[:,None]),("current",dj)):
            total=float(np.vdot(x,x).real)
            ev=float(np.linalg.norm(vo.conj().T@x)**2)
            ef=float(np.linalg.norm(q.conj().T@x)**2)
            rows.append(dict(task_id=str(tid),method=m["method"],amplitude=m["amplitude"],noise=m["noise"],
                             quantity=quantity,total_energy=total,receiver_strong_energy=ev,
                             second_fold_energy=ef,remainder_energy=total-ev-ef))
    with (R/"actual_tsom_chart.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    ag=[]
    for method in ("native_trf","reciprocal_trf"):
        for noise in ("clean","noisy"):
            for a in (.4,2.5,12.):
                for quantity in ("material","current"):
                    rr=[r for r in rows if r["method"]==method and r["noise"]==noise and float(r["amplitude"])==a and r["quantity"]==quantity]
                    total=sum(r["total_energy"]for r in rr)
                    ag.append(dict(method=method,noise=noise,amplitude=a,quantity=quantity,
                                   receiver_strong_pct=100*sum(r["receiver_strong_energy"]for r in rr)/total,
                                   second_fold_pct=100*sum(r["second_fold_energy"]for r in rr)/total,
                                   remainder_pct=100*sum(r["remainder_energy"]for r in rr)/total))
    receipt=dict(classification="DIRECT_LEGACY_TSOM_CHART_FROM_SAVED_ARRAYS",tau=.1,
                 receiver_rank=vo.shape[1],second_fold_rank=q.shape[1],
                 remainder_rank=vo.shape[0]-vo.shape[1]-q.shape[1],
                 wall_seconds=time.perf_counter()-start,new_LU=0,new_RHS=0,
                 scope="Orthogonal chart groups; second fold jointly projects into receiver weak+complement with shared coefficients, not union of separately parameterized spans",
                 aggregate=ag)
    (R/"actual_tsom_chart_summary.json").write_text(json.dumps(receipt,indent=2))
    print(json.dumps({k:v for k,v in receipt.items()if k!="aggregate"}))
    for r in ag:
        if r["method"]=="native_trf" and r["noise"]=="clean":print(r)
if __name__=="__main__":
    with threadpool_limits(limits=1):main()


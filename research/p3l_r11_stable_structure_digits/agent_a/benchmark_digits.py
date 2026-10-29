"""Explicitly released digit occupancy kernel inputs; no online/truth reads."""
import os
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
import csv,json,time
import numpy as np
from threadpoolctl import threadpool_limits
from benchmark import HERE,Ledger,Full256,load_host,SIGMA,write,summarize
from solver import Solver


def run():
    contract=dict(id="R11_S3_A_DIGIT_KERNEL_V1",frozen_utc=time.time(),
        input="data/kernel_digit_shapes.npz occupancy,dataset explicitly released by main; not available to inversion",
        shape_count=18,delta=[[.8,.03],[1.8,.18]],workloads=["forward_only","complete_occupancy_and_property_J"],
        ordering="alternate support/Full256 by shape, node, workload",threads=1,
        planned_RHS_all_orders=3456,lifetime_cap=20000,truth_files_access=False,
        scope="known supplied occupancy kernel benchmark, not inference or blind image quality")
    if (HERE/"CONTRACT_DIGIT_KERNEL.json").exists():raise RuntimeError("contract exists")
    write(HERE/"CONTRACT_DIGIT_KERNEL.json",contract)
    book=Ledger("digit_kernel");rows=[]
    try:
        s=load_host();support=Solver(s["C"],s["O"],s["ER"],s["Fh"],SIGMA);full=Full256(s["C"],s["O"],s["ER"],s["Fh"])
        with np.load(HERE.parent/"data/kernel_digit_shapes.npz",allow_pickle=False) as z:
            occupancy=z["occupancy"].copy();datasets=z["dataset"].copy()
        for ishape,(a,ds) in enumerate(zip(occupancy,datasets)):
            ns=int(np.count_nonzero(a>0));identity=f"kernel_{ishape:02d}_{ds}"
            for node,(dr,di) in enumerate(contract["delta"]):
                delta=complex(dr,di)
                for need_jac in (False,True):
                    methods=("support","Full256") if (ishape+node+int(need_jac))%2 else ("Full256","support");result={}
                    for method in methods:
                        order=ns if method=="support" else 256;nrhs=40 if need_jac else 8
                        book.reserve(order,nrhs,1,f"{identity}:node{node}:jac{need_jac}:{method}")
                        result[method]=(support if method=="support" else full).evaluate(a,delta,need_jac)
                    f=result["Full256"];su=result["support"]
                    row=dict(shape=identity,dataset=str(ds),node=node,need_jac=need_jac,support_size=ns,
                        delta_real=dr,delta_imag=di,Full256_seconds=f["diagnostics"]["wall_seconds"],support_seconds=su["diagnostics"]["wall_seconds"],
                        white_F_error=float(np.sqrt(2)*np.linalg.norm(f["F"]-su["F"])/SIGMA),
                        relative_J_occ_error=float(np.linalg.norm(f["J_occ"]-su["J_occ"])/np.linalg.norm(f["J_occ"])) if need_jac else 0.,
                        relative_J_delta_error=float(np.linalg.norm(f["J_delta"]-su["J_delta"])/np.linalg.norm(f["J_delta"])) if need_jac else 0.,
                        zero_columns_nonzero=int(np.count_nonzero(np.linalg.norm(su["J_occ"][:,a==0],axis=0)>0)) if need_jac else 0,
                        zero_columns_total=int(np.count_nonzero(a==0)))
                    rows.append(row);book.event("PAIRED_NODE_COMPLETE",**row)
        with (HERE/"DIGIT_KERNEL_NODES.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        summary=dict(status="COMPLETE",shapes=18,support_size_min=min(r["support_size"] for r in rows),support_size_max=max(r["support_size"] for r in rows),
            all_datasets=summarize(rows),by_dataset={str(ds):summarize([r for r in rows if r["dataset"]==ds]) for ds in ("N","P")},
            qualification="released designed digit occupancy, fixed two legal properties; no inference/truth scoring; one controlled kernel window")
        write(HERE/"DIGIT_KERNEL_SUMMARY.json",summary);summary["receipt"]=book.close("COMPLETE");print(json.dumps(summary,indent=2))
    except BaseException as e:book.close("FAILED",repr(e));raise


if __name__=="__main__":
    with threadpool_limits(limits=1):run()

"""Audited R10 physical kernel benchmark; strict <=20,000 RHS total.

Only named occupancy/delta arrays and fixed-host ROI cache are read.  No ground
truth files or full576 solver are opened.  Set all physics BLAS threads to one.
"""
import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

from pathlib import Path
import argparse
import csv
import json
import time
import resource
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from threadpoolctl import threadpool_limits, threadpool_info
from solver import Solver

HERE = Path(__file__).resolve().parent
B1 = HERE.parents[3]
HIST = B1 / "b1-p3l-r10-bilateral-structure/research/p3l_r10_bilateral_structure/followup_shape_reuse"
HOST = B1 / "b1-p3l-r9-stabilized-subspace/research/p3l_r9_stabilized_subspace/results/stable_roi/SETUP.npz"
SIGMA = 0.00119211798514519
CAP = 20000


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")


class Ledger:
    def __init__(self, run):
        self.f = (HERE / "COST_LEDGER.jsonl").open("a")
        self.rhs = 0
        self.lu = {}
        self.by_order = {}
        self.sequence = 0
        self.run = run
        self.t0 = time.perf_counter()
        # The budget is lifetime across every run of this benchmark.
        for line in (HERE / "COST_LEDGER.jsonl").read_text().splitlines():
            event = json.loads(line)
            if event["event"] == "RESERVE_SOLVE":
                self.rhs += event["RHS"]
                order = str(event["order"])
                self.by_order[order] = self.by_order.get(order, 0) + event["RHS"]
                self.lu[order] = self.lu.get(order, 0) + event.get("LU", 0)
            self.sequence = max(self.sequence, event["sequence"])
        self.event("RUN_START", cap=CAP, lifetime_RHS_before=self.rhs)

    def event(self, name, **fields):
        self.sequence += 1
        self.f.write(json.dumps(dict(sequence=self.sequence, utc=time.time(), run=self.run, event=name, **fields)) + "\n")
        self.f.flush()

    def reserve(self, order, nrhs, lu, label):
        if self.rhs + nrhs > CAP:
            raise RuntimeError("lifetime kernel RHS cap exhausted")
        self.rhs += nrhs
        key = str(order)
        self.by_order[key] = self.by_order.get(key, 0) + nrhs
        self.lu[key] = self.lu.get(key, 0) + lu
        self.event("RESERVE_SOLVE", order=int(order), RHS=int(nrhs), LU=int(lu), label=label, lifetime_RHS=self.rhs)

    def close(self, status, error=None):
        receipt = dict(status=status, error=error, RHS_all_orders=self.rhs, RHS_by_order=self.by_order,
                       LU_by_order=self.lu, new_RHS576=0, GT_access=False,
                       wall_seconds=time.perf_counter()-self.t0,
                       peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                       BLAS_threads=threadpool_info())
        self.event("RUN_END", **receipt)
        self.f.close()
        write(HERE / f"RECEIPT_{self.run}.json", receipt)
        return receipt


class Full256:
    """Independent full-size current solve with implicit output sensitivities."""
    def __init__(self, C, O, E, Fh):
        self.C, self.O, self.E, self.Fh = C, O, E, Fh

    def evaluate(self, a, delta, need_jac=False):
        t = time.perf_counter()
        n = len(a)
        A = np.eye(n, dtype=complex) - delta * a[:, None] * self.C
        lu = lu_factor(A, check_finite=False)
        h = lu_solve(lu, delta * a[:, None] * self.E, check_finite=False)
        F = self.Fh + self.O @ h
        out = dict(F=F)
        if need_jac:
            e = self.E + self.C @ h
            Z = lu_solve(lu, self.O.T, trans=1, check_finite=False).T
            out["J_occ"] = np.einsum("rj,js->srj", delta*Z, e, optimize=True).reshape(self.E.shape[1]*self.O.shape[0], n)
            out["J_delta"] = (Z @ (a[:, None]*e)).T.reshape(-1)
        out["diagnostics"] = dict(wall_seconds=time.perf_counter()-t)
        return out


def load_host():
    with np.load(HOST, allow_pickle=False) as z:
        return {k: z[k].copy() for k in ("C", "O", "ER", "Fh", "ROI")}


def inputs_r10(roi):
    for p in sorted(HIST.glob("R10BI_FRESH_*_ANCHOR.npz")):
        with np.load(p, allow_pickle=False) as z:
            a, ds = z["fraction"][roi].copy(), z["delta"].copy()
        yield p.stem, a, ds
        for i in range(3):
            pp = p.with_name(p.name.replace("_ANCHOR", f"_t{i}"))
            with np.load(pp, allow_pickle=False) as z:
                a, ds = z["fraction"][roi].copy(), z["delta"].copy()
            yield pp.stem, a, ds


def summarize(rows):
    by = {}
    for j in (False, True):
        rr = [r for r in rows if r["need_jac"] == j]
        full = sum(r["Full256_seconds"] for r in rr)
        support = sum(r["support_seconds"] for r in rr)
        by["full_jacobian" if j else "forward_only"] = dict(
            nodes=len(rr), Full256_seconds=full, support_seconds=support,
            support_saving=1-support/full,
            median_paired_speed_ratio=float(np.median([r["Full256_seconds"]/r["support_seconds"] for r in rr])),
            faster_nodes=sum(r["support_seconds"] < r["Full256_seconds"] for r in rr),
            max_white_F_error=max(r["white_F_error"] for r in rr),
            max_relative_J_occ_error=max(r["relative_J_occ_error"] for r in rr),
            max_relative_J_delta_error=max(r["relative_J_delta_error"] for r in rr),
        )
    return by


def run_r10(optimized=False):
    tag = "R10_OPTIMIZED" if optimized else "R10"
    contract = dict(id="R11_S3_A_EXACT_SUPPORT_R10_OPT_V2" if optimized else "R11_S3_A_EXACT_SUPPORT_R10_V1", frozen_utc=time.time(),
        inputs="8 R10 anchors plus 24 saved path occupancies; fraction/delta keys only; fixed ROI cache",
        support="strict a>0, no threshold; include every zero-cell occupancy derivative",
        delta_indices=[0] if optimized else [0,2,4,6], workloads=["forward_only", "complete_occupancy_and_property_J"],
        execution="alternate paired ordering by node/workload; no hidden warmup; single BLAS thread",
        cap_RHS_all_orders_lifetime=CAP, planned_R10_RHS=3072 if optimized else 12288, new_RHS576=0,
        GT_access=False, scientific_scope="kernel workload; no end-to-end inversion or material claims")
    if (HERE / f"CONTRACT_{tag}.json").exists():
        raise RuntimeError("R10 kernel contract already exists; preserve append-only history")
    write(HERE / f"CONTRACT_{tag}.json", contract)
    book = Ledger("r10_optimized" if optimized else "r10")
    rows = []
    try:
        st = load_host()
        support = Solver(st["C"], st["O"], st["ER"], st["Fh"], SIGMA)
        full = Full256(st["C"], st["O"], st["ER"], st["Fh"])
        for ishape, (identity, a, ds) in enumerate(inputs_r10(st["ROI"])):
            ns = int(np.count_nonzero(a > 0))
            for inode, node in enumerate(contract["delta_indices"]):
                delta = complex(ds[node])
                for need_jac in (False, True):
                    order = ("support", "Full256") if (ishape+inode+int(need_jac))%2 else ("Full256", "support")
                    result = {}
                    for method in order:
                        n = ns if method == "support" else 256
                        nrhs = (40 if need_jac else 8) if n else 0
                        book.reserve(n, nrhs, int(n>0), f"{identity}:node{node}:jac{int(need_jac)}:{method}")
                        solver = support if method == "support" else full
                        result[method] = solver.evaluate(a, delta, need_jac=need_jac)
                    f, s = result["Full256"], result["support"]
                    row = dict(shape=identity, node=node, need_jac=need_jac, support_size=ns,
                        delta_real=delta.real, delta_imag=delta.imag,
                        Full256_seconds=f["diagnostics"]["wall_seconds"], support_seconds=s["diagnostics"]["wall_seconds"],
                        white_F_error=float(np.sqrt(2)*np.linalg.norm(f["F"]-s["F"])/SIGMA),
                        relative_J_occ_error=float(np.linalg.norm(f["J_occ"]-s["J_occ"])/np.linalg.norm(f["J_occ"])) if need_jac else 0.,
                        relative_J_delta_error=float(np.linalg.norm(f["J_delta"]-s["J_delta"])/np.linalg.norm(f["J_delta"])) if need_jac else 0.,
                        zero_columns_nonzero=int(np.count_nonzero(np.linalg.norm(s["J_occ"][:,a==0],axis=0)>0)) if need_jac else 0,
                        zero_columns_total=int(np.count_nonzero(a==0)))
                    rows.append(row)
                    book.event("PAIRED_NODE_COMPLETE", **row)
        with (HERE / f"{tag}_KERNEL_NODES.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        summary=dict(status="COMPLETE", shapes=32, nodes=len(rows),
            support_size_min=min(r["support_size"] for r in rows),support_size_max=max(r["support_size"] for r in rows),
            costs=summarize(rows), qualification="single controlled physical-kernel workload; no population speedup, no inversion claim")
        write(HERE/f"{tag}_KERNEL_SUMMARY.json",summary)
        summary["receipt"] = book.close("COMPLETE")
        print(json.dumps(summary,indent=2))
    except BaseException as e:
        book.close("FAILED",repr(e))
        raise


if __name__ == "__main__":
    parser=argparse.ArgumentParser();parser.add_argument("action",choices=["r10","r10_optimized"])
    args=parser.parse_args()
    with threadpool_limits(limits=1):
        if args.action=="r10":run_r10()
        elif args.action=="r10_optimized":run_r10(optimized=True)

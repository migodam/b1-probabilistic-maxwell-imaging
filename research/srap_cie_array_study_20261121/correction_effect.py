"""Match every light output to its own paid seed in already-scored 2D trajectories."""
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"


def read(name):
    with (RESULTS/name).open() as f:
        return list(csv.DictReader(f))


def main():
    out=[]
    for cohort in ("confirmation","followup","development","legacy","legacy_mac"):
        if not (RESULTS/f"{cohort}_trajectories.csv").exists():
            continue
        trajectories=defaultdict(list)
        for r in read(f"{cohort}_trajectories.csv"):
            trajectories[(r["case"],r["method"])].append(r)
        for r in read(f"{cohort}_metrics.csv"):
            if r["view"]!="light" or not r.get("rel_L2"):
                continue
            budget=int(r["budget"])
            candidates=[t for t in trajectories[(r["case"],r["method"])] if float(t["seconds"])<=.8*budget]
            if not candidates:
                continue
            seed=min(candidates,key=lambda q:float(q["phi"]))
            row={k:r[k] for k in ("case","array","method","budget","host")}
            row.update(cohort=cohort,seed_step=seed["step"],seed_seconds=seed["seconds"],
                       audited_correction_steps=r["iterations_corrector"])
            for key in ("phi","rel_L2","half_IoU","residual"):
                row["seed_"+key]=seed[key]
                row["final_"+key]=r[key]
                row["change_"+key]=float(r[key])-float(seed[key])
            out.append(row)
    with (RESULTS/"corrector_seed_changes.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)


if __name__=="__main__":
    main()

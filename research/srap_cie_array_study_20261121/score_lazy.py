from score_study import score_case,write_csv
import run_study as run
from concurrent.futures import ProcessPoolExecutor
arms=("dual","srap","srap_lazy")
directory=run.ROOT/"results"/"lazy_fresnel"
def score(c):return score_case(c,arms,directory)[0]
if __name__=="__main__":
    cs=[c for c in run.cases() if c["suite"]=="LA_FRESNEL"]
    missing=[(c["case"],m) for c in cs for m in arms if not (directory/f"{c['case']}__{m}.npz").is_file()]
    if missing:raise SystemExit(str(missing))
    with ProcessPoolExecutor(max_workers=2) as pool: rows=[x for group in pool.map(score,cs) for x in group]
    write_csv(run.ROOT/"results"/"fresnel_lazy_metrics.csv",rows)

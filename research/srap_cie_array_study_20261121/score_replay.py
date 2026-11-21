"""Score CPU crossover only when all twelve public synthetic cases are present."""
from score_study import score_case,write_csv
import run_study as run
from concurrent.futures import ProcessPoolExecutor
methods=("dual","srap","u11")
folder=run.ROOT/"results"/"legacy_mac"
def score(c):return score_case(c,methods,folder)
if __name__=="__main__":
    cases=[c for c in run.legacy_cases() if c["suite"]=="LA_SYNTH"]
    if not all((folder/f"{c['case']}__{m}.npz").is_file() for c in cases for m in methods):raise SystemExit("Incomplete crossover cohort")
    with ProcessPoolExecutor(max_workers=2) as ex: output=list(ex.map(score,cases))
    write_csv(run.ROOT/"results"/"legacy_mac_metrics.csv",[r for s,_ in output for r in s])
    write_csv(run.ROOT/"results"/"legacy_mac_trajectories.csv",[r for _,s in output for r in s])

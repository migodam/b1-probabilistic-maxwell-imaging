"""Re-anchor the broad fixed patch at the sealed 140-trial checkpoint."""
import json
from time import perf_counter, process_time

from continuation import HERE, Model, run_path, np, threadpool_limits, write_json, write_arrays
from replay import archived_reduction

CASE = "n32_c01"
METHOD = "ball_min_margin"


def main():
    start = perf_counter(), process_time()
    out = HERE/f"{CASE}_paths_broad.json"
    document = json.loads(out.read_text())
    if document.get("checkpoint_restart"):
        raise ValueError("The one registered endpoint restart has already been run")
    original = document["paths"][METHOD]
    with np.load(out.with_suffix(".npz"), allow_pickle=False) as a:
        arrays = dict(a)
    m = Model.read(CASE)
    r = archived_reduction(CASE,"_broad",m)
    patch = arrays["patch"]
    is_patch = np.isin(r.indices,patch)
    previous = arrays[METHOD+"_materials"]
    seed = previous[-1]
    lo = float(seed[patch][0])
    if not np.all(seed[patch]==lo):
        raise ValueError("Sealed patch endpoint cells differ")
    phase, materials, knots = run_path(r,is_patch,arrays["roi_mask"],METHOD,2.,2.**-10,
        document["cell_area_m2"],seed=seed,start_target=lo,start_step=.03125,attempt_budget=24)
    if not np.array_equal(materials[0],seed):
        raise ValueError("Endpoint restart has a material discontinuity")
    merge = dict(**phase)
    merge["attempts"] = original["attempts"]+phase["attempts"]
    for k in ("wall_s","cpu_s","strict_segments"):
        merge[k] = original[k]+phase[k]
    merge["accepted_ball_upper_over_r"] = max(q for q in
        (original["accepted_ball_upper_over_r"],phase["accepted_ball_upper_over_r"]) if q is not None)
    for k in ("certified_patch_mean_interval","certified_registered_region_mean_interval"):
        merge[k][0]=original[k][0]
    replay=m.forward(materials[-1],fresh=True)
    merge["terminal_independent_original_replay"]=dict(
        residual_over_r=replay["residual_ratio"],
        reduced_vs_original_over_r=float(np.linalg.norm(replay["pred"]-knots[-1].prediction)/m.radius),
        state_defect=replay["state_relative_defect"])
    arrays[METHOD+"_materials"]=np.concatenate((previous,materials[1:]))
    if len(knots)>1:
        pairs = [("_predictions",np.stack([k.prediction for k in knots[1:]])),
                 ("_support_currents",np.stack([k.current for k in knots[1:]]))]
        for label, raw in pairs:
            arrays[METHOD+label]=np.concatenate((arrays[METHOD+label],raw))
    arrays[METHOD+"_terminal_original_current"]=replay["j"]
    arrays[METHOD+"_terminal_original_prediction"]=replay["pred"]
    document["paths"][METHOD]=merge
    document["checkpoint_restart"]=dict(start_at_uniform_patch=lo,
        already_paid_attempts=len(original["attempts"]),first_phase_stop=original["stop"],
        charged_attempt_cap=24,attempts=len(phase["attempts"]),strict_new_segments=phase["strict_segments"],
        new_native_counters=dict(m.cost),wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1])
    write_arrays(out.with_suffix(".npz"),**arrays)
    write_json(out,document)


if __name__=="__main__":
    with threadpool_limits(1):
        main()

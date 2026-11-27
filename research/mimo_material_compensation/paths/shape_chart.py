"""Strictly constant data with a tangent chosen by a regional-mean task."""
from pathlib import Path
from time import perf_counter, process_time

from continuation import HERE, Model, exact_chart, np, threadpool_limits, write_arrays, write_json
from replay import archived_reduction


def main():
    for case, suffix in (("n32_c01","_whole_region"),("n96_c01","_ball")):
        start = perf_counter(), process_time()
        model = Model.read(case)
        r = archived_reduction(case,suffix,model)
        source = HERE/f"{case}_paths{suffix}.npz"
        with np.load(source,allow_pickle=False) as a:
            target = a["roi_mask"].nonzero()[0] if suffix=="_whole_region" else a["patch"]
        is_task = np.isin(r.indices,target)
        rec,raw = exact_chart(r,is_task,free_shape=True)
        rec.update(case=case,source=source.name,target_cells=len(target),target_area_m2=(2/model.side)**2*len(target),
                   target_definition="whole registered ROI mean" if suffix=="_whole_region" else "mean on the declared connected subregion",
                   cost_wall_s=perf_counter()-start[0],cost_cpu_s=process_time()-start[1])
        base=HERE/f"{case}_exact_shape"
        write_json(base.with_suffix(".json"),rec)
        write_arrays(base.with_suffix(".npz"),**raw,target=target,support=r.indices,m0=r.x0)
        print(case,rec["strict"],rec.get("patch_mean_interval"),rec["direction_on_task_min"],
              rec["direction_on_task_max"],flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()

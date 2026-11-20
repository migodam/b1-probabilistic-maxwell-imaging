"""Compare the same proposed current/voxel steps with an independent full Maxwell solve."""

import json
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import lu_solve
from threadpoolctl import threadpool_limits

from cohorts import load_cohort
from current_cycle import (
    OUT,
    compatible_blend,
    norm,
    physical,
    prepare,
    regularizer,
    schur_current,
)


def trial(model, case, m, bank):
    j, y, fac = physical(model, m)
    field = model["E0"] + model["G"] @ j
    s = model["S"]
    sensitivity = lu_solve(fac, s.T, trans=1, check_finite=False).T
    observed = case["Y"]
    den = norm(observed) ** 2
    rows = []

    def compare(t, mp, jp, ep, defect):
        actual = physical(model, mp)[1]
        born = y + sensitivity @ ((mp - m)[:, None] * field)
        transported = y + sensitivity @ ((mp - m)[:, None] * ep)
        mix = compatible_blend(j, field, jp, ep, mp)
        blended = born + mix * (transported - born)
        forward_step = max(norm(actual - y), 1.0e-20)
        prior = regularizer(mp, model["shape"])[0]

        def value(p):
            return 0.5 * norm(p - observed) ** 2 / den + prior

        rows.append(
            dict(
                t=t,
                true_value=value(actual),
                born_value=value(born),
                transport_value=value(transported),
                blend_value=value(blended),
                blend_rel=norm(actual - blended) / forward_step,
                mix=mix,
                born_rel=norm(actual - born) / forward_step,
                transport_rel=norm(actual - transported) / forward_step,
                step_norm=forward_step,
                defect=norm(defect),
            )
        )

    threshold = case["sigma"] * np.sqrt(len(bank["sigma"]) * observed.shape[1])
    schur_current(
        model, j, m, observed, bank, "shared_cubic", threshold, fac, observe=compare
    )
    choices = {}
    for key in ("true", "born", "transport", "blend"):
        index = min(range(len(rows)), key=lambda i: rows[i][key + "_value"])
        choices[key] = index
    return dict(
        candidates=rows,
        selection=choices,
        start_value=0.5 * norm(y - observed) ** 2 / den
        + regularizer(m, model["shape"])[0],
    )


def main():
    rows = []
    t0, c0 = perf_counter(), process_time()
    for cohort in ("r10_twod", "r12_confirm24"):
        model, cases = load_cohort(cohort)
        bank = prepare(model)
        with np.load(OUT / f"current_cycle_equal_{cohort}.npz") as a:
            estimates = {
                cid: im[-1]
                for (_, cid, arm), im in zip(a["key"], a["images"])
                if arm == "voxel_only"
            }
        for case in cases:
            for stage, m in (
                ("background", np.full(len(model["G"]), 0.1)),
                ("voxel_fit", estimates[case["id"]]),
            ):
                rows.append(
                    dict(
                        cohort=cohort,
                        case=case["id"],
                        stage=stage,
                        **trial(model, case, m, bank),
                    )
                )
    (OUT / "receiver_predictions.json").write_text(
        json.dumps(
            dict(rows=rows, wall_s=perf_counter() - t0, cpu_s=process_time() - c0),
            indent=1,
        )
        + "\n"
    )


if __name__ == "__main__":
    with threadpool_limits(1):
        main()

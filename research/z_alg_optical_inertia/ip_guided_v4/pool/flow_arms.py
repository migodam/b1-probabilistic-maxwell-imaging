"""Compose measured A13 stage costs along identical candidate/refinement paths."""

from collections import Counter
import json

from flow_audit import OUTPUT
from pool_common import write_json


def read(name):
    return json.loads((OUTPUT/f"{name}.json").read_text())


def union_cost(*items):
    return {key: sum(v[key] for v in items) for key in ("wall_seconds", "cpu_seconds")}


def merge_counts(*items):
    merged = Counter()
    for item in items:
        merged.update(item)
    return merged


def report():
    first = read("online")
    second = read("continuation")
    scores = read("continuation_offline")
    checks = read("strict_point_witnesses")
    native_l4 = read("native_l4")
    proposals = read("screen_proposals")
    neumann = read("native_neumann")
    ideal = read("native_neumann_ideal")
    cases = {}
    for case, orig in first["cases"].items():
        extended = second["current_pass"][case]
        valid = {int(c) for c, ci in checks["cases"][case]["boxes"].items() if ci["strictly_feasible"]}
        excluded = {int(c) for c, ci in native_l4["cases"][case]["boxes"].items() if ci["status"]=="STRICTLY_EXCLUDED"}
        if valid & excluded:
            raise ArithmeticError("A certified feasible endpoint fell in a pruned box.")
        branches = orig["pipelines"]["all_cluster"]["visit"]

        proof_costs = [union_cost(*[item["costs"][stage] for stage in
                                   ("numeric_center", "exact_center_transport", "quadratic_completion")])
                       for item in native_l4["cases"][case]["boxes"].values()
                       if item["status"]=="STRICTLY_EXCLUDED"]
        shifted_search = union_cost(*proposals["cases"][case]["boxes"].values(), *proof_costs)
        neumann_time = union_cost(*[v["total"] for v in neumann["cases"][case]["branches"].values()])
        ineffectual = ideal["cases"][case]["rows"]
        if len(ineffectual)!=len(branches) or not all(row["classical_inapplicable"] or
               max(row["ideal_standard_screen_upper"], row["ideal_projected_max_upper"]) <= orig["radius"]
               for row in ineffectual.values()):
            raise ArithmeticError("Neumann control retained a possible decisive row.")

        first_fit_position = next(i for i,c in enumerate(branches) if c in valid)
        generated_and_ranked = union_cost(orig["setup"], orig["point_cost"])
        plan = dict(all_cluster=(branches, union_cost()), top_1=(branches[:1], union_cost()),
                    top_3=(branches[:3], union_cost()), top_6=(branches[:6], union_cost()),
                    top_early=(branches[:first_fit_position+1], union_cost()),
                    standard_neumann=(branches, neumann_time),
                    shifted=([c for c in branches if c not in excluded], shifted_search))
        methods = {}
        for method, (visited, screen_cost) in plan.items():
            first_pass = [orig["branch_cost"][c] for c in visited]
            second_pass = [extended["refined_box_results"][str(c)] for c in visited]
            selected = min(visited, key=lambda c: extended["refined_box_results"][str(c)]["residual"])
            inherited = method=="top_early"
            if inherited:
                selected = visited[-1]
            point_receipt = checks["cases"][case]["boxes"].get(str(selected))
            extra = union_cost() if point_receipt is None else point_receipt["total_full_wave_and_audit"]
            paid = merge_counts(dict(lu=len(orig["point_residuals"]),
                                     forward_rhs=2*len(orig["point_residuals"])),
                                *(d["counts"] for d in first_pass+second_pass))
            if extra["wall_seconds"]:
                paid.update(dict(lu=1,forward_rhs=2))
            methods[method] = dict(visited=visited, winner=selected,
                selected_physical_residual=extended["refined_box_results"][str(selected)]["residual"],
                selected_strict_feasible=(selected in valid),
                selected_strict_residual_upper=point_receipt.get("strict_residual_upper") if point_receipt else None,
                measured_stage_time=union_cost(generated_and_ranked,screen_cost,*first_pass,*second_pass,extra),
                screening_time=screen_cost, optimizer_and_point_counts=dict(paid),
                screening_lu=len(branches) if method=="standard_neumann" else len(excluded) if method=="shifted" else 0)
        best_box = scores[case]["best_box"]
        cases[case] = dict(radius=orig["radius"], original_true_material_contained_boxes=0,
            strict_material_feasible_boxes=sorted(valid), strict_excluded_boxes=sorted(excluded),
            unresolved_boxes=sorted(set(branches)-valid-excluded),
            first_strict_fit_point_rank=first_fit_position+1, best_full_replay_box=best_box,
            best_box_scores=scores[case]["scores"]["continued_boxes"], methods=methods)
    write_json(OUTPUT/"native_method_paths.json", dict(cases=cases, notes=(
        "Actual sampled, solved, and verified phases charged by selection path, "
        "independent aggregate repeat runs timed only in separate 24² study. "
        "Each selected start refines with the same 9 + 30 objective evaluations budget. "
        "All methods use 64 paid point evaluations for selecting the initial start per Ward cluster. "
        "Shared checkpoint/geometry setup is separate. Shifted screening phase costs above exclude optional Neumann ideal diagnostics.")))
    for key, row in cases.items():
        print(key,"cert-fit",len(row["strict_material_feasible_boxes"]),
              "prune",len(row["strict_excluded_boxes"]),"top until",row["first_strict_fit_point_rank"],
              [(k, round(v["measured_stage_time"]["wall_seconds"],3)) for k,v in row["methods"].items()])


if __name__ == "__main__":
    report()

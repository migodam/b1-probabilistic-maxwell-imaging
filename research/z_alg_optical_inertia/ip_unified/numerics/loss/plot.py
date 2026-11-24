"""Render outward endpoint evidence."""
import json

import matplotlib.pyplot as plt

from run import HERE


def load(name):
    return json.loads((HERE/f"{name}.json").read_text())


def main():
    finite = load("finite_0_verified")
    scalar = load("finite_0_pure_shift_verified")
    direct = load("finite_0_directional")
    state = load("finite_0_state_verified")
    native = load("native_verified")
    control = load("exact_control")
    first = next(r for block in finite["projected"] if block["cap"]==1. and block["tangent"]==.2
                 for r in block["rows"] if r["id"]==scalar["id"])
    r = first["rho"]
    directional = next(row for group in direct["profiles"] if group["cap"]==1. and group["tangent"]==.2
                       for row in group["rows"] if row["id"]==first["id"])
    shifted = next(row for row in state["rows"]
                   if row["id"]==first["id"] and row["tangent"]==.2 and row["cap"]==1.25)
    fig, ax = plt.subplots(1,3,figsize=(12.5,3.8),layout="constrained",gridspec_kw={"width_ratios":[1.15,1.2,1]})
    lower, upper = "#156751", "#875925"
    entries = [
        ("Receiver ellipsoid",first["dual_lower"]/r,lower),
        ("Right-singular directional prices",directional["methods"]["right_singular"]["verified"]["distance_lower"]/r,lower),
        ("Any loss-only scalar shift",scalar["all_tau_distance_upper"]/r,upper),
        ("Scalar current cap",first["scalar"]["scalar_primal_upper"]/r,upper),
        ("Sector + complete shift, $u=1.25$",shifted["state"]["shift"]["distance_lower"]/r,lower),
    ]
    for i,(_,distance,color) in enumerate(entries):
        ax[0].barh(i,distance,color=color,height=.58)
        ax[0].text(distance+.02,i,f"{distance:.3f}",va="center",fontsize=9)
    ax[0].set_yticks(range(len(entries)),[e[0] for e in entries],fontsize=9)
    ax[0].axvline(1,color="black",ls="--",lw=1)
    ax[0].invert_yaxis()
    ax[0].set(xlim=(0,1.48),xlabel="Distance / radius",title="Analytic 12 × 12 TM, $t=0.2$")

    for index,tangent in enumerate((.01,.05,.2)):
        entry = next((v for v in native["projected"] if v["cap"]==.75 and v["tangent"]==tangent),None)
        if entry is None:
            continue
        for row in entry["rows"]:
            low_peak = row["id"].split("__")[0].endswith("_low")
            ratio = row["primal"]["primal_upper"]/row["rho"]
            ax[1].plot(ratio,index,marker="x" if low_peak else "o",markersize=5,color=upper)
    ax[1].axvline(1,color="black",ls="--",lw=1)
    ax[1].set(yticks=range(3),yticklabels=("0.01","0.05","0.20"),xlim=(0,1.1),
              ylabel="Loss tangent $t$",xlabel="Upper distance / radius",title="Native 96 × 96 TM, $u=0.75$")

    names = ["Old rank tail","All-negative scalar cap","Mixed with one free","Complete shift"]
    values = [control["unshifted_distance"],control["full_negative_common_matrix_cap_distance"],
              control["mixed_distance_lower"],control["full_shift_distance"]]
    for i,distance in enumerate(values):
        ax[2].barh(i,distance,color=lower,height=.58)
        ax[2].text(distance+.02,i,f"{distance:.3f}",va="center",fontsize=9)
    ax[2].set_yticks(range(4),names,fontsize=9)
    ax[2].invert_yaxis()
    ax[2].axvline(.5,color="black",ls="--",lw=1)
    ax[2].set(xlim=(0,1.04),xlabel="Distance",title="Exact reciprocal control")
    fig.savefig(HERE/"loss_summary.png",dpi=170)


if __name__ == "__main__":
    main()

"""Physical ablations and full-state checks for both certified point-fiber roots."""

from time import perf_counter, process_time
import json

from flint import acb, acb_mat, arb, ctx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Patch
import numpy as np

from certify_fold import (
    HALFWIDTH, PRECISION, complex_norm, diagonal, float_interval, full_matrix, interval,
    null_from_minor, pack_rows, pulse_operators, real_column, reduce_host, response,
    sigma_min_enclosure, unit,
)
from host_model import RESULTS, full_forward


def to_numpy(A):
    return np.array([
        [complex(float(A[i, j].real.mid()), float(A[i, j].imag.mid()))
         for j in range(A.ncols())] for i in range(A.nrows())
    ])


def singular_values(J):
    approximate = to_numpy(J)
    _, values, Vh = np.linalg.svd(approximate)
    # The floating SVD only proposes the threshold and upper-bound witness.
    # The enclosure is established from the interval matrix itself.
    threshold = format(0.999 * values[-1], ".8g")
    _, upper, pivots = sigma_min_enclosure(acb_mat(J), threshold, Vh[-1].conjugate())
    return {
        "singular_values_at_midpoint": values.tolist(),
        "certified_sigma_min_strictly_exceeds": threshold,
        "certified_sigma_min_upper": float_interval(upper)["hi"],
        "minimum_shifted_gram_ldl_pivot_lo": min(float_interval(p)["lo"] for p in pivots),
    }


def receiver_nodes(T):
    intervals = [float_interval(T[0, j].abs_lower().union(T[0, j].abs_upper())) for j in range(T.ncols())]
    return {"each_modulus_interval": intervals, "minimum_modulus_lower": min(x["lo"] for x in intervals)}


def curvature(K, T, E, left, right):
    Dv = diagonal(right)
    r = len(right)

    def contribution(matrix):
        h = real_column(pack_rows((2 * T * Dv * matrix * Dv * E).transpose()))
        return sum((a * b for a, b in zip(left, h)), arb(0))

    local_return = diagonal([K[i, i] for i in range(r)])
    between_slots = acb_mat([[K[i, j] if i != j else 0 for j in range(r)] for i in range(r)])
    onsite = contribution(local_return)
    interslot = contribution(between_slots)
    total = contribution(K)
    return {
        "on_site_return": float_interval(onsite),
        "between_distinct_slots": float_interval(interslot),
        "sum_of_two": float_interval(onsite + interslot),
        "full": float_interval(total),
        "interslot_absolute_fraction_of_total": float_interval(abs(interslot / total)),
    }


def verify_root_boxes(numeric, fibers, H, B, R, y0, G):
    verified = []
    for side, row in fibers.items():
        for index, root in enumerate(row["physical_roots"]):
            z = [arb(float(lo)).union(arb(float(hi))) for lo, hi in root["material_intervals"]]
            A = full_matrix(z, G)
            _, values, Vh = np.linalg.svd(to_numpy(A))
            _, upper, pivots = sigma_min_enclosure(A, "0.25", Vh[-1].conjugate())
            y, _, metrics = full_forward(
                np.asarray(root["material_midpoint"]), *[numeric[key] for key in ("G", "E0", "S", "m0")]
            )
            _, _, _, _, T = response(z, H, B, R, y0)
            verified.append({
                "side": side, "root_index_0based": index,
                "material_intervals": root["material_intervals"],
                "material_midpoint": root["material_midpoint"],
                "sigma_min_at_midpoint": float(values[-1]),
                "certified_sigma_min_strictly_exceeds": "0.25",
                "certified_sigma_min_upper": float_interval(upper)["hi"],
                "minimum_shifted_gram_ldl_pivot_lo": min(float_interval(p)["lo"] for p in pivots),
                "full_forward_midpoint_residual_l2": float(np.linalg.norm(y - numeric["y_" + side])),
                "full_forward_midpoint_state_residual": metrics["state_relative_residual"],
                "receiver_nodes": receiver_nodes(T),
            })
    return verified


def geometry_figure(numeric):
    radius = float(numeric["radius"])
    centers, sources, receivers = (numeric[key] for key in ("centers", "sources", "receivers"))
    m0, nactive = numeric["m0"], len(numeric["active"])
    ink, muted = "#172939", "#587082"
    active_color, host_weak, host_strong = "#F5A26F", "#BBD9DD", "#548D98"
    source_color, receiver_color = "#714B8A", "#B25446"

    with plt.rc_context({
        "font.family": "DejaVu Sans", "font.size": 10, "text.color": ink,
        "axes.edgecolor": "#8B9FAA", "axes.labelcolor": ink, "axes.linewidth": 0.9,
        "xtick.color": muted, "ytick.color": muted, "svg.fonttype": "none", "pdf.fonttype": 42,
    }):
        fig, ax = plt.subplots(figsize=(6.8, 7.4))
        fig.subplots_adjust(left=.12, right=.97, top=.91, bottom=.235)
        ax.set_aspect("equal")
        ax.set(xlim=(-.78, .73), ylim=(-.78, .96),
               xlabel=r"$x/\lambda_0$", ylabel=r"$y/\lambda_0$")
        ax.set_xticks(np.arange(-.6, .61, .2))
        ax.set_yticks(np.arange(-.6, .81, .2))
        ax.grid(color="#D9E3E8", linewidth=.65, alpha=.8)
        ax.set_axisbelow(True)
        ax.set_title("Known host with six unknown disks", loc="left",
                     fontsize=15, fontweight="bold", pad=16)

        for i, (x, y) in enumerate(centers):
            face = active_color if i < nactive else host_weak if m0[i] < 2 else host_strong
            ax.add_patch(Circle((x, y), radius, facecolor=face, edgecolor=ink, lw=.8, zorder=3))
            if i < nactive:
                ax.text(x, y, str(i + 1), ha="center", va="center", fontsize=10,
                        fontweight="bold", zorder=4)

        ax.scatter(*sources.T, marker="^", s=56, color=source_color, edgecolor=ink, lw=.55, zorder=5)
        ax.scatter(*receivers.T, marker="D", s=48, color=receiver_color, edgecolor=ink, lw=.55, zorder=5)
        for i, (x, y) in enumerate(sources):
            offset = (-.036, .043) if i == 1 else (.035, .037)
            ax.text(x + offset[0], y + offset[1], f"Tx {i + 1}", ha="right" if i == 1 else "left",
                    color=source_color, weight="bold", fontsize=9)
        for x, y in receivers:
            ax.text(x + .052, y + .018, "Rx", color=receiver_color, weight="bold", fontsize=9)

        items = [
            Patch(facecolor=active_color, edgecolor=ink, label="Six unknown sites"),
            Patch(facecolor=host_weak, edgecolor=ink, label=r"Known host, $\varepsilon_r=2.4$"),
            Patch(facecolor=host_strong, edgecolor=ink, label=r"Known host, $\varepsilon_r=3.6$"),
            Line2D([0], [0], marker="^", color=source_color, markeredgecolor=ink,
                   linestyle="None", label="Line source (Tx)"),
            Line2D([0], [0], marker="D", color=receiver_color, markeredgecolor=ink,
                   linestyle="None", label="Receiver (Rx)"),
        ]
        fig.legend(handles=items, loc="lower center", bbox_to_anchor=(.52, .07), ncol=2,
                   frameon=False, columnspacing=1.8, labelspacing=.8, fontsize=9)
        fig.text(.52, .035, r"Disk radius $a=0.04\lambda_0$", ha="center",
                 color=muted, fontsize=9)
        figures = RESULTS.parent / "figures"
        figures.mkdir(exist_ok=True)
        fig.savefig(figures / "host_geometry.pdf", metadata={"Title": "Known-host six-slot geometry"})
        fig.savefig(figures / "host_geometry.png", dpi=210)
        plt.close(fig)


def main():
    ctx.prec = PRECISION
    start_wall, start_cpu = perf_counter(), process_time()
    numeric = np.load(RESULTS / "host_operators.npz")
    fibers = json.loads((RESULTS / "point_fibers.json").read_text())["fibers"]
    G, E0, S = pulse_operators()
    H, B, R, y0, *_ = reduce_host(G, E0, S)
    z = [arb(float(t)) for t in numeric["zc"]]
    axis, h = int(numeric["axis"]), arb(HALFWIDTH)
    z[axis] = (z[axis] - h).union(z[axis] + h)
    r = len(z)
    born_J = pack_rows(acb_mat([
        [R[0, j] * B[j, s] for j in range(r)] for s in range(B.ncols())
    ]))
    born_stats = singular_values(born_J)
    D = diagonal(z)
    diagonal_H = diagonal([H[i, i] for i in range(r)])
    no_interslot_y, no_interslot_J, *_ = response(z, diagonal_H, B, R, y0)
    y, J, C, E, T = response(z, H, B, R, y0)
    left, right, *_ = null_from_minor(J, int(np.argmax(abs(numeric["normal"]))),
                                     int(np.argmax(abs(numeric["null"]))))

    radius, pi = arb(str(float(numeric["radius"]))), arb.pi()
    ka = 2 * pi * radius
    c = pi * ka * acb(ka).bessel_j(1).real / 2
    if not c > 0:
        raise ArithmeticError("Outgoing Green prefactor sign not verified.")
    result = {
        "scope": "One 36-disk known-host layout, Arb enclosures from exact geometry. Ablations retain the same dressed B, R and y0.",
        "precision_bits": ctx.prec,
        "green_imaginary_gram": {
            "equation": "Im_H G = c J0(k d), c = pi k a J1(k a)/2, including the self disk d=0.",
            "strictly_positive_prefactor_exact_interval": interval(c),
            "strictly_positive_prefactor_display": float_interval(c),
        },
        "dressed_born": {
            "equation": "y = y0 + R diag(z) B. The material Jacobian is constant.",
            "global_nonsingularity": True,
            **born_stats,
            "jacobian_determinant": float_interval(born_J.det()),
            "fold_segment_full_vs_born_response_l2_difference": float_interval(
                complex_norm(y - y0 - R * D * B)
            ),
        },
        "pointwise_suppress_between_slot_feedback": {
            "equation": "H is replaced by diag(H), with all loaded susceptibilities and dressed B, R, y0 unchanged.",
            "at_certified_fold_material_segment": singular_values(no_interslot_J),
            "fold_segment_full_vs_ablation_response_l2_difference": float_interval(complex_norm(y - no_interslot_y)),
        },
        "full_at_fold": {
            "determinant_contains_zero": J.det().contains(0),
            "midpoint_response_j_singular_values": np.linalg.svd(to_numpy(J), compute_uv=False).tolist(),
            "receiver_nodes": receiver_nodes(T),
            "curvature_decomposition": curvature(H * C, T, E, unit(left), unit(right)),
        },
        "host_reference_receiver_nodes": receiver_nodes(R),
        "full_state_root_box_checks": verify_root_boxes(numeric, fibers, H, B, R, y0, G),
    }
    result["numeric_and_certificate_wall_seconds"] = perf_counter() - start_wall
    result["numeric_and_certificate_cpu_seconds"] = process_time() - start_cpu
    (RESULTS / "host_mechanism.json").write_text(json.dumps(result, indent=2) + "\n")
    geometry_figure(numeric)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

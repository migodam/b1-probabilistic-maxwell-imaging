"""Display two fixed recorded cases using actual coefficients and saved P3 DOFs.

C4_1 materials and strong_center fields are different predeclared cases.  This
script performs display sampling, P3 indexing and an analytic incident-field
evaluation only: no FEM object, operator assembly, factorization or solve.
Run through tools/run_job.py; the wrapper measures this plotting CPU once.
"""
from pathlib import Path
import argparse
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import tri as mtri
from matplotlib.colors import BoundaryNorm, ListedColormap, Normalize
from scipy.special import hankel1


ROOT = Path(__file__).resolve().parents[1]
MATERIAL_BOX = (-2., 2., .25, 3.)
FREQUENCY = 1.
POLYGON_N = 384
SOURCE0 = np.array([-1.5, 0.])
MAGNITUDE_CMAP = "viridis"
DIFFERENCE_CMAP = "RdBu_r"
PHASE_CMAP = "twilight"


def rel(path):
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def not_run(name, inputs, reason):
    return dict(name=name, status="NOT_RUN", inputs=[rel(path) for path in inputs], reason=reason)


def save_plot(fig, output, stem):
    files = []
    for extension in ("png", "pdf"):
        path = output / f"{stem}.{extension}"
        fig.savefig(path, bbox_inches="tight", **({"dpi": 180} if extension == "png" else {}))
        files.append(rel(path))
    plt.close(fig)
    return files


def limits(values):
    finite = np.asarray(values)[np.isfinite(values)]
    return dict(minimum=float(finite.min()), maximum=float(finite.max())) if finite.size else None


def format_axes(ax):
    ax.set(xlim=MATERIAL_BOX[:2], ylim=MATERIAL_BOX[2:], xlabel=r"$x/\lambda_0$", ylabel=r"$y/\lambda_0$")
    ax.set_aspect("equal")


def materials(output):
    source = ROOT / "results/g0_stress/C4_1.npz"
    stem = "C4_1_ACTUAL_MATERIALS"
    if not source.exists():
        return not_run(stem, [source], "Saved C4_1 endpoint pose is absent; no material display was created")
    with np.load(source, allow_pickle=False) as saved:
        if "endpoint_pose" not in saved.files:
            return not_run(stem, [source], "endpoint_pose was not recorded; no substitute endpoint is used")
        endpoint_pose = np.asarray(saved["endpoint_pose"], dtype=float).copy()
    if endpoint_pose.shape != (2,) or not np.all(np.isfinite(endpoint_pose)):
        raise ValueError("C4_1 endpoint_pose must be a finite two-coordinate array")
    # Use the exact campaign World coefficient, which imports only public
    # geometry and has no observation/truth loader or physical solver.
    sys.path.insert(0, str(ROOT / "vendor"))
    from wslam.world import World
    module_path = Path(sys.modules[World.__module__].__file__).resolve()
    if module_path != (ROOT / "vendor/wslam/world.py").resolve():
        raise RuntimeError("Material display must use the frozen campaign vendor World")
    baseline = World(np.array([[1, 0], [1, 1]], bool), 3.75, .125, (.2, 2.),
                     radius=.12, contrast=3., target_loss=.02)
    endpoint = World(np.array([[1, 0], [1, 1]], bool), 3.85, .125, endpoint_pose,
                     radius=.12, contrast=3., target_loss=.02)
    same_pose = bool(np.array_equal(endpoint_pose, baseline.pose))
    pair_caption = ("This pair changes wall material only; no new displaced-target multi-solution is established."
                    if same_pose else "Recorded wall/pose endpoint display; no imaging recovery or target-compensation outcome is asserted.")
    pair_caption_zh = ("本pair仅wall材料变化；未建立新目标位移多解" if same_pose
                       else "仅显示已记录wall与pose端点；不宣称成像恢复或目标补偿成功")
    # 0.01 lambda0 is a display-sampling interval, never a new PDE mesh.
    x = np.linspace(MATERIAL_BOX[0], MATERIAL_BOX[1], 401)
    y = np.linspace(MATERIAL_BOX[2], MATERIAL_BOX[3], 276)
    xx, yy = np.meshgrid(x, y)
    locations = np.column_stack((xx.ravel(), yy.ravel()))
    first = baseline.coefficient(FREQUENCY, POLYGON_N)(locations).reshape(xx.shape)
    second = endpoint.coefficient(FREQUENCY, POLYGON_N)(locations).reshape(xx.shape)
    difference = second - first
    if not (np.isfinite(first).all() and np.isfinite(second).all()):
        raise ValueError("Nonfinite coefficient in the actual material display")
    data_path = output / f"{stem}_DATA.npz"
    np.savez_compressed(data_path, x=x, y=y, epsilon_baseline=first, epsilon_endpoint=second,
                        epsilon_endpoint_minus_baseline=difference, baseline_pose=baseline.pose,
                        endpoint_pose=endpoint_pose, support=baseline.support,
                        frequency=np.array(FREQUENCY), polygon_n=np.array(POLYGON_N),
                        material_box=np.array(MATERIAL_BOX))
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 8.1), constrained_layout=True)
    for row, component, label in ((0, np.real, r"Re $\epsilon_r$"), (1, np.imag, r"Im $\epsilon_r$")):
        left, right, delta = component(first), component(second), component(difference)
        norm = Normalize(vmin=min(float(left.min()), float(right.min())),
                         vmax=max(float(left.max()), float(right.max())))
        image = None
        for col, values in ((0, left), (1, right)):
            image = axes[row, col].pcolormesh(x, y, values, shading="nearest", cmap=MAGNITUDE_CMAP,
                                             norm=norm, rasterized=True)
        fig.colorbar(image, ax=list(axes[row, :2]), label=label, shrink=.82)
        maximum = float(np.max(np.abs(delta)))
        if maximum == 0:
            # A single zero bin avoids inventing a nonzero difference range.
            cmap = ListedColormap([plt.get_cmap(DIFFERENCE_CMAP)(.5)])
            delta_image = axes[row, 2].pcolormesh(x, y, delta, shading="nearest", cmap=cmap,
                norm=BoundaryNorm([-.5, .5], 1), rasterized=True)
            fig.colorbar(delta_image, ax=axes[row, 2], ticks=[0], label=f"Delta {label}: all samples exactly zero", shrink=.82)
        else:
            delta_image = axes[row, 2].pcolormesh(x, y, delta, shading="nearest", cmap=DIFFERENCE_CMAP,
                norm=Normalize(-maximum, maximum), rasterized=True)
            fig.colorbar(delta_image, ax=axes[row, 2], label=f"Delta {label}", shrink=.82)
        for col in range(3):
            format_axes(axes[row, col])
    axes[0, 0].set_title(r"Baseline $\epsilon_\infty=3.75$, target $(.2,2)$")
    axes[0, 1].set_title(f"Endpoint epsinf=3.85\ntarget ({endpoint_pose[0]:.5g}, {endpoint_pose[1]:.5g})")
    axes[0, 2].set_title("Endpoint minus baseline")
    fig.suptitle("C4_1: actual prescribed materials, Re and Im epsilon\n"
                 "Known-support narrow TX0/RX8 mechanism fixture; coefficient display sampling only\n"
                 + pair_caption, fontsize=11)
    figures = save_plot(fig, output, stem)
    return dict(name=stem, status="CREATED", figures=figures, data_npz=rel(data_path), input=rel(source),
                coefficient_source=rel(module_path), frequency=FREQUENCY, polygon_n=POLYGON_N,
                display_grid_shape=list(xx.shape), display_grid_step_lambda0=.01,
                baseline_world=baseline.metadata(), endpoint_world=endpoint.metadata(),
                same_target_pose=same_pose, pair_caption=pair_caption, pair_caption_zh=pair_caption_zh,
                real_difference_range=limits(difference.real), imaginary_difference_range=limits(difference.imag),
                unknown_wall_image_claim=False, full_current_ambiguity_claim=False,
                target_compensation_success_claim=False, new_physical_solve=False)


def stored_field(output):
    field_path = ROOT / "results/g0_qualification/strong_center.npz"
    mesh_path = ROOT / "results/preflight/FIXED_MESH.npz"
    stem = "STRONG_CENTER_SOURCE0_ACTUAL_TOTAL_FIELD"
    inputs = [field_path, mesh_path]
    if any(not path.exists() for path in inputs):
        return not_run(stem, inputs, "The fixed field or original fixed mesh is absent; no field display was created")
    with np.load(field_path, allow_pickle=False) as saved:
        if "fixed_field" not in saved.files:
            return not_run(stem, inputs, "fixed_field was not stored")
        free_field = saved["fixed_field"].copy()
    with np.load(mesh_path, allow_pickle=False) as saved:
        if not {"points", "triangles"}.issubset(saved.files):
            return not_run(stem, inputs, "Original mesh coordinates/connectivity were not stored")
        points, triangles = saved["points"].copy(), saved["triangles"].copy()
    if (points.ndim != 2 or points.shape[0] != 2 or triangles.ndim != 2 or triangles.shape[0] != 3
            or not np.isfinite(points).all() or triangles.dtype.kind not in "iu"):
        raise ValueError("Saved mesh coordinates/connectivity have an invalid shape or type")
    if (free_field.ndim != 2 or free_field.shape[1] != 2 or not np.iscomplexobj(free_field)
            or not np.isfinite(free_field).all()):
        raise ValueError("The stored fixed_field must contain finite complex free-DOF-by-two-source values")
    from skfem import Basis, ElementTriP3, MeshTri
    mesh = MeshTri(points, triangles)
    if not (np.array_equal(mesh.p, points) and np.array_equal(mesh.t, triangles)):
        raise ValueError("Reconstructed mesh changed its recorded ordering; saved field cannot be mapped safely")
    # Reproduce only the original indexing operations.  No forms/asm/FEM/LU.
    basis = Basis(mesh, ElementTriP3(), intorder=8)
    boundary = basis.get_dofs().all()
    free = np.setdiff1d(np.arange(basis.N), boundary)
    if len(free) != free_field.shape[0]:
        raise ValueError("Reconstructed P3 free-DOF count does not match stored fixed_field")
    scattered = np.zeros((basis.N, 2), dtype=np.complex128)
    scattered[free] = free_field
    dof_points = basis.doflocs.T.copy()
    xl, xr, yb, yt = MATERIAL_BOX
    inside = ((dof_points[:, 0] >= xl) & (dof_points[:, 0] <= xr)
              & (dof_points[:, 1] >= yb) & (dof_points[:, 1] <= yt))
    material_indices = np.flatnonzero(inside)
    display_points = dof_points[inside]
    distance = np.linalg.norm(display_points - SOURCE0, axis=1)
    if not len(distance) or np.any(distance == 0):
        raise ValueError("Material-box field nodes must exclude the singular line source")
    # Exact formula in frozen vendor line_green, evaluated only in y >= .25.
    incident = .25j * hankel1(0, 2*np.pi*FREQUENCY*distance)
    total = incident + scattered[inside, 0]
    if not np.isfinite(total).all():
        raise ValueError("The stored-field plus analytic-incident sum is nonfinite")
    amplitude = np.abs(total)
    phase = np.angle(total)
    undefined_phase = amplitude == 0
    phase[undefined_phase] = np.nan
    # Triangulate the actual P3 node locations for rendering only.  This is
    # neither an FE space change nor a new P1 solve/interpolated PDE solution.
    display_tri = mtri.Triangulation(display_points[:, 0], display_points[:, 1])
    phase_tri = mtri.Triangulation(display_points[:, 0], display_points[:, 1], triangles=display_tri.triangles.copy())
    if np.any(undefined_phase):
        phase_tri.set_mask(np.any(undefined_phase[phase_tri.triangles], axis=1))
    data_path = output / f"{stem}_DATA.npz"
    np.savez_compressed(data_path, original_mesh_points=points, original_mesh_triangles=triangles,
        global_p3_dof_points=dof_points, free_dof_indices=free, boundary_dof_indices=boundary,
        saved_free_scattered_field=free_field, mapped_global_scattered_dofs=scattered,
        material_global_dof_indices=material_indices, display_points=display_points,
        display_triangles=display_tri.triangles, source0=SOURCE0, source0_incident=incident,
        source0_scattered=scattered[inside, 0], source0_total=total, source0_total_magnitude=amplitude,
        source0_total_phase_rad=phase, undefined_phase=undefined_phase,
        frequency=np.array(FREQUENCY), material_box=np.array(MATERIAL_BOX))
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.5), constrained_layout=True)
    magnitude_image = axes[0].tripcolor(display_tri, amplitude, shading="gouraud", cmap=MAGNITUDE_CMAP,
                                       rasterized=True)
    phase_image = axes[1].tripcolor(phase_tri, np.where(undefined_phase, 0, phase), shading="gouraud",
                                   cmap=PHASE_CMAP, vmin=-np.pi, vmax=np.pi, rasterized=True)
    fig.colorbar(magnitude_image, ax=axes[0], label=r"$|E_{total}|$ (unit line source)", shrink=.88)
    phase_bar = fig.colorbar(phase_image, ax=axes[1], label=r"arg $E_{total}$ / rad", shrink=.88,
                             ticks=[-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
    phase_bar.ax.set_yticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"])
    for ax in axes:
        format_axes(ax)
    axes[0].set_title("Source 0: total-field magnitude")
    axes[1].set_title("Source 0: total-field principal phase")
    fig.suptitle("strong_center: actual saved P3 scattered DOFs + analytic incident\n"
                 "Full 2x2 wall; source (-1.5, 0). Material-box nodes only; triangulation is display only.", fontsize=12)
    figures = save_plot(fig, output, stem)
    return dict(name=stem, status="CREATED", figures=figures, data_npz=rel(data_path),
        inputs=[rel(path) for path in inputs], full_dofs=int(basis.N), free_dofs=len(free),
        material_box_nodes=len(display_points), minimum_source_distance_lambda0=float(distance.min()),
        source0=SOURCE0.tolist(), frequency=FREQUENCY, magnitude_range=limits(amplitude), phase_range_rad=limits(phase),
        undefined_phase_nodes=int(undefined_phase.sum()), incident_formula="0.25j*hankel1(0, 2*pi*f*distance)",
        field_identity="saved P3 scattered nodal field plus analytic unit line incident; source 0",
        display="Gouraud colors on a Matplotlib triangulation of actual P3 node locations; no new P1 solve",
        new_physical_solve=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/g0_material_fields")
    parser.add_argument("--summary-dir", type=Path, default=ROOT / "research/delegated/results_support")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.summary_dir.mkdir(parents=True, exist_ok=True)
    records = [materials(args.output_dir), stored_field(args.output_dir)]
    summary = dict(status="ACTUAL_FIXED_CASE_DISPLAY_SNAPSHOT", records=records,
        cases_fixed_before_display=["C4_1 materials", "strong_center source0 stored field"],
        cases_are_different=True, coordinate_units="lambda0", new_full_fe=0, new_rhs=0,
        physical_solves_executed=False, unknown_wall_image_claim=False, gates_modified=False,
        cpu_accounting="existing serial wrapper completion, once; display sampling is not a physical solve",
        current_artifact_policy="Only CREATED entries describe this invocation; NOT_RUN does not refresh any preexisting plot.")
    destination = args.summary_dir / "MATERIAL_FIELD_FIGURE_SUMMARY.json"
    destination.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(summary=rel(destination), created_figures=sum(row["status"] == "CREATED" for row in records),
                          not_run=sum(row["status"] == "NOT_RUN" for row in records), new_full_fe=0, new_rhs=0)))


if __name__ == "__main__":
    main()

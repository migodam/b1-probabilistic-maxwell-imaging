"""Pack the redistributable physical revision at its original relative paths."""

import argparse
import json
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
ROOT = STUDY.parents[1]
OUTPUT = ROOT / ".artifacts/z_alg_physical_revision/physical_revision_evidence.zip"

REVISION_FILES = {
    ".": "REPRODUCE.md HANDOFF.md package_release.py",
    "array": """
        ARRAY_SCAN_EN.md ARRAY_SCAN_ZH.md FREEZE.md array_scaling.pdf array_scaling.png
        array_scan.csv bare_inertia.json bare_inertia.py forward_audit.json
        gap_refinement_p0.json gap_refinement_p20.json gap_refinement_p32.json
        high_grid_p20.json high_grid_p20_recovery.json high_grid_p32.json high_grid_p32_recovery.json
        native_p20.json native_p32.json prepare.py refine_gaps.py screen.json screen.py
        screen_scalar_grid.npz source_spectra.csv source_spectra.pdf source_spectra.png
        summarize.py summary.json verify_native.py
    """,
    "born": """
        REPORT.md born_bounds.csv cost.json directions.npz extinction.csv extinction.py
        extinction_geometry.json fixed_port_weights.npz run.py
    """,
    "born/controls": "FROZEN_DESIGN.md RESULTS.md arrays.npz replay.py results.json",
    "born/ports": "REPORT_ZH.md replay.py results.json",
    "electrical": """
        ANALYSIS.md FIGURE_NOTES.md angular.py angular_interpolation.csv chamber_by_frequency.csv
        electrical_scaling_en.pdf electrical_scaling_en.png electrical_scaling_zh.pdf electrical_scaling_zh.png
        geometry.csv manitoba_large_nylon_cap_en.pdf manitoba_large_nylon_cap_en.png
        manitoba_large_nylon_cap_zh.pdf manitoba_large_nylon_cap_zh.png
        plot_electrical.py primary_metadata.json scale.py state_counts.csv
    """,
    "inversion": """
        COMPACT_FOLLOWUP.md PROTOCOL.md REGISTRATION_ZH.md RESULTS_EN.md RESULTS_ZH.md
        SYNTHETIC_DESIGN.md compact.py compact_followup.py fit.py gate.py
        original_gate.json original_offline_scores.json plot.py score_original.py
    """,
    "inversion/compact_synthetic": "gate.json offline_scores.json",
    "inversion/compact_synthetic/offline_truth": "generation.json linear_diagnostic.json reference.npz",
    "modes": """
        REPORT_EN.md REPORT_ZH.md analyze.py fourier_dark_subspaces.pdf fourier_dark_subspaces.png
        native_contrast_maps.pdf native_contrast_maps.png negative_state_currents.pdf
        negative_state_currents.png plot.py plot_data.npz plot_timing.json results.json
        retained_source_spectra.pdf retained_source_spectra.png
    """,
}

PRIOR_FILES = {
    "core": "certificate.py inertia.py",
    "fresnel": "physical.py REPORT_ZH.md",
    "results/fresnel": "ERROR_LEDGER.csv SENSITIVITY.csv",
    "results/theory": "fine_feasibility.json fine_verified.json",
    "submission/curve": """
        curve_inputs.py envelope.py fixed_blocks.py global.5.json global.75.json
        local.75.json native.py noise_curve.csv noise_curve_summary.json observation_audit.json search.py
    """,
    "submission/single": "coherent_general_profiles.json coherent_general_verified.json",
    "synthetic": "FREEZE.md README.md nonvacuity.py",
    "validation/manitoba": "calibration_margins.csv physical.py physical_audit.json",
    "validation/srap": "online.py",
}

PRIOR_R8 = Path("research/delegated/max_struct_r1/r8_worktree/research")
PACKED_R8 = ROOT / ".artifacts/z_alg_optical_release/dependencies" / PRIOR_R8
R8_FILES = {
    "research_fbg_r3_experiments": "__init__.py cie_tsom_baseline.py physics.py",
    "research_max_struct_r2": "cie_global.py",
    "research_max_struct_r4": "cie_onehop.py",
    "research_max_struct_r6": "fft_tsom_ops.py",
    "research_max_struct_r8": "fresnel_adapter.py paper_physics.py",
    "research_z_alg_r4_diversity": "__init__.py physics.py solver.py",
}


def named(base, selection):
    return [base / dirname / file for dirname, names in selection.items() for file in names.split()]


def observations(base):
    rows = json.loads((base / "manifest.json").read_text(encoding="utf-8"))["observations"]
    return [base / "manifest.json", *(base / Path(row["file"]).name for row in rows)]


def members():
    files = named(HERE, REVISION_FILES) + named(STUDY, PRIOR_FILES)
    files += observations(HERE / "array/observations")
    files.append(HERE / "array/observations/arrays.npz")
    files += observations(STUDY / "submission/curve/observations")
    # This entire upstream synthetic cohort has already been released.
    for subdir in ("public", "evaluation"):
        files += sorted((STUDY / "synthetic" / subdir).glob("*.npz"))
        files += sorted((STUDY / "synthetic" / subdir).glob("*.json"))
    inversion = HERE / "inversion"
    for subdir, patterns in (
        ("figures", ("*.pdf", "*.png")),
        ("original_fits", ("optfresh_*.json", "optfresh_*.npz")),
        ("compact_synthetic/fits", ("compact_tm_*.json", "compact_tm_*.npz")),
        ("compact_synthetic/exposed_followup", ("compact_tm_*.json", "compact_tm_*.npz", "offline_scores.json")),
    ):
        for pattern in patterns:
            files += sorted((inversion / subdir).glob(pattern))
    compact = inversion / "compact_synthetic/public"
    ids = json.loads((compact/"manifest.json").read_text(encoding="utf-8"))["observations"]
    files.extend([compact / "manifest.json", compact / "geometry.npz"])
    files.extend(compact / (r["id"]+".npz") for r in ids)
    files.extend([
        ROOT / "research/z_alg_secant_tasks/highdim/public_vie_inputs.npz",
        ROOT / "research/z_alg_xudong_cie/minpeak_native/search.py",
        ROOT / ("research/delegated/max_struct_r1/r14_worktree/research/"
                "research_fbg_r3_experiments/src/cie_tsom_baseline.py"),
    ])
    sources = {f.relative_to(ROOT): f for f in files}
    # Use the R8 source snapshot from the original public evidence package.
    for project, items in R8_FILES.items():
        for name in items.split():
            path = Path(project) / "src" / name
            archived = PACKED_R8 / path
            sources[PRIOR_R8 / path] = archived if archived.is_file() else ROOT / PRIOR_R8 / path
    return sources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    sources = members()
    missing = [str(rel) for rel, path in sources.items() if not path.is_file()]
    if missing:
        parser.error("Missing release source files:\n" + "\n".join(missing))
    output = args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".zip.tmp")
    with ZipFile(temporary, "w", compression=ZIP_DEFLATED, compresslevel=6) as archive:
        for rel, source in sorted(sources.items()):
            archive.write(source, rel.as_posix())
    temporary.replace(output)
    print(f"{output.resolve()} · {len(sources)} files · {output.stat().st_size:,} bytes")


if __name__ == "__main__":
    main()

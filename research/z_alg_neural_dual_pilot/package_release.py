#!/usr/bin/env python3
"""Assemble the standalone Stage 2 replay ZIP from selected public inputs."""
from __future__ import annotations

import argparse
import ast
from contextlib import nullcontext
import json
from pathlib import Path, PurePosixPath
from tempfile import NamedTemporaryFile
from zipfile import ZIP_DEFLATED, ZipFile


PILOT = Path(__file__).resolve().parent
ROOT = PILOT.parents[1]
PILOT_REL = PILOT.relative_to(ROOT).as_posix()
PAPER = "research/z_alg_optical_inertia"
ARCHIVE_NAME = "Z_ALG_Optical_Inertia_Validation.zip"
BASE_ARCHIVE = ROOT / ".artifacts/z_alg_optical_validation_release" / ARCHIVE_NAME
BASE_RELEASE = (
    "https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/"
    "z-alg-optical-inertia-validation-v1"
)
OUTPUT = ROOT / ".artifacts/z_alg_neural_dual_pilot_release/Z_ALG_Neural_Dual_Pilot.zip"

PAIR_DIR = f"{PAPER}/submission/curve/observations"
PAIR_MANIFEST = f"{PAIR_DIR}/manifest.json"
SRAP_RECORD = f"{PAPER}/validation/srap/online.json"
PAPER_MEMBERS = (
    "core/certificate.py",
    "core/inertia.py",
    "build_site.py",
    "research.css",
    "synthetic/nonvacuity.py",
    "synthetic/public/coarse_green.npz",
    "synthetic/public/fine_geometry.npz",
    "synthetic/public/pw360_10x20.npz",
    "submission/curve/curve_inputs.py",
    "submission/curve/envelope.py",
    "submission/curve/fixed_blocks.py",
    "submission/curve/native.py",
    "submission/curve/search.py",
    "submission/curve/fixed_blocks.75.json",
    "submission/curve/fixed_blocks_higher_caps.json",
    "submission/curve/global.5.json",
    "submission/curve/global.75.json",
    "submission/highcaps/certificates.json",
    "submission/single/coherent_general_profiles.json",
    "submission/single/coherent_general_verified.json",
    "validation/srap/online.py",
    "validation/srap/offline_feasible.py",
    "validation/srap/online.json",
    "validation/srap/public_observations.npz",
)
SRAP = "research/srap_cie_array_study_20261121"
R8_SNAPSHOT = "research/delegated/max_struct_r1/r8_worktree"
GENERATOR_SOURCES = (
    f"{SRAP}/legacy_la.py",
    f"{SRAP}/data2d.py",
    f"{SRAP}/data_followup.py",
)
R8_IMPORTS = (
    "research/research_fbg_r3_experiments/src/__init__.py",
    "research/research_fbg_r3_experiments/src/cie_tsom_baseline.py",
    "research/research_fbg_r3_experiments/src/physics.py",
    "research/research_max_struct_r2/src/cie_global.py",
    "research/research_max_struct_r4/src/cie_onehop.py",
    "research/research_max_struct_r6/src/fft_tsom_ops.py",
    "research/research_max_struct_r8/src/fresnel_adapter.py",
    "research/research_max_struct_r8/src/paper_physics.py",
    "research/research_z_alg_r4_diversity/src/__init__.py",
    "research/research_z_alg_r4_diversity/src/physics.py",
    "research/research_z_alg_r4_diversity/src/solver.py",
)

PILOT_REQUIRED = (
    "REPRODUCE.md",
    "METHODS_TRAINING_ZH.md",
    "PILOT_DESIGN.md",
    "COARSE_METHOD.md",
    "NATIVE_VERIFICATION.md",
    "SUPPLEMENT_DESIGN.md",
    "package_release.py",
    "data/generation_receipt.json",
    "data/synthetic_observations.npz",
    "results/coarse_gradient_check.json",
    "results/fixed_baselines.csv",
    "results/fixed_baselines.json",
    "results/holdout_open.json",
    "results/native_restriction_training_check.json",
    "results/native_training_check.json",
    "results/native_verifier_check.json",
    "results/neural_dual_training.json",
    "results/neural_dual_weights.pt",
    "results/no_feature_training.json",
    "results/no_feature_weights.pt",
    "results/training_checkpoint_audit.json",
    "results/training_recipe_frozen.json",
    "src/chart.py",
    "src/check_coarse_gradients.py",
    "src/check_native_restriction.py",
    "src/coarse.py",
    "src/fixed_baselines.py",
    "src/generate.py",
    "src/generate_srap_bands.py",
    "src/holdout.py",
    "src/native.py",
    "src/native_check.py",
    "src/network.py",
    "src/projected_upper.py",
    "src/recover_unknown.py",
    "src/propose.py",
    "src/score_batch.py",
    "src/train.py",
    "src/validation_swaps.py",
)
PILOT_FINAL = (
    "README.md",
    "RESEARCH_REPORT_ZH.md",
    "FOLLOWUP_RESULTS_ZH.md",
    "GEOMETRY_INERTIA_ZH.md",
    "EXPERIMENTAL_RESULTS.md",
    "MATHEMATICAL_NOTES_ZH.md",
    "NEXT_PRO_HANDOFF.md",
    "ADVISOR_ZH.md",
    "data/srap_supplement.npz",
    "data/srap_supplement.json",
    "data/followup_observations.npz",
    "data/followup_external8.npz",
    "data/transfer16.npz",
    "results/native_recovery_training_check.json",
    "results/validation_swaps.json",
    "results/native_recovery_bordered_training_check.json",
    "results/cold_latency.json",
    "results/geometry_inertia.json",
    "results/paired_timing.csv",
    "src/cold_latency.py",
    "src/diagnose_inertia.py",
    "src/paired_tables.py",
    "src/figures.py",
    "src/render_report.py",
    "src/summarize.py",
    "results/cases.csv",
    "results/summary.json",
    "results/supplement/fixed_baselines.json",
    "results/supplement/fixed_baselines.csv",
)
FOLLOWUP_REQUIRED = (
    "followup/DESIGN_ZH.md",
    "followup/THEORY_ZH.md",
    "followup/NATIVE_JENSEN_REPLAY.md",
    "followup/NESTED_SEGMENT_COVER.md",
    "followup/MASKED_PORTS_REPLAY.md",
    "followup/MASKED_SEGMENT.md",
    "followup/GALERKIN_COVER_ZH.md",
    "followup/TRANSFERRED_WEIGHT_ZH.md",
    "followup/SEGMENT_COVER_FROZEN.json",
    "followup/MASKED_PORTS_FROZEN.json",
    "followup/GALERKIN_FROZEN.json",
    "followup/TRANSFER16_DESIGN.json",
    "followup/TRANSFERRED_WEIGHT_FROZEN.json",
    "followup/generation.json",
    "followup/external8_generation.json",
    "followup/offline_scores.json",
    "followup/external8_offline_scores.json",
    "followup/transfer16_generation.json",
    "followup/transfer16_offline.json",
    "followup/generate_fields.py",
    "followup/generate_transfer16.py",
    "followup/jensen_port.py",
    "followup/segment_cover.py",
    "followup/masked_ports.py",
    "followup/masked_segment.py",
    "followup/galerkin_cover.py",
    "followup/plot_results.py",
    "followup/transferred_weight.py",
    "followup/results/port_enclosures.json",
    "followup/results/segment_enclosures.json",
    "followup/results/galerkin_enclosures.json.gz",
    "followup/results/archived_endpoint_crosscheck.json",
    "followup/results/masked_geometry.json",
    "followup/results/masked_continuous_inertia.json",
    "followup/results/masked_direct_checks.json",
    "followup/results/masked_segment.json",
    "followup/results/galerkin_paired_summary.json",
    "followup/results/transferred_weight_common.json",
    "followup/results/transferred_weight_offline_audit.json",
    *(f"followup/results/transferred_weight_{cohort}.json"
      for cohort in ("donor", "historical", "prior32", "external8", "transfer16")),
    *(f"followup/results/{cohort}_jensen.json" for cohort in ("validation", "fresh", "historical")),
    *(f"followup/results/{cohort}_{method}.json"
      for cohort in ("validation", "fresh", "historical", "external8")
      for method in ("segment_cover", "masked_ports")),
    *(f"followup/results/galerkin_{cohort}.json"
      for cohort in ("validation", "fresh", "historical", "external8", "validation_focus", "historical_focus")),
)
TRACE_DIRS = (
    ("proposals", ".jsonl"),
    ("optimization_traces", ".jsonl.gz"),
    ("native", ".jsonl"),
)
COUNTS = {"paired_objects": 6, "transfer_objects": 22, "public_observations": 40}


def public_dependencies(read):
    manifest = json.loads(read(PAIR_MANIFEST))
    paired = [row for row in manifest["observations"] if row["noise_pct"] in (1, 2, 10)]
    paired_objects = {row["object_id"] for row in paired}
    srap = json.loads(read(SRAP_RECORD))["objects"]
    transfer_objects = {row["object_id"] for row in srap}
    if (len(paired_objects) != COUNTS["paired_objects"]
            or len(transfer_objects) != COUNTS["transfer_objects"]
            or len(srap) != len(transfer_objects)
            or len(paired) + len(srap) != COUNTS["public_observations"]
            or paired_objects & transfer_objects):
        raise ValueError("The public dependencies do not match the frozen 6+22 pilot cohort.")
    for row in paired:
        name = row["file"]
        if PurePosixPath(name).name != name or not name.endswith(".npz"):
            raise ValueError(f"Invalid paired-input filename {name!r}.")
    paired_paths = [f"{PAIR_DIR}/{row['file']}" for row in paired]
    files = [*GENERATOR_SOURCES, *(f"{PAPER}/{name}" for name in PAPER_MEMBERS), *paired_paths]
    for name in files:
        yield name, read(name)
    for name in R8_IMPORTS:
        yield name, read(name, archived_name=f"{R8_SNAPSHOT}/{name}")
    trimmed = dict(manifest, observations=paired)
    yield PAIR_MANIFEST, (json.dumps(trimmed, indent=2, ensure_ascii=False) + "\n").encode()


def rendered_pages():
    # Keep the archive routes synchronized with the renderer without importing
    # its templating dependencies or executing it during packaging.
    program = ast.parse((PILOT / "src/render_report.py").read_text())
    route_lists = {}
    for node in program.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id in ("MAIN", "DETAILS"):
                route_lists[target.id] = ast.literal_eval(node.value)
    if set(route_lists) != {"MAIN", "DETAILS"}:
        raise ValueError("HTML route lists are missing from the research renderer.")
    pages = [page for group in ("MAIN", "DETAILS") for _, page, _ in route_lists[group]]
    if len(pages) != len(set(pages)) or any(PurePosixPath(p).name != p or not p.endswith(".html") for p in pages):
        raise ValueError("HTML report route names are ambiguous.")
    return pages


def pilot_paths(objects, supplement_objects, preview):
    ready = list(PILOT_REQUIRED)
    expected = [*PILOT_FINAL, *FOLLOWUP_REQUIRED, *rendered_pages()]
    for folder, cohort in (("results", objects), ("results/supplement", supplement_objects)):
        expected.extend(f"{folder}/{kind}/{oid}{suffix}"
                        for oid in sorted(cohort) for kind, suffix in TRACE_DIRS)

    if preview:
        ready.extend(name for name in expected if (PILOT / name).is_file())
    else:
        missing = [name for name in expected if not (PILOT / name).is_file()]
        if missing:
            raise FileNotFoundError("Missing final Stage 2 files: " + ", ".join(missing))
        ready.extend(expected)
        worker_patterns = ("proposals/setup_*.json", "proposals/worker_*.json", "native/worker-*.json", "recovery/*.jsonl", "recovery/worker-*.json")
        for folder in ("results", "results/supplement"):
            marker = f"{folder}/proposals_done.json"
            if (PILOT / marker).is_file():
                ready.append(marker)
            for pattern in worker_patterns:
                ready.extend(p.relative_to(PILOT).as_posix()
                             for p in sorted(PILOT.glob(f"{folder}/{pattern}")) if p.is_file())
    for pattern in ("*.md", "*.py", "*.json", "results/**/*.json", "results/**/*.jsonl",
                    "results/**/*.jsonl.gz", "results/**/*.json.gz", "results/**/*.npz"):
        ready.extend(p.relative_to(PILOT).as_posix()
                     for p in sorted((PILOT / "followup").glob(pattern)) if p.is_file())
    for location in ("figures", "followup/figures"):
        for suffix in (".png", ".pdf", ".svg"):
            ready.extend(p.relative_to(PILOT).as_posix()
                         for p in sorted(PILOT.glob(f"{location}/*{suffix}")) if p.is_file())
    for name in dict.fromkeys(ready):
        path = PILOT / name
        if not path.is_file() or path.is_symlink():
            raise FileNotFoundError(f"Required regular pilot output is missing: {path}")
        yield name, path


def package(output, base_archive, preview=False, check_inputs=False):
    if base_archive is None and BASE_ARCHIVE.is_file():
        base_archive = BASE_ARCHIVE
    if base_archive is None and not (ROOT / "PACKAGE_CONTENTS.json").is_file():
        raise FileNotFoundError(f"Locate the local public validation archive with --base-archive: {BASE_ARCHIVE}")
    source_context = ZipFile(base_archive) if base_archive is not None else nullcontext()
    with source_context as archived:
        def read(name, *, archived_name=None):
            if archived is not None and name != f"{SRAP}/data_followup.py":
                return archived.read(archived_name or name)
            return (ROOT / name).read_bytes()

        external = list(public_dependencies(read))
        paired = json.loads(dict(external)[PAIR_MANIFEST])["observations"]
        srap = json.loads(dict(external)[SRAP_RECORD])["objects"]
        objects = {row["object_id"] for row in paired + srap}
        supplement = {row["object_id"] for row in srap}
        main = list(pilot_paths(objects, supplement, preview))
        if check_inputs:
            print(f"Inputs ready: {len(main)} pilot files and {len(external)} public dependencies.")
            return 0
        output.parent.mkdir(parents=True, exist_ok=True)
        with NamedTemporaryFile(dir=output.parent, prefix="." + output.name + ".", delete=False) as temp:
            pending = Path(temp.name)
        try:
            with ZipFile(pending, "w", compression=ZIP_DEFLATED, compresslevel=6) as z:
                for name, data in external:
                    z.writestr(name, data)
                for name, path in main:
                    z.write(path, f"{PILOT_REL}/{name}")
                z.writestr("README.md", "# Z-ALG 双乘子复现包\n\n"
                           f"从[研究网页]({PILOT_REL}/index.html)或[复现说明]({PILOT_REL}/REPRODUCE.md)进入。\n")
                listed = [dict(path=info.filename, bytes=info.file_size, compressed_bytes=info.compress_size)
                          for info in z.infolist()]
                z.writestr("PACKAGE_CONTENTS.json", json.dumps({
                    "status": "local_preview" if preview else "frozen_replay_packet",
                    "base_public_release": BASE_RELEASE,
                    "base_archive": ARCHIVE_NAME,
                    "copied_dependency_count": len(external),
                    "open_inputs": COUNTS,
                    "supplement_observations_expected": 2 * len(supplement),
                    "followup_inputs": {"fresh_objects": 32, "fresh_views": 96,
                                        "external_objects": 8, "external_views": 24,
                                        "transfer16_objects": 16, "transfer16_views": 48},
                    "r8_source_mapping": "Eleven archived Python modules installed at their canonical research import paths.",
                    "paired_manifest_selection": "Preserved the public rows for 1%, 2%, and 10%.",
                    "members": listed,
                }, indent=2, ensure_ascii=False) + "\n")
            pending.replace(output)
        finally:
            pending.unlink(missing_ok=True)
    output_bytes = output.stat().st_size
    print(f"{output}  {output_bytes:,} bytes, {len(listed)+1} entries "
          f"({len(external)} archived dependencies)")
    return output_bytes


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--base-archive", type=Path,
                     help="Archived public validation ZIP; an extracted packet can use its embedded copies.")
    cli.add_argument("--output", type=Path, default=OUTPUT)
    cli.add_argument("--preview", action="store_true",
                     help="Use only for local dependency inspection while workers are still writing results.")
    cli.add_argument("--check-inputs", action="store_true",
                     help="Check the input allowlist without creating an archive.")
    args = cli.parse_args()
    package(args.output, args.base_archive, preview=args.preview, check_inputs=args.check_inputs)

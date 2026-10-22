#!/usr/bin/env python3
"""Build the P3-L Chinese report without xeCJK or network access.

The input Markdown is read into a snapshot and NEVER modified.  Requires
pandoc, XeLaTeX, the standard TeX packages listed in p3l_header.tex, and the
two supplied Noto Serif CJK font files.  All layout transformations affect
only generated TeX.  Usage from the project root:

  python build/build_report.py
  python build/build_report.py --output output/Research_P3L_Package/Research_P3L_Report_ZH.pdf
  python typesetting/build_report.py --source Research_P3L_Report_ZH.md --font-dir /path/to/fonts --output Research_P3L_Report_ZH.pdf
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def first_existing(candidates: list[Path]) -> Path:
    return next((path for path in candidates if path.exists()), candidates[0])


def run(cmd: list[str], *, cwd: Path, log: Path) -> str:
    result = subprocess.run(cmd, cwd=cwd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, check=False)
    log.write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        print("\n".join(result.stdout.splitlines()[-65:]), file=sys.stderr)
        raise SystemExit(f"Build failed; full log: {log}")
    return result.stdout


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "--source", dest="input", type=Path,
                        default=first_existing([
                            ROOT / "Research_P3L_Report_ZH.md",
                            ROOT / "output/Research_P3L_Package/Research_P3L_Report_ZH.md"]))
    parser.add_argument("--output", type=Path,
                        default=HERE / "Research_P3L_draft.pdf")
    parser.add_argument("--font-dir", type=Path,
                        default=first_existing([ROOT / "fonts", ROOT / "tmp/fonts"]))
    parser.add_argument("--passes", type=int, default=3)
    parser.add_argument("--no-fit-math", action="store_true",
                        help="Keep raw displays; useful for diagnosing natural formula widths.")
    args = parser.parse_args()
    source = args.input.resolve()
    target = args.output.resolve()
    fonts = args.font_dir.resolve()
    for exe in ("pandoc", "xelatex"):
        if not shutil.which(exe):
            raise SystemExit(f"Missing required executable: {exe}")
    for filename in ("NotoSerifCJKsc-Regular.otf", "NotoSerifCJKsc-Bold.otf"):
        if not (fonts / filename).is_file():
            raise SystemExit(f"Missing font: {fonts / filename}")
    # The generated TeX embeds local font paths, not user-controlled commands.
    if any(c in str(fonts) for c in "{}%\\\n"):
        raise SystemExit("Font path contains TeX-special characters; use a plain local directory.")
    if args.passes < 2:
        raise SystemExit("Use at least two passes so the table of contents is populated.")
    target.parent.mkdir(parents=True, exist_ok=True)
    work = HERE / "latex"
    work.mkdir(parents=True, exist_ok=True)
    raw = source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    snapshot = work / "report_snapshot.md"
    snapshot.write_bytes(raw)
    font_header = work / "font_path.tex"
    font_header.write_text("\\newcommand{\\PThreeFontPath}{" + fonts.as_posix() + "/}\n",
                           encoding="utf-8")
    tex_path = work / "report.tex"
    start = time.monotonic()
    pandoc_cmd = ["pandoc", str(snapshot), "--from=markdown+tex_math_dollars+raw_tex",
                  "--to=latex", "--standalone", "--no-highlight", "--wrap=none",
                  "--template", str(HERE / "p3l_template.tex"),
                  "--include-in-header", str(font_header),
                  "--include-in-header", str(HERE / "p3l_header.tex"),
                  "--resource-path", os.pathsep.join([str(source.parent), str(ROOT)]),
                  "--output", str(tex_path)]
    run(pandoc_cmd, cwd=ROOT, log=work / "pandoc.log")
    tex = tex_path.read_text(encoding="utf-8")
    resolved_images = []
    def resolve_image(match: re.Match) -> str:
        options, name = match.group(1) or "", match.group(2)
        candidates = [Path(name)] if Path(name).is_absolute() else [source.parent / name, ROOT / name]
        chosen = next((p.resolve() for p in candidates if p.is_file()), None)
        if chosen is None:
            raise SystemExit(f"Missing image referenced by the report: {name}")
        if any(c in str(chosen) for c in "{}%\\\n"):
            raise SystemExit(f"Image path contains TeX-special characters: {chosen}")
        resolved_images.append(str(chosen))
        return "\\includegraphics" + options + "{" + chosen.as_posix() + "}"
    tex = re.sub(r"\\includegraphics(\[[^\]]*\])?\{([^}]+)\}", resolve_image, tex)
    # File identifiers stay identical while being allowed to wrap after underscores.
    preamble, document = tex.split("\\begin{document}", 1)
    document = document.replace("\\_", "\\_\\allowbreak{}")
    tex = preamble + "\\begin{document}" + document
    fitted = 0
    if not args.no_fit_math:
        # Displays with explicit equation numbering cannot be enclosed in a box.
        # The current report has no such environments; leave any future ones alone.
        def fit(match: re.Match) -> str:
            nonlocal fitted
            body = match.group(1)
            if re.search(r"\\(?:tag|label)\b|\\begin\{(?:equation|align|gather|multline)", body):
                return match.group(0)
            fitted += 1
            return "\\PThreeDisplay{" + body + "}"
        preamble, document = tex.split("\\begin{document}", 1)
        document = re.sub(r"\\\[(.*?)\\\]", fit, document, flags=re.S)
        tex = preamble + "\\begin{document}" + document
    tex_path.write_text(tex, encoding="utf-8")
    for index in range(1, args.passes + 1):
        run(["xelatex", "-halt-on-error", "-interaction=nonstopmode", "-file-line-error",
             "-no-shell-escape", "report.tex"], cwd=work,
            log=work / f"xelatex_pass{index}.stdout.log")
    pdf = work / "report.pdf"
    if not pdf.is_file():
        raise SystemExit("XeLaTeX returned without a PDF.")
    shutil.copyfile(pdf, target)
    log = (work / "report.log").read_text(encoding="utf-8", errors="replace")
    warnings = []
    for line in log.splitlines():
        if ("Overfull " in line or "Missing character:" in line or
            ("Warning:" in line and not "rerunfilecheck" in line)):
            warnings.append(line)
    pages_match = re.search(r"Output written on .*?\((\d+) pages?", log)
    math_fits = [
        {"display": int(n), "natural_width_pt": float(w), "available_width_pt": float(a),
         "scale": round(float(a) / float(w), 4)}
        for n, w, a in re.findall(r"P3L-MATH-FIT\s+(\d+): natural=([0-9.]+)pt; available=([0-9.]+)pt", log)
    ]
    summary = {
        "input": str(source), "input_sha256": digest,
        "input_bytes": len(raw), "output": str(target),
        "pages": int(pages_match.group(1)) if pages_match else None,
        "elapsed_seconds": round(time.monotonic() - start, 2),
        "display_math_guard_count": fitted, "warnings": warnings,
        "scaled_displays": math_fits,
        "resolved_images": resolved_images,
        "snapshot": str(snapshot), "tex": str(tex_path),
        "source_unchanged_during_build": hashlib.sha256(source.read_bytes()).hexdigest() == digest,
        "font_directory": str(fonts),
    }
    summary_path = HERE / "build_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Built {target}: {summary['pages']} pages in {summary['elapsed_seconds']} s")
    print(f"Read-only source snapshot: {digest[:16]} ({len(raw)} bytes)")
    print(f"Math width guards: {fitted}; layout/font warnings: {len(warnings)}")
    if math_fits:
        print("Display math scaled only when needed: " + ", ".join(
              f"#{x['display']} at {x['scale']:.3f}x" for x in math_fits))
    for warning in warnings[:18]:
        print("  " + warning)
    if len(warnings) > 18:
        print(f"  ... full warning list: {summary_path}")
    if not summary["source_unchanged_during_build"]:
        print("Source changed during this draft build; rerun once the report is complete.")
    print(f"Build summary: {summary_path}")


if __name__ == "__main__":
    main()

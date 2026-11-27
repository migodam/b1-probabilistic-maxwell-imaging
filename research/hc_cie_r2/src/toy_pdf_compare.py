#!/usr/bin/env python3
"""Compare the two locally supplied reports using Poppler."""
import argparse
import json
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory


def command(*args):
    return subprocess.check_output(args, stderr=subprocess.PIPE)


def compare(original, zip_copy):
    original_size, zip_size = original.stat().st_size, zip_copy.stat().st_size
    with TemporaryDirectory() as tmp:
        first, second = [Path(tmp)/part for part in ("standalone", "zip")]
        first.mkdir()
        second.mkdir()
        for source, out in [(original, first), (zip_copy, second)]:
            command("pdftoppm", "-r", "75", "-png", source, out/"page")
        files = sorted(first.glob("*.png"))
        counterparts = sorted(second.glob("*.png"))
        equal = len(files) == len(counterparts) and all(
            a.read_bytes() == b.read_bytes() for a, b in zip(files, counterparts))
    attachment_lists = [command("pdfdetach", "-list", x).decode().splitlines()
                        for x in (original, zip_copy)]
    return dict(
        poppler=subprocess.run(["pdfinfo", "-v"], capture_output=True, check=True,
                               text=True).stderr.splitlines()[0],
        bytes=dict(standalone=original_size, zip=zip_size, identical=(
            original.read_bytes() == zip_copy.read_bytes())),
        text=pdftotext_equal(original, zip_copy),
        pages=dict(standalone=len(files), zip=len(counterparts), dpi=75, identical=equal),
        external_urls_equal=command("pdfinfo", "-url", original)==command("pdfinfo", "-url", zip_copy),
        named_destinations_equal=command("pdfinfo", "-dests", original)==command("pdfinfo", "-dests", zip_copy),
        embedded_files=dict(standalone=attachment_lists[0], zip=attachment_lists[1]),
        standalone_c2pa_filespec=(b"/AFRelationship/C2PA_Manifest" in original.read_bytes()),
        credential_signature_verified=False,
    )


def pdftotext_equal(a, b):
    return command("pdftotext", "-layout", a, "-")==command("pdftotext", "-layout", b, "-")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--standalone", type=Path, required=True)
    p.add_argument("--packaged-pdf", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    result = compare(args.standalone, args.packaged_pdf)
    args.out.write_text(json.dumps(result, indent=2)+"\n")


if __name__ == "__main__":
    main()

"""Rebuild local chamber inputs from the providers' original downloads."""

import argparse
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory
from zipfile import ZipFile

from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path[:0] = [str(HERE/"chambers"), str(HERE.parent/"srap_cie_array_study_20261121")]
import training
import legacy_la


def prepare_fresnel(source):
    required = {f"{case.target}.exp" for case in training.CASES.values() if case.acquisition=="Fresnel"}
    target = training.LEGACY
    for folder in ("public","heldout"):
        (target/folder).mkdir(exist_ok=True,parents=True)
    with TemporaryDirectory(dir=target,prefix="fresnel-replay-") as temp_name:
        temp = Path(temp_name)
        if source.is_dir():
            raw = source
        else:
            raw = temp/"inputs"
            raw.mkdir()
            with ZipFile(source) as archive:
                entries = {Path(name).name:name for name in archive.namelist() if Path(name).name in required}
                for basename,member in entries.items():
                    (raw/basename).write_bytes(archive.read(member))
        missing = sorted(name for name in required if not (raw/name).is_file())
        if missing:
            raise FileNotFoundError("Missing official Fresnel 2005 files: "+", ".join(missing))
        legacy_la.FRESNEL_DIR = raw
        legacy_la.RUNTIME = temp
        for folder in ("public","heldout"):
            (temp/folder).mkdir()
        requested = {case.id for case in training.CASES.values() if case.acquisition=="Fresnel"}
        for case in legacy_la.all_specs():
            name = case["case"]
            if name not in requested:
                continue
            if all((target/folder/f"{name}.npz").is_file() for folder in ("public","heldout")):
                print(name,"already prepared")
                continue
            legacy_la.prepare_one(case)
            for folder in ("public","heldout"):
                path = target/folder/f"{name}.npz"
                if not path.is_file():
                    shutil.copy2(temp/folder/path.name,path)


def prepare_manitoba(source):
    source = source.resolve(strict=True)
    target = training.MANITOBA
    target.mkdir(parents=True,exist_ok=True)
    archive = target/"Tomographic_Dataset_1.zip"
    if not archive.exists():
        shutil.copy2(source,archive)
    elif not archive.samefile(source):
        raise FileExistsError("An original archive is already in the replay directory. "
                              "Use its path explicitly to reuse it.")
    prepared = target/"prepared_private.npz"
    if not prepared.exists():
        physical = training.load_module(training.OPTICAL/"validation/manitoba/physical.py",
                                        "manitoba_rehydration")
        physical.main(archive,target/"physical_replay_local.json",prepared)
    if not training.MANITOBA_TRAIN.exists():
        training.prepare_manitoba_training()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fresnel",type=Path,help="Original supp.zip or a directory of official .exp files")
    parser.add_argument("--manitoba",type=Path,help="Official Tomographic_Dataset_1.zip")
    options=parser.parse_args()
    if options.fresnel is None and options.manitoba is None:
        parser.error("Supply at least one official download.")
    with threadpool_limits(1):
        if options.fresnel is not None:
            prepare_fresnel(options.fresnel)
        if options.manitoba is not None:
            prepare_manitoba(options.manitoba)


if __name__=="__main__":
    main()

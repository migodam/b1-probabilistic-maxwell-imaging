"""Recreate the 18 legacy inputs without rerunning legacy solvers."""
import os
from pathlib import Path
from run_study import ROOT, la

if __name__=="__main__":
    target=Path(os.environ.get("SRAP_HIST",ROOT/"legacy_inputs"))
    la.RUNTIME=target/".runs"
    la.FRESNEL_DIR=Path(os.environ.get("FRESNEL_2005_DIR",la.FRESNEL_DIR))
    for part in ("public","heldout"):(la.RUNTIME/part).mkdir(parents=True,exist_ok=True)
    for c in la.all_specs(): la.prepare_one(c)

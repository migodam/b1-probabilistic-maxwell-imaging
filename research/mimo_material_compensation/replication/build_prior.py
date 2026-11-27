"""Select the exact source and input subset of the 2026-11-27 public archive."""
from argparse import ArgumentParser
from io import BytesIO
from pathlib import Path
from time import perf_counter, process_time
from zipfile import ZipFile, ZIP_STORED, ZIP_DEFLATED
import json
import subprocess

import numpy as np

HERE=Path(__file__).resolve().parent
PREFIX="research/z_alg_optical_inertia/blind_uq_mimo_theory"
TAG="z-alg-blind-uq-mimo-20261127"
COMMIT="5f806c9ac2cd568d28ea23bea02f6ccd6ee32c93"
REPO="https://github.com/migodam/b1-probabilistic-maxwell-imaging"
STUDY="blind-uq-mimo-study.zip"
OPERATORS="blind-uq-mimo-operators.zip"
SCORING="blind-uq-mimo-scoring-truths.zip"
CASES=("n32_c01","n32_c02","n32_c03","n32_c04","n32_c05","n32_c06","n32_c07","n32_c08","n96_c01","n96_c02")
WITNESSES=tuple(x for x in CASES if x.startswith("n32") and x!="n32_c07")
SRC=(
    "requirements.txt","ports/model.py",
    "diagnostics/verified_wave.py","diagnostics/native_replay.py",
    "diagnostics/aperture_control.py","diagnostics/aperture_control_freeze.json",
    "diagnostics/aperture_joint.py","diagnostics/aperture_joint_freeze.json",
    "imaging/shared.py","imaging/FREEZE.json",
)
APERTURE_KEYS=("zero","patch_candidate","roi","patch","A_one_way_Y","A_one_way_r")
OPERATORS_USED=("n32_f400_a08","n32_f460_a08_c048","n96_f400_a08","n96_f460_a08_c048")

def main():
    p=ArgumentParser()
    p.add_argument("checkout",type=Path)
    p.add_argument("original_asset_dir",type=Path)
    args=p.parse_args()
    start=perf_counter(),process_time()
    git=lambda *a: subprocess.check_output(["git","-C",str(args.checkout),*a])
    assert git("rev-parse","HEAD").decode().strip()==COMMIT
    committed=git("ls-tree","--name-only","-r",COMMIT,"--",PREFIX).decode().splitlines()
    vendored=[x.removeprefix(PREFIX+"/") for x in committed if x.startswith(PREFIX+"/vendor_r8/") and x.endswith(".py")]
    inputs={f"ports/public/cases/{c}.npz":None for c in CASES}
    for c in WITNESSES:
        more=(("B_interleaved_Y","B_interleaved_r","B_interleaved_refit_zero","B_interleaved_refit_patch")
              if c in ("n32_c01","n32_c02") else ())
        inputs[f"results/aperture_control/{c}.npz"]=APERTURE_KEYS+more
    inputs["results/aperture_control/geometry.npz"]=("E0_both","S_both")
    inputs["results/aperture_joint/n32_c01.npz"]=("roi_zero_x","connected_patch_floor_2_x")
    inputs["results/aperture_joint/n32_c02.npz"]=("roi_zero_x",)
    for stem in ("stable_profiles/n32_c01/roi_01","gn_profiles/n96_c01/roi_02","gn_profiles/n96_c02/roi_02"):
        inputs[f"results/{stem}.npz"]=("material","roi")
    for stem in ("stable_profiles/n32_c01/roi_01_inner","gn_profiles/n96_c01/roi_02_inner"):
        inputs[f"results/{stem}.json"]=None
    for c in ("n32_c01","n32_c02","n96_c01","n96_c02"):
        inputs[f"results/extended/{c}.npz"]=("roi","patch","material")
        inputs[f"results/extended/{c}.json"]=None
    rows=[]
    def add(bundle,rel,data,source,member,keys=None):
        fmt=ZIP_STORED if rel.endswith(".npz") else ZIP_DEFLATED
        bundle.writestr(PREFIX+"/"+rel,data,compress_type=fmt)
        rows.append(dict(path=PREFIX+"/"+rel,source=source,member=member,
                         retained_npz_keys=list(keys) if keys else None,bytes=len(data),bundle=bundle.filename.rsplit("/",1)[-1]))
    def select(archive,rel,keys=None):
        name="blind_uq_mimo_theory/"+rel
        data=archive.read(name)
        if keys is not None:
            buf=BytesIO()
            with np.load(BytesIO(data),allow_pickle=False) as z:
                np.savez_compressed(buf,**{k:z[k] for k in keys})
            data=buf.getvalue()
        return name,data
    with (ZipFile(HERE/"prior-core.zip","w") as f,
          ZipFile(args.original_asset_dir/STUDY) as old,
          ZipFile(args.original_asset_dir/OPERATORS) as ops):
        for rel in (*SRC,*vendored):
            member=PREFIX+"/"+rel
            data=git("show",COMMIT+":"+member)
            alt_member,alt_data=select(old,rel)
            if alt_data!=data:
                raise ValueError(f"Git and published asset differ for {rel}")
            add(f,rel,data,"git:"+COMMIT,member)
            rows[-1]["alternate_original_asset"]=STUDY
            rows[-1]["alternate_original_member"]=alt_member
        geom=old.read("blind_uq_mimo_theory/results/aperture_control/geometry.npz")
        operator=ops.read("blind_uq_mimo_theory/ports/public/operators/n32_f400_a08.npz")
        with np.load(BytesIO(geom),allow_pickle=False) as a,np.load(BytesIO(operator),allow_pickle=False) as b:
            if not np.array_equal(a["G"],b["G"]):
                raise ValueError("Original geometry's dense Green differs from the released 400 MHz operator")
        for rel,keys in inputs.items():
            member,data=select(old,rel,keys)
            add(f,rel,data,STUDY,member,keys)
        for op in OPERATORS_USED:
            rel=f"ports/public/operators/{op}.npz"
            member,data=select(ops,rel)
            add(f,rel,data,OPERATORS,member)
    with ZipFile(HERE/"prior-scoring-opt-in.zip","w") as f, ZipFile(args.original_asset_dir/SCORING) as old:
        for c in ("n32_c02","n32_c07"):
            rel=f"ports/private/{c}.npz"
            member,data=select(old,rel,("material",))
            add(f,rel,data,SCORING,member,("material",))
    (HERE/"prior-manifest.json").write_text(json.dumps(dict(
        source_repository=REPO,source_tag=TAG,source_commit=COMMIT,
        original_assets={key:REPO+"/releases/download/"+TAG+"/"+key for key in (STUDY,OPERATORS,SCORING)},
        archive_paths="workspace-relative; historical study release member paths begin blind_uq_mimo_theory/",
        geometry_restore="For compact geometry.npz the sole missing Joint key G is copied from the exact archived n32_f400_a08 operator by restore_prior.py.",
        files=rows,build=dict(wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1])),
        indent=2,ensure_ascii=False)+"\n")


if __name__=="__main__":
    main()

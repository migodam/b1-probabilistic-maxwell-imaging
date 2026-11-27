"""Restore only the archived paths used by the new work, without overwriting files."""
from argparse import ArgumentParser
from io import BytesIO
from pathlib import Path
from time import perf_counter, process_time
from zipfile import ZipFile
from contextlib import ExitStack
import json

import numpy as np

HERE=Path(__file__).resolve().parent
GD="results/aperture_control/geometry.npz"
OP="ports/public/operators/n32_f400_a08.npz"
ROOT="research/z_alg_optical_inertia/blind_uq_mimo_theory"


def compatible(existing,data):
    if existing.read_bytes()==data:
        return True
    if existing.suffix!=".npz":
        return False
    with np.load(existing,allow_pickle=False) as old,np.load(BytesIO(data),allow_pickle=False) as new:
        return set(new.files)<=set(old.files) and all(
            old[key].dtype==new[key].dtype and np.array_equal(old[key],new[key]) for key in new.files)


def safe_install(base,rel,data):
    target=(base/rel).resolve()
    if not target.is_relative_to(base) or not target.is_relative_to(base/ROOT):
        raise ValueError(f"Unsafe archived path {rel}")
    if target.exists():
        if not compatible(target,data):
            raise FileExistsError(f"Existing path differs: {target}")
        return False
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open("xb") as f:
        f.write(data)
    return True


def restore(base,include_scoring=False,original=None):
    start=perf_counter(),process_time()
    base=Path(base).resolve()
    core=HERE/"prior-core.zip"
    scored=HERE/"prior-scoring-opt-in.zip"
    manifest=json.loads((HERE/"prior-manifest.json").read_text())
    restored,existing=0,0
    with ExitStack() as ctx:
        open_zip=lambda x:ctx.enter_context(ZipFile(x))
        mini=open_zip(core) if original is None else None
        full={}
        if original is not None:
            full={name:open_zip(original/name)
                  for name in manifest["original_assets"]
                  if not name.endswith("truths.zip") or include_scoring}
        scoring=open_zip(scored) if include_scoring and original is None else None
        for record in manifest["files"]:
            rel=record["path"]
            sealed=record["bundle"]==scored.name
            if sealed and not include_scoring:
                continue
            if original is None:
                data=(scoring if sealed else mini).read(rel)
            else:
                source=record.get("alternate_original_asset",record["source"])
                member=record.get("alternate_original_member",record["member"])
                data=full[source].read(member)
            if rel==ROOT+"/"+GD and original is None:
                with np.load(BytesIO(data),allow_pickle=False) as geom:
                    fields={k:geom[k] for k in geom.files}
                with np.load(BytesIO(mini.read(ROOT+"/"+OP)),allow_pickle=False) as op:
                    fields["G"]=op["G"]
                z=BytesIO()
                np.savez_compressed(z,**fields)
                data=z.getvalue()
            if safe_install(base,rel,data):
                restored+=1
            else:
                existing+=1
    return dict(new_files=restored,matching_files=existing,
                sealed_scoring=include_scoring,source_mode="original_assets" if original else "portable",
                wall_s=perf_counter()-start[0],cpu_s=process_time()-start[1])


if __name__=="__main__":
    p=ArgumentParser()
    p.add_argument("workspace",type=Path)
    p.add_argument("--include-sealed-scoring",action="store_true")
    p.add_argument("--original-asset-dir",type=Path)
    a=p.parse_args()
    print(json.dumps(restore(a.workspace,a.include_sealed_scoring,a.original_asset_dir),indent=2))

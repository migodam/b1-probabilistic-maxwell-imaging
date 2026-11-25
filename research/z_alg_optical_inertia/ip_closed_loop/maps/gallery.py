"""Prepare or reveal a label-free replay of the four-morphology library."""
import argparse
import json
import shutil

from common import FREEZE, HERE, save_json


SHAPES=("solid","ring","double","austria")
DIR=HERE/"public"/"gallery"
SEALED=HERE/"results"/"gallery_sealed"
MAPPING=HERE/"offline"/"gallery_reveal.json"


def prepare():
    DIR.mkdir(parents=True,exist_ok=True)
    roster=json.loads((HERE/"public"/"acquisition.json").read_text())["observations"]
    ids=sorted(k for k,v in roster.items() if v["side"]==32 and v["scene"]=="weak")
    # Frozen permutation for the replay. The opaque IDs carry no class.
    p=(7,2,13,8,0,9,1,12,4,14,10,5,15,11,3,6)
    aliases={f"sample-{i+1:02d}":ids[q] for i,q in enumerate(p)}
    for name,oid in aliases.items():
        shutil.copyfile(HERE/"public"/f"{oid}.npz",DIR/f"{name}.npz")
    save_json(DIR/"acquisition.json",
              {"observations":{name:{"side":32} for name in aliases}})
    save_json(MAPPING,aliases)


def reveal():
    aliases=json.loads(MAPPING.read_text())
    b=FREEZE["model"]["background_upper"]
    seen=[]
    for shape in SHAPES:
        record=json.loads((SEALED/f"sector_32_{shape}_library_all.json").read_text())
        assert record["opaque_public_ids"]
        assert set(record["summary"])==set(aliases)
        for alias,oid in aliases.items():
            s=SEALED/f"{shape}__{alias}__b{b:g}__library_all_sector.npz"
            d=HERE/"results"/f"{shape}__{oid}__b{b:g}__library_all_sector.npz"
            assert s.exists()
            seen.append((s,d))
        record["summary"]={aliases[k]:v for k,v in record["summary"].items()}
        record["reveal_source"]=str(MAPPING.relative_to(HERE))
        save_json(HERE/"results"/f"sector_32_{shape}_library_all.json",record)
    # Complete and validate the four sets before writing any labelled arrays.
    for source,target in seen:
        shutil.copyfile(source,target)


if __name__=="__main__":
    a=argparse.ArgumentParser()
    a.add_argument("stage",choices=("prepare","reveal"))
    if a.parse_args().stage=="prepare":
        prepare()
    else:
        reveal()

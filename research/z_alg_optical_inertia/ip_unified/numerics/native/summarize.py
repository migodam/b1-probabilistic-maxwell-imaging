"""Collect source-spectrum upper/lower controls and measured run cost."""
import json
from pathlib import Path

from flint import arb, ctx
import numpy as np

HERE=Path(__file__).resolve().parent


def read(name):
    return json.loads((HERE/f"{name}.json").read_text())


def main():
    result={}
    ctx.prec=112
    for case in ("u075","u100_best","u100_stable","u125"):
        run=read(case)
        nu=run["verified_state"]["exact_nu"]
        rec=[]
        for row in run["rows"]:
            src=row["outward_baseline"]
            radius=arb(row["padded_radius"])
            hi,lo=[],[]
            for c in range(nu+1):
                for key, array, nextdir in (("source_eigenvalues_upper",hi,np.inf),
                                              ("source_eigenvalues_lower",lo,-np.inf)):
                    tail=sum((arb(max(0.,v)) for v in src[key][nu-c:]),arb(0))
                    bound=(tail.sqrt()/radius)
                    value=bound.upper() if nextdir==np.inf else bound.lower()
                    array.append(max(0.,float(np.nextafter(float(value),nextdir))))
            rec.append(dict(id=row["id"],verified_old_ratio=hi[0],
                optimistic_center_free_ratio_upper_by_c=hi,
                optimistic_center_free_ratio_lower_by_c=lo,
                minimum_c_needed_center_free=next((i for i,v in enumerate(hi) if v>1),None)))
        result[case]=dict(negative_count=nu,rows=rec)
    (HERE/"OUTWARD_OPTIMISTIC_CONTROL.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    main()

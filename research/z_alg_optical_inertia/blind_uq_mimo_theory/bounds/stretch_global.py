"""Finer whole-field cap scans after the frozen 0.25 report; no material labels."""
import argparse
import json
import time

import numpy as np
from scipy.linalg import norm
from threadpoolctl import threadpool_limits

from mimo import Transfer
from sector import HERE, PUBLIC, Green
from verify import ExactUniform, ctx

CAPS=(0.25,0.3125,0.375,0.4375,0.5,0.625)
RATIOS=(0.,0.25,1.,4.,16.,64.,256.)
SCALES=(0.012,0.06,0.3,1.5)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--side",type=int,choices=(32,96))
    args=parser.parse_args()
    ctx.prec,ctx.threads=112,2
    start=time.perf_counter()
    grouped={}
    for obs in sorted((PUBLIC/"cases").glob("n*_c*.npz")):
        if obs.stem.endswith("_mimo") or (args.side and not obs.stem.startswith(f"n{args.side}_")):
            continue
        data=np.load(obs)
        grouped.setdefault(str(data["operator_id"]),[]).append((obs.stem,data["Y"],float(data["radius"])))
    output=HERE/"results"/"stretch_global.json"
    old=json.loads(output.read_text()) if args.side and output.exists() else {"rows":[],"wall_seconds":0.}
    out=[r for r in old["rows"] if args.side and not r["case"].startswith(f"n{args.side}_")]
    skipped=[]
    for label, rows in grouped.items():
        geom=PUBLIC/"operators"/f"{label}.npz"
        files=np.load(geom)
        green=Green.read(geom,int(label.split("_")[0][1:]))
        s,b=files["S"],files["E0"]
        scale=norm(s,2)**2
        candidates={}
        for cap in CAPS:
            comp=[]
            for ratio in RATIOS:
                kept=np.count_nonzero(1.-cap*green.spectrum.real+ratio*green.spectrum.imag<.5)
                if kept>850:
                    skipped.append(dict(operator=label,cap=cap,ratio=ratio,retained_fourier_modes=int(kept)))
                    continue
                comp.extend(Transfer(green,s,b,cap,ratio,w*scale) for w in SCALES)
            candidates[cap]=comp
        for case,y,r in rows:
            curve=[]
            for cap,comp in candidates.items():
                dist,plan=max(((t.distance(y),t) for t in comp),key=lambda a:a[0])
                curve.append(dict(cap=cap,min_surrogate_nu=min(t.nu for t in comp),
                                  proposal=dict(ratio=plan.ratio,weight=plan.w,
                                                distance=dist,margin=dist-r)))
            proofs=[]
            for item in reversed(curve):
                m=item["proposal"]
                if m["margin"]<=0:
                    continue
                cert=ExactUniform(green,s,b,y).proof(item["cap"],m["ratio"],m["weight"],r)
                proofs.append(cert)
                if cert["proved"]:
                    break
            out.append(dict(case=case,operator=label,noise_radius=r,curve=curve,proofs=proofs))
            print(case,curve[-2]["min_surrogate_nu"],[v["upper"] for v in proofs if v["proved"]],flush=True)
        output.write_text(json.dumps(dict(rows=out,skipped=skipped,wall_seconds=time.perf_counter()-start,
                                           previously_completed_wall_seconds=old["wall_seconds"]),
                                                                          indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    with threadpool_limits(1):
        main()

"""Geometry-only DPCS noise amplification. Uses public acquisition and no material."""
import numpy as np
import run_study as R
from score_study import write_csv
def main():
    cases=[c for c in R.cases() if c["cohort"]=="confirmation"]
    rows=[];used=set()
    for c in cases:
        if c["array"] in used:continue
        used.add(c["array"]);p=R.load_public(c);S=p["S"]
        sig=np.linalg.svd(S,compute_uv=False)
        keep=sig[sig**2>1e-12*sig[0]**2][:15]
        for cap in (False,True):
            kk=keep[keep>=.1*sig[0]] if cap else keep
            rows.append(dict(array=c["array"],cap=cap,rank=len(kk),weakest_ratio=kk[-1]/sig[0],
                             per_port_rms_dimensionless_noise_gain=np.sqrt(np.sum((sig[0]/kk)**2)/len(sig)),
                             inverse_trace=np.sum(1/kk**2),receiver_rows=len(sig)))
    write_csv(R.ROOT/"results"/"geometry_noise.csv",rows)
if __name__=="__main__":main()

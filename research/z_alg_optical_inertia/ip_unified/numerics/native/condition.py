"""Recompute row probes after orthogonalizing the stored complex Ritz span."""
import argparse
from time import perf_counter,process_time

import numpy as np
from threadpoolctl import threadpool_limits

from replay import (HERE, CASES, configure, frame_weight_upper, global_screen, normalized_modes,
                    projected_gram, save_json, select_rows, times)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case",choices=CASES)
    args=parser.parse_args()
    start=perf_counter(),process_time()
    geo,native,_,_=configure(args.case)
    base=HERE/args.case
    data=np.load(base.with_suffix(".npz"),allow_pickle=False)
    raw=data["negative_vec"]
    raw_min=float(np.linalg.eigvalsh(raw.conj().T@raw)[0])
    eig,vec=normalized_modes(native,raw)
    defect=float(np.linalg.norm(native.numeric_action(vec)-vec*eig))
    ortho=float(np.linalg.norm(vec.conj().T@vec-np.eye(len(eig))))
    n=vec/np.sqrt(-eig)
    gn=geo.numerical_green(n)
    wc=frame_weight_upper(native,vec,eig)
    probe,xx,ev,basis=global_screen(native,n,gn,wc)
    weak=basis[:,probe["contraction_dimension"]:]
    if not weak.shape[1]:
        weak=basis[:,-1:]
    pix=select_rows(n,gn,native.w.cap,weak)
    print(args.case,raw_min,defect,ortho,"new global proposal",probe,flush=True)
    _,trial,gram,cost=projected_gram(native,vec,pix)
    out=HERE/f"{args.case}_conditioned"
    np.savez_compressed(out.with_suffix(".npz"),negative_vec=vec,negative_eig=eig,negative_gn=gn,
                        pixel_index=pix,endpoint_trial=trial,endpoint_gram=gram,homogeneous_x=xx,
                        homogeneous_spectra=ev,homogeneous_subspace=basis,frame_weight_upper_numeric=wc)
    save_json(out.with_suffix(".json"),dict(case=args.case,raw_ritz_gram_min=raw_min,
        checked_ritz_defect=defect,checked_ritz_orthogonality_defect=ortho,global_probe=probe,
        projected_gram_columns=len(cost),projected_gram_solves=cost,
        operations=dict(geo.cost),**times(*start)))
    print(args.case,"done",times(*start),flush=True)


if __name__=="__main__":
    with threadpool_limits(1):
        main()

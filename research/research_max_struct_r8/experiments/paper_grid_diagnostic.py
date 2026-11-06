"""Paid offline physical adequacy on already-sealed public eps3.5 benchmark.

Fixed analytic published object, no fit/tuning. Compare independent square grids
to150 and200 using the identical equal-area kernel/acquisition and strict solve.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time
import numpy as np
from research.research_max_struct_r8.src.paper_physics import acquisition,austria,fft_forward
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,write

def run():
    case='austria_eps3p5_400MHz';assert (OUT/'original64'/case/'cie_single/PREDICTION_SEAL.json').exists()
    observations=np.load(OUT/case/'PUBLIC_DATA.npz')['Y'];states={};rows=[]
    for n in [32,48,64,96,150,200]:
        p=acquisition(n,400e6,dense=False);m=austria(p['centers'],3.5);y,cost=fft_forward(p,m)
        states[n]=y;rows.append(dict(side=n,physical_area_m2=float(np.count_nonzero(m)*4/n**2),observed_noisy_residual=float(np.linalg.norm(y-observations)/np.linalg.norm(observations)),cost=cost))
    for row in rows:
        y=states[row['side']];row['clean_relative_difference_to150']=float(np.linalg.norm(y-states[150])/np.linalg.norm(states[150]));row['clean_relative_difference_to200']=float(np.linalg.norm(y-states[200])/np.linalg.norm(states[200]))
    write(OUT/'OFFLINE_GRID_ADEQUACY.json',dict(case=case,analytic_area_m2=float(np.pi*.35),purpose='offline physical discretization after original64 prediction seal; no online fitting',rows=rows));print(rows,flush=True)

if __name__=='__main__':run()

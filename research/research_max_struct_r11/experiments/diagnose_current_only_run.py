"""Process-local dispatch for the isolated GreenCurrentOrdinaryMaterial DEV control."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import sys,json
ROOT=Path(__file__).absolute().parents[1];sys.path.insert(0,str(ROOT.parents[1]))
from research.research_max_struct_r11.experiments import twod_run
from research.research_max_struct_r8.src import benchmark_solver
from research.research_max_struct_r11.src.diagnostic_green_current_only import run_cie_fft

def main():
    # This wrapper never mutates source files or dispatch in another process.
    benchmark_solver.run_cie_fft=run_cie_fft
    return twod_run.main()
if __name__=='__main__':raise SystemExit(main())

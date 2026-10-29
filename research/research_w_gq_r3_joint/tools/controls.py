"""Single entry point for frozen ordinary controls and isolated known-wall diagnostic."""
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[_thread_key]='1'
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
def main():
    if '--method' not in sys.argv:raise ValueError('--method PRIOR_ONLY|WALL_FIRST_FIXED|KNOWN_WALL_ORACLE required')
    i=sys.argv.index('--method');method=sys.argv[i+1]
    if method=='KNOWN_WALL_ORACLE':
        del sys.argv[i:i+2]
        from oracle_control import main as run
    else:
        from ordinary_controls import main as run
    run()
if __name__=='__main__':main()

"""Record complete non-physics training/analysis subprocess time in new budget."""
from pathlib import Path
import argparse,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import Book,write
def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',required=True);p.add_argument('--job',required=True);p.add_argument('--max-wall',type=float,default=14400);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    cmd=a.command[1:] if a.command[:1]==['--'] else a.command
    b=Book(a.stage,a.job,max_wall=a.max_wall,analysis_only=True);status='COMPLETE';error=None
    write(b.jobdir/'EXECUTION.json',dict(command=cmd,no_physics_in_child=True,CPU_wall_includes_subprocess=True))
    try:
        code=subprocess.call(cmd,cwd=ROOT)
        if code:raise RuntimeError('child exit '+str(code))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)
if __name__=='__main__':main()

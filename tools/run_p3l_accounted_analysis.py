"""Bill ordinary analysis wall while a separate approved process reads GT."""
from pathlib import Path
import sys,argparse,subprocess
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import Book

def main():
    a=argparse.ArgumentParser();a.add_argument('analysis',choices=['gaussian_images','joint_gaussian_images','locked_imaging']);args=a.parse_args()
    scripts={'gaussian_images':('B','evaluate_gaussian_reference_images','evaluate_p3l_gaussian_images.py'),
             'joint_gaussian_images':('B','evaluate_joint_gaussian_images','evaluate_p3l_joint_gaussian.py'),
             'locked_imaging':('D','evaluate_locked_imaging','evaluate_p3l_locked_imaging.py')}
    stage,job,script=scripts[args.analysis];b=Book(stage,job,max_wall=1800,analysis_only=True);status='COMPLETE';error=None
    try:subprocess.run([sys.executable,str(ROOT/'tools'/script)],check=True)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()

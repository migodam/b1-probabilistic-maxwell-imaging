"""Serialize final exposed DEV projection validation after C pilot closes."""
from pathlib import Path
import json
import sys
import subprocess
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write


def main():
    folder=OUT/'schedules/legal_reference_projection';folder.mkdir(parents=True,exist_ok=True)
    dependency=OUT/'schedules/surrogate_force_hmc/CLOSED.json'
    if (folder/'PROTOCOL.json').exists():raise FileExistsError('sequence already started, no implicit replay')
    write(folder/'PROTOCOL.json',dict(dependency=str(dependency.relative_to(OUT)),
        parents=[f'DEV_{i:03d}' for i in range(32)],maximum_RHS=256,
        locked_methods_changed=False,GT_online=False,
        purpose='EXPOSED_ONLY_DIAGNOSTIC_NOT_FULL_POSTERIOR_BAYES_ACTION'))
    while not dependency.exists():time.sleep(30)
    if json.loads(dependency.read_text())['status']!='COMPLETE':raise RuntimeError('C pilot interrupted; preserve all receipts')
    with (folder/'validation.log').open('x') as log:
        subprocess.run([sys.executable,str(ROOT/'tools/validate_p3l_legal_reference_projection.py')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    book=Book('B','score_legal_reference_projection_all32',max_wall=600,analysis_only=True)
    status,error='ERROR',None
    try:
        with (folder/'scoring.log').open('x') as log:
            subprocess.run([sys.executable,str(ROOT/'tools/evaluate_p3l_legal_reference_projection.py')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:book.close(status,error)
    write(folder/'CLOSED.json',dict(status='COMPLETE',new_RHS_upper=256))


if __name__=='__main__':main()

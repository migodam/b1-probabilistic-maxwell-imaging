"""Continue unpaid projection validation after the retained parser recovery."""
from pathlib import Path
import sys,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,write


def main():
    folder=OUT/'schedules/legal_reference_projection_fix1';folder.mkdir(parents=True,exist_ok=True)
    dependency=OUT/'schedules/surrogate_force_hmc_fix1/CLOSED.json'
    if (folder/'PROTOCOL.json').exists():raise FileExistsError('no replay of continuation')
    if (OUT/'B/jobs/validate_legal_reference_projection_all32').exists():raise RuntimeError('physical validation already started')
    write(folder/'PROTOCOL.json',dict(dependency=str(dependency.relative_to(OUT)),
        retained_previous='schedules/legal_reference_projection',maximum_RHS=256,
        unchanged_paid_job_IDs=True,ledger_reset=False,locked_methods_changed=False,
        study='EXPOSED_DEV_ONLY_NOT_FULL_POSTERIOR_BAYES_ACTION'))
    while not dependency.exists():time.sleep(30)
    if json.loads(dependency.read_text())['status']!='COMPLETE':raise RuntimeError('preserve C failure before validation')
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
    write(folder/'CLOSED.json',dict(status='COMPLETE',maximum_new_RHS=256))


if __name__=='__main__':main()

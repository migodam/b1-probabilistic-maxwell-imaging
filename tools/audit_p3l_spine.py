"""Audit the actually completed four-chain full-target spine pilot."""
from pathlib import Path
import json
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target
from p3l.reference_audit import audit_reference

def main():
    paths=[OUT/f'B/jobs/spine_hmc_DEV030_pilot_c{c}/CHAIN.npz' for c in range(4)]
    for p in paths:
        if json.loads((p.parent/'CLOSED.json').read_text())['status']!='COMPLETE':
            raise RuntimeError('all four completed pilot chains required')
    b=Book('B','audit_spine_DEV030_pilot',max_wall=600)
    status,error='COMPLETE',None
    try:
        stats=audit_reference(paths,Target(b).prior,512,OUT/'B/references/DEV_030/spine_hmc_pilot')
        print(json.dumps(dict(status=stats['status'],MCSE_RMS=stats['full_grid_mean_rms_mcse'],summaries=stats['summaries'])))
    except BaseException as exc:
        status,error='ERROR',repr(exc)
        raise
    finally:b.close(status,error)

if __name__=='__main__':main()

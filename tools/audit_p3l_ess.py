"""Full native Gaussian64/uniform2 blocked ESS; completed four-chain audit."""
from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target
from p3l.reference_audit import audit_reference

def main():
    paths=[OUT/f'A/jobs/reference_DEV_030_ess_mode_c{c}/CHAIN.npz' for c in range(4)]
    for p in paths:
        if json.loads((p.parent/'CLOSED.json').read_text())['status']!='COMPLETE':raise RuntimeError('all four closures required')
    b=Book('A','audit_DEV_030_ess_mode',max_wall=600);status,error='COMPLETE',None
    try:
        stats=audit_reference(paths,Target(b).prior,1024,OUT/'A/references/DEV_030/ess_mode')
        print(json.dumps(dict(status=stats['status'],MCSE_RMS=stats['full_grid_mean_rms_mcse'],summaries=stats['summaries'])))
    except BaseException as exc:status,error='ERROR',repr(exc);raise
    finally:b.close(status,error)

if __name__=='__main__':main()

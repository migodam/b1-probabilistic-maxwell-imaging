from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write
from p3l.reference_audit import audit_reference

def main():
    parents=sys.argv[1:] or json.loads((OUT/'A/DEEP_REFERENCE_PARENTS.json').read_text())['parents']
    for parent in parents:
        paths=[OUT/f'A/jobs/reference_{parent}_mala_c{c}/CHAIN.npz' for c in range(4)]
        if not all(p.exists() and (p.parent/'CLOSED.json').exists() for p in paths):continue
        dest=OUT/f'A/references/{parent}/mala'
        if (dest/'REFERENCE_STATS.json').exists():continue
        b=Book('A',f'audit_{parent}_mala',max_wall=600);status='COMPLETE';error=None
        try:
            stats=audit_reference(paths,Target(b).prior,1024,dest)
            print(json.dumps(dict(parent=parent,status=stats['status'],MCSE_RMS=stats['full_grid_mean_rms_mcse'],summaries=stats['summaries'])),flush=True)
        except BaseException as e:status='ERROR';error=repr(e);raise
        finally:b.close(status,error)

if __name__=='__main__':main()

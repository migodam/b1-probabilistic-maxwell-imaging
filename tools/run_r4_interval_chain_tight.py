"""Revision1: same fixed-kernel audit, tighter deterministic inverse norm bound.

Original closed v1 is never edited/replayed. This uses prospective development
parents32..39. All probabilities, beta, refinement cap and persistence stay the
same; only the valid operator-norm upper bound is sharpened.
"""
from pathlib import Path
import argparse,json,subprocess,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'tools')]
import run_r4_interval_chain as engine
from p3r4.decision_refresh_tight import DecisionRefreshTight

engine.EXPERIMENT='interval_chain_tight_r1'
engine.PARENTS=tuple(f'R4DEV_{k:05d}' for k in range(32,40))
engine.SEED=20261025101
engine.RESULT=engine.OUT/'results/II'/engine.EXPERIMENT
engine.DecisionRefresh=DecisionRefreshTight

def main():
    ap=argparse.ArgumentParser();ap.add_argument('arm',nargs='?',default='pair',choices=['pair','interval','full','compare'])
    ap.add_argument('--release',type=Path,required=True);a=ap.parse_args()
    if a.arm=='compare':return engine.compare()
    if a.arm!='pair':return engine.run(a.arm,a.release)
    release=json.loads(a.release.read_text())
    if not release.get('released') or release['experiment']!=engine.EXPERIMENT:raise ValueError('explicit revision release required')
    marker=engine.OUT/'results/II/interval_chain_tight_pair'
    if marker.exists() or engine.RESULT.exists():raise FileExistsError('no paid replay')
    marker.mkdir(parents=True);records=[];start=time.monotonic()
    for arm in ['interval','full','compare']:
        p=subprocess.run([sys.executable,'-u',str(Path(__file__).resolve()),arm,'--release',str(a.release.resolve())],cwd=ROOT)
        records.append(dict(arm=arm,returncode=p.returncode))
        if p.returncode:
            (marker/'FAILED.json').write_text(json.dumps(dict(status='FAILED',records=records),indent=2)+'\n')
            raise SystemExit(p.returncode)
    comparison=json.loads((engine.RESULT/'COMPARISON.json').read_text())
    (marker/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',records=records,
        comparison=str(engine.RESULT/'COMPARISON.json'),comparison_status=comparison['status'],
        coordinator_wall_not_double_counted=time.monotonic()-start),indent=2)+'\n')
if __name__=='__main__':main()

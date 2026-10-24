"""Wait for the current acquisition block to finish, then run one paired revision."""
from pathlib import Path
import sys,time,subprocess,json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'

def main():
    release=OUT/'configs/INTERVAL_CHAIN_TIGHT_RELEASE.json'
    if release.exists():raise FileExistsError('revision already registered')
    release.write_text(json.dumps(dict(released=True,campaign_id='P3LR4-20261025',experiment='interval_chain_tight_r1',
        hypothesis='min(Frobenius,sqrt(one*infinity)) inverse norm upper bound reduces conservative fallback and end-to-end cost',
        old_result='v1 identical paths but355.79s interval versus229.78s Full, not a speedup',
        change='bound only; no new learned error, target, noise, proposal or acceptance threshold',
        parents=[f'R4DEV_{k:05d}' for k in range(32,40)],steps=1024,beta=.002,seed=20261025101,
        per_arm_wall_cap=1200,combined_RHS_cap=200000,
        success='all common-random-number paths and potential enclosures agree; at least20percent end-to-end wall decrease vs new matched Full',
        floating_not_directed_interval=True,posterior_mixing_NOT_QUALIFIED=True),indent=2)+'\n')
    while not (OUT/'results/CONTEXT_EXPANSION_COMPLETE.json').exists():time.sleep(10)
    p=subprocess.run([sys.executable,'-u','tools/run_r4_interval_chain_tight.py','--release',str(release)],cwd=ROOT)
    raise SystemExit(p.returncode)
if __name__=='__main__':main()

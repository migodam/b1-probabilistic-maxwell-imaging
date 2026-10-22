"""Finite explicit sequence; no daemon, no extra authorization or retries."""
from pathlib import Path
import sys,subprocess,json,time,fcntl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,write
P=sys.executable
def run(name,args):
    log=OUT/f'{name}.log';status=OUT/'CHECKPOINT.json'
    write(status,dict(current=name,command=args,status='RUNNING',utc_unix=time.time()))
    with log.open('x') as f:r=subprocess.run([P,*args],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    if r.returncode:
        write(status,dict(current=name,status='STOPPED_ERROR',exit_code=r.returncode,log=str(log)));raise SystemExit(r.returncode)
    print(name+' COMPLETE',flush=True)
def main():
    # Waiting for an already authorized paid worker does not create a second
    # solver; each child still obtains the same exclusive Book lock itself.
    lock=(OUT/'paid_worker.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX);fcntl.flock(lock,fcntl.LOCK_UN)
    old='../b1-p3l-local-posterior/research/p3l_local_posterior_imaging/A/anchors/DEV_030.npz'
    run('B_full_r0',['tools/run_p3r2_posterior.py','--mode','full','--job','B_full_nonlinear_DEV030_r0','--parent','DEV_030','--anchor',old])
    train=[f'DEV_{i:03d}' for i in range(8)]
    run('C_train_r0',['tools/run_p3r2_geometry.py','--job','C_train_r0','--parents',*train,'--arms','frozen_fixed','dynamic_fixed','dynamic_curvature','dynamic_bandit','ordinary_random','--split','train','--steps','8','--max-wall','1200'])
    run('C_policy_bandit',['tools/train_p3r2_policy.py','--job','C_policy_bandit','--trajectory-jobs',str(OUT/'C/jobs/C_train_r0')])
    run('C_eval_r0',['tools/run_p3r2_geometry.py','--job','C_exposed_eval_r0','--parents',*[f'DEV_{i:03d}' for i in range(24,32)],'--arms','frozen_fixed','dynamic_fixed','dynamic_curvature','dynamic_bandit','ordinary_random','--split','exposed_eval','--policy',str(OUT/'C/policy/FROZEN_POLICY.json'),'--steps','8','--max-wall','1200'])
    run('A_anchor_train',['tools/run_p3r2_anchors.py','train','--count','512','--native','48','--probit','48'])
    run('A_anchor_validation',['tools/run_p3r2_anchors.py','validation','--count','128','--native','48','--probit','48'])
    run('A_pack_train',['tools/pack_p3r2_anchors.py','train','--count','512'])
    run('A_pack_validation',['tools/pack_p3r2_anchors.py','validation','--count','128'])
    write(OUT/'CHECKPOINT.json',dict(status='INITIAL_PHYSICS_SEQUENCE_COMPLETE'))
if __name__=='__main__':main()

"""No neural contention during optimized interval/Full timing, then DINO controls."""
from pathlib import Path
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def call(args):
    p=subprocess.run([sys.executable,*args],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)
def main():
    release=OUT/'configs/INTERVAL_CHAIN_FAST_RELEASE.json'
    if release.exists():raise FileExistsError('fast revision release already exists')
    release.write_text(json.dumps(dict(released=True,experiment='interval_chain_fast_v1',campaign_id='P3LR4-20261025',
        classification='common lossless I/O plus exact row-norm/cache implementation optimization; not a new physics theorem',
        unchanged_target_and_kernel=True,common_codec='stored',endpoint_storage='one canonical latest numerical state per same material; exact repeat-state references, append-only refinement proof scalars',
        preserve_all_previous_pairs=True,new_parents=list(range(40,48)),steps=1024,beta=.002,seed=20261025131,
        combined_RHS_cap=200000,wall_perarm=1200,success='zeroflags/paths/allproofmisses and>=20percent wholeprocesswall decrease vs commonoptimizedFull'),indent=2)+'\n')
    for marker in ['ANCHOR_NEURAL_COMPLETE.json','RISK_FUSION_COMPLETE.json','REP_DIAGNOSTICS_COMPLETE.json']:
        while not (OUT/'results'/marker).exists():time.sleep(10)
    call(['tools/run_r4_interval_chain_fast.py','--release',str(release)])
    plan=[dict(mode=m,seed=s,updates=1000,coordinates=66) for m in ['value','jvp'] for s in [2026102501,2026102502,2026102503]]
    config=OUT/'configs/DINO_CONTROL_ACTUAL_RUNS.json'
    config.write_text(json.dumps(dict(released=True,plan=plan,teacher_owner_RHS=459712,new_RHS=0,
        target='saved exact nonlinear whitened finite response remainder; derivative-informed control uses two JVPs',
        independent_parent_unit=True,all_operator_test_development=True,not_posterior_or_certification=True),indent=2)+'\n')
    for row in plan:
        job=f"forward_{row['mode']}_66_s{row['seed']}"
        call(['tools/train_forward_jvp.py','--train','--mode',row['mode'],'--coordinates','66','--seed',str(row['seed']),
            '--updates','1000','--job',job,'--max-wall','900'])
        call(['tools/train_forward_jvp.py','--evaluate','--checkpoint',str(OUT/'results/II/forward_jvp_training'/job/'BEST.pt'),
            '--job',job+'_evaluation','--max-wall','300'])
    (OUT/'results/FAST_AND_DINO_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=len(plan)),indent=2)+'\n')
if __name__=='__main__':main()

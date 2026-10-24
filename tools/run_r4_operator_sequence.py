"""Explicit release of three small empirical operator-action fits, no new physics."""
from pathlib import Path
import sys,json,time,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def main():
    release=OUT/'configs/OPERATOR_TRAIN_RELEASE.json'
    if release.exists():raise FileExistsError('operator fits already registered')
    seeds=[2026102501,2026102502,2026102503]
    release.write_text(json.dumps(dict(released=True,seeds=seeds,dimension=66,
        hypothesis='predict material GN-action error to select among paid-state cheap Jacobian updates',
        TRAIN=16,validation=8,operator_test=8,all_exposed_development=True,
        teacher_cost_RHS=459712,new_RHS=0,wall_perfit=900,scope='empirical ranker, no certified or posterior claim'),indent=2)+'\n')
    while not (OUT/'results/EXTENDED_NEURAL_COMPLETE.json').exists():time.sleep(10)
    for seed in seeds:
        job=f'operator_ranker66_s{seed}'
        p=subprocess.run([sys.executable,'tools/train_operator.py','--train','--coordinates','66',
            '--seed',str(seed),'--job',job,'--max-wall','900'],cwd=ROOT)
        if p.returncode:raise SystemExit(p.returncode)
        checkpoint=OUT/'results/II/operator_training'/job/'BEST.pt'
        p=subprocess.run([sys.executable,'tools/train_operator.py','--evaluate-operator-test',
            '--checkpoint',str(checkpoint),'--job',job+'_evaluation','--max-wall','300'],cwd=ROOT)
        if p.returncode:raise SystemExit(p.returncode)
    (OUT/'results/OPERATOR_TRAIN_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',fits=3),indent=2)+'\n')
if __name__=='__main__':main()

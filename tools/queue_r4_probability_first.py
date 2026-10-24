"""Predeclared fixed-index Full Maxwell checks for the first full66 generators."""
from pathlib import Path
import json,sys,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'

def call(args):
    p=subprocess.run([sys.executable,*args],cwd=ROOT)
    if p.returncode:raise SystemExit(p.returncode)

def main():
    folder=OUT/'results/probability_physics_first'
    if folder.exists():raise FileExistsError('do not replay paid probability checks')
    folder.mkdir(parents=True);generated=OUT/'results/development_density_first/online'
    parents=[f'R4DEV_{k:05d}' for k in range(128)]
    groups=[]
    for name,density in [('gaussian_s1','explicit'),('spline_s1','explicit'),('triangular_s1','explicit'),('cfm_s1__nfe32','none')]:
        groups.append(dict(id=name,directory=str(generated/name),online_gate=str(generated/'ONLINE_COMPLETE.json'),
            observations=str(OUT/'inputs/online'),parents=parents,extended_parents=parents[:16],density=density))
    spec=dict(generated=groups,likelihood_samples=64,importance_weights=True,
        cohort='NEW_R4_DEVELOPMENT_ALREADY_EXPOSED',fixed_sample=0,no_residual_selection=True,
        hypothesis='full66 likelihood residual law and data-dependent ranks; full66 NLL alone does not establish Maxwell posterior',
        required_Maxwell_RHS=36352,new_posterior_teacher=False)
    specpath=folder/'SPEC.json';frozen=folder/'FROZEN.json'
    specpath.write_text(json.dumps(spec,indent=2)+'\n')
    while not (generated/'ONLINE_COMPLETE.json').exists():time.sleep(10)
    call(['tools/run_r4_job.py','--stage','I','--job','freeze_probability_physics_first','--max-wall','60','--',
          sys.executable,'tools/run_r4_probability_audit.py','freeze','--spec',str(specpath),'--out',str(frozen)])
    # The physical acquisition and matched timing experiment have a fixed queue.
    while not (OUT/'results/II/interval_chain_tight_pair/COMPLETE.json').exists():
        if (OUT/'results/II/interval_chain_tight_pair/FAILED.json').exists():raise RuntimeError('reconcile earlier physical failure before more queries')
        time.sleep(10)
    call(['tools/run_r4_probability_audit.py','run','--freeze',str(frozen),'--out',str(folder/'physics'),
          '--stage','I','--job','probability_physics_first','--max-rhs','36352','--max-wall','1200'])
    call(['tools/run_r4_job.py','--stage','I','--job','score_probability_physics_first','--max-wall','300','--',
          sys.executable,'tools/run_r4_probability_audit.py','score','--freeze',str(frozen),'--online',str(folder/'physics'),
          '--truth',str(OUT/'inputs/offline_only'),'--out',str(folder/'scores')])
    (folder/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',no_GT_in_physics=True,all_fixed_indices=True),indent=2)+'\n')
if __name__=='__main__':main()

"""Evaluate the final phase-head initialization hypothesis, without new physics."""
from pathlib import Path
import sys,json,subprocess,time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'research/p3l_r4_conditional_transport'
def main():
    folder=OUT/'results/development_phase_revision2'
    if folder.exists():raise FileExistsError('no repeated phase evaluation')
    while not (OUT/'results/ANCHORBIN_REVISION2_COMPLETE.json').exists():time.sleep(10)
    folder.mkdir(parents=True);methods=[]
    for family,prefix,loader in [('direct_cnn','risk_local_cnn','r4'),('anchorbin8','risk_anchorbin_phase8','anchorbin'),('anchorbin16','risk_anchorbin_phase16','anchorbin')]:
        for seed in [2026102501,2026102502,2026102503]:
            methods.append(dict(id=f'{family}_s{seed}',family=family,loader=loader,seed=seed,checkpoint=str(OUT/'models/risk'/f'{prefix}_s{seed}'/'BEST.pt')))
    spec=dict(parents=[f'R4DEV_{i:05d}' for i in range(128)],methods=methods,fallback='anchor',cohort='DEVELOPMENT_EXPOSED',
        comparisons=[dict(id=f'{family}_vs_reoptimized_R2class',candidate_family=family,reference_family='direct_cnn',action='clipped') for family in ['anchorbin8','anchorbin16']],
        scope='last phase-family mechanism revision; do not tune eps or expand family after this result')
    path=folder/'METHODS.json';frozen=folder/'FROZEN.json';path.write_text(json.dumps(spec,indent=2)+'\n')
    commands=[['freeze','--methods',str(path),'--out',str(frozen)],
        ['infer','--freeze',str(frozen),'--contexts',str(OUT/'results/contexts/dev'),'--out',str(folder/'online')],
        ['score','--freeze',str(frozen),'--online',str(folder/'online'),'--truth',str(OUT/'inputs/offline_only'),'--out',str(folder/'scores')]]
    for args in commands:
        p=subprocess.run([sys.executable,'tools/evaluate_r4_risk.py',*args],cwd=ROOT)
        if p.returncode:raise SystemExit(p.returncode)
    (folder/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',phase_family_expansion='STOP',new_physics_RHS=0),indent=2)+'\n')
if __name__=='__main__':main()

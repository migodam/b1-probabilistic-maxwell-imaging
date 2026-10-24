"""Fixed exposed-development test of paid-anchor full66 models; no physics.

Inference and offline scoring use distinct child processes. This does not
select locked methods, retrain, sample from a corrector or read locked truth.
"""
from pathlib import Path
import json, subprocess, sys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r4_conditional_transport'

def main():
    folder=OUT/'results/development_anchor_density_first'
    if folder.exists():raise FileExistsError('preserve existing evaluation')
    export=OUT/'inputs/anchor_conditions/dev'
    complete=json.loads((export/'EXPORT_COMPLETE.json').read_text())
    if complete['status']!='COMPLETE' or complete['truth_read'] or complete['parents']!=128:
        raise ValueError('proven observation-only anchor export required')
    folder.mkdir()
    methods=[]
    for kind in ['spline','diag_gaussian']:
        for coordinate in ['none','gn','random']:
            name=f'anchor_{kind}_{coordinate}_s2026102501'
            methods.append(dict(id=name,family=kind+'_'+coordinate,seed=2026102501,
                checkpoint=str(OUT/'models/anchor_density'/name/'BEST.pt')))
    name='anchor_cfm_none_s2026102501'
    methods.append(dict(id=name,family='cfm_none',seed=2026102501,
        checkpoint=str(OUT/'models/anchor_density'/name/'BEST.pt')))
    comparisons=[dict(candidate_family=c,reference_family=r,metric='full66_nll')
        for c,r in [('spline_none','diag_gaussian_none'),('spline_gn','spline_none'),('spline_gn','spline_random')]]
    spec=dict(methods=methods,parents=[f'R4DEV_{i:05d}' for i in range(128)],comparisons=comparisons,
        sample_seed=20261025031,K=128,cohort='DEVELOPMENT_EXPOSED',proposal_selection='none',
        paid_context='common native96+probit64 plus receiverJ; cold costs must include this')
    path=folder/'METHODS.json';path.write_text(json.dumps(spec,indent=2)+'\n')
    frozen=folder/'FROZEN.json'
    commands=[['freeze','--methods',str(path),'--out',str(frozen)],
        ['infer','--freeze',str(frozen),'--contexts',str(export),'--out',str(folder/'online'),'--device','cpu'],
        ['score','--freeze',str(frozen),'--online',str(folder/'online'),'--truth',str(OUT/'inputs/offline_only'),
         '--out',str(folder/'scores'),'--device','cpu']]
    for command in commands:
        result=subprocess.run([sys.executable,'tools/evaluate_r4_density.py',*command],cwd=ROOT)
        if result.returncode:raise SystemExit(result.returncode)
    (folder/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',new_physics_RHS=0,
        parents=128,models=len(methods),reference_qualification='NOT_ESTABLISHED',test_role='development only'))+'\n')

if __name__=='__main__':main()

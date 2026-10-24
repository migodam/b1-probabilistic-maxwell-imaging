"""One-off development evaluation, with inference and truth scoring in separate PIDs.

This does not open locked cohorts or launch Maxwell queries. The caller profiles
the complete sequence through run_r4_job. All model/parent choices are registered
before the first evaluation and existing outputs are never silently replayed.
"""
from pathlib import Path
import argparse,json,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'research/p3l_r4_conditional_transport'
OLD=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('kind',choices=['risk','density']);a=ap.parse_args()
    folder=OUT/f'results/development_{a.kind}_first'
    if folder.exists():raise FileExistsError('new evaluation destination required')
    folder.mkdir(parents=True)
    parents=[f'R4DEV_{k:05d}' for k in range(128)]
    methods=[];comparisons=[]
    if a.kind=='density':
        for kind in ['gaussian','spline','triangular','cfm']:
            job=f'density_{kind}_32768_s2026102501'
            methods.append(dict(id=kind+'_s1',family=kind,seed=2026102501,
                                checkpoint=str(OUT/'models/density'/job/'BEST.pt')))
        comparisons=[dict(candidate_family='spline',reference_family='gaussian',metric='full66_nll'),
                     dict(candidate_family='triangular',reference_family='spline',metric='material_energy_score'),
                     dict(candidate_family='cfm__nfe32',reference_family='spline',metric='raw_sample_mean_full_mse')]
        spec=dict(methods=methods,parents=parents,comparisons=comparisons,sample_seed=20261025031,
                  cohort='NEW_R4_DEVELOPMENT_EXPOSED_AFTER_SCORING',K=128,proposal_selection='none')
        contexts=OUT/'inputs/online'
    else:
        methods=[dict(id='probit_anchor',family='anchor',loader='anchor',seed=0),
                 dict(id='r2_frozen',family='r2',loader='r2',seed=33,
                      checkpoint=str(OLD/'A/residual_training/n512/local_geometry_cnn_s33/BEST.pt'))]
        for family,pattern in [('local_cnn','risk_local_cnn'),('no_geometry_cnn','risk_no_geometry_cnn'),
                               ('phase8','risk_phase8'),('phase16','risk_phase16'),
                               ('bin_phase8','risk_bin_phase8'),('bin_phase16','risk_bin_phase16')]:
            for seed in [2026102501,2026102502,2026102503]:
                checkpoint=OUT/'models/risk'/f'{pattern}_s{seed}'/'BEST.pt'
                if not checkpoint.exists():raise FileNotFoundError(checkpoint)
                methods.append(dict(id=f'{family}_s{seed}',family=family,seed=seed,loader='r4',checkpoint=str(checkpoint)))
        for kind in ['affine','sharpen','spatial']:
            methods.append(dict(id=kind,family=kind,seed=2026102501,loader='r4',
                checkpoint=str(OUT/'models/risk'/f'risk_{kind}_s2026102501'/'BEST.pt')))
        comparisons=[dict(candidate_family=c,reference_family=r,action='clipped')
                     for c,r in [('local_cnn','r2'),('bin_phase8','local_cnn'),('bin_phase16','local_cnn')]]
        spec=dict(methods=methods,parents=parents,comparisons=comparisons,fallback='anchor',
                  cohort='NEW_R4_DEVELOPMENT_EXPOSED_AFTER_SCORING',GT_selection=False)
        contexts=OUT/'results/contexts/dev'
    methods_path=folder/'METHODS.json';methods_path.write_text(json.dumps(spec,indent=2)+'\n')
    frozen=folder/'FROZEN.json';online=folder/'online';scores=folder/'scores'
    tool=f'tools/evaluate_r4_{a.kind}.py'
    tasks=[['freeze','--methods',str(methods_path),'--out',str(frozen)],
           ['infer','--freeze',str(frozen),'--contexts',str(contexts),'--out',str(online),'--device','cpu'],
           ['score','--freeze',str(frozen),'--online',str(online),'--truth',str(OUT/'inputs/offline_only'),'--out',str(scores),'--device','cpu']]
    for args in tasks:
        process=subprocess.run([sys.executable,tool,*args],cwd=ROOT)
        if process.returncode:raise SystemExit(process.returncode)
    (folder/'COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',parents=len(parents),methods=len(methods),
        scientific_role='DEVELOPMENT_ONLY',physics_calls=0,locked_GT_read=False),indent=2)+'\n')

if __name__=='__main__':main()

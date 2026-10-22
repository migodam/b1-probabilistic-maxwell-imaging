"""One pre-observation method freeze for256 genuinely fresh P1 parents."""
from pathlib import Path
import json, subprocess, sys, time

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r2.core import OUT,Book,write


def main():
    marker=OUT/'METHODS_FROZEN_FOR_LOCKED.json'
    if marker.exists() or any((OUT/'inputs/online').glob('R2_locked_*')):
        raise FileExistsError('Methods cannot be selected after fresh observations exist')
    revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    # Freeze all deployed/model code and both prespecified evaluators. Later
    # report/plot scripts may be added without changing a frozen online method.
    source_paths=['src/p3r2','src/p1_material',*['tools/'+name for name in (
        'generate_p3r2_data.py','run_p3r2_anchors.py','run_p3r2_geometry.py',
        'run_p3r2_fresh_evaluation.py','run_p3r2_accounted.py',
        'infer_p3r2_locked.py','evaluate_p3r2_locked.py',
        'integrate_p3r2_fresh_trajectories.py','evaluate_p3r2_fresh_trajectories.py',
        'evaluate_p3r2_decision_forward.py','freeze_p3r2_methods.py')]]
    subprocess.run(['git','diff','--exit-code',revision,'--',*source_paths],cwd=ROOT,check=True)
    if subprocess.check_output(['git','ls-files','--others','--',*source_paths],cwd=ROOT,text=True).strip():
        raise RuntimeError('Commit the exact new source before freezing')
    groups=[('direct',OUT/'A/training'),('direct',OUT/'A/training32768'),
            ('residual',OUT/'A/residual_training'),('residual',OUT/'A/residual_controls'),
            ('residual',OUT/'A/residual_local_controls')]
    all_runs=[];ablations=[]
    for kind,folder in groups:
        if not (folder/'SELECTED.json').exists():raise FileNotFoundError(folder/'SELECTED.json')
        result=json.loads((folder/'RESULTS.json').read_text())
        for row in result['results']:
            path=Path(row['path'])/'BEST.pt'
            if not path.is_absolute():path=ROOT/path
            size=int(row['train_size']);seed=int(row['seed'])
            name=f"{row['kind']}_n{size}_s{seed}"
            if any(v['name']==name for v in ablations):raise ValueError('Duplicate frozen ablation')
            ablations.append(dict(name=name,kind=kind,checkpoint=str(path.resolve())))
            all_runs.append(dict(kind=kind,checkpoint=str(path.resolve()),result=row))
    direct=min((r for r in all_runs if r['kind']=='direct'),key=lambda r:r['result']['validation']['complex_roi_mse'])
    residual=min((r for r in all_runs if r['kind']=='residual'),key=lambda r:r['result']['validation']['complex_roi_mse'])
    density=json.loads((OUT/'A/density_training/SELECTED.json').read_text())
    weak_density=json.loads((OUT/'A/weak_density_training/SELECTED.json').read_text())
    fusion_path=OUT/'D/fusion/FROZEN.json';fusion=json.loads(fusion_path.read_text())
    parents=[f'R2_locked_{i:05d}' for i in range(256)]
    jobs=[]
    for start in range(0,256,32):
        jobs.append(dict(job=f'C_fresh_curvature_{start:03d}',parents=parents[start:start+32],
            arms=['dynamic_curvature'],seed=2026102400+start))
    jobs.append(dict(job='C_fresh_five_controls_000',parents=parents[:32],
        arms=['frozen_fixed','dynamic_fixed','dynamic_bandit','ordinary_random','dynamic_mbrl'],seed=2026102400))
    artifacts=[Path(r['checkpoint']) for r in all_runs]+[
        Path(density['checkpoint']),Path(weak_density['checkpoint']),fusion_path,
        Path(fusion['selected']['checkpoint']),OUT/'D/classical/FROZEN.json',
        OUT/'D/trajectory_integration/FROZEN.json',OUT/'C/policy/FROZEN_POLICY.json',
        OUT/'C/policy/FROZEN_MODEL_POLICY.json',OUT/'configs/TARGET.json',
        OUT/'configs/AUTHORIZATION.json',OUT/'A/FIXED_RISK_FRAME.npz',
        ROOT/'research/p1_material_target_20261018/inputs/geometry.json']
    records=[]
    for path in sorted(set(p.resolve() for p in artifacts)):
        stat=path.stat();records.append(dict(path=str(path),size_bytes=stat.st_size,mtime_ns=stat.st_mtime_ns))
    closures=['L/jobs/locked_inference/CLOSED.json','L/jobs/locked_trajectory_integration/CLOSED.json',
        'L/jobs/locked_decision_forward/CLOSED.json']+[f"C/jobs/{j['job']}/CLOSED.json" for j in jobs]
    b=Book('L','freeze_fresh_protocol',max_wall=120,analysis_only=True);status='COMPLETE';error=None
    try:
        write(marker,dict(schema_version='P3LR2-fresh-method-freeze-v1',status='FROZEN',
            source_revision=revision,frozen_unix=time.time(),frozen_source_paths=source_paths,
            target_id=json.loads((OUT/'configs/TARGET.json').read_text())['target_id'],
            master_seed=2026102202,locked_seed_group=2,expected_count=256,parents=parents,
            direct_checkpoint=direct['checkpoint'],residual_checkpoint=residual['checkpoint'],
            density_checkpoint=density['checkpoint'],weak_density_checkpoint=weak_density['checkpoint'],
            fusion_json=str(fusion_path),classical_action_json=str(OUT/'D/classical/FROZEN.json'),
            trajectory_weights_json=str(OUT/'D/trajectory_integration/FROZEN.json'),
            ablation_checkpoints=ablations,frozen_files=records,trajectory_jobs=jobs,
            include_matched_point176=True,
            generator_closed_json='L/jobs/data_locked_0_256/CLOSED.json',
            context_closed_jsons=['L/jobs/anchors_locked_0_256_96_64_m16/CLOSED.json'],
            additional_online_closures=closures,
            forward_decision_methods=['fused','residual','projected_direct','classical_postprocess','dynamic_curvature__arc_length_trapezoid'],
            primary_metric='mean independent-parent full-grid complex material MSE',
            primary_comparisons=['validation-selected supervised fusion versus native96+probit64',
                'dynamic-curvature arc-length material risk action versus native96+probit64'],
            familywise_intervals='Two comparisons: nominal95% and Bonferroni two-sided97.5% bootstrap intervals both reported',
            bootstrap=dict(draws=2000,seed=2026102292,unit='paired independent parent'),
            positive_research_target=dict(relative_MSE_reduction=.05,adjusted_CI_excludes_zero=True),
            secondary_design='All completed individual training seeds are scored; no locked selection. Same-seed capacity/permutation/decoder controls and first32 trajectory arms are secondary.',
            training_context='512train and128validation native48+probit48',
            fresh_context='256 originalprior parents, native96+probit64 and16extraordinary-costcontrol; context shift is disclosed',
            truth_roles='new simulator writes isolated truth; no estimator opens truth; scoring follows all256 completed online outcomes and closures',
            fallback='network failure->same paid anchor; trajectory incomplete/invalid->same paid anchor; no parent exclusion',
            probability_scope='Supervised material-risk actions; two fixed10 coordinate densities. FULL reference remains unqualified; no residual-softmax probabilities.',
            posterior_sampling_claim=False,exploration_exact_fiber_claim=False))
        print(json.dumps(dict(status='FROZEN',source_revision=revision,parents=256,checkpoints=len(ablations))),flush=True)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)


if __name__=='__main__':main()

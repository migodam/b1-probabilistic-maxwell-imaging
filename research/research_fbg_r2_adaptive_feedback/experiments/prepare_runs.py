"""Declared finite run manifests; all initial states are public and truth-free."""
from pathlib import Path
import argparse,datetime,hashlib,importlib.util,json
ROOT=Path(__file__).resolve().parents[1]
OLD='D:/AI/ZALG_R4_DIVERSITY_20261030/project'
REMOTE='D:/AI/FBG_R2_20261030'
METHODS=['native_trf','reciprocal_trf','data_lbfgs','gd_all','feedback_beta:0.2','feedback_qbeta']

def initial_module():
    spec=importlib.util.spec_from_file_location('fbg_public_starts',ROOT/'src/multistart.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def common(stage,tasks):
    return {'stage':stage,'status':'DEVELOPMENT_DISCLOSED','declared_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'lanes':4,'vendor_path':OLD+'/vendor','operators_path':OLD+'/public_data/operators.npz',
            'output_dir':REMOTE+'/results/'+stage,'tasks':tasks,'all_failures_counted':True}

def cfg(wall,divisor=1):
    return {'wall_seconds':wall,'max_lu':2000//divisor,'max_rhs':60000//divisor,'lower':.02,'upper':20.,
            'stages':[8,12,24],'stage_fractions':[.2,.25,.55],'maxiter':2000,'maxfun':4000,
            'smooth_weight':1e-4,'tv_weight':1e-4,'tv_epsilon':.01}

def stage2():
    cases=json.loads((ROOT/'data_contract/stage1_cases.json').read_text())['cases']
    cases=[c for c in cases if c['case_id'] in ('dev_moderate_ring','dev_high_digit8','dev_high_ring')]
    m=initial_module();starts=m.make_starts();tasks=[]
    for regime,K,wall in [('per_start',8,45.),('fixed_total',4,45./4),('fixed_total',8,45./8)]:
        for si in range(K):
            for ci,c in enumerate(cases):
                shift=(si+ci)%len(METHODS);order=METHODS[shift:]+METHODS[:shift]
                for method in order:
                    divisor=K if regime=='fixed_total' else 1
                    tasks.append({'task_id':c['case_id']+'__'+method.replace(':','_')+'__'+regime+'_K'+str(K)+'_s'+str(si),
                                  'case_id':c['case_id'],'method':method,'regime':regime,'K_run':K,'start_index':si,'start_name':m.START_NAMES[si],
                                  'initial_pixels':starts[si].tolist(),'public_path':OLD+'/public_data/'+Path(c['public_path']).name,
                                  'config':cfg(wall,divisor)})
    manifest=common('stage2',tasks)
    manifest.update({'protocol':'perstart K1/4/8 are nested prefixes of one K8 run; fixedtotalB45 K4/8 use shorter independent allocations',
                     'automatic_beta':'public uniform .1 median exact gradientflip threshold, upper half d, clip[.05,.4], empty0',
                     'start_seed_public':20261030,'selection':'top1 original full physical + commonregularizer; bestK offlineonly'})
    return manifest

def stage3():
    source=json.loads((ROOT/'data_contract/stage3_dev/manifest.json').read_text());tasks=[]
    for ci,c in enumerate(source['records']):
        for method in ('native_trf','gd_all','feedback_beta:0.2','feedback_qbeta'):
            tasks.append({'task_id':c['observation_id']+'__'+method.replace(':','_'),'case_id':c['case_id'],
                          'method':method,'p':c['p'],'layout':c['design'],'energy_protocol':'per_source',
                          'public_path':REMOTE+'/stage3_dev/'+c['public_path'],
                          'operators_path':REMOTE+'/stage3_dev/'+c['operators_path'],
                          'initial_material':.1,'config':cfg(45.)})
    manifest=common('stage3',tasks)
    manifest['energy_control']='fixedtotal clean is analytic global amplitude equivalence; no duplicateinverse; fixedabsolutenoise NOT_RUN'
    return manifest

def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['stage2','stage3']);a=ap.parse_args()
    manifest=globals()[a.stage]()
    (ROOT/'data_contract'/f'{a.stage}_remote.json').write_text(json.dumps(manifest,indent=2))
    files=list((ROOT/'src').glob('*.py'))
    (ROOT/'audit'/f'{a.stage.upper()}_FREEZE.json').write_text(json.dumps({'utc':manifest['declared_utc'],'files':[
        {'path':str(f.relative_to(ROOT)),'blake2b128':hashlib.blake2b(f.read_bytes(),digest_size=16).hexdigest()} for f in files],
        'tasks':len(manifest['tasks']),'fresh_truth':'NOT_GENERATED'},indent=2))
    print(json.dumps({'stage':a.stage,'tasks':len(manifest['tasks']),'maximum_sum_taskwall_s':sum(t['config']['wall_seconds'] for t in manifest['tasks'])}))
if __name__=='__main__':main()

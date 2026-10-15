"""Explicit resumable phases; never runs experiments on import."""
import argparse,json,sys,time,ast,os
from pathlib import Path
import numpy as np
from b1_joint.accounting import dump
from b1_joint.material import from_geometry,decode,geometry,material_distance,prior_event_masses
from .budget import Budget
from .target import Target,ROOT,OUT,HISTORY,TAU,to_p,profile
from .search import basic_reconstruction,endpoint,search

def e0():
    target=Target(Budget(OUT));data=np.load(HISTORY/'COMMON_RESPONSE_BANK.npz')
    latent=np.array([[0,0,2.5,.1],[1,1,2.5,.1],[-.5,.5,2.5,.1]])
    assert all(material_distance(latent[0],x)==0 for x in latent)
    assert np.allclose(prior_event_masses().sum(),1)
    response=target.proxy.predict(latent)
    assert np.max(abs(response-response[0]))==0
    f=data['fields_1.6'][0,0,2,2];assert np.max(abs(f-response[0]))<1e-12
    quals=json.loads((HISTORY/'PHYSICS_QUALIFICATION.json').read_text())
    assert np.array_equal(target.sigma[:32],np.repeat(quals['noise_sigma_per_frequency'][0],32))
    assert f.shape==(32,)
    # Scan only online implementation import/name contracts. No teacher/truth
    # object is accepted by either search or correction functions.
    import inspect
    from .search import search as s,correct
    assert all('truth' not in str(inspect.signature(fn)) for fn in [s,correct])
    online_tree=ast.parse((ROOT/'src/research_a/search.py').read_text())
    imports=[n.module or '' for n in ast.walk(online_tree) if isinstance(n,ast.ImportFrom)]
    assert not any('offline' in m for m in imports)
    receipt=dict(status='PASS',new_full_physics=0,inactive_latent_material_dedup=True,bank_response_and_packing_replay=True,original_noise_replay=True,prior_event_partition=True,online_signature_and_import_truth_isolation=True,scope='reused vector-Mie bank interface; package8dipolealgebra separate. New full backend checks belong to E1 audit.',historical_cost='P2 bank1250 prior material-frequency solves; notfreefirstsetup')
    dump(OUT/'results/E0_INTERFACE.json',receipt);print(receipt)

def prepare():
    from . import offline
    data=json.loads((HISTORY/'DEVELOPMENT_CASES.json').read_text())['cases']
    wanted=['solid','large_cavity','annular_cavity','combined'];selected=[next(d for d in data if d['id']==i) for i in wanted]
    dump(OUT/'results/DEV_OBSERVATIONS.json',[dict(id=d['id'],y=d['y'][:32],sigma=d['sigma'][:32]) for d in selected])
    dump(OUT/'results/offline_truth/DEV_TRUTH.json',[dict(id=d['id'],p=to_p(d['latent'])) for d in selected])

def audit(target):
    rng=np.random.default_rng(261015112);rows=[]
    for p in [np.array([0,0,2.5,.1]),np.array([.25,0,2.5,.1]),np.array([0,.08,2.5,.1])]:
        f,_,_,meta=target.forward(p,kind='INTERFACE_FULL_AUDIT')
        g,_,_,mi=target.forward(p,kind='INDEPENDENT_RECURSIVE_AUDIT',implementation='recursive')
        higher=target.forward(p,kind='MULTIPOLE_AUDIT',nmax=18)[0]
        rows.append(dict(p=p,recursive_relative=float(np.linalg.norm(f-g)/np.linalg.norm(f)),n14_n18_whitened=float(np.linalg.norm((f-higher)/target.sigma[:32])),boundary_backward_residual=meta.get('boundary_backward_residual'),Qabs=meta.get('Qabs')))
    assert max(r['recursive_relative'] for r in rows)<1e-10
    dump(OUT/'results/E1_BACKEND_AUDIT.json',dict(status='PASS',checks=rows,new_attempts=9,scope='numericalMietargetcheck; notcurrentVIESchurreplay'))

def campaign(target,name,ncalls,revision):
    observations=json.loads((OUT/f'results/{"DEV" if name.startswith("DEV") else "CONFIRM"}_OBSERVATIONS.json').read_text())
    path=OUT/f'results/{name}_CAMPAIGN.json'
    saved=json.loads(path.read_text()) if path.exists() else [];done={d['id'] for d in saved}
    c1=Path(os.environ.get('C1_P2',str(ROOT.parents[1]/'C1/c1-joint-p2/p2')));sys.path.insert(0,str(c1))
    from b1_joint.proposal_correction import load_conditioned_proposal,measurement_condition
    model=c1/'results/models/full_distill-26101321.pt'
    for i,obs in enumerate(observations):
        if obs['id'] in done:continue
        rng=np.random.default_rng(261015101+1000*i+(0 if name.startswith('DEV') else 100000))
        y=np.array(obs['y']);sigma=np.array(obs['sigma']);target.budget.context=dict(parent=obs['id'],role='coarse_reconstruction')
        p,receipt=basic_reconstruction(target,y,sigma,rng)
        base=endpoint(target,p,y,sigma,'SHARED_COARSE_ENDPOINT',dict(method='four proxy bounded conventional multistarts from best96prior',corrector=receipt))
        methods=[]
        parameters,mdnreceipt,density=load_conditioned_proposal(model,measurement_condition(y,sigma))
        arms=['adaptive'] if name=='DEV3' else (['finite_static','adaptive'] if name=='DEV2' else ['finite_static','multistart','population','mdn','adaptive'])
        for j,arm in enumerate(arms):
            armrng=np.random.default_rng(261015101+1000*i+j+(0 if name.startswith('DEV') else 100000))
            draws=density.sample_parameters(parameters,ncalls,armrng)['latent'][0] if arm=='mdn' else None
            target.budget.context=dict(parent=obs['id'],arm=arm)
            result=search(target,y,sigma,base,arm,ncalls,armrng,revision=revision,mdn_draws=draws)
            methods.append(result)
        saved.append(dict(id=obs['id'],base=base,methods=methods,observation=y,sigma=sigma,coarse_cost='one shared new endpoint; allarmsreceiveit; baseproposalwall separately accounted',mdn_receipt=mdnreceipt,revision=revision))
        dump(path,saved);print(json.dumps(dict(parent=obs['id'],worlds={r['arm']:len(r['archive']) for r in methods},physics=target.budget.count(),wall=target.budget.used_wall())),flush=True)
    return saved

def score(name):
    from .offline import score
    c=json.loads((OUT/f'results/{name}_CAMPAIGN.json').read_text())
    t=json.loads((OUT/f'results/offline_truth/{"DEV" if name.startswith("DEV") else "CONFIRM"}_TRUTH.json').read_text())
    threshold=None
    frozen=OUT/'configs/CONFIRM_FREEZE.json'
    if frozen.exists():threshold=json.loads(frozen.read_text())['coverage_normalized_mae']
    rows=score(c,t,threshold)
    dump(OUT/f'results/{name}_METRICS.json',rows)
    arms=sorted({r['arm'] for r in rows});print(json.dumps({a:dict(Best4=float(np.mean([r['by_k']['4']['best']['normalized_mae'] for r in rows if r['arm']==a and r['by_k']['4']])),worlds=float(np.mean([r['distinct_archive'] for r in rows if r['arm']==a]))) for a in arms},indent=2))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--phase',required=True,choices=['E0','prepare','audit','dev1','dev2','dev3','truth','confirm','score-dev1','score-dev2','score-dev3','score-confirm']);args=parser.parse_args()
    if args.phase=='E0':e0();return
    if args.phase=='prepare':prepare();return
    if args.phase.startswith('score'):
        score({'score-dev1':'DEV1','score-dev2':'DEV2','score-dev3':'DEV3','score-confirm':'CONFIRM'}[args.phase]);return
    ledger=Budget(OUT)
    with ledger.job(args.phase):
        target=Target(ledger)
        if args.phase=='audit':audit(target)
        elif args.phase=='truth':
            if not (OUT/'configs/CONFIRM_FREEZE.json').exists():raise RuntimeError('freeze before freshtruthgeneration')
            from .offline import generate
            generate(target)
        elif args.phase=='dev1':campaign(target,'DEV1',8,1)
        elif args.phase=='dev2':campaign(target,'DEV2',8,2)
        elif args.phase=='dev3':campaign(target,'DEV3',8,3)
        else:campaign(target,'CONFIRM',12,3)
    dump(OUT/'results/RESOURCE_SUMMARY_E1.json',ledger.summary())
    print(ledger.summary())

if __name__=='__main__':main()

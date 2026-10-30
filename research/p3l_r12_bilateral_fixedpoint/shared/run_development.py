"""Matched-RHS, multi-anchor first-order material exploration; never reads GT."""
from .physics import *
from .budget import BudgetStop
from ..agent_b.engine import Config,run_arm,LocalBudgetStop
from ..agent_a.modes import proposal_basis
import argparse,traceback,hashlib

ARMS={
 'H0_RML':('H0','full_prior',16),
 'R9_STYLE_FIXED8':('H1','full_prior',8),
 'H1_FULL16':('H1','full_prior',16),
 'H2_FULL16':('H2','full_prior',16),
 'HA_FULL16':('HA','full_prior',16),
 'H1_RECEIVER16':('H1','receiver_prior',16),
 'H1_PAIR16':('H1','pair_prior',16),
 'H1_FULL32':('H1','full_prior',32),
}

def anchors(parent):
    p=ROOT/'results/P1'/parent/'COMMON';a=dict(np.load(p/'ANCHORS.npz'));meta=json.loads((p/'SEALED.json').read_text())
    states=[]
    for i in range(len(a['xs'])):
        states.append(dict(x=a['xs'][i],material=a['materials'][i],F=a['Fs'][i],phi=float(a['phi'][i]),residual=white(a['Fs'][i]-np.load(ROOT/'data/online'/f'{parent}.npz')['observation']),jac=a['jac'][i],dmaterial=a['dmaterial'][i],L=a['L'][i],R=a['R'][i],prior_legal=True,ticket='COMMON_CACHED_'+str(i)))
    for s in states:
        pr=Problem();s['occupancy'],s['delta'],_,_=pr.decode(s['x'])
    eligible=[i for i in range(1,len(states)) if states[i]['phi']<=282.8737921622988]
    second=max(eligible,key=lambda i: material_rms(states[i]['material'],states[0]['material'])) if eligible else int(np.argmin(a['phi'][1:])+1)
    return [states[0],states[second]],meta

def merge_cloud(out,results,states,ph,common,parent,arm,started):
    cloud=[];qualified=[];layers=[];branches=[];traces=[]
    for b,res in enumerate(results):
        cloud+=res['states'];qualified+=res['qualified'];layers+=res['layer'];branches += [b]*len(res['states'])
        traces.extend([dict(branch=b,**r) for r in res['trace']])
    if not cloud:cloud=states;qualified=[s['phi']<=282.8737921622988 for s in cloud];layers=[0]*len(cloud);branches=[0]*len(cloud)
    available=[i for i,q in enumerate(qualified) if q]
    selected=min(available,key=lambda i: cloud[i]['phi']+.5*float(cloud[i]['x']@cloud[i]['x'])) if available else 0
    # Independent final complete nonlinear validation, included in arm cap.
    validation='NOT_RUN_LOCAL_CAP'
    try:
        final=ph.evaluate(cloud[selected]['x'],False);validation='EXACT_ROI_VALIDATED'
        if abs(final['phi']-cloud[selected]['phi'])>1e-6:raise RuntimeError('Final state response mismatch')
    except BudgetStop:pass
    noise_qualified=qualified.copy();unique=[];qualified=[]
    for i,(s,ok) in enumerate(zip(cloud,noise_qualified)):
        distinct=bool(ok and all(material_rms(s['material'],cloud[j]['material'])>.002 for j in unique))
        qualified.append(distinct)
        if distinct:unique.append(i)
    nonanchor=np.array(qualified)&np.array([material_rms(s['material'],states[0]['material'])>.002 for s in cloud])
    np.savez_compressed(out/'CLOUD.npz',xs=np.stack([s['x'] for s in cloud]),materials=np.stack([s['material'] for s in cloud]),Fs=np.stack([s['F'] for s in cloud]),phi=np.array([s['phi'] for s in cloud]),qualified=np.array(qualified),noise_qualified=np.array(noise_qualified),nonanchor_qualified=nonanchor,legal=np.ones(len(cloud),bool),layer=np.array(layers),branch=np.array(branches),prior_cost=np.array([.5*s['x']@s['x'] for s in cloud]),online_selected_index=selected,occupancyROI=np.stack([s['occupancy'] for s in cloud]),delta=np.array([s['delta'] for s in cloud]))
    atomic_json(out/'TRACE.json',traces)
    summary=dict(parent=parent,arm=arm,GT_used=False,utc=time.time(),selected=selected,final_validation=validation,new_RHS=ph.rhs,new_LU=ph.lu,common_new_RHS=common['new_RHS'],common_new_LU=common['new_LU'],standalone_RHS=ph.rhs+common['new_RHS'],standalone_LU=ph.lu+common['new_LU'],shared_historical_anchor='R9 read-only; notnewcharges',wall_seconds=time.perf_counter()-started,cloud_size=len(cloud),qualified=int(sum(qualified)),statuses=[r['status'] for r in results],refreshes=sum(r['refreshes'] for r in results),branches=2,prediction_blake2b=hashlib.blake2b((out/'CLOUD.npz').read_bytes(),digest_size=32).hexdigest())
    atomic_json(out/'SEALED.json',summary)
    print(parent,arm,'cloud',len(cloud),'qual',sum(qualified),'RHS',summary['standalone_RHS'],'refresh',summary['refreshes'],'wall',round(summary['wall_seconds'],2),flush=True)

def run(parent,arm):
    out=ROOT/'results/P1'/parent/arm
    if (out/'SEALED.json').exists():return
    out.mkdir(parents=True,exist_ok=True);ts=time.perf_counter()
    states,common=anchors(parent);y=np.load(ROOT/'data/online'/f'{parent}.npz')['observation'];problem=Problem('P1')
    cap=5400-common['new_RHS'];lucap=240-common['new_LU'];ph=Physics(problem,y,'P1/'+parent+'/'+arm,out/'states',local_rhs_cap=cap,local_lu_cap=lucap)
    atype,kind,rank=ARMS[arm];results=[]
    for b,s in enumerate(states):
        cfg=Config(arm=atype,basis_kind=kind,weak_rank=rank,layers=3,proposals_per_layer=2,corrector_states=5,local_rhs_cap=cap//2-8,local_lu_cap=(lucap-1)//2,seed=20261030+int(parent[-3:])*100+b,noise_phi_max=282.8737921622988,layer_radii=(.02,.05,.10),max_latent_radius=4.,adaptive_max_updates=2,distinct_rms=.002)
        rng=np.random.default_rng(np.random.SeedSequence([20261030,120012,201,int(parent[-3:]),b]))
        draws=rng.normal(size=(6,66));prior_targets=rng.normal(size=(6,66));targets=rng.normal(size=(6,512))
        factory=lambda st,cf: proposal_basis(st['jac'],st['dmaterial'],kind=cf.basis_kind,rank=cf.weak_rank,L=st['L'],R=st['R'])
        try:
            bp=out/f'branch{b}'
            if (bp/'SUMMARY.json').exists():
                old=json.loads((bp/'SUMMARY.json').read_text());arr=dict(np.load(bp/'CLOUD.npz'))
                recovered=[]
                for n in range(len(arr['xs'])):
                    recovered.append(dict(x=arr['xs'][n],material=arr['materials'][n],F=arr['Fs'][n],phi=float(arr['phi'][n]),occupancy=arr['occupancyROI'][n],delta=arr['delta'][n],prior_legal=True))
                result=dict(states=recovered,qualified=arr['qualified'].tolist(),layer=arr['layer'].tolist(),trace=[json.loads(z) for z in (bp/'TRACE.jsonl').read_text().splitlines()],status=old['status'],refreshes=old['refreshes'])
            else:
                result=run_arm(ph.evaluate,s['x'],cfg,initial_state=s,decoder=problem.decode,basis_factory=factory,random_draws=draws,prior_targets=prior_targets,residual_targets=targets,output=bp)
            results.append(result)
        except (BudgetStop,LocalBudgetStop) as e:
            atomic_json(out/f'BRANCH_{b}_BUDGET_STOP.json',dict(error=str(e),utc=time.time(),newRHS=ph.rhs));break
        except BaseException as e:
            atomic_json(out/f'BRANCH_{b}_FAILURE.json',dict(error=repr(e),traceback=traceback.format_exc(),utc=time.time(),newRHS=ph.rhs));raise
    merge_cloud(out,results,states,ph,common,parent,arm,ts)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--parents',default='0,1,2,3');ap.add_argument('--arms',default=','.join(ARMS));a=ap.parse_args()
    for i in map(int,a.parents.split(',')):
        for arm in a.arms.split(','):run(f'R9V_FRESH_{i:03d}',arm)

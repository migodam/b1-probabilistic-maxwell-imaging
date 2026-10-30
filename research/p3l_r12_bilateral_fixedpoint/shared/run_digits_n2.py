"""R11 N observation-only material exploration, original unchanged TV objective."""
from .physics import *
from .budget import BudgetStop
from ..agent_b.engine import Config,run_arm
from ..agent_a.modes import material_canonical_basis
from scipy.linalg import svd,cholesky
import argparse,hashlib,traceback

ARMS={'H0_GN':('H0',16),'H0_DATA_RML':('H0',16),'R9_STYLE_FIXED8':('H1',8),'H1_FULL16':('H1',16),'H2_FULL16':('H2',16),'HA_FULL16':('HA',16)}

def common(oid):
    path=ROOT/'results/N'/oid/'COMMON';path.mkdir(parents=True,exist_ok=True)
    y=np.load(ROOT/'data/online'/f'{oid}.npz')['y'];problem=Problem('NATIVE')
    x=np.load(R11/'results/N/FULL256'/oid/'PREDICTION.npz')['x']
    ph=Physics(problem,y,'N/'+oid+'/COMMON',path/'states',local_rhs_cap=40,local_lu_cap=1)
    state=ph.evaluate(x,True)
    if not (path/'ANCHOR.npz').exists():
        np.savez_compressed(path/'ANCHOR.npz',x=x,material=state['material'],F=state['F'],phi=state['phi'],occupancy=np.r_[np.zeros(0)],delta=state['delta'],dmaterial=state['dmaterial'])
        atomic_json(path/'SEALED.json',dict(utc=time.time(),Phi=state['phi'],new_RHS=ph.rhs,new_LU=ph.lu,GT_used=False,reused='R11FULL256blindinverseanchor',historical_R11_anchor_RHS='separatelocalhistoryreference, notnew'))
    return problem,state,y

def run(oid,arm):
    out=ROOT/'results/N2'/oid/arm
    if (out/'SEALED.json').exists():return
    out.mkdir(parents=True,exist_ok=True);ts=time.perf_counter();problem,anchor,y=common(oid)
    ph=Physics(problem,y,'N2/'+oid+'/'+arm,out/'states',local_rhs_cap=2400,local_lu_cap=140)
    # One frozen prior curvature for all exploration frames of this object;
    # frame refresh therefore cannot be credited with adapting a new image prior.
    factor_prior=cholesky(anchor['reg_hessian'],lower=False,check_finite=False)
    def factory(state,config):
        T=material_canonical_basis(state['dmaterial'])
        stacked=np.vstack((state['jac']@T,factor_prior@T))
        _,ss,vh=svd(stacked,full_matrices=False,check_finite=False)
        W=T@vh[-min(config.weak_rank,len(ss)):].T
        return dict(W=W,info=dict(kind='currentFullJbilateral+frozenoriginalTVpriorcurvature',material_rank=T.shape[1],selected_rank=W.shape[1],data_or_prior_information_added=False))
    atype,rank=ARMS[arm]
    config=Config(arm=atype,weak_rank=rank,layers=3,proposals_per_layer=2,corrector_states=5,local_rhs_cap=2392,local_lu_cap=139,layer_radii=(.02,.05,.10),noise_phi_max=282.8737921622988,seed=20261030+int(oid[-3:]),full_coordinate_corrector=True,distinct_rms=.002,native_box_retraction=True,active_set_iterations=259)
    rng=np.random.default_rng(np.random.SeedSequence([20261030,120012,301,int(oid[-3:])]))
    draws=rng.normal(size=(6,258));targets=rng.normal(size=(6,512))
    if arm=='H0_GN':targets=np.zeros((6,512))
    result=run_arm(ph.evaluate,anchor['x'],config,initial_state=anchor,decoder=problem.decode,basis_factory=factory,random_draws=draws,residual_targets=targets,output=out/'raw')
    states=result['states'];raw=np.array(result['qualified'],bool);kept=[];qualified=[]
    for i,(s,ok) in enumerate(zip(states,raw)):
        q=bool(ok and all(material_rms(s['material'],states[j]['material'])>.002 for j in kept));qualified.append(q)
        if q:kept.append(i)
    # Public original objective, never oracle image quality or digit topology.
    feasible=np.flatnonzero(raw)
    selected=min(feasible,key=lambda i: states[i]['phi']+problem.regularization(states[i]['x'])[0]) if len(feasible) else 0
    validation='EXACT_NONLINEAR_STATE_RECEIPT'
    try:ph.evaluate(states[selected]['x'],False)
    except BudgetStop:validation='ALREADY_EXACT_PAID_STATE_NO_EXTRA_FINAL_CAP'
    occupancies=np.zeros((len(states),576))
    occupancies[:,problem.idx]=np.stack([s['occupancy'] for s in states])
    np.savez_compressed(out/'CLOUD.npz',xs=np.stack([s['x'] for s in states]),materials=np.stack([s['material'] for s in states]),Fs=np.stack([s['F'] for s in states]),phi=np.array([s['phi'] for s in states]),qualified=np.array(qualified),noise_qualified=raw,legal=np.ones(len(states),bool),layer=np.array(result['layer']),occupancy=occupancies,delta=np.array([s['delta'] for s in states]),prior_cost=np.array([problem.regularization(s['x'])[0] for s in states]),online_selected_index=selected)
    atomic_json(out/'TRACE.json',result['trace'])
    atomic_json(out/'SEALED.json',dict(utc=time.time(),parent=oid,arm=arm,scope='N2box-feasibleDEVcorrection, notP1; N1boundaryfailurepreserved',GT_used=False,label_used=False,selected_index=int(selected),selection='originalunrandomizedPhi+R11TVregularizer',noise_qualified=int(raw.sum()),distinct_qualified=int(sum(qualified)),new_RHS=ph.rhs,new_LU=ph.lu,shared_common_RHS=40,shared_common_LU=1,refreshes=result['refreshes'],wall_seconds=time.perf_counter()-ts,status=result['status'],final_validation=validation,prediction_blake2b=hashlib.blake2b((out/'CLOUD.npz').read_bytes(),digest_size=32).hexdigest()))
    print(oid,arm,'qual',sum(qualified),'RHS',ph.rhs,'select',selected,'refresh',result['refreshes'],'wall',round(time.perf_counter()-ts,2),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--objects',default=','.join(map(str,range(9))));ap.add_argument('--arms',default=','.join(ARMS));ns=ap.parse_args()
    for i in map(int,ns.objects.split(',')):
        for arm in ns.arms.split(','):run(f'N_{i:03d}',arm)

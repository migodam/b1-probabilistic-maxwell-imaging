"""Frozen one-to-three-block numerical active sensing; no Maxwell imports.

Acquisition functions accept current observations and candidate predictions.
The Simulator is the only reader of future outcomes and truth. All comparisons
use parent/action/block-indexed common potential outcomes and persistent gain.
"""
from __future__ import annotations
import argparse,json,sys,time,csv
from pathlib import Path
import numpy as np
from scipy.stats import chi2
REPO=Path(__file__).resolve().parents[3];sys.path.insert(0,str(REPO))
from src.b1_joint.material import from_geometry,geometry,event_ids,void_volume
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.likelihood import conditional_gain
from acquisition_core import effective_information
from run_e0 import physical_jacobian
from sensing import Action,normalized,resample_indices,expected_utility,conditional_log_increment,field_l1
from e0_candidate_controls import profile
from cost import Cost,ROOT,dump

FREQS=[1.2,1.6,2.,2.4,2.8,3.2]
IDS=[f'frequency_{k:g}' if k!=1.6 else 'repeat_1.6' for k in FREQS]
POLICIES=['none','fixed_1.2','fixed_1.6','fixed_2','fixed_2.4','fixed_2.8','fixed_3.2','random','uniform',
    'raw','classical_Fisher','joint_projected_Fisher','event_EIG','expected_material_risk','robust','conflict_diagnostic',
    'calibration_only','calibrate_then_risk','geometry_only']

def menu(F,sigma,calibration_sd):
    return [Action(a,F[:,j],np.full(32,sigma[j]),kind='repeat' if k==1.6 else 'frequency') for j,(a,k) in enumerate(zip(IDS,FREQS))]+[
        Action('reference_calibration',np.ones((len(F),1)),np.array([calibration_sd]),kind='reference_calibration')]

def box_pair_min(A,B,C,lo1,hi1,lo2,hi2):
    """Minimum of A*x²+B*y²-2C*x*y on independent feasible gain intervals."""
    values=[]
    for x in [lo1,hi1]:
        y=np.clip(C*x/np.maximum(B,1e-300),lo2,hi2);values.append(A*x*x+B*y*y-2*C*x*y)
    for y in [lo2,hi2]:
        x=np.clip(C*y/np.maximum(A,1e-300),lo1,hi1);values.append(A*x*x+B*y*y-2*C*x*y)
    # The homogeneous unconstrained optimum is (0,0), excluded by positive gain.
    return np.maximum(np.min(values,axis=0),0.)

def candidate_scores(x,Fhistory,y,sigma,Fa,sa,extra=None):
    """Candidate-only finite composite screen; no normalized material weights."""
    f=Fhistory/sigma;yy=y/sigma;aa=np.sum(f*f,axis=1);bb=f@yy;cc=float(yy@yy)
    gm=bb/aa;g=np.clip(gm,.7,1.3);q=cc-2*bb*g+aa*g*g
    if extra is not None:q=q+extra
    cutoff=float(q.min()+9.)
    keep=np.argsort(q,kind='stable')[:48];keep=keep[q[keep]<=cutoff]
    if len(keep)<2:return dict(raw=0.,robust=0.,conflict=0.,active_candidates=len(keep),pair_count=0)
    # Bounded current-data compatibility; this is a set, not a posterior region.
    b0=cc-bb[keep]**2/aa[keep]+(0 if extra is None else extra[keep]);delta=np.sqrt(np.maximum((cutoff-b0)/aa[keep],0))
    lo=np.maximum(.7,gm[keep]-delta);hi=np.minimum(1.3,gm[keep]+delta)
    d=field_l1(x[keep],x[keep]);i,j=np.where(np.triu(d>=.02,1))
    if not len(i):return dict(raw=0.,robust=0.,conflict=0.,active_candidates=len(keep),pair_count=0)
    a=Fa[keep]/sa;A=np.sum(a*a,axis=1);cross=a@a.T
    robust=box_pair_min(A[i],A[j],cross[i,j],lo[i],hi[i],lo[j],hi[j])
    raw=A[i]+A[j]-2*cross[i,j]
    # Relative-gain secant diagnostic, retained as E0 negative algorithm.
    old=f[keep];oldA=np.sum(old*old,axis=1);oldcross=old@old.T
    c0=np.clip(oldcross[i,j]/oldA[j],.7/1.3,1.3/.7)
    ca=np.clip(cross[i,j]/np.maximum(A[j],1e-300),.7/1.3,1.3/.7)
    cj=np.clip((oldcross[i,j]+cross[i,j])/(oldA[j]+A[j]),.7/1.3,1.3/.7)
    fit=lambda c,ax,by,xy: ax+c*c*by-2*c*xy
    conflict=fit(cj,oldA[i]+A[i],oldA[j]+A[j],oldcross[i,j]+cross[i,j])-fit(c0,oldA[i],oldA[j],oldcross[i,j])-fit(ca,A[i],A[j],cross[i,j])
    return dict(raw=float(np.mean(raw)),robust=float(np.min(robust)),
        conflict=float(np.mean(np.maximum(conflict,0)*d[i,j])),active_candidates=len(keep),pair_count=len(i))

def local_fisher(response,base,history,obs,sd,actions,sigma,calibration_sd):
    jac=physical_jacobian(response,np.atleast_2d(base),FREQS)[0].reshape(6,32,4)
    B=[];N=[]
    for aid in history:
        if aid=='reference_calibration': B.append(np.zeros((1,4)));N.append(np.ones((1,1)))
        else:
            j=IDS.index(aid);B.append(jac[j]);N.append(response.predict(base,(FREQS[j],))[0,:,None])
    # Reconstruct the same full current history at the online base material.
    f=np.concatenate([v[:,0] for v in N]);m,v=conditional_gain(f[None],obs,sd);gain=float(np.clip(m[0],.7,1.3))
    B=np.vstack(B)*gain/sd[:,None];N=np.vstack(N)/sd[:,None]
    old=effective_information(B,N,np.array([[1/.15**2]]));out={}
    for a in actions:
        if a.id=='reference_calibration':ba=np.zeros((1,4));na=np.ones((1,1))/calibration_sd
        else:
            j=IDS.index(a.id);ba=gain*jac[j]/sigma[j];na=a.fields[0,:,None]/sigma[j]
        new=effective_information(np.vstack((B,ba)),np.vstack((N,na)),np.array([[1/.15**2]]))
        out[a.id]={'classical_Fisher':float(np.linalg.slogdet(np.eye(4)+B.T@B+ba.T@ba)[1]-np.linalg.slogdet(np.eye(4)+B.T@B)[1]),
            'joint_projected_Fisher':float(np.linalg.slogdet(np.eye(4)+new)[1]-np.linalg.slogdet(np.eye(4)+old)[1])}
    return out

def recommend(x,w,Fhist,y,sd,F,sigma,decisions,cx,CFhist,CF,c_extra,history,response,calibration_sd,seed,full=True,bundles=False):
    start=time.perf_counter();actions=menu(F,sigma,calibration_sd);cmenu=menu(CF,sigma,calibration_sd)
    cs={a.id:candidate_scores(cx,CFhist,y,sd,a.fields,a.sigma,c_extra) for a in cmenu}
    current_q,_=profile(CFhist,y,sd);current_q+=c_extra;base=cx[int(np.argmin(current_q))]
    baseF=response.predict(base,FREQS).reshape(1,6,32);bm=menu(baseF,sigma,calibration_sd)
    fisher=local_fisher(response,base,history,y,sd,bm,sigma,calibration_sd)
    utility={}
    if full:
        for a in actions:utility[a.id]=expected_utility(x,w,Fhist,y,sd,a,decisions,seed=seed,outer=128)
    picks={k:max(cs,key=lambda a:cs[a][k]) for k in ['raw','robust','conflict']}
    for k in ['classical_Fisher','joint_projected_Fisher']:picks[k]=max(fisher,key=lambda a:fisher[a][k])
    if full:
        picks['event_EIG']=max(utility,key=lambda a:utility[a]['event_EIG'])
        picks['expected_material_risk']=max(utility,key=lambda a:utility[a]['risk_reduction'])
    bundle_utilities={}
    if bundles:
        for a in actions:
            if a.id=='reference_calibration':continue
            b=Action('calibration+'+a.id,np.column_stack((np.ones(len(x)),a.fields)),np.r_[calibration_sd,a.sigma],blocks=2.,kind='bundle')
            u=expected_utility(x,w,Fhist,y,sd,b,decisions,seed=seed,outer=128)
            bundle_utilities[a.id]={k:v for k,v in u.items() if k!='risk_draws'}
        picks['two_step_bundle']=max(bundle_utilities,key=lambda a:bundle_utilities[a]['risk_reduction'])
    return dict(picks=picks,set_scores=cs,Fisher=fisher,utilities={a:{k:v for k,v in u.items() if k!='risk_draws'} for a,u in utility.items()},
        bundle_utilities=bundle_utilities,compute_seconds=time.perf_counter()-start)

def policy_action(policy,step,recommendation,parent,noise):
    if policy=='none':return None
    if policy.startswith('fixed_'):
        k=float(policy[6:]);return 'repeat_1.6' if k==1.6 else f'frequency_{k:g}'
    if policy=='random':return IDS[int(np.random.default_rng(26101565+parent*101+noise*17+step).integers(6))]
    if policy=='uniform':return [IDS[0],IDS[3],IDS[5]][step]
    if policy=='calibration_only':return 'reference_calibration'
    if policy=='geometry_only':return 'external_geometry'
    if policy=='calibrate_then_risk':return 'reference_calibration' if step==0 else recommendation['picks']['expected_material_risk']
    if policy=='two_step_bundle':
        return 'reference_calibration' if step==0 else recommendation['picks']['two_step_bundle'] if step==1 else None
    if policy=='risk_cost_stop':
        u=max(recommendation['utilities'].values(),key=lambda u:u['risk_reduction'])
        return None if u['risk_reduction']+2.58*u['mc_se']<.002 else recommendation['picks']['expected_material_risk']
    if policy=='robust_abstain':
        return None if max(v['robust'] for v in recommendation['set_scores'].values())<.25 else recommendation['picks']['robust']
    key='conflict' if policy=='conflict_diagnostic' else policy
    return recommendation['picks'][key]

class Simulator:
    def __init__(self,path):self.data=np.load(path,allow_pickle=False)
    def initial(self,p,r):return self.data['outcomes'][p,r,0,FREQS.index(1.6)].copy()
    def measure(self,p,r,step,aid):
        if aid=='reference_calibration':return np.atleast_1d(self.data['calibration'][p,r,step])
        if aid=='external_geometry':return np.atleast_1d(self.data['context'][p,r,step])
        return self.data['outcomes'][p,r,step+1,IDS.index(aid)].copy()
    def score(self,p,selected,decisions,event_decision):
        truth=self.data['truth'][p]
        return dict(selected_MAE=float(field_l1(selected,truth)[0,0]),best_of_K_MAE=float(field_l1(decisions,truth).min()),
            Coverage_at_K_005=bool(field_l1(decisions,truth).min()<=.05),
            selected_structure_correct=bool(event_ids(selected)==event_ids(truth)),event_decision_correct=bool(event_decision==event_ids(truth)))

def run(cost,blocks,stage,job):
    response=Response(ROOT/'results/E1_RESPONSE_BANK.npz');ref=Reference(response);bank=ref.bank(13);cfg=json.loads((ROOT/'configs/E1_FROZEN.json').read_text())
    anchors=np.load(ROOT/'results/E1_NOISE_ANCHORS.npz');sigma=anchors['sigma'];calibration_sd=sigma[1]/np.linalg.norm(anchors['F'][1])
    mesh=np.stack([v.ravel() for v in np.meshgrid(*[np.load(ROOT/'results/E1_RESPONSE_BANK.npz')[f'axis_{j}'] for j in range(4)],indexing='ij')],axis=1)
    cx=np.array([from_geometry(*v) for v in mesh]);raw=np.load(ROOT/'results/E1_RESPONSE_BANK.npz');CF=np.stack([raw[f'fields_{k:g}'].reshape(-1,32) for k in FREQS],axis=1)
    sim=Simulator(ROOT/'results/E1_FRESH_SIMULATOR.npz');rows=[];rankings=[];initial_reference_audit=[]
    policies=POLICIES if stage=='E1' else POLICIES+['risk_cost_stop','robust_abstain','two_step_bundle']
    dest=ROOT/'results'/f'{stage}_POLICY_RESULTS.json';rdest=ROOT/'results'/f'{stage}_RANKINGS.json'
    for parent in range(32):
        for noise in range(2):
            cost.check_time();y0=sim.initial(parent,noise);sd0=np.full(32,sigma[1]);rr=ref.solve(y0,sd0,order=13)
            # This comparison uses original prior weights, never equal cubature-node weights.
            if stage=='E1':
                fine=ref.solve(y0,sd0,order=17);initial_reference_audit.append(dict(parent=parent,noise=noise,
                    event_mass_max=float(abs(rr['event_mass']-fine['event_mass']).max()),mean_max=float(abs(rr['mean']-fine['mean']).max())))
            ids=resample_indices(rr['mass'],2048,26101530+parent*17+noise);x=bank['x'][ids];w0=np.full(2048,1/2048);F=response.predict(x,FREQS).reshape(-1,6,32)
            did=resample_indices(rr['mass'],96,26101531);decisions0=bank['x'][did]
            initial=recommend(x,w0,F[:,1],y0,sd0,F,sigma,decisions0,cx,CF[:,1],CF,np.zeros(len(cx)),[IDS[1]],response,calibration_sd,
                seed=26101570+parent*7+noise,bundles=stage=='E2')
            rankings.append(dict(parent=parent,noise=noise,step=0,**initial))
            candidate_loss=field_l1(cx,x);initial_loss=field_l1(decisions0,x)
            for policy in policies:
                w=w0.copy();histF=F[:,1].copy();chy=CF[:,1].copy();yy=y0.copy();sd=sd0.copy();history=[IDS[1]];cextra=np.zeros(len(cx));decisions=decisions0.copy()
                recommendation=initial;compute=0.;stop=False;geometry_count=0
                for step in range(blocks+1):
                    if step>0:
                        aid=policy_action(policy,step-1,recommendation,parent,noise)
                        if aid is None:
                            rows[-1]['stop_status']='PRESERVE_AMBIGUITY_FINITE_SET' if policy=='robust_abstain' else 'STOP_SINGLE_STEP_NOT_WORTH_DECLARED_BLOCK_COST' if policy=='risk_cost_stop' else 'STOP_TWO_BLOCK_BUNDLE_COMPLETE'
                            break
                        new=sim.measure(parent,noise,step-1,aid)
                        if aid=='external_geometry':
                            lw=-.5*((void_volume(x)-new[0])/.01)**2;w=normalized(np.log(np.maximum(w,1e-300))+lw)
                            cextra+=((void_volume(cx)-new[0])/.01)**2
                            geometry_count+=1
                        else:
                            j=None if aid=='reference_calibration' else IDS.index(aid)
                            ff=np.ones((len(x),1)) if j is None else F[:,j];cf=np.ones((len(cx),1)) if j is None else CF[:,j]
                            ns=np.array([calibration_sd]) if j is None else np.full(32,sigma[j])
                            increment=conditional_log_increment(histF,yy,sd,ff,new,ns)[0];w=normalized(np.log(np.maximum(w,1e-300))+increment)
                            histF=np.column_stack((histF,ff));chy=np.column_stack((chy,cf));yy=np.r_[yy,new];sd=np.r_[sd,ns];history.append(aid)
                        # An explicit reference proposal refresh expands the output catalogue.
                        # This is finite proxy-assisted material proposal, not new A discovery.
                        if blocks>1:
                            dd=resample_indices(w,64,26101580+step);decisions=np.concatenate((decisions,x[dd]))
                    loss=initial_loss if blocks==1 or step==0 else field_l1(decisions,x);di=int(np.argmin(loss@w));ev=int(np.argmax(np.bincount(event_ids(x),weights=w,minlength=12)))
                    cq,cg=profile(chy,yy,sd);cq+=cextra;ci=int(np.argmin(cq));prediction=float(cq[ci]);sentinel=prediction>chi2.ppf(.99,max(len(yy)+geometry_count-1,1))
                    reference_score=sim.score(parent,decisions[di],decisions,ev);candidate_score=sim.score(parent,cx[ci],cx,int(event_ids(cx[ci])))
                    # Offline reference risk of actual selected candidate, not a fabricated candidate posterior.
                    crisk=float(candidate_loss[ci]@w);oracle=float(np.min(candidate_loss@w))
                    rows.append(dict(parent=parent,noise=noise,policy=policy,step=step,selected_action=None if step==0 else aid,
                        mode1_candidate_only=candidate_score,mode2_reference_selection=reference_score,
                        candidate_profile_chi2=prediction,residual_sentinel=bool(sentinel),candidate_reference_regret=crisk-oracle,
                        finite_reference_risk=float((loss@w)[di]),posterior_ESS=float(1/np.sum(w*w)),reference_K=len(decisions),candidate_only_K=len(cx),
                        measurement_blocks=step,compute_acquisition_seconds=compute,history=history.copy(),
                        acquisition_input='reference-proxy' if policy in ['event_EIG','expected_material_risk','calibrate_then_risk','risk_cost_stop','two_step_bundle'] else 'candidate-only or fixed independent menu',
                        residual_sentinel_status='heuristic only; finite-candidate selection and bounded gain invalidate exact chi-square calibration',geometry_observations=geometry_count,
                        candidate_weights_kind='none',reference_status='finite conditional samples of continuous saved-response proxy',
                        exposure='method families frozen before fresh generation; completed implementation before outcome evaluation; numerical simulation' if blocks==1 else 'same E1 parents reused for explicitly exposed closed-loop followup',
                        generator_refresh=bool(blocks>1 and step>0),refresh_engine='finite reference proposal catalogue expansion' if blocks>1 else 'none'))
                    if step==blocks or policy=='none':break
                    if step==0:
                        recommendation=initial
                        if policy in ['raw','robust','conflict_diagnostic','classical_Fisher','joint_projected_Fisher','event_EIG','expected_material_risk','risk_cost_stop','robust_abstain','two_step_bundle']:
                            compute+=initial['compute_seconds']
                    elif policy in ['raw','robust','conflict_diagnostic','classical_Fisher','joint_projected_Fisher','event_EIG','expected_material_risk','calibrate_then_risk','risk_cost_stop','robust_abstain']:
                        recommendation=recommend(x,w,histF,yy,sd,F,sigma,decisions,cx,chy,CF,cextra,history,response,calibration_sd,
                            seed=26101570+parent*7+noise+step*1009,full=policy in ['event_EIG','expected_material_risk','calibrate_then_risk','risk_cost_stop'])
                        compute+=recommendation['compute_seconds']
                        rankings.append(dict(parent=parent,noise=noise,step=step,policy=policy,**recommendation))
            if noise==1: print(stage,'parent',parent,'rows',len(rows),'physics',cost.count(),flush=True)
        dump(dest,rows);dump(rdest,rankings)
        if stage=='E1':dump(ROOT/'results/E1_INITIAL_REFERENCE_AUDIT.json',initial_reference_audit)
    dump(ROOT/'results'/f'{stage}_STATUS.json',dict(status='EXECUTED_NUMERICAL_MAXWELL_MEASUREMENT_SIMULATION',parents=32,noise_per_parent=2,
        max_new_blocks=blocks,policies=policies,new_hardware_data=False,mode1='full-Mie625 candidate set and shared bounded-gain profile selection',
        mode2='continuous-proxy-derived2048 particles; uncertainty and support limits disclosed',A_real_package='WAITING_FOR_A_REAL_PACKAGE'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['E1','E2'],required=True);parser.add_argument('--job',required=True);args=parser.parse_args()
    with Cost().job(args.job,args.stage) as cost:run(cost,1 if args.stage=='E1' else 3,args.stage,args.job)

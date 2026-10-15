"""Real saved-response E0, zero solver imports/calls. Simulation, not hardware."""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from scipy.stats import truncnorm

REPO=Path(__file__).resolve().parents[3];sys.path.insert(0,str(REPO))
from src.b1_joint.material import geometry, from_geometry, RMAX, event_ids
from src.b1_joint.response import Response
from src.b1_joint.reference import Reference
from src.b1_joint.likelihood import conditional_gain, log_likelihood
from acquisition_core import effective_information
from cost import Cost,ROOT,dump
from sensing import Action,field_l1,resample_indices,expected_utility,conditional_log_increment,normalized,set_scores

SOURCE=REPO.parent/'b1-joint-p2/research/joint_p2/results'

def actions(f1,f2,sigma):
    pol=np.r_[np.arange(0,16,2),np.arange(16,32,2)]
    front=np.r_[np.arange(0,8),np.arange(16,24)]
    return [Action('repeat_1.6',f1,sigma[:32],kind='repeat'),
        Action('frequency_2.4',f2,sigma[32:]),
        Action('frequency_2.4_pol0',f2[:,pol],sigma[32:][pol]),
        Action('frequency_2.4_front',f2[:,front],sigma[32:][front]),
        Action('reference_calibration',np.ones((len(f1),1)),np.array([.03]),kind='reference_calibration')]

def physical_jacobian(response,x,frequencies=(1.6,2.4)):
    p=np.column_stack((geometry(x),x[:,2:]));lo=np.array([0,0,2.2,.04]);hi=np.r_[RMAX,2.8,.16];span=hi-lo
    result=[]
    for j in range(4):
        h=span[j]*1e-4;pa=p.copy();pb=p.copy();pa[:,j]=np.minimum(hi[j],pa[:,j]+h);pb[:,j]=np.maximum(lo[j],pb[:,j]-h)
        xa=np.array([from_geometry(*v) for v in pa]);xb=np.array([from_geometry(*v) for v in pb])
        result.append((response.predict(xa,frequencies)-response.predict(xb,frequencies))/(pa[:,j]-pb[:,j])[:,None]*span[j])
    return np.stack(result,axis=-1)

def fisher_scores(response,x,f0,y0,s0,menu):
    J=physical_jacobian(response,x);m,v=conditional_gain(f0,y0,s0);g=np.clip(m,.7,1.3);B0=J[:,:32]*g[:,None,None]/s0[None,:,None];N0=f0/s0
    vals={}
    pol=np.r_[np.arange(0,16,2),np.arange(16,32,2)];front=np.r_[np.arange(0,8),np.arange(16,24)]
    for a in menu:
        if a.id=='reference_calibration':Ja=np.zeros((len(x),1,4))
        elif a.id=='repeat_1.6':Ja=J[:,:32]
        elif a.id=='frequency_2.4':Ja=J[:,32:]
        elif a.id.endswith('pol0'):Ja=J[:,32:][:,pol]
        else:Ja=J[:,32:][:,front]
        classical=[];projected=[]
        for i in range(len(x)):
            b=B0[i];n=N0[i,:,None];ba=Ja[i]*g[i]/a.sigma[:,None];na=a.fields[i,:,None]/a.sigma[:,None]
            old=effective_information(b,n,np.array([[1/.15**2]]))
            new=effective_information(np.vstack((b,ba)),np.vstack((n,na)),np.array([[1/.15**2]]))
            classical.append(np.linalg.slogdet(np.eye(4)+b.T@b+ba.T@ba)[1]-np.linalg.slogdet(np.eye(4)+b.T@b)[1])
            projected.append(np.linalg.slogdet(np.eye(4)+new)[1]-np.linalg.slogdet(np.eye(4)+old)[1])
        vals[a.id]=dict(classical_Fisher=float(np.mean(classical)),joint_projected_Fisher=float(np.mean(projected)))
    return vals

def choose(scores, utilities, seed):
    ids=list(utilities);rng=np.random.default_rng(seed)
    picks={'no_measurement':'no_measurement','random':ids[int(rng.integers(len(ids)))],
        'predetermined_second_frequency':'frequency_2.4','uniform_one_block':'frequency_2.4'}
    for criterion in ['raw','robust','conflict','classical_Fisher','joint_projected_Fisher']:
        picks[criterion]=max(ids,key=lambda k:scores[k][criterion])
    picks['event_EIG']=max(ids,key=lambda k:utilities[k]['event_EIG'])
    picks['expected_material_risk']=max(ids,key=lambda k:utilities[k]['risk_reduction'])
    return picks

def score_realized(case,x,w,f0,sigma,decisions,picks):
    """Offline evaluator: truth and new independent noise enter only here."""
    F=np.asarray(case['F']);gain=case['gain'];rng=np.random.default_rng(26101552)
    menu={a.id:a for a in actions(f0,np.zeros_like(f0),sigma)}
    truth_menu={a.id:a for a in actions(F[None,:32],F[None,32:],sigma)}
    full=response.predict(x,(1.6,2.4))
    pred={a.id:a for a in actions(full[:,:32],full[:,32:],sigma)}
    initial=int(np.argmin(field_l1(decisions,x)@w));truth=np.asarray(case['latent']);best=float(field_l1(decisions,truth).min());out=[]
    by_action={}
    for aid,a in pred.items():
        obs=gain*truth_menu[aid].fields[0]+rng.normal(size=len(a.sigma))*a.sigma
        pp=normalized(np.log(np.maximum(w,1e-300))+conditional_log_increment(f0,np.array(case['y'])[:32],sigma[:32],a.fields,obs,a.sigma)[0])
        d=int(np.argmin(field_l1(decisions,x)@pp));event=int(np.argmax(np.bincount(event_ids(x),weights=pp,minlength=12)))
        by_action[aid]=(d,event,pp)
    y=np.asarray(case['y'])[:32]
    for policy,aid in picks.items():
        d,ev,pp=(initial,int(np.argmax(np.bincount(event_ids(x),weights=w,minlength=12))),w) if aid=='no_measurement' else by_action[aid]
        out.append(dict(case=case['id'],policy=policy,action=aid,selected_material_MAE=float(field_l1(decisions[d],truth)[0,0]),
            initial_selected_MAE=float(field_l1(decisions[initial],truth)[0,0]),best_of_K_MAE=best,K=len(decisions),
            structure_correct=bool(ev==int(event_ids(truth))),measurement_blocks=int(aid!='no_measurement'),
            exposure='existing exposed development case, new independent simulated noise',generator_refresh=False))
    return out

def run(cost):
    global response
    cfg=json.loads((ROOT/'configs/E0_FROZEN.json').read_text());response=Response(SOURCE/'COMMON_RESPONSE_BANK.npz');ref=Reference(response)
    cases=json.loads((SOURCE/'DEVELOPMENT_CASES.json').read_text())['cases'];inventory=[]
    for a in ['repeat_1.6','frequency_2.4','frequency_2.4_pol0','frequency_2.4_front','reference_calibration','no_measurement']:
        inventory.append(dict(action=a,status='ANALYTIC_SIMULATION' if a=='reference_calibration' else 'AVAILABLE',
            cache='COMMON_RESPONSE_BANK.npz and package full64 endpoints',new_Maxwell=0,independent_new_noise=a!='no_measurement'))
    for k in [1.2,2.,2.8,3.2]:inventory.append(dict(action=f'frequency_{k}',status='NOT_AVAILABLE',reason='not in E0 paid cache',new_Maxwell=0))
    dump(ROOT/'E0_ACTION_INVENTORY.json',dict(actions=inventory,already_conditioned_channels='1.6 Re16/Im16; not new evidence',retrospective_subset='NOT_RUN',hardware_status='SIMULATION_ONLY'))
    allrows=[];realized=[];ranking={};audit=[]
    for case in cases:
        name=case['id'];y=np.asarray(case['y'])[:32];sigma=np.asarray(case['sigma'])
        a=ref.solve(y,sigma[:32],order=13);fine=ref.solve(y,sigma[:32],order=17)
        audit.append(dict(case=name,event_mass_13_17=float(abs(a['event_mass']-fine['event_mass']).max()),mean_abs_13_17=float(abs(a['mean']-fine['mean']).max()),
            logZ_13_17=abs(a['log_evidence']-fine['log_evidence']),numerical_target='continuous saved-response proxy, not full posterior'))
        ids=resample_indices(fine['mass'],cfg['inner_particles'],cfg['seed']);bank=fine['bank'];x=bank['x'][ids];f0=bank['F1'][ids];f2=bank['F2'][ids];w=np.full(len(ids),1/len(ids))
        di=resample_indices(fine['mass'],cfg['decisions'],cfg['seed']+1);decisions=bank['x'][di]
        package=json.loads((SOURCE/'packages'/f'{name}.json').read_text());h=package['hypotheses'];cx=np.array([v['representative_latent'] for v in h]);cf=np.array([v['forward_prediction']['unscaled_field'] for v in h]);cmenu=actions(cf[:,:32],cf[:,32:],sigma)
        c_scores={a.id:set_scores(cx,cf[:,:32],y,sigma[:32],a) for a in cmenu}
        fish=fisher_scores(response,cx,cf[:,:32],y,sigma[:32],cmenu)
        scores={k:{**v,**fish[k]} for k,v in c_scores.items()};utilities={}
        for action in actions(f0,f2,sigma):
            cost.check_time();u=expected_utility(x,w,f0,y,sigma[:32],action,decisions,seed=cfg['seed'],outer=cfg['outer_predictive_draws'])
            utilities[action.id]=u
            row=dict(case=name,action=action.id,**{k:v for k,v in u.items() if k!='risk_draws'},**scores[action.id],new_Maxwell=0,
                probability_status='finite conditional approximation to continuous proxy reference',candidate_only_status='full-Mie endpoint set, no posterior weights')
            allrows.append(row)
        picks=choose(scores,utilities,cfg['seed']);ranking[name]=dict(selected_actions=picks,reference_event_mass=fine['event_mass'],candidate_count=len(cx),
            candidate_probability_kind='none',rank_by_material_risk=sorted(utilities,key=lambda k:-utilities[k]['risk_reduction']))
        realized+=score_realized(case,x,w,f0,sigma,decisions,picks)
        # Same physical reference; check whether representative-only approximation changes ranking.
        cp=np.array([fine['event_mass'][v['event_id']] for v in h]);rem=1-cp.sum()
        # Remainder is preserved as particles; no top-K mass renormalization.
        selected_events=np.array([v['event_id'] for v in h]);rest=~np.isin(event_ids(x),selected_events)
        xx=np.concatenate((cx,x[rest]));ww=np.r_[cp,np.full(rest.sum(),rem/max(rest.sum(),1))]
        if len(xx) and abs(ww.sum()-1)<1e-8:
            ff=response.predict(xx,(1.6,2.4));med={}
            for action in actions(ff[:,:32],ff[:,32:],sigma):
                z=expected_utility(xx,ww,ff[:,:32],y,sigma[:32],action,decisions,seed=cfg['seed'],outer=cfg['outer_predictive_draws'])
                med[action.id]=z['risk_reduction']
            ranking[name]['medoids_plus_remainder_risk_reduction']=med
        dump(ROOT/'results/E0_ACTION_RANKINGS.json',allrows);dump(ROOT/'results/E0_POLICY_REALIZED.json',realized);dump(ROOT/'results/E0_CRITERION_COMPARISON.json',ranking);dump(ROOT/'results/E0_REFERENCE_AUDIT.json',audit)
        np.savez_compressed(ROOT/'results'/f'E0_{name}_conditional.npz',latent=x,weights=w,F0=f0,F2=f2,decisions=decisions,y=y,sigma=sigma)
        print(name,'risk ranking',ranking[name]['rank_by_material_risk'],'candidate conflict',picks['conflict'],flush=True)
    import csv
    with (ROOT/'results/E0_POLICY_REALIZED.csv').open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(realized[0]));wr.writeheader();wr.writerows(realized)
    decision=dict(status='E0_EXECUTED',new_Maxwell=cost.count(),reference='saved-response continuous proxy; finite-particle predictive integration',
        criterion_disagreements=sum(len(set(v['selected_actions'].values())-{'no_measurement'})>1 for v in ranking.values()),
        signal_for_new_physics='Compare new-frequency material risk versus shared-gain calibration and repeat; raw/Fisher/set criteria differ, negative cases retained',
        not_claimed=['optimal physical frequency','exact full posterior utility','hardware measurement benefit','new candidate discovery'])
    dump(ROOT/'results/E0_DECISION.json',decision)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--job',default='e0-v1');args=parser.parse_args()
    cost=Cost()
    with cost.job(args.job,'E0'):run(cost)

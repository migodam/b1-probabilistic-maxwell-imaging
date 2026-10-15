"""Full-validated candidate-only selection and missing-family challenge."""
import sys,json
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from run_e0 import SOURCE,actions
from sensing import field_l1
from src.b1_joint.material import geometry,event_ids
from cost import ROOT,Cost,dump

def profile(F,y,sigma):
    a=F/sigma;b=np.asarray(y)/sigma
    g=np.clip(a@b/np.sum(a*a,axis=1),.7,1.3)
    return np.sum((g[:,None]*a-b)**2,axis=1),g

def run():
    rank=json.loads((ROOT/'results/E0_CRITERION_COMPARISON.json').read_text());cases=json.loads((SOURCE/'DEVELOPMENT_CASES.json').read_text())['cases'];rows=[];stress=[]
    for c in cases:
        p=json.loads((SOURCE/'packages'/f"{c['id']}.json").read_text());h=p['hypotheses'];x=np.array([a['representative_latent'] for a in h]);F=np.array([a['forward_prediction']['unscaled_field'] for a in h]);sigma=np.array(c['sigma']);y=np.array(c['y'])[:32]
        pred={a.id:a for a in actions(F[:,:32],F[:,32:],sigma)};truth={a.id:a for a in actions(np.array(c['F'])[None,:32],np.array(c['F'])[None,32:],sigma)}
        rng=np.random.default_rng(26101552);new={k:c['gain']*truth[k].fields[0]+rng.normal(size=len(v.sigma))*v.sigma for k,v in pred.items()}
        for deletion,keep in [('none',np.ones(len(x),bool)),('central_family_deleted',geometry(x)[:,0]==0),('annular_family_deleted',geometry(x)[:,1]==0)]:
            if not keep.any(): continue
            xx=x[keep];old=F[keep,:32];rr,_=profile(old,y,sigma[:32]);d0=int(np.argmin(rr));best=float(field_l1(xx,c['latent']).min())
            for policy,aid in rank[c['id']]['selected_actions'].items():
                ff=old;yy=y;sd=sigma[:32]
                if aid!='no_measurement':
                    a=pred[aid];ff=np.column_stack((ff,a.fields[keep]));yy=np.r_[yy,new[aid]];sd=np.r_[sd,a.sigma]
                residual,g=profile(ff,yy,sd);d=int(np.argmin(residual));mae=float(field_l1(xx[d],c['latent'])[0,0]);df=len(yy)-1
                row=dict(case=c['id'],policy=policy,action=aid,candidate_deletion=deletion,K=len(xx),
                    candidate_only_selected_MAE=mae,initial_candidate_only_MAE=float(field_l1(xx[d0],c['latent'])[0,0]),
                    best_of_K_MAE=best,selected_event=int(event_ids(xx[d])),true_event=int(event_ids(c['latent'])),
                    selected_structure_correct=bool(event_ids(xx[d])==event_ids(c['latent'])),
                    min_joint_profile_chi2=float(residual[d]),residual_sentinel=bool(residual[d]>chi2.ppf(.99,df)),
                    candidate_weights_kind='none',likelihood_used_for_ranking_only=True,posterior_probability=None,
                    missing_support_mass=None,full_Mie_candidates=True,generator_refresh=False)
                (rows if deletion=='none' else stress).append(row)
    dump(ROOT/'results/E0_CANDIDATE_ONLY_REALIZED.json',rows);dump(ROOT/'results/E0_COVERAGE_STRESS.json',stress)
    print('Full candidate-only rows',len(rows),'deletion challenge rows',len(stress),'new Maxwell 0')

if __name__=='__main__':
    with Cost().job('e0-candidate-controls-v1','E0'):run()

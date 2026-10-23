"""Observation-only integration and branch-preserving response compression.

All input F arrays are paid full-Maxwell PUBLIC prior nodes. No oracle labels
are loaded. Conditional intervals use newly integrated truncated intervals.
"""
import json,time
import numpy as np
from scipy.special import logsumexp
from .accounting import ROOT,write
from .posterior_math import log_marginal_likelihood

def posterior(rule,ll):
    ll=np.asarray(ll,float);prior=np.asarray(rule['weights'],float)
    if ll.shape!=prior.shape or not np.any(np.isfinite(ll)) or np.any(np.isnan(ll)) or np.any(ll==np.inf):raise ValueError('invalid likelihood array')
    # Shift likelihood BEFORE adding log volume; otherwise an enormous common
    # offset can erase the prior factors before logsumexp sees the array.
    anchor=float(np.max(ll));relative=ll-anchor+np.log(prior)
    relative_logz=float(logsumexp(relative));logz=anchor+relative_logz
    w=np.exp(relative-relative_logz);w/=w.sum()
    x=rule['x'];mean=float(w@x)
    branch_logs=[float(logsumexp(relative[rule['branches']==b])-relative_logz) for b in (0,1)]
    log_odds=branch_logs[1]-branch_logs[0]
    odds=float(np.exp(log_odds)) if np.isfinite(log_odds) and log_odds<=np.log(np.finfo(float).max) else None
    return dict(logZ=logz,weights=w,branch1=float(w[rule['branches']==1].sum()),
        x_mean=mean,x_sd=float(np.sqrt(w@((x-mean)**2))),contrast_mean=float(w@rule['contrast']),
        posterior_log_odds=log_odds,posterior_odds=odds,
        posterior_odds_status='FINITE' if odds is not None else 'FLOAT_RANGE_EXCEEDED_USE_LOG_ODDS')

def cdf_coherence(values,tolerance=1e-8):
    a=np.asarray(values,float)
    return bool(np.isfinite(a).all() and np.all(a>=-tolerance)
        and np.all(a<=1+tolerance) and np.all(np.diff(a)>=-tolerance))

def marginal(f,y):
    previous=None;history=[];selected=None
    for n in [11,21,41,65]:
        answer=log_marginal_likelihood(f,y,n=n)
        values=np.asarray(answer['logZ'],float)
        error=None if previous is None else float(np.max(abs(values-previous)))
        history.append(dict(order=n,max_logZ_change=error))
        selected=answer
        if error is not None and error<=2e-4:break
        previous=values
    return selected,history,bool(history[-1]['max_logZ_change'] is not None and history[-1]['max_logZ_change']<=2e-4)

def lagrange(nodes,points):
    nodes=np.asarray(nodes);points=np.asarray(points);n=len(nodes)
    delta=nodes[:,None]-nodes[None,:];np.fill_diagonal(delta,1.)
    weights=1/np.prod(delta,axis=1);out=np.empty((len(points),n))
    for i,p in enumerate(points):
        hit=np.where(abs(p-nodes)<1e-14)[0]
        if len(hit):out[i]=0.;out[i,hit[0]]=1.
        else:
            coeff=weights/(p-nodes);out[i]=coeff/coeff.sum()
    return out

def interpolate(rule3,b,u,c):
    result=np.empty((len(b),5,4,12),complex)
    for branch in [0,1]:
        ids=np.asarray(b)==branch
        if not np.any(ids):continue
        use=rule3['branches']==branch;un=np.unique(rule3['u'][use]);cn=np.unique(rule3['contrast'][use])
        cube=rule3['F'][use].reshape(len(un),len(cn),5,4,12)
        lu=lagrange(un,np.asarray(u)[ids]);lc=lagrange(cn,np.asarray(c)[ids])
        result[ids]=np.einsum('bi,bj,ijfrs->bfrs',lu,lc,cube)
    return result

def main():
    start=time.process_time();protocol=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text())
    directory=ROOT/'results/P1A/rules'
    if not (directory/'COMPLETE.json').exists():raise RuntimeError('full paid quadrature bank incomplete')
    rules={p.stem:dict(np.load(p)) for p in directory.glob('*.npz')}
    for key,r in rules.items():
        if not np.all(np.isfinite(r['F'])) or not np.all(r['weights']>0):raise RuntimeError('invalid probability quadrature arrays: '+key)
        if key.startswith('full_') and abs(float(sum(r['weights']))-1)>1e-12:raise RuntimeError('full prior quadrature is not normalized')
    obs=np.load(ROOT/'results/P1A/data/OBSERVATIONS.npz');ys=obs['y']
    out=ROOT/'results/P1A/reference';out.mkdir(parents=True,exist_ok=True)
    rows=[];weights=[];likelihoods=[]
    for index,y in enumerate(ys):
        records={};lik={};nuisance={};qualified=True;compressed_ok=True;raw3_ok=True
        for order in protocol['full_orders']:
            key=f'full_{order}';a,h,ok=marginal(rules[key]['F'],y)
            nuisance[key]=h
            if order in protocol['reference_orders']:qualified &=ok
            else:raw3_ok &=ok
            lik[key]=np.asarray(a['logZ'],float)
            records[key]=posterior(rules[key],lik[key])
        ref=records['full_9'];coarse=records['full_5'];cdf={};compressed_cdf={};direct3_cdf={}
        full3=rules['full_3'];fine=rules['full_9']
        approximation=interpolate(full3,fine['branches'],fine['u'],fine['contrast'])
        comp_l,h,ok=marginal(approximation,y);compressed_ok &=ok;nuisance['compressed']=h
        comp=posterior(fine,np.asarray(comp_l['logZ'],float))
        for i,cut in enumerate(protocol['cdf_checkpoints']):
            values=[]
            for order in protocol['reference_orders']:
                key=f'cdf_{i}_{order}';r=rules[key]
                a,h,ok=marginal(r['F'],y);qualified &=ok;nuisance[key]=h
                ll=np.asarray(a['logZ'],float)
                values.append(float(np.exp(logsumexp(np.log(r['weights'])+ll)-records[f'full_{order}']['logZ'])))
            cdf[str(cut)]=dict(order5=values[0],order9=values[1],error=abs(values[1]-values[0]))
            key3=f'cdf_{i}_3'
            if key3 in rules:
                r3=rules[key3];a3,h3,ok3=marginal(r3['F'],y);raw3_ok &=ok3;nuisance[key3]=h3
                direct3_cdf[str(cut)]=float(np.exp(logsumexp(np.log(r3['weights'])+np.asarray(a3['logZ'],float))-records['full_3']['logZ']))
            r=rules[f'cdf_{i}_9'];af=interpolate(full3,r['branches'],r['u'],r['contrast'])
            al,ch,ok=marginal(af,y);compressed_ok &=ok;nuisance[f'compressed_cdf_{i}']=ch
            compressed_cdf[str(cut)]=float(np.exp(logsumexp(np.log(r['weights'])+np.asarray(al['logZ'],float))-comp['logZ']))
        event_error=abs(ref['branch1']-coarse['branch1']);cdf_error=max(x['error'] for x in cdf.values())
        reference_coherence=all(cdf_coherence([cdf[str(cut)][f'order{n}'] for cut in protocol['cdf_checkpoints']]) for n in protocol['reference_orders'])
        compressed_coherence=cdf_coherence([compressed_cdf[str(cut)] for cut in protocol['cdf_checkpoints']])
        qualified &=event_error<=.02 and cdf_error<=.02 and reference_coherence
        compressed_ok &=compressed_coherence
        raw3=records['full_3']
        response_error=float(np.max(np.linalg.norm((approximation-fine['F']).reshape(len(approximation),-1),axis=1)/.01))
        row=dict(observation=index,parent=int(obs['parent_ids'][index]),
            target='CONTROLLED_W_FIXED_SUPPORT_JOINT_V1_CONTINUOUS',reference_status='QUALIFIED_EMPIRICAL_QUADRATURE' if qualified else 'REFERENCE_UNQUALIFIED',
            branch1_posterior_mass=ref['branch1'],prior_branch1_mass=.5,
            posterior_odds=ref['posterior_odds'],posterior_log_odds=ref['posterior_log_odds'],posterior_odds_status=ref['posterior_odds_status'],
            x_mean=ref['x_mean'],x_sd=ref['x_sd'],contrast_mean=ref['contrast_mean'],log_evidence=ref['logZ'],
            quadrature_event_error=event_error,quadrature_CDF_max_error=cdf_error,CDF=cdf,
            quadrature_x_mean_error=abs(ref['x_mean']-coarse['x_mean']),quadrature_x_SD_error=abs(ref['x_sd']-coarse['x_sd']),
            reference_CDF_coherence=reference_coherence,compressed_CDF_coherence=compressed_coherence,
            compressed_nuisance_qualified=compressed_ok,direct_order3_nuisance_qualified=raw3_ok,
            compressed_branch_error=abs(comp['branch1']-ref['branch1']),
            compressed_CDF_max_error=max(abs(compressed_cdf[k]-v['order9']) for k,v in cdf.items()),
            compressed_x_mean_error=abs(comp['x_mean']-ref['x_mean']),compressed_response_max_error_white=response_error,
            direct_order3_branch_error=abs(raw3['branch1']-ref['branch1']),nuisance_convergence=nuisance,
            direct_order3_CDF=direct3_cdf,
            direct_order3_CDF_coherence=cdf_coherence(list(direct3_cdf.values())) if direct3_cdf else None,
            direct_order3_CDF_max_error=max(abs(direct3_cdf[k]-v['order9']) for k,v in cdf.items()) if direct3_cdf else None,
            outside_model_support_mass=None,uncovered_other_wall_families='UNQUANTIFIED')
        weights.append(ref['weights']);likelihoods.append(lik['full_9']);rows.append(row)
        write(out/'RESULTS.json',rows);print(json.dumps({k:row[k] for k in ['observation','reference_status','branch1_posterior_mass','quadrature_event_error','quadrature_CDF_max_error','compressed_branch_error']}),flush=True)
    r=rules['full_9'];np.savez_compressed(out/'POSTERIOR_NODES.npz',posterior_weights=np.array(weights),log_marginal_likelihood=np.array(likelihoods),F=r['F'],branches=r['branches'],u=r['u'],contrast=r['contrast'],x=r['x'],prior_weights=r['weights'])
    qualified=all(r['reference_status']=='QUALIFIED_EMPIRICAL_QUADRATURE' for r in rows)
    gate=dict(status='P1_A_POSTERIOR_VALID_COMPRESSION_NOT_YET_COST_QUALIFIED' if qualified else 'REFERENCE_UNQUALIFIED',
        probability_reference_qualified=qualified,observations=len(rows),independent_parents=len(np.unique(obs['parent_ids'])),
        compressed_numerically_qualified=all(r['compressed_nuisance_qualified'] and r['compressed_branch_error']<=.02 and r['compressed_CDF_max_error']<=.02 for r in rows),
        max_event_error=max(r['quadrature_event_error'] for r in rows),max_CDF_error=max(r['quadrature_CDF_max_error'] for r in rows),
        compressed_event_error=max(r['compressed_branch_error'] for r in rows),compressed_CDF_error=max(r['compressed_CDF_max_error'] for r in rows),
        original_wall_prior_unresolved=True,probability_scope='explicit newly frozen finite-family prior only',
        posterior_mean_not_a_new_physical_world=True,cpu_seconds=time.process_time()-start,
        new_physics_rhs=0,cost_note='All full-reference/source-bank physics was paid in P1A; shared teacher setup is not free deployment.')
    write(out/'GATE.json',gate);print(json.dumps(gate),flush=True)

if __name__=='__main__':main()

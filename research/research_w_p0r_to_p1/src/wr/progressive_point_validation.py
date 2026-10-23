"""Pay for exact full-Maxwell checks of the two observation-only point arms.

No generation module or offline labels are imported. Interpolation failure
falls back to a globally enumerated shared-gain profile over the SAME paid
full9 rule, independently of truth. This is a point decision, not posterior
profiling or an alteration of the marginalized probability target.
"""
import json
from pathlib import Path
import numpy as np
from .accounting import ROOT, Book, write
from .probability_family import Family


def profile_points(F,y):
    f=np.asarray(F).reshape(len(F),-1);v=np.asarray(y).ravel()
    cross=f@v.conj();norm=np.sum(abs(f)**2,axis=1);ynorm=float(np.sum(abs(v)**2))
    angles=np.column_stack([np.full(len(f),-.01),np.full(len(f),.01),np.clip(-np.angle(cross),-.01,.01)])
    real=np.real(np.exp(1j*angles)*cross[:,None])
    amps=np.clip(real/norm[:,None],.99,1.01)
    phi=(amps**2*norm[:,None]-2*amps*real+ynorm)/(2*.01**2)
    jj=np.argmin(phi,axis=1);ii=np.arange(len(f))
    return phi[ii,jj],amps[ii,jj],angles[ii,jj]


def key(b,u,c):return (int(b),float(u).hex(),float(c).hex())


def main():
    gate=json.loads((ROOT/'results/P1A/reference/GATE.json').read_text())
    if gate.get('probability_reference_qualified') is not True:raise RuntimeError('P1A reference gate closed')
    source=ROOT/'results/P1B/inference'
    if not (source/'GATE.json').exists() and not (source/'COMPLETE.json').exists():
        raise RuntimeError('Freeze inference output before point validation')
    rows=json.loads((source/'RESULTS.json').read_text())
    observations=np.load(ROOT/'results/P1B/data/OBSERVATIONS.npz')
    if len(rows)!=len(observations['y']):raise RuntimeError('Missing point arms/observations; do not omit failed cases')
    protocol=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text())
    limits=protocol['P1B']['point_full_validation']
    family=Family.load();book=Book('P1B');out=ROOT/'results/P1B/point_validation';out.mkdir(parents=True,exist_ok=True)
    rule=dict(np.load(ROOT/'results/P1A/rules/full_9.npz'))
    cache={key(b,u,c):dict(F=f,source='P1A/full_9',source_rhs=20,incremental_rhs=0) for b,u,c,f in zip(rule['branches'],rule['u'],rule['contrast'],rule['F'])}
    physical_rule=ROOT/'results/P1B/full_material_rule/FULL_RULE.npz'
    if physical_rule.exists():
        full=dict(np.load(physical_rule))
        for b,u,c,f in zip(full['branches'],full['u'],full['contrast'],full['F']):
            cache[key(b,u,c)]=dict(F=f,source='P1B/full_material_rule actual paid Full-Maxwell node',source_rhs=20,incremental_rhs=0)
    results=[];method_files={'direct joint point':'direct_joint_point','progressive block point':'progressive_block_point'}
    for row in rows:
        i=row['observation'];parent=int(observations['parent_ids'][i]);y=observations['y'][i]
        dest=out/f'observation_{i:03d}';dest.mkdir(exist_ok=True);record=dest/'RESULT.json'
        if record.exists():
            old=json.loads(record.read_text())
            if old.get('inference_input_snapshot')!=row['methods']:raise RuntimeError('Point resume inference differs')
            results.append(old);continue
        result=dict(observation=i,parent=parent,methods={},inference_input_snapshot=row['methods'],truth_access=False)
        for name,stem in method_files.items():
            point=row['methods'][name];approx=np.load(source/f'observation_{i:03d}'/row['attempt']/(stem+'.npz'))
            b,u,c,a,p=(float(v) for v in approx['params']);b=int(b)
            if not np.allclose([b,u,c,a,p],[point[k] for k in ('branch','u','contrast','amplitude','phase')],rtol=0,atol=0):raise RuntimeError('Point row and arrays differ')
            k=key(b,u,c)
            if k not in cache:
                folder=dest/stem;attempt=len(list(dest.glob(stem+'_attempt_*')))
                folder=dest/f'{stem}_attempt_{attempt:03d}'
                f=family.point(b,u,c,book,str(folder.relative_to(ROOT/'results')))
                cache[k]=dict(F=f,source=str(folder.relative_to(ROOT)),source_rhs=20,incremental_rhs=20)
                # Persist before comparing anything; a failed check is still a paid endpoint.
                np.savez_compressed(dest/(stem+'_VALIDATION_F.npz'),F=f,params=np.array([b,u,c,a,p]))
                increment=20
            else:
                f=cache[k]['F'];increment=0
            f=np.asarray(f);gamma=a*np.exp(1j*p)
            exact_phi=float(.5*np.sum(abs(gamma*f-y)**2)/.01**2)
            response_error=float(np.linalg.norm((f-approx['F_approx']).ravel())/.01)
            potential_error=abs(exact_phi-float(approx['Phi']))
            fallback=response_error>limits['whitened_response_error_max'] or potential_error>limits['potential_error_max']
            r=dict(branch=b,u=u,contrast=c,amplitude=a,phase=p,x=family.coordinates(b,u,c)['center'][0],
                Phi=exact_phi,approximate_Phi=float(approx['Phi']),response_error_white=response_error,
                potential_error=potential_error,fallback_used=bool(fallback),new_rhs=increment,
                physical_source=cache[k]['source'],source_endpoint_rhs=cache[k]['source_rhs'],
                full_validation_status='FULL_CHECK_PASSED' if not fallback else 'INTERPOLATION_FAILED_FALLBACK_TO_FULL9',
                independent_grid_validation=False,
                full_validation_meaning='Exact response at the common qualified production solver; existing paid full nodes are reused, not relabeled as new independent-grid solves.',
                original_point=point,weight_semantics='POINT_ESTIMATE_NO_POSTERIOR_MASS')
            if fallback:
                phis,amps,phases=profile_points(rule['F'],y);j=int(np.argmin(phis))
                b=int(rule['branches'][j]);u=float(rule['u'][j]);c=float(rule['contrast'][j]);a=float(amps[j]);p=float(phases[j]);f=rule['F'][j]
                exact_phi=float(.5*np.sum(abs(a*np.exp(1j*p)*f-y)**2)/.01**2)
                r.update(branch=b,u=u,contrast=c,amplitude=a,phase=p,x=float(rule['x'][j]),Phi=exact_phi,
                    physical_source='P1A/full_9 node '+str(j),fallback_rule='Shared-gain global box profile over all paid full9 nodes',
                    fallback_source_cost_rhs=20,fallback_incremental_rhs=0)
            np.savez_compressed(dest/(stem+'_FINAL.npz'),F=f,params=np.array([b,u,c,a,p]),Phi=exact_phi)
            result['methods'][name]=r
        write(record,result);results.append(result);write(out/'RESULTS.json',results)
        print(json.dumps(dict(observation=i,methods={k:dict(Phi=v['Phi'],fallback=v['fallback_used'],white_error=v['response_error_white']) for k,v in result['methods'].items()})),flush=True)
    write(out/'COMPLETE.json',dict(status='FULL_POINT_VALIDATION_COMPLETE',observations=len(results),
        inference_frozen=True,truth_access=False,fallbacks=sum(v['fallback_used'] for r in results for v in r['methods'].values()),
        physical_budget=book.totals()['P1B'],cache_cost='P1A bank setup disclosed separately; exact reuse is zero incremental RHS, not free cold deployment'))


if __name__=='__main__':main()

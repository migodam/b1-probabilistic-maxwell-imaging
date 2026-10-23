#!/usr/bin/env python3
"""Read-only P1A probability array replay. No Book, physics, labels or main imports.

Reads exactly full5/full9, five directly truncated CDF5/9 rule archives,
observation-only y/parent_ids and the saved P1A reference RESULTS.json.
Writes only a fresh independent AUDIT/array_replay run directory.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
import resource
import sys
import time
import numpy as np
from scipy.special import logsumexp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
# This module is algebra only: it imports no Book, backend or file reader.
from wr.posterior_math import log_marginal_likelihood,nuisance_quadrature

CUTS=(-.909,-.9,-.85,-.8,-.791)
ORDERS=(5,9)
NUISANCE=(11,21,41,65)
REPLAY_ATOL=2e-9


def write(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def expected_rule(order,cut=None):
    z,w=np.polynomial.legendre.leggauss(order);b,u,c,x,m=[],[],[],[],[]
    for branch,center in enumerate((-.9,-.8)):
        lo=-.01;hi=.01 if cut is None else min(.01,cut-center)
        if hi<=lo:continue
        for i,offset in enumerate((hi+lo)/2+(hi-lo)*z/2):
            for j,contrast in enumerate(3+.06*z):
                b.append(branch);u.append(offset);c.append(contrast);x.append(center+offset)
                m.append(.5*(w[i]/2*(hi-lo)/.02)*(w[j]/2))
    return dict(branches=np.array(b),u=np.array(u),contrast=np.array(c),x=np.array(x),weights=np.array(m))


def load_rule(path,order,cut):
    with np.load(path,allow_pickle=False) as archive:
        rule={k:archive[k].copy() for k in ('F','branches','u','contrast','x','weights')}
    expected=expected_rule(order,cut)
    for k,v in expected.items():
        if rule[k].shape!=v.shape or not np.allclose(rule[k],v,rtol=0,atol=2e-14):
            raise ValueError('Frozen direct tensor/truncated prior rule mismatch: '+str(path)+':'+k)
    n=len(rule['weights'])
    if rule['F'].shape!=(n,5,4,12) or not np.iscomplexobj(rule['F']) or not np.isfinite(rule['F']).all():raise ValueError('Invalid Full response array')
    if not np.all(rule['weights']>0):raise ValueError('Nonpositive prior integration weights')
    mass=float(sum(rule['weights']));target=1. if cut is None else sum(.5*np.clip((cut-center+.01)/.02,0,1) for center in (-.9,-.8))
    if abs(mass-target)>2e-14:raise ValueError('Prior/truncated numerator mass mismatch')
    return rule,dict(rule=path.name,node_count=n,prior_mass=mass,analytic_prior_mass=float(target),
                     truncated_weights_renormalized=False,F_bytes=int(rule['F'].nbytes))


def integrate(rule,y,batch_size=16):
    """Proper normalized uniform gain integral ONCE after all 240 complex channels.

    Each batch calls algebra-only long-double sufficient statistics + normalized
    amplitude/phase GL tensor integration, including the 480-real Gaussian
    density constant. Truncated material weights retain ORIGINAL prior density.
    Compare the SAME global max log-likelihood changes as the saved reference.
    """
    previous=None;history=[];selected=None;gain_norm=[]
    for n in NUISANCE:
        chunks=[]
        for first in range(0,len(rule['F']),batch_size):
            out=log_marginal_likelihood(rule['F'][first:first+batch_size],y,n=n)
            chunks.append(np.asarray(out['logZ'],float))
            gain_norm.append(abs(float(out['prior_quadrature_weight_sum'])-1.))
        ll=np.concatenate(chunks)
        change=None if previous is None else float(np.max(abs(ll-previous)))
        history.append(dict(order=n,max_logZ_change=change));selected=ll
        if change is not None and change<=2e-4:break
        previous=ll
    # Shift before adding prior weights so large common log densities do not
    # erase unequal material integration volumes. Never node-count normalize.
    anchor=float(np.max(selected));relative=selected-anchor+np.log(rule['weights'])
    normalizer=float(logsumexp(relative));weights=np.exp(relative-normalizer);weights/=weights.sum()
    mean=float(weights@rule['x'])
    return dict(log_evidence=anchor+normalizer,branch1=float(weights[rule['branches']==1].sum()),
        x_mean=mean,x_sd=float(np.sqrt(weights@((rule['x']-mean)**2))),
        contrast_mean=float(weights@rule['contrast']),posterior_weight_sum=float(weights.sum()),
        gain_weight_normalization_error=max(gain_norm),nuisance_history=history,
        nuisance_qualified=bool(history[-1]['max_logZ_change'] is not None and history[-1]['max_logZ_change']<=2e-4))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch-size',type=int,default=16)
    args=parser.parse_args(argv)
    if not 1<=args.batch_size<=32:raise ValueError('batch-size must be 1..32 to bound temporary arrays')
    cpu,wall=time.process_time(),time.perf_counter()
    dest=ROOT/'AUDIT/array_replay'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    dest.mkdir(parents=True,exist_ok=False)
    allowed=[f'results/P1A/rules/full_{n}.npz' for n in ORDERS]
    allowed += [f'results/P1A/rules/cdf_{i}_{n}.npz' for i in range(5) for n in ORDERS]
    allowed += ['results/P1A/data/OBSERVATIONS.npz','results/P1A/reference/RESULTS.json']
    manifest=dict(read_whitelist=allowed,truth_access=False,new_rhs=0,Book_constructed=False,
        physics_executed=False,original_outputs_modified=False,shared_gamma_integrations=1,
        amplitude_prior=[.99,1.01],phase_prior=[-.01,.01],sigma_per_real_quadrature=.01,
        Gaussian_real_dimensions=480,material_reference_orders=list(ORDERS),
        cdf_checkpoints=list(CUTS),nuisance_orders=list(NUISANCE),batch_size=args.batch_size,
        saved_value_absolute_tolerance=REPLAY_ATOL,
        scientific_GO_evaluated=False,source=str(Path(__file__).relative_to(ROOT)))
    write(dest/'MANIFEST.json',manifest)
    results=[];checks=[]
    try:
        with np.load(ROOT/allowed[-2],allow_pickle=False) as obs:
            ys=obs['y'].copy();parents=obs['parent_ids'].copy()
        saved=json.loads((ROOT/allowed[-1]).read_text())
        if ys.shape!=(32,5,4,12) or not np.isfinite(ys).all() or len(saved)!=32:raise ValueError('Expected 32 saved P1A observations')
        lookup={r['observation']:r for r in saved}
        if set(lookup)!=set(range(32)):raise ValueError('Saved reference IDs are not exactly 0..31')
        _,counts=np.unique(parents,return_counts=True)
        if len(counts)!=16 or not np.all(counts==2):raise ValueError('Expected 16 parents with two noise replicas each')
        # Only the twelve small response rules are retained; no field/J/LU arrays.
        rules={}
        for n in ORDERS:
            name=f'full_{n}';rules[name],check=load_rule(ROOT/f'results/P1A/rules/{name}.npz',n,None);checks.append(check)
        for i,cut in enumerate(CUTS):
            for n in ORDERS:
                name=f'cdf_{i}_{n}';rules[name],check=load_rule(ROOT/f'results/P1A/rules/{name}.npz',n,cut);checks.append(check)
        write(dest/'RULE_NORMALIZATION.json',checks)
        for index,y in enumerate(ys):
            s=lookup[index]
            if s['parent']!=int(parents[index]):raise ValueError('Parent identity mismatch')
            calculated={name:integrate(rule,y,args.batch_size) for name,rule in rules.items()}
            full={n:calculated[f'full_{n}'] for n in ORDERS}
            cdf={str(cut):{f'order{n}':float(np.exp(calculated[f'cdf_{i}_{n}']['log_evidence']-full[n]['log_evidence'])) for n in ORDERS} for i,cut in enumerate(CUTS)}
            for entry in cdf.values():entry['error']=abs(entry['order9']-entry['order5'])
            errors={k:abs(full[9][v]-float(s[k])) for k,v in (
                ('branch1_posterior_mass','branch1'),('log_evidence','log_evidence'),
                ('x_mean','x_mean'),('x_sd','x_sd'),('contrast_mean','contrast_mean'))}
            event_error=abs(full[9]['branch1']-full[5]['branch1']);cdf_error=max(r['error'] for r in cdf.values())
            errors['quadrature_event_error']=abs(event_error-s['quadrature_event_error'])
            errors['quadrature_CDF_max_error']=abs(cdf_error-s['quadrature_CDF_max_error'])
            for cut,entry in cdf.items():
                for k in ('order5','order9','error'):errors[f'CDF:{cut}:{k}']=abs(entry[k]-s['CDF'][cut][k])
            coherence=all(all(-1e-8<=x<=1+1e-8 for x in vals) and np.all(np.diff(vals)>=-1e-8)
                          for vals in ([cdf[str(cut)][f'order{n}'] for cut in CUTS] for n in ORDERS))
            qualified=all(v['nuisance_qualified'] for v in calculated.values()) and event_error<=.02 and cdf_error<=.02 and coherence
            weight_error=max(abs(v['posterior_weight_sum']-1.) for v in calculated.values())
            gain_error=max(v['gain_weight_normalization_error'] for v in calculated.values())
            pass_=max(errors.values())<=REPLAY_ATOL and weight_error<=2e-14 and gain_error<=2e-14 and qualified==(s['reference_status']=='QUALIFIED_EMPIRICAL_QUADRATURE')
            row=dict(observation=index,parent=int(parents[index]),numeric_readback_pass=bool(pass_),
                full=full,CDF=cdf,event_order_error=event_error,CDF_order_max_error=cdf_error,
                posterior_normalization_max_error=weight_error,gain_prior_normalization_max_error=gain_error,
                saved_value_errors=errors,saved_value_max_error=max(errors.values()),
                replay_reference_qualified=bool(qualified),saved_reference_status=s['reference_status'],
                CDF_coherent=bool(coherence),nuisance={k:v['nuisance_history'] for k,v in calculated.items()})
            results.append(row);write(dest/'RESULTS.json',results)
            print(json.dumps(dict(observation=index,numeric_readback_pass=bool(pass_),max_error=max(errors.values()),new_rhs=0)),flush=True)
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss_bytes=int(rss if sys.platform=='darwin' else rss*1024)
        ok=all(r['numeric_readback_pass'] for r in results) and rss_bytes<1_000_000_000
        summary=dict(status='ARRAY_NUMERIC_READBACK_PASS' if ok else 'ARRAY_NUMERIC_READBACK_FAIL',
            observations=len(results),rules=len(rules),all_saved_values_match=all(r['numeric_readback_pass'] for r in results),
            max_saved_value_absolute_error=max(r['saved_value_max_error'] for r in results),
            peak_RSS_bytes=rss_bytes,memory_below_1GB=rss_bytes<1_000_000_000,
            cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
            new_rhs=0,truth_access=False,Book_constructed=False,physics_executed=False,
            scope='Saved numerical readback only; not scientific GO or independent continuous PDE certification',
            audit_directory=str(dest.relative_to(ROOT)))
        write(dest/'SUMMARY.json',summary);print(json.dumps(summary),flush=True)
        return 0 if ok else 2
    except BaseException as exc:
        write(dest/'FAILURE.json',dict(error_type=type(exc).__name__,error=str(exc),observations_completed=len(results),
            new_rhs=0,truth_access=False,cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
        raise


if __name__=='__main__':raise SystemExit(main())

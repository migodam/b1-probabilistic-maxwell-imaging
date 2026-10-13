"""Single resumable owner of P2 physics, reference, and new truth generation."""
from __future__ import annotations
import os
import argparse,json,time,sys
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import truncnorm
from .accounting import Ledger,dump
from .material import BOUNDS,RMAX,WIDTH,decode,from_geometry,geometry,prior_sample,prior_event_masses,event_ids,material_distance,void_volume
from .forward import MaterialForward
from .response import nodes,save_bank,Response
from .likelihood import log_likelihood,context_loglik
from .reference import Reference

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/joint_p2';RESULTS=OUT/'results'

def cases():
    return [dict(id=name,latent=from_geometry(r,h,er,ei),gain=g) for name,r,h,er,ei,g in [
        ('solid',0,0,2.5,.1,1.),('small_cavity',.08,0,2.5,.1,1.),
        ('large_cavity',.25,0,2.5,.1,1.),('annular_cavity',0,.08,2.5,.1,1.),
        ('combined',.2,.06,2.5,.1,1.),('prior_boundary',.29,.11,2.21,.159,1.),
        ('gain_negative_control',.2,.06,2.5,.1,1.2)]]

def physics_bank(ledger,f):
    ledger.stage='mechanism';entries=[];rows=[];rng=np.random.default_rng(261013201)
    anchor=[]
    for k in [1.6,2.4]:anchor.append(f.predict(from_geometry(0,0),k,kind='FIXED_NOISE_ANCHOR')[0])
    sigma=np.array([.12*np.sqrt(np.mean(v*v)) for v in anchor])
    for item in cases():
        fields=[]
        for i,k in enumerate([1.6,2.4]):
            primary=f.predict(item['latent'],k,kind='STRESS_TRUTH_AND_MULTIPOLE')[0]
            lower=f.predict(item['latent'],k,kind='MULTIPOLE10',nmax=10)[0]
            higher=f.predict(item['latent'],k,kind='MULTIPOLE18',nmax=18)[0]
            recursive=f.predict(item['latent'],k,kind='INDEPENDENT_ADMITTANCE',implementation='recursive')[0]
            _,_,_,info=f.predict(item['latent'],k)
            rows.append(dict(observation=item['id'],k=k,n10_n14_relative=float(np.linalg.norm(lower-primary)/np.linalg.norm(primary)),n14_n18_relative=float(np.linalg.norm(higher-primary)/np.linalg.norm(primary)),n14_n18_whitened=float(np.linalg.norm(higher-primary)/sigma[i]),recursive_relative=float(np.linalg.norm(recursive-primary)/np.linalg.norm(primary)),boundary_residual=info['boundary_backward_residual'],Qabs=info['Qabs']))
            fields.append(primary)
        y=np.concatenate([item['gain']*v+sigma[i]*rng.normal(size=32) for i,v in enumerate(fields)])
        entries.append({**item,'y':y,'F':np.concatenate(fields),'sigma':np.repeat(sigma,32),'target_exposure':'preregistered development stress case'})
    dump(RESULTS/'PHYSICS_QUALIFICATION.json',dict(checks=rows,physical_dimension='3D vector concentric layered Maxwell',noise_sigma_per_frequency=sigma,new_budget=ledger.summary(),no_noise_tuning=True))
    dump(RESULTS/'DEVELOPMENT_CASES.json',dict(cases=entries,seed=261013201,ordering='frozen before posterior inspection'))
    make_bank(ledger,f)

def make_bank(ledger,f):
    ledger.stage='reference';axes=nodes(5);fields=np.empty((2,5,5,5,5,32))
    for a,k in enumerate([1.6,2.4]):
        for idx in np.ndindex((5,)*4):
            r,h,er,ei=[axes[j][idx[j]] for j in range(4)]
            x=from_geometry(r,h,er,ei)
            fields[(a,*idx)]=f.predict(x,k,kind='COMMON_RESPONSE_BANK')[0]
        print(json.dumps(dict(phase='bank',frequency=k,physics=ledger.count())),flush=True)
    save_bank(RESULTS/'COMMON_RESPONSE_BANK.npz',axes,fields)
    dump(RESULTS/'BANK_METADATA.json',dict(target='full-vector-Mie-matching-n14',surrogate_target='tensor-Chebyshev-five-per-physical-axis-radial-v1',axes=[a.tolist() for a in axes],physical_nodes_per_frequency=625,material_prior='coefficient-uniform; interpolation nodes are not discrete prior',freqs=[1.6,2.4],source='existing B1 Maxwell solver',actual_physics=ledger.summary(),uniform_bound=False))

def mechanism(ledger,f):
    ledger.stage='mechanism';response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');data=json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text());sigma=np.array(data['cases'][0]['sigma']);xa=from_geometry(0,0)
    fulla=np.concatenate([f.predict(xa,k)[0] for k in [1.6,2.4]])
    records=[]
    # Only geometry varies between worlds; material compensation is explicit.
    for item in cases()[1:6]:
        xb=item['latent'];raw=np.concatenate([f.predict(xb,k)[0] for k in [1.6,2.4]])
        for frequencies in [(1.6,),(1.6,2.4)]:
            n=32*len(frequencies);fa=fulla[:n];sd=sigma[:n]
            for permit in ['epsilon','epsilon_gain']:
                def prediction(t):
                    x=xb.copy();x[2:]=t[:2];fb=response.predict(x,frequencies)[0]
                    gain=t[2] if len(t)==3 else 1.
                    return gain*fb
                def residual(t):return (prediction(t)-fa)/sd
                lo=[2.2,.04]+([.7] if permit=='epsilon_gain' else []);hi=[2.8,.16]+([1.3] if permit=='epsilon_gain' else [])
                opt=least_squares(residual,[2.5,.1]+([1.] if permit=='epsilon_gain' else []),bounds=(lo,hi),xtol=1e-12,ftol=1e-12,gtol=1e-12,max_nfev=80)
                xc=xb.copy();xc[2:]=opt.x[:2];gain=float(opt.x[2]) if len(opt.x)==3 else 1.
                fb=np.concatenate([f.predict(xc,k,kind='COMPENSATION_SELECTED_FULL_ENDPOINT')[0] for k in frequencies])
                d0=float(np.linalg.norm((raw[:n]-fa)/sd));d=float(np.linalg.norm((gain*fb-fa)/sd))
                records.append(dict(geometry=item['id'],frequencies=frequencies,compensation=permit,raw_distance=d0,full_compensated_distance=d,absorbed_fraction=1-d/d0 if d0 else 0.,optimized_epsilon=opt.x[:2],gain=gain,material_distance=material_distance(xa,xc),gain_prior_z=(gain-1)/.15,minimum_status='feasible best-found upper bound; surrogate-guided optimizer, selected full-Mie endpoint',surrogate_predicted_distance=float(np.linalg.norm(opt.fun)),proxy_evaluations=opt.nfev,epsilon_boundary_fraction=np.minimum((opt.x[:2]-BOUNDS[2:,0])/WIDTH[2:],(BOUNDS[2:,1]-opt.x[:2])/WIDTH[2:])))
    # Local white Jacobian at one physical anchor. Gain derivative is analytic.
    J=[]
    for j,h in [(2,1e-4),(3,1e-5)]:
        xp=xa.copy();xm=xa.copy();xp[j]+=h;xm[j]-=h
        fp=np.concatenate([f.predict(xp,k,kind='EPSILON_CENTRAL_FD')[0] for k in [1.6,2.4]])
        fm=np.concatenate([f.predict(xm,k,kind='EPSILON_CENTRAL_FD')[0] for k in [1.6,2.4]])
        J.append((fp-fm)/(2*h)/sigma)
    A=np.column_stack([*J,fulla/sigma]);q=np.linalg.svd(A,full_matrices=False)[0];local=[]
    for r in [.02,.04,.08,.16,.25]:
        x=from_geometry(r,0);fb=np.concatenate([f.predict(x,k,kind='SMALL_CAVITY_FULL_SECANT')[0] for k in [1.6,2.4]])
        s=(fb-fulla)/sigma
        local.append(dict(radius=r,raw=float(np.linalg.norm(s)),projected=float(np.linalg.norm(s-q@(q.T@s))),per_single_frequency=float(np.linalg.norm(s[:32]))))
    dump(RESULTS/'COMPENSATION.json',dict(records=records,local_projection=local,anchor=xa,shared_calibration=True,prior='bounded proper gain [.7,1.3]; no both-world gain shrink-to-zero',profile_is_not_evidence=True,physical_dimension='3D vector Mie',new_budget=ledger.summary()))

def independent_forward(ledger,f):
    ledger.stage='fresh';response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');rng=np.random.default_rng(261013202)
    x=prior_sample(rng,112);full=[];pred=[];rows=[]
    data=json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text());sigma=np.array(data['cases'][0]['sigma'])
    for i,xx in enumerate(x):
        ff=np.concatenate([f.predict(xx,k,kind='INDEPENDENT_PRIOR_RESPONSE_CHECK',cache=False)[0] for k in [1.6,2.4]])
        pp=response.predict(xx,[1.6,2.4])[0];full.append(ff);pred.append(pp)
        rows.append(dict(point=i,relative_field_error=float(np.linalg.norm(ff-pp)/np.linalg.norm(ff)),whitened_error=float(np.linalg.norm((ff-pp)/sigma))))
    full=np.array(full);pred=np.array(pred);changes={}
    for c in data['cases']:
        y=np.array(c['y']);sd=np.array(c['sigma']);changes[c['id']]={}
        for condition,n in [('C',32),('D',64)]:
            a=log_likelihood(full[:,:n],y[:n],sd[:n]);b=log_likelihood(pred[:,:n],y[:n],sd[:n]);delta=a-b
            changes[c['id']][condition]=dict(max_abs_loglik=float(abs(delta).max()),median_abs_loglik=float(np.median(abs(delta))),weighted_mean_delta=float(np.exp(a-np.logaddexp.reduce(a))@delta))
    np.savez_compressed(RESULTS/'INDEPENDENT_FORWARD_ARRAYS.npz',latent=x,F_full=full,F_proxy=pred,sigma=sigma)
    dump(RESULTS/'INDEPENDENT_FORWARD_AUDIT.json',dict(seed=261013202,points=112,new_material_frequency_evaluations=224,checks=rows,max_relative=float(max(r['relative_field_error'] for r in rows)),max_whitened=float(max(r['whitened_error'] for r in rows)),likelihood_changes=changes,uniform_certificate=False,prior_probes_not_confirmation_of_all_modes=True))

def pilot_reference(ledger,f):
    response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);data=json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text());rows={}
    widths=np.array([*RMAX,.6,.12]);rng=np.random.default_rng(261013203)
    for item in data['cases']:
        ledger.check_time();y=np.array(item['y']);sigma=np.array(item['sigma']);name=item['id'];case={}
        for condition in ['A','B','C','D']:
            eps=np.array(item['latent'])[2:] if condition=='A' else None
            yy=y if condition=='D' else y[:32];sd=sigma if condition=='D' else sigma[:32]
            results=[ref.solve(yy,sd,order=q,condition=condition,known_eps=eps) for q in [9,13,17]]
            compact=[{k:r[k] for k in ['event_mass','log_event_mass','mean','SD','covariance','conditional','log_evidence','order']} for r in results]
            a,b=results[-2:];err=dict(event_mass=float(abs(a['event_mass']-b['event_mass']).max()),mean_fraction=float((abs(a['mean']-b['mean'])/widths).max()),SD_relative=float((abs(a['SD']-b['SD'])/np.maximum(b['SD'],widths*1e-6)).max()))
            case[condition]=dict(orders=compact,refinement=err,numerical_target='tensor-Chebyshev-five-per-physical-axis; full-Mie fidelity separate')
            if condition=='C':
                # Only the primary shared target needs the costly full CDF audit.
                s13=ref.continuous_summary(a,grid_size=9)
                s17=ref.continuous_summary(b,grid_size=9,conditional_events=np.flatnonzero(b['event_mass']>=.05))
                err['CDF']=float(abs(np.array(s13[0]['CDF'])-s17[0]['CDF']).max());err['quantile_fraction']=float((abs(np.array(s13[0]['quantiles_05_50_95'])-s17[0]['quantiles_05_50_95'])/widths[:,None]).max())
                case[condition]['continuous']=s17
                adaptive=ref.adaptive(b,max_cells=512);case[condition]['independent_adaptive']=adaptive
                case[condition]['adaptive_difference']=dict(event=float(abs(adaptive['event_mass']-b['event_mass']).max()),mean_fraction=float((abs(adaptive['mean']-b['mean'])/widths).max()),SD_relative=float((abs(adaptive['SD']-b['SD'])/np.maximum(b['SD'],widths*1e-6)).max()))
                particles=ref.particles(b,8192,rng)
                logp=log_likelihood(response.predict(particles),yy,sd)-np.log(WIDTH).sum()-b['log_evidence']
                np.savez_compressed(RESULTS/f'{name}_REFERENCE_PARTICLES.npz',latent=particles,weights=np.full(len(particles),1/len(particles)),reference_log_density=logp,y=yy,sigma=sd,event_mass=b['event_mass'],target_identity=np.array('joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1'))
        rows[name]=case;dump(RESULTS/'PILOT_REFERENCE.json',dict(cases=rows,complete=False,empirical=True,uniform_certificate=False))
        print(json.dumps(dict(case=name,C_event_mass=case['C']['orders'][-1]['event_mass'].tolist(),C_refinement=case['C']['refinement'],adaptive=case['C']['adaptive_difference'])),flush=True)
    dump(RESULTS/'PILOT_REFERENCE.json',dict(cases=rows,complete=True,empirical=True,uniform_certificate=False,polynomial_evaluations=response.evaluations))

def split_observations(ledger,f,split,n,seed):
    ledger.stage='fresh';rng=np.random.default_rng(seed);x=prior_sample(rng,n);meta=json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text());sigma=np.array(meta['cases'][0]['sigma']);fields=[];y=[];truth=[];parents=[];gains=[]
    for i,xx in enumerate(x):
        ff=np.concatenate([f.predict(xx,k,kind=f'{split.upper()}_PARENT',cache=False)[0] for k in [1.6,2.4]]);fields.append(ff)
        for noise in range(2):
            gain=float(truncnorm.ppf(rng.random(),-2,2,loc=1,scale=.15));yy=gain*ff+sigma*rng.normal(size=64)
            y.append(yy);truth.append(xx);parents.append(i);gains.append(gain)
    p=RESULTS/f'{split.upper()}_OBSERVATIONS.npz'
    np.savez_compressed(p,y=np.array(y),truth=np.array(truth),parent=np.array(parents),gain=np.array(gains),F_parent=np.array(fields),parent_latent=x,sigma=sigma,seed=np.array(seed))
    dump(RESULTS/f'{split.upper()}_MANIFEST.json',dict(split=split,parents=n,noise_per_parent=2,seed=seed,new_physics=2*n,prior='joint-p2-radial-v1',truth='newly simulated in current authorization; not historical locked',exposure='training/validation development' if split!='fresh_final' else 'created after weights frozen; no model selection'))
    return dict(y=np.array(y),truth=np.array(truth),parent=np.array(parents),sigma=sigma)

def training_data(ledger,f):
    response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);rng=np.random.default_rng(261013204);arrays={}
    for split,n,seed,order in [('train',128,261013205,9),('validation',32,261013206,13)]:
        data=split_observations(ledger,f,split,n,seed);sigma=data['sigma'][:32];teacher=[];checks=[]
        for i,yy in enumerate(data['y']):
            r=ref.solve(yy[:32],sigma,order=order);teacher.append(ref.particles(r,512,rng))
            if i<8:
                fine=ref.solve(yy[:32],sigma,order=17);checks.append(dict(observation=i,event_mass_max_difference=float(abs(r['event_mass']-fine['event_mass']).max()),mean_max_fraction=float((abs(r['mean']-fine['mean'])/np.array([*RMAX,.6,.12])).max())))
        arrays[f'{split}_y']=np.column_stack((data['y'][:,:32],np.full(len(data['y']),np.log(sigma[0])),np.full(len(data['y']),.15)))
        arrays[f'{split}_teacher']=np.array(teacher);arrays[f'{split}_parent']=data['parent']
        if split=='train':arrays['train_x']=data['truth']
        dump(RESULTS/f'{split.upper()}_TEACHER_AUDIT.json',dict(integration_order=order,fine_checks=checks,teacher_draws=512,teacher_type='randomized stratified quadrature resampling; inactive latent uniform conditional prior',teacher_is_additional_supervision=True,full_physics_setup='shared bank disclosed, no forward at training'))
    arrays['target_identity']=np.array('joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1')
    np.savez_compressed(RESULTS/'C1_TRAINING_DATA.npz',**arrays)
    print(json.dumps(dict(training_dataset=str(RESULTS/'C1_TRAINING_DATA.npz'),physics=ledger.count(),proxy_evaluations=response.evaluations)),flush=True)

def adaptive_refinement(ledger,f):
    response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);cases_data=json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text())['cases'];rows={}
    old=json.loads((RESULTS/'PILOT_REFERENCE.json').read_text())
    for item in cases_data:
        r=ref.solve(np.array(item['y'])[:32],np.array(item['sigma'])[:32],order=17)
        a=ref.adaptive_physical(r,max_cells=2048)
        widths=np.array([*RMAX,.6,.12])
        a['tensor17_difference']=dict(event=float(abs(a['event_mass']-r['event_mass']).max()),mean_fraction=float((abs(a['mean']-r['mean'])/widths).max()),SD_relative=float((abs(a['SD']-r['SD'])/np.maximum(r['SD'],widths*1e-6)).max()))
        a['original_coefficient_adaptive_preserved']=dict(converged=old['cases'][item['id']]['C']['independent_adaptive']['converged'],estimated_relative_error=old['cases'][item['id']]['C']['independent_adaptive']['estimated_relative_error'])
        rows[item['id']]=a;dump(RESULTS/'ADAPTIVE_PHYSICAL_REFERENCE.json',dict(cases=rows,complete=False,unchanged_prior_and_target=True))
        print(json.dumps(dict(case=item['id'],converged=a['converged'],cells=a['cells'],relative=a['estimated_relative_error'],difference=a['tensor17_difference'])),flush=True)
    dump(RESULTS/'ADAPTIVE_PHYSICAL_REFERENCE.json',dict(cases=rows,complete=True,unchanged_prior_and_target=True,polynomial_evaluations=response.evaluations))

def teacher_refinement(ledger,f):
    response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);rng=np.random.default_rng(261013207)
    old=np.load(RESULTS/'C1_TRAINING_DATA.npz');arrays={k:old[k] for k in old.files};data=np.load(RESULTS/'TRAIN_OBSERVATIONS.npz');sigma=data['sigma'][:32];teacher=[];checks=[]
    for i,y in enumerate(data['y']):
        r=ref.solve(y[:32],sigma,order=13);teacher.append(ref.particles(r,512,rng))
        if i<8:
            fine=ref.solve(y[:32],sigma,order=17);checks.append(dict(observation=i,event_mass_max_difference=float(abs(r['event_mass']-fine['event_mass']).max()),mean_max_fraction=float((abs(r['mean']-fine['mean'])/np.array([*RMAX,.6,.12])).max())))
    arrays['train_teacher']=np.array(teacher);np.savez_compressed(RESULTS/'C1_TRAINING_DATA_REF13.npz',**arrays)
    dump(RESULTS/'TRAIN_TEACHER_AUDIT_REF13.json',dict(integration_order=13,fine_checks=checks,change='refine teacher integration only; unchanged posterior target/prior/observations; no new physical truths',previous='TRAIN_TEACHER_AUDIT.json preserved',new_physics=0))
    print(json.dumps(dict(max_event=max(x['event_mass_max_difference'] for x in checks),max_mean_fraction=max(x['mean_max_fraction'] for x in checks),dataset=str(RESULTS/'C1_TRAINING_DATA_REF13.npz'))),flush=True)

def corrected_demo(ledger,f):
    from .correction import corrected_reference
    ledger.stage='reserve';response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');ref=Reference(response)
    adaptive=json.loads((RESULTS/'ADAPTIVE_PHYSICAL_REFERENCE.json').read_text())['cases'];data=json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text())['cases'];rows={}
    for name,n,seed in [('large_cavity',32,261013209),('solid',16,261013210)]:
        item=next(x for x in data if x['id']==name);r=ref.solve(np.array(item['y'])[:32],np.array(item['sigma'])[:32],order=17)
        summary,arrays=corrected_reference(r,adaptive[name],response,f,n,np.random.default_rng(seed))
        rows[name]=summary;np.savez_compressed(RESULTS/f'{name}_FULL_CORRECTION.npz',**arrays)
        dump(RESULTS/'FULL_MIE_CORRECTED_REFERENCE.json',dict(cases=rows,complete=False))
        print(json.dumps(dict(case=name,event_mass=summary['event_mass'].tolist(),max_physical_correction=summary['max_event_change'],max_MCse=float(max(summary['event_mass_MC_se'])))),flush=True)
    dump(RESULTS/'FULL_MIE_CORRECTED_REFERENCE.json',dict(cases=rows,complete=True,uniform_error_bound=False))

def generic_smc(ledger,f):
    from .smc import run_smc
    ledger.stage='reference';response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');item=next(x for x in json.loads((RESULTS/'DEVELOPMENT_CASES.json').read_text())['cases'] if x['id']=='large_cavity');rows=[]
    for seed in [261013211,261013212]:
        r,a=run_smc(response,f,np.array(item['y'])[:32],np.array(item['sigma'])[:32],seed,n=64,full_cap=340)
        rows.append(r);np.savez_compressed(RESULTS/f'SMC_{seed}.npz',**a)
        dump(RESULTS/'GENERIC_SMC.json',dict(observation_id='large_cavity',runs=rows,comparison_scope='same full-Mie numerical target as corrected reference; 64 particles, not asymptotic guarantee'))
        print(json.dumps(dict(seed=seed,status=r['status'],beta=r['beta'],full_calls=r['actual_full_evaluations'],event_mass=r['event_mass'].tolist())),flush=True)

def fresh_final(ledger,f):
    c1=Path(os.environ.get('C1_P2',ROOT.parent.parent/'C1/c1-joint-p2'));selected=c1/'p2/results/models/TRAINING_SUMMARY.json'
    if not selected.exists() or json.loads(selected.read_text())['status']!='MODELS_FROZEN_BEFORE_FRESH_TEST':raise RuntimeError('fresh final awaits frozen models')
    data=split_observations(ledger,f,'fresh_final',64,261013213);response=Response(RESULTS/'COMMON_RESPONSE_BANK.npz');ref=Reference(response);rng=np.random.default_rng(261013214);samples=[];logps=[];masses=[];checks=[]
    for i,y in enumerate(data['y']):
        r=ref.solve(y[:32],data['sigma'][:32],order=13);x=ref.particles(r,4096,rng);samples.append(x)
        logps.append(log_likelihood(response.predict(x),y[:32],data['sigma'][:32])-np.log(WIDTH).sum()-r['log_evidence']);masses.append(r['event_mass'])
        if i<8:
            fine=ref.solve(y[:32],data['sigma'][:32],order=17)
            checks.append(dict(observation=i,max_event_refinement=float(abs(r['event_mass']-fine['event_mass']).max()),mean_fraction=float((abs(r['mean']-fine['mean'])/np.array([*RMAX,.6,.12])).max())))
    yy=np.column_stack((data['y'][:,:32],np.full(len(data['y']),np.log(data['sigma'][0])),np.full(len(data['y']),.15)))
    np.savez_compressed(RESULTS/'C1_FINAL_DATA.npz',y=yy,truth=data['truth'],reference_samples=np.array(samples),reference_log_density=np.array(logps),reference_event_mass=np.array(masses),parent_ids=data['parent'],target_identity=np.array('joint-p2-radial-v1:C:tensor-Chebyshev-five-per-physical-axis-radial-v1'))
    dump(RESULTS/'FRESH_REFERENCE_AUDIT.json',dict(models_frozen_before_generation=True,parents=64,observations=128,new_physics=128,reference_target='same declared continuous proxy density; full-Mie residual correction separately tested on demos',independent_refinement_checks=checks,seed=261013213))
    print(json.dumps(dict(final_dataset=str(RESULTS/'C1_FINAL_DATA.npz'),physics=ledger.count())),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['physics_bank','make_bank','mechanism','independent_forward','pilot_reference','training_data','adaptive_refinement','teacher_refinement','corrected_demo','generic_smc','fresh_final']);parser.add_argument('--job-id');a=parser.parse_args()
    ledger=Ledger(RESULTS);f=MaterialForward(ROOT,ledger)
    with ledger.job(a.job_id or 'p2-'+a.phase):globals()[a.phase](ledger,f)
    dump(RESULTS/'RESOURCE_SUMMARY.json',ledger.summary());print(json.dumps(ledger.summary()),flush=True)

if __name__=='__main__':main()

"""Round2: density-explicit forward interpolation, separately audited.
No full-grid gate reclassification; alternative reference is empirical.
"""
from __future__ import annotations
import json,sys,time,subprocess
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import legvander2d
from .physics import Forward,topology_layers,angular_data
from .budget import Ledger
from .inference import infer,refinement,quadrature_nodes
from .pilot import dump,compact,save_case


def fit(bank,order):
    H=bank['H'];coeff=[];lo=bank['bounds'][:,0];hi=bank['bounds'][:,1]
    for h in range(3):
        z=2*(bank['theta'][H==h]-lo)/(hi-lo)-1
        V=legvander2d(z[:,0],z[:,1],[order-1,order-1]).reshape(-1,order*order)
        coeff.append(np.linalg.solve(V,bank['F'][H==h]))
    return np.asarray(coeff)


def predict(coeff,H,theta,bounds):
    theta=np.atleast_2d(theta);z=2*(theta-bounds[:,0])/(bounds[:,1]-bounds[:,0])-1
    order=int(np.sqrt(coeff.shape[1]));V=legvander2d(z[:,0],z[:,1],[order-1,order-1]).reshape(-1,order*order)
    return V@coeff[int(H)]


def interpolate_bank(coeff,order,c,bounds):
    theta,weights=quadrature_nodes(order,c);F=np.concatenate([predict(coeff,h,theta,bounds) for h in range(3)])
    return dict(H=np.repeat(np.arange(3),len(theta)),theta=np.tile(theta,(3,1)),weights=np.tile(weights,3),F=F,bounds=bounds)


def run(root,c):
    root=Path(root);out=root/'results/refinement';out.mkdir(exist_ok=True)
    if (out/'SUMMARY.json').exists():raise RuntimeError('round2 completed; do not rerun physical validation')
    led=Ledger(root/'results',c['physics_new_cap']);f=Forward(c,led,root/'results/forward_cache');start=time.perf_counter();before=led.count()
    b5=dict(np.load(root/'results/banks/order5.npz'));b7=dict(np.load(root/'results/banks/order7.npz'));c5=fit(b5,5);c7=fit(b7,7);bounds=b7['bounds']
    np.savez_compressed(out/'FORWARD_POLYNOMIAL.npz',coeff5=c5,coeff7=c7,bounds=bounds)
    lo=bounds[:,0];hi=bounds[:,1];rng=np.random.default_rng(c['seed']+700);design=[]
    for h in range(3):
        for i,j in [(0,0),(0,1),(1,0),(1,1)]:design.append(dict(H=h,theta=[float([lo[0],hi[0]][i]),float([lo[1],hi[1]][j])],type='prior_box_corner'))
        for _ in range(4):design.append(dict(H=h,theta=(lo+rng.uniform(size=2)*(hi-lo)).tolist(),type='fresh_uniform_development'))
    dump(out/'VALIDATION_MANIFEST.json',dict(seed=c['seed']+700,design=design,adaptive_development=True,new_physical_budget=28))
    checks=[];fullrecords=[]
    for item in design:
        h=item['H'];theta=item['theta'];full=f.run(h,*theta,kind='ROUND2_FORWARD_INTERPOLATION_AUDIT',cache=False)[0]
        pred=predict(c7,h,theta,bounds)[0];coarse=predict(c5,h,theta,bounds)[0]
        rel=float(np.linalg.norm(pred-full)/np.linalg.norm(full));ab=float(np.max(abs(pred-full)))
        checks.append(dict(**item,relative_field_error=rel,max_component_abs_error=ab,order5_relative_field_error=float(np.linalg.norm(coarse-full)/np.linalg.norm(full)),passed=rel<1e-5))
        fullrecords.append(dict(H=h,theta=theta,F=full.tolist(),F_interp7=pred.tolist()))
    # Four independent existing treams evaluations; reuse measured main outputs.
    # Even a subprocess/import failure has already been reserved and remains paid.
    independent=[]
    for idx in [0,8,16,19]:
        item=design[idx];h=item['H'];er,ei=item['theta'];rr,ee=topology_layers(h,complex(er,ei),c)
        params=dict(radii=rr.tolist(),eps_layers=[[z.real,z.imag] for z in ee],k=c['wavenumber'],nmax=14)
        def oracle():
            r=subprocess.run([sys.executable,str(root/'tools/treams_one.py'),json.dumps(params)],text=True,capture_output=True)
            if r.returncode:raise RuntimeError(r.stderr[-1500:])
            d=json.loads(r.stdout);a=np.array([complex(*v) for v in d['a']]);b=np.array([complex(*v) for v in d['b']]);return angular_data(a,b,c['angles_degrees'],c['wavenumber']),d
        ref,receipt=led.evaluate('EXISTING_A15_TREAMS_INDEPENDENT',dict(H=h,theta=item['theta'],**params),oracle)
        full=np.array(fullrecords[idx]['F']);err=float(np.linalg.norm(full-ref)/np.linalg.norm(ref))
        independent.append(dict(**item,relative_field_error=err,passed=err<1e-10,environment=receipt['environment'],oracle_coefficients_a=receipt['a'],oracle_coefficients_b=receipt['b'],oracle_data=ref.tolist()))
    dump(out/'FULL_FORWARD_CHECKS.json',dict(checks=checks,independent_treams=independent,full_outputs=fullrecords,
              empirical_field_error_not_uniform_certificate=True,actual_new_calls=led.count()-before))
    cb517=interpolate_bank(c5,17,c,bounds);cb717=interpolate_bank(c7,17,c,bounds);cb725=interpolate_bank(c7,25,c,bounds);reports={}
    for folder in sorted((root/'results/cases').iterdir()):
        if not folder.is_dir():continue
        p=json.loads((folder/'probabilities.json').read_text());obs=np.load(folder/'posterior.npz');y=obs['data'];sigma=float(obs['noise_real_sd']);pi=obs['topology_prior'];sd=.15 if folder.name=='gain_marginalized' else 0.
        r517=infer(cb517,y,sigma,pi,sd);r717=infer(cb717,y,sigma,pi,sd);r725=infer(cb725,y,sigma,pi,sd)
        quad=refinement(r717,r725,c);model=refinement(r517,r717,c)
        ok=all(q['passed'] for q in checks) and all(q['passed'] for q in independent) and quad['status']=='NUMERICAL_REFERENCE_PASS' and model['status']=='NUMERICAL_REFERENCE_PASS'
        status='EMPIRICALLY_REFINED_REFERENCE' if ok else 'ALTERNATIVE_REFERENCE_UNRESOLVED'
        reference=dict(status=status,integration_17_vs_25=quad,forward_interpolation_order5_vs7=model,independent_forward_tests=28,
                       certified_uniform_remainder=False,target='approximation to original continuous full-Mie posterior, 7x7 polynomial forward; interpolation error empirically audited, not certified')
        name=folder.name
        save_case(root,'refined_'+name,c,cb725,y,sigma,pi,r725,reference,gain_sd=sd)
        # At independent full points, report actual likelihood change; no bound extrapolation.
        local_defects=[]
        for rec in fullrecords:
            ll,_=__import__('src.b1.inference',fromlist=['log_likelihood']).log_likelihood(np.array(rec['F'])[None],y,sigma,sd)
            la,_=__import__('src.b1.inference',fromlist=['log_likelihood']).log_likelihood(np.array(rec['F_interp7'])[None],y,sigma,sd)
            local_defects.append(float(abs(ll[0]-la[0])))
        reports[name]=dict(result=compact(r725),alternative_reference=reference,original_full_grid_reference=p['reference'],
                            topology_change_from_original_order7=float(max(abs(np.array(p['result']['topology_probabilities'])-r725['topology_probabilities']))),
                            tested_log_likelihood_error_max=max(local_defects),unseen_tail_error_certified=False)
    summary=dict(cases=reports,new_physics=led.count()-before,resources=led.summary(),wall_seconds=time.perf_counter()-start,
                 max_empirical_relative_field_error=max(q['relative_field_error'] for q in checks),max_independent_treams_relative_field_error=max(q['relative_field_error'] for q in independent),
                 integration_response_evaluations=dict(poly5_order17=3*17**2,poly7_order17=3*17**2,poly7_order25=3*25**2),
                 status='empirical reference refinement; original full-grid failures retained',reference_bias_remaining='no uniform interpolation remainder or out-of-template mass certification')
    dump(out/'SUMMARY.json',summary);dump(root/'results/RESOURCE_SUMMARY.json',led.summary())
    print(json.dumps(dict(cases={n:r['result']['topology_probabilities'] for n,r in reports.items()},statuses={n:r['alternative_reference']['status'] for n,r in reports.items()},resources=led.summary(),max_error=summary['max_empirical_relative_field_error'])))

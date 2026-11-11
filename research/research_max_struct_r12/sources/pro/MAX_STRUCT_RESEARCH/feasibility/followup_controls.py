"""Posthoc mechanism controls after the initial four-case toy, not new test data.
1. Start DR at the direct solver's actual terminal point, two declared failures.
2. Two fixed smooth-random direct initializations on all four cases.
No GT is used for online choices. All controls and their costs are retained.
"""
from pathlib import Path
import json,time
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
from scipy.stats import chi2
from consensus_reflection_toy import setup,Physics,Projectors,stats

OUT=Path(__file__).parent
xy,G,E,S=setup()
names=['ring_high','split_shell_high','bridge_high','disk_low']
results={'status':'POSTHOC_MECHANISM_CONTROL','random_start_protocol':'two fixed seeds, smooth iid Uniform(0,12), smoothing sigma .8, each max600 direct steps; min physical data loss selection', 'cases':[]}
for name in names:
    data=np.load(OUT/(name+'.npz'))
    truth,Y,yclean=data['truth'],data['Y'],data['yclean']
    sigma=.01*np.linalg.norm(yclean)/np.sqrt(yclean.size)
    delta=sigma*np.sqrt(chi2.ppf(.99,2*Y.size)/2)
    entry={'case':name,'random_starts':[]}
    saved={}
    for seed in [9401,9402]:
        m0=gaussian_filter(np.random.default_rng(seed).uniform(0,12,(10,10)),.8).ravel()
        p=Physics(G,E,S,Y);tt=time.perf_counter()
        sol=minimize(lambda m:p.forward(m,True),m0,jac=True,method='L-BFGS-B',bounds=[(0,12)]*len(m0),
             options={'maxiter':600,'maxls':30,'ftol':1e-14,'gtol':1e-10})
        f,_,_=p.forward(sol.x)
        entry['random_starts'].append({'seed':seed,**stats(sol.x,truth,f),'seconds':time.perf_counter()-tt,
              'lu_count':p.lu_count,'rhs_count':p.rhs_count,'iterations':int(sol.nit)})
        saved['random_'+str(seed)]=sol.x
    entry['selected_random_seed']=min(entry['random_starts'],key=lambda row:row['relative_data_l2'])['seed']
    if name in names[:2]:
        m0=data['direct']
        p=Physics(G,E,S,Y);tt=time.perf_counter();pr=Projectors(G,E,S,Y,delta)
        f,J,Q=p.forward(m0);X=J/pr.scale
        best=(f,m0.copy(),0);hist=[]
        for it in range(1,3001):
            cx,cq,cm=pr.material(X,Q)
            ax,aq=pr.affine_ball(2*cx-X,2*cq-Q)
            gap=float(np.hypot(np.linalg.norm(ax-cx),np.linalg.norm(aq-cq)))
            X+=ax-cx;Q+=aq-cq
            if it%50==0:
                ff,_,_=p.forward(cm)
                if ff<best[0]:best=(ff,cm.copy(),it)
                hist.append([it,ff,gap])
        ff,mm,at=best
        entry['dr_from_direct_terminal']={**stats(mm,truth,ff),'seconds':time.perf_counter()-tt,
              'lu_count':p.lu_count,'rhs_count':p.rhs_count,'selected_iteration_by_data':at,'last_gap':hist[-1][2],
              'fixed_background_cholesky':1,'fixed_background_rhs':24+3000*6}
        saved['dr_from_direct_terminal']=mm
        saved['dr_terminal_history']=np.asarray(hist)
        # Endpoint chord diagnostic uses algorithmically obtained endpoints, no GT.
        ts=np.array([0,.001,.01,.05,.1,.25,.5,.75,1.])
        path=[]
        f0,g0=p.forward(m0,True)
        for t in ts:
            ft,_,_=p.forward(m0+t*(mm-m0))
            path.append([float(t),float(ft)])
        entry['endpoint_chord']={'path':path,'start_directional_derivative':float(g0@(mm-m0)),
              'min_gradient_at_direct_terminal':float(g0.min()),'diagnostic_extra_lu':10,'diagnostic_extra_rhs':66}
        saved['endpoint_chord']=np.asarray(path)
    np.savez_compressed(OUT/(name+'_controls.npz'),**saved)
    results['cases'].append(entry)
    (OUT/'toy_controls.json').write_text(json.dumps(results,indent=2))
    print(json.dumps(entry),flush=True)

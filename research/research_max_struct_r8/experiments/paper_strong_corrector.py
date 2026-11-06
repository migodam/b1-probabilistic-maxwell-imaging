"""Paid long physical continuation; preserve early80 results unchanged.

No GT loaded. All candidates selected by same original Phi. Full incumbent
materials and times at already-paid evaluations permit offline budget curves.
"""
import os
for k in ('OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r8.experiments.paper_benchmark import ROOT,OUT,CASES,ARMS,write
from research.research_max_struct_r8.src.paper_physics import acquisition,make_model

def run():
    jobs=[(case['id'],'physical_lbfgs') for case in CASES]+[('austria_eps4_400MHz',a) for a in ('nn_cie_dual','nn_a_plus','cie_dual','a_plus','a_plus_sourcewise','aug_nn_cie_dual','aug_nn_a_plus','aug_nn_a_plus_sourcewise')]
    for case,arm in jobs:
        source=OUT/case/arm;dest=OUT/'strong_corrector'/case/arm
        if (dest/'PREDICTION_SEAL.json').exists():continue
        dest.mkdir(parents=True,exist_ok=True);d=np.load(OUT/case/'PUBLIC_DATA.npz');side=int(d['inverse_side']);Y=d['Y']
        model=make_model(acquisition(side,float(d['frequency'])),prior=not arm.startswith('paper_'));m=np.load(source/'PREDICTION.npz')['m'];best=[m.copy(),model.physical_score(m,Y)]
        c,t=time.process_time(),time.perf_counter();history=[];materials=[];accepted=[]
        def fg(z):
            v,g=model.data_value_gradient(z,Y)
            if v<best[1]:best[:]=[z.copy(),v]
            history.append(dict(evaluation=len(history)+1,Phi=v,best_Phi=best[1],cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,LU=model.cost.get('LU',0),complex_RHS=model.cost.get('actual_complex_RHS',0),gradient_inf=float(np.linalg.norm(g,np.inf))))
            if len(history)==1 or len(history)%10==0:materials.append((len(history),best[0].astype(np.float32).copy()))
            return 1e6*v,1e6*g
        def callback(z):accepted.append(dict(iteration=len(accepted)+1,evaluation=len(history),Phi=history[-1]['Phi']))
        opt=minimize(fg,m,jac=True,method='L-BFGS-B',bounds=[(0,20)]*model.n,callback=callback,options=dict(maxiter=1000,maxfun=4000,gtol=.01,ftol=1e-12,maxls=30))
        final=best[0];st=model.state(final);_,gradient=model.data_value_gradient(final,Y)
        pg=np.where(((final<=0)&(gradient>0))|((final>=20)&(gradient<0)),0,gradient)
        np.savez_compressed(dest/'PREDICTION.npz',m=final,Y_prediction=st['Y'],budget_material=np.array([x[1] for x in materials]),budget_evaluations=np.array([x[0] for x in materials]))
        write(dest/'HISTORY.json',dict(evaluations=history,accepted=accepted));write(dest/'RECEIPT.json',dict(status='COMPLETED',case=case,arm=arm,additional_accepted_iterations=int(opt.nit),additional_evaluations=int(opt.nfev),success=bool(opt.success),stop=str(opt.message),cpu_seconds=time.process_time()-c,wall_seconds=time.perf_counter()-t,physical_cost=dict(model.cost),final_Phi=best[1],relative_data_residual=float(np.linalg.norm(st['Y']-Y)/np.linalg.norm(Y)),terminal_projected_gradient_inf=float(np.linalg.norm(pg,np.inf)),GT_loaded=False,qualification='1000-step continuation from sealed early80 endpoint; same full pixels/box/prior; originalR7scale/ftol/gtol'))
        write(dest/'PREDICTION_SEAL.json',dict(status='SEALED_BEFORE_OFFLINE_GT',time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),GT_loaded=False));print('STRONG',case,arm,round(time.perf_counter()-t,2),opt.nit,str(opt.message),flush=True)

if __name__=='__main__':run()

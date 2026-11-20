"""Design material witnesses from old weak material modes, then test nonlinear fine-grid fields."""
import json
from time import perf_counter, process_time
import numpy as np
from scipy.linalg import eigh, solve_triangular
from threadpoolctl import threadpool_limits
from material_witness import OUT, ARMS, PAIRS, FREQUENCIES, defect_patch, physical_patch, gram, gap
from designs import SIGMA_FRACTION


def modes(g_old,g_new,penalty):
    chol=np.linalg.cholesky(penalty)
    w=solve_triangular(chol.T,np.eye(len(chol)),lower=False)
    eig,U=eigh(w.T@g_old@w);base=w@U
    # Evaluate successive weak subspaces. Each direction maximizes the other schedule's
    # local signal within that subspace under the same spatial roughness penalty.
    sizes=sorted(set([1,2,4,8,12,18,24,32,40,50,64,80,len(eig)]))
    for k in sizes:
        if k>len(eig):continue
        b=base[:,:k];z=b@eigh(b.T@g_new@b,subset_by_index=(k-1,k-1))[1][:,-1]
        yield k,z/np.max(abs(z))


def main():
    wall,cpu=perf_counter(),process_time();ids,w,pen=defect_patch()
    with np.load(OUT/'acquisition.npz') as arr:sd=float(arr['reference_amplitude'])*SIGMA_FRACTION
    terms={n:{f:physical_patch(n,f,ids,w) for f in FREQUENCIES} for n in (24,48)}
    results=[]
    for pair in PAIRS:
        for hard,easy in (pair,pair[::-1]):
            da,db=ARMS[hard],ARMS[easy];used=set(da.freqs+db.freqs)
            pick24={f:terms[24][f] for f in used};pick48={f:terms[48][f] for f in used}
            ga,gb=gram(pick24,da,sd),gram(pick24,db,sd)
            records=[];directions=[]
            for k,z in modes(ga,gb,pen):
                a=min(.625,.8/max(2*np.sqrt(max(0,z@ga@z)),1.e-12))
                for _ in range(3):
                    pred,_=gap(pick24,False,z,w,a,sd);meas=pred(da)['gap_over_sigma']
                    if a==.625 and meas<=.8:break
                    a=min(.625,a*.8/max(meas,1.e-12))
                p,_=gap(pick24,False,z,w,a,sd);ha,eb=p(da),p(db)
                size=float(np.linalg.norm(2*a*w*z)/np.sqrt(len(z)))
                records.append(dict(mode_count=k,amplitude=float(a),patch_rms=size,hard=ha,easy=eb));directions.append(z)
            eligible=[i for i,v in enumerate(records) if v['hard']['equal_prior_error']>=.20 and v['easy']['equal_prior_error']<=.05]
            if eligible:chosen=max(eligible,key=lambda i:records[i]['patch_rms'])
            else:chosen=max(range(len(records)),key=lambda i:records[i]['hard']['equal_prior_error']-records[i]['easy']['equal_prior_error'])
            z=directions[chosen];v=records[chosen];a=v['amplitude'];p,err=gap(pick48,True,z,w,a,sd)
            g=2*a*w*z
            row=dict(hard_design=hard,easy_design=easy,selected_on_24=chosen,met_coarse_targets=bool(eligible),
                     patch_rms=v['patch_rms'],image_rms=float(np.linalg.norm(g)/24),max_cell_difference=float(np.max(abs(g))),
                     coarse_candidates=records,fine_hard=p(da),fine_easy=p(db),identity_rel=err,
                     pixel_ids=ids.tolist(),minus=(-.625*w-.5*g).tolist(),plus=(-.625*w+.5*g).tolist());results.append(row)
            print(hard,easy,'rms',round(v['patch_rms'],3),'peak',round(row['max_cell_difference'],3),'risks',row['fine_hard']['equal_prior_error'],row['fine_easy']['equal_prior_error'],flush=True)
    (OUT/'material_weak_modes.json').write_text(json.dumps(dict(cases=results,cpu_s=process_time()-cpu,wall_s=perf_counter()-wall,
                 acceptance='Coarse hard pair Pe >= .20, easy pair Pe <= .05; maximize local RMS among proposals; retain all tested directions and failures.'),indent=2)+'\n')

if __name__=='__main__':
    with threadpool_limits(1):main()

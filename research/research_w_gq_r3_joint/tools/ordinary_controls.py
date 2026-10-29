"""Measurement-only ordinary controls, with no private input path or sampler import."""
from pathlib import Path
import argparse,json,time,sys,resource
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[_thread_key]='1'
import numpy as np
from scipy.stats import truncnorm
ROOT=Path(__file__).resolve().parents[1];R2=ROOT.parent/'research_w_gq_r2'
sys.path.insert(0,str(ROOT/'frozen/src' if (ROOT/'frozen/src').exists() else ROOT/'src'))
for p in [ROOT/'src',R2/'src',R2/'reuse/previous_gq/src',R2/'reuse/fixed_vendor',ROOT.parents[2]/'research-w-continuation/.deps']:sys.path.append(str(p))
from r3.physics import JointModel,default_config
from r3.accounting import BudgetBook,plain
from r3.inverse import optimize,objective
from r3.io import Trace,save_endpoint

def target_control(model,y,out,base,method,options=None):
    options={} if options is None else options;trace=Trace(out);nwall=len(base)-7;seeds=[];candidates=[]
    for x in [-.6,0.,.6]:
        for yy in [1.8,2.1,2.4]:
            for rx,ry,ang in [(.12,.12,0.),(.15,.09,.35)]:
                z=base.copy();z[-7:]=[x,yy,rx,ry,3.,.06,ang];p=model.forward(z,label=method+' public coarse seed');f=objective(model,y,z,p);cid=trace.record(z,p,f,stage='coarse_initialization',kind='coarse');candidates.append((f,z,p,cid))
    candidates.sort(key=lambda r:r[0]);starts=options.get('starts',4)
    for row in candidates:
        if not seeds or all(np.linalg.norm((row[1][-7:-5]-s[1][-7:-5])/np.array([.5,.35]))>.35 for s in seeds):seeds.append(row)
        if len(seeds)>=starts:break
    if len(seeds)<starts:
        for row in candidates:
            if all(row[3]!=s[3] for s in seeds):seeds.append(row)
            if len(seeds)>=starts:break
    np.savez_compressed(out/'INITIALIZATION.npz',wall_z=base,seeds=np.asarray([s[1] for s in seeds]),coarse_objectives=np.asarray([s[0] for s in seeds]));endpoints=[]
    for k,row in enumerate(seeds):
        z,p,rec=optimize(model,y,row[1],trace,f'joint_start_{k}',active=np.arange(nwall,len(base)),max_evaluations=options.get('max_evaluations',120),max_accepted=options.get('max_accepted',100))
        endpoints.append((rec['objective'],z,p,rec));np.savez_compressed(out/f'START_{k}_ENDPOINT.npz',z=z,data=p.data)
    endpoints.sort(key=lambda r:r[0]);f,z,p,rec=endpoints[0]
    summary=dict(method=method,objective=f,selected_start=rec['stage'],starts=[e[3] for e in endpoints],data_chi2=float(np.sum(abs(p.data-y)**2)/model.config['sigma']**2),selection='minimum full measurement objective including frozen public prior',converged=rec['stop'] in ['PROJECTED_GRADIENT','RELATIVE_OBJECTIVE'] and rec.get('final_exact_projected_gradient',float('inf'))<.005,wall_frozen_during_target_fit=True,stationarity_scope='seven active target coordinates conditional on fixed wall; not joint 19-coordinate stationarity')
    save_endpoint(model,z,p,out,summary);return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('--method',choices=['PRIOR_ONLY','WALL_FIRST_FIXED'],required=True);p.add_argument('--observation',required=True,type=Path);p.add_argument('--config',type=Path);p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    cfg=default_config();cfg.update(json.loads(a.config.read_text()) if a.config else {});book=BudgetBook(root=ROOT);start=time.process_time();startwall=time.time();before=book.snapshot();a.out.mkdir(parents=True,exist_ok=True)
    if (a.out/'endpoint.npz').exists() or (a.out/'TRACE.jsonl').exists():raise FileExistsError('Control output already contains scientific artifacts')
    model=JointModel(cfg,book);y=np.load(a.observation)['y'];model.observation=y;z=np.asarray(model.spec['initial'],float);nwall=len(z)-7
    if a.method=='PRIOR_ONLY':
        n=nwall//2;z[:n]=truncnorm.mean((1-3)/1.5,(6-3)/1.5,loc=3,scale=1.5);z[n:nwall]=truncnorm.mean((0-.06)/.15,(.25-.06)/.15,loc=.06,scale=.15);lo=np.asarray(model.spec['lower']);hi=np.asarray(model.spec['upper']);z[nwall:]=(lo[nwall:]+hi[nwall:])/2
        pred=model.forward(z,label='prior parameter-mean representative');f=objective(model,y,z,pred);Trace(a.out).record(z,pred,f,stage='prior_parameter_mean',kind='prior')
        save_endpoint(model,z,pred,a.out,dict(method=a.method,objective=f,data_chi2=float(np.sum(abs(pred.data-y)**2)/cfg['sigma']**2),converged=None,selection='independent parameter prior means; no measurement optimization',limitation='parameter mean representative, not the ensemble mean material image or a predictive posterior'))
    else:
        trace=Trace(a.out);zw=z.copy();zw[-3:-1]=0.;wallz,_,rec=optimize(model,y,zw,trace,'wall_only_initialization',active=np.arange(nwall),max_evaluations=35,max_accepted=25,gtol=1e-3)
        # Avoid candidate filename collisions between initialization and target stage.
        init=a.out/'WALL_INITIALIZATION';init.mkdir()
        phasefiles=list(a.out.glob('*.npz'))+[a.out/'TRACE.jsonl',a.out/'STEPS.jsonl']
        for f in phasefiles:
            if f.exists():f.rename(init/f.name)
        summary=target_control(model,y,a.out,wallz,a.method,cfg.get('inverse',{}));summary['wall_initialization']=rec
        (a.out/'CONTROL_WALL_INIT.json').write_text(json.dumps(rec,default=plain,indent=2)+'\n')
    after=book.snapshot();receipt=dict(pid=os.getpid(),source_root=str(ROOT/'frozen/src' if (ROOT/'frozen/src').exists() else ROOT/'src'),peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),process_cpu_seconds=time.process_time()-start,wall_seconds=time.time()-startwall,physical_cost={k:after[k]-before[k] for k in before},method=a.method)
    (a.out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n');book.event('R3','control_job_completion',**{k:v for k,v in receipt.items() if k!='pid'})
if __name__=='__main__':main()

"""R3 serial entry point; thread policy is set before numerical imports."""
from pathlib import Path
import os,sys,json,time,resource,argparse
sys.dont_write_bytecode=True
os.environ['PYTHONDONTWRITEBYTECODE']='1'
ROOT=Path(__file__).resolve().parents[1];R2=ROOT.parent/'research_w_gq_r2'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'results/.mpl'))
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:os.environ[k]='1'
source_root=Path(os.environ.get('R3_SOURCE_ROOT',str(ROOT/'src')))
for p in reversed([source_root,R2/'src',R2/'reuse/previous_gq/src',R2/'reuse/fixed_vendor',ROOT.parents[2]/'research-w-continuation/.deps']):sys.path.insert(0,str(p))
import numpy as np
from r3.physics import JointModel,default_config,parameter_spec
from r3.accounting import BudgetBook,plain

def main():
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['benchmark','chart_benchmark','inverse','run']);ap.add_argument('--config');ap.add_argument('--observation','--measurement',dest='observation');ap.add_argument('--out','--output',dest='out');ap.add_argument('--method',default='FULL_GN');a=ap.parse_args()
    c=default_config();c.update(json.loads(Path(a.config).read_text()) if a.config else {});book=BudgetBook(root=ROOT);book.check(new_job=True);cpu=time.process_time();wall=time.time();out=Path(a.out) if a.out else ROOT/'results'/a.action;out.mkdir(parents=True,exist_ok=True)
    book.event('R3','job_start',action=a.action,config=c,output=str(out))
    try:
        model=JointModel(c,book)
        if a.action=='benchmark':
            z=np.asarray(parameter_spec(c)['initial']);z[-7:]=[.13,2.17,.137,.093,3.7,.075,.27];p=model.forward(z,jacobian=True,label='R3 necessary derivative audit base');rng=np.random.default_rng(331029);v=rng.normal(size=len(z))*np.asarray(model.spec['scale']);v/=np.linalg.norm(v/np.asarray(model.spec['scale']));h=1e-5
            pp=model.forward(z+h*v,label='R3 gradient audit plus');pm=model.forward(z-h*v,label='R3 gradient audit minus');fd=(pp.data-pm.data).ravel()/(2*h);jv=p.jacobian@v
            result=dict(relative_directional_error=float(np.linalg.norm(fd-jv)/np.linalg.norm(fd)),absolute_error=float(np.linalg.norm(fd-jv)),jv_norm=float(np.linalg.norm(jv)),fd_norm=float(np.linalg.norm(fd)),config=c,z=z,diagnostics=p.diagnostics,physical_cost=book.snapshot())
            (out/'BENCHMARK.json').write_text(json.dumps(result,default=plain,indent=2));print(json.dumps(result,default=plain),flush=True)
        elif a.action=='chart_benchmark':
            from r3.accounting import fem_lock
            z=np.asarray(model.spec['initial']);z[-7:]=[.13,2.17,.137,.093,3.7,.075,.27];model.initialize_charts(rank=8);pf=model.forward(z,jacobian=True,label='R3 chart necessary exact reference');pc=model.forward(z,jacobian='FIXED_CURRENT_CHART',label='R3 fixed rank8 audit');rng=np.random.default_rng(331029);v=rng.normal(size=len(z))*np.asarray(model.spec['scale']);v/=np.linalg.norm(v/np.asarray(model.spec['scale']));events,fallback=model.audit_chart_step(z,v,pc);
            with fem_lock(book.root):
                for ev,ch in zip(model.evaluators,model.charts):
                    ch.rebuild(ev._receiver_adjoint(ev._last,'R3 rank16 exactness audit adjoint'),rank=16,label='necessary current rank16 exactness audit')
            p16=model.forward(z,jacobian='FIXED_CURRENT_CHART',label='R3 rank16 projected Jacobian audit');result=dict(rank8_J_relative_error=float(np.linalg.norm(pc.jacobian-pf.jacobian)/np.linalg.norm(pf.jacobian)),rank16_J_relative_error=float(np.linalg.norm(p16.jacobian-pf.jacobian)/np.linalg.norm(pf.jacobian)),adaptive_audit=events,adaptive_fallback=fallback is not None,metric='D=diag of consistent FE mass on WALL_BOX union TARGET_BOX; no row-sum lumping',physical_cost=book.snapshot());(out/'CHART_BENCHMARK.json').write_text(json.dumps(result,default=plain,indent=2));print(json.dumps(result,default=plain),flush=True)
        else:
            from r3.inverse import invert
            obs=np.load(a.observation);y=obs['y'];invert(model,y,out,method=a.method,options=c.get('inverse',{}))
    except Exception as e:
        book.failure('R3',a.action,e);raise
    finally:
        artifact_alias_error=None
        if (out/'ONLINE_SEALED.json').exists():
            try:
                from output_artifacts import complete_output_aliases
                complete_output_aliases(out)
            except Exception as error:
                artifact_alias_error=str(error);book.failure('R3','postoutput_alias',error)
        receipt=dict(action=a.action,pid=os.getpid(),source_root=str(source_root),python_version=sys.version,blas_threads=1,artifact_alias_error=artifact_alias_error,process_cpu_seconds=time.process_time()-cpu,wall_seconds=time.time()-wall,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024),physical_cost=book.snapshot(),job_physical_cost=book.snapshot(pid=os.getpid()));book.event('R3','job_completion',**{k:v for k,v in receipt.items() if k!='pid'});(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2));(out/'runtime_cost.json').write_text(json.dumps(receipt,indent=2))

if __name__=='__main__':main()

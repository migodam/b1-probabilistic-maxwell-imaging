"""Same complete buffered trace policy, exact equations, stronger frozen-J control."""
import json,time
import numpy as np
from accounting import P,Book,write,identity
from physics import model_for,quality
from charts import MaterialChart,CurrentChart,partition

class TraceBook(Book):
    def __init__(self,run):super().__init__(run);self.buffer=[];self.snap=0
    def save(self,action,**arrays):
        self.buffer.append((str(action),{k:np.asarray(v).copy() for k,v in arrays.items()}))
        self.append(dict(event='PAID_ARRAYS_BUFFERED',action=action,array_bytes=sum(np.asarray(v).nbytes for v in arrays.values())))
    def algebra_snapshot(self,**arrays):
        self.snap+=1;self.buffer.append(('alg'+str(self.snap),{k:np.asarray(v).copy() for k,v in arrays.items()}))
    def flush(self,path):
        data={}
        for label,arrays in self.buffer:
            for key,a in arrays.items():data['s'+label+'_'+key]=a
        np.savez_compressed(path,**data)
        self.append(dict(event='COMPLETE_TRACE_BATCH_SAVED',path=path.relative_to(P).as_posix(),snapshots=len(self.buffer),bytes=path.stat().st_size,blake2b=identity(path)))
        self.buffer=[]

class TraceCurrent(CurrentChart):
    def residual_jac(self,z,free_values,need_jac=True):
        r,D,states=super().residual_jac(z,free_values,need_jac)
        m,j,E=states
        self.book.algebra_snapshot(z=z,free=free_values,m=m,j=j,E=E,r=r,D=np.array([]) if D is None else D)
        return r,D,states

class FrozenProxy:
    def __init__(self,model,J):self.model=model;self.J=J;self.book=model.book
    def __getattr__(self,name):return getattr(self.model,name)
    def jacobian(self,state=None):
        self.book.append(dict(event='FROZEN_J_CACHE_HIT',actual_complex_RHS=0));return self.J

def main():
    c=json.loads((P/'EXPERIMENT_CONTRACT.json').read_text());out=P/'results/COST_CONTROL';out.mkdir(parents=True,exist_ok=False)
    b=TraceBook('OPTIMIZED_EQUAL_TRACE_COST_CONTROL');rows=[]
    offsets=[np.array([a,d]) for a in (-.15,0.,.15) for d in (-.15,0.,.15)]
    methods=['DYNAMIC_FULLJ','FROZEN_J','CURRENT_FULL_KERNEL']
    try:
        for parent in range(3):
            for repeat in range(3):
                model=model_for(b);mref=np.array(c['E1']['public_reference_materials'][parent]);ref=model.forward(mref,jac=True,label='QA_COMMON_REFERENCE');dep,free=partition(ref['J']);target=ref['y']
                b.flush(out/f'COMMON_P{parent}_R{repeat}.npz')
                order=methods[(parent+repeat)%3:]+methods[:(parent+repeat)%3]
                for method in order:
                    tsetup=time.perf_counter();before=b.run_used
                    if method=='CURRENT_FULL_KERNEL':chart=TraceCurrent(model,target,ref,dep,free);initial=chart.encode(ref['j'],mref)
                    else:chart=MaterialChart(FrozenProxy(model,ref['J']) if method=='FROZEN_J' else model,target,mref,dep,free);initial=mref
                    setup=time.perf_counter()-tsetup;setup_rhs=b.run_used-before
                    for idx,offset in enumerate(offsets):
                        before=b.run_used;t=time.perf_counter()
                        if np.all(offset==0):
                            r=dict(success=True,m=mref,iterations=0,history=[]);b.append(dict(event='EXACT_REFERENCE_CACHE_HIT',actual_complex_RHS=0))
                        else:r=chart.correct(initial,mref[free]+offset,'QA_'+method)
                        solve_seconds=time.perf_counter()-t
                        t=time.perf_counter();validation=model.forward(r['m'],independent=True,label='QA_COMMON_FULL_VALIDATION');q=quality(model,validation,target,mref)
                        validation_seconds=time.perf_counter()-t
                        filename=f'P{parent}_R{repeat}_{method}_{idx}.npz';t=time.perf_counter();b.flush(out/filename);persist=time.perf_counter()-t
                        rows.append(dict(parent=parent,repeat=repeat,method=method,target=idx,qualified=q['qualified'],corrector_success=r['success'],
                            pair_relative=q['pair_relative'],material=r['m'].tolist(),RHS=b.run_used-before,
                            solver_seconds=solve_seconds,validation_seconds=validation_seconds,persistence_seconds=persist,total_seconds=solve_seconds+validation_seconds+persist,
                            setup_seconds=setup if idx==0 else 0.,setup_RHS=setup_rhs if idx==0 else 0.,trace=filename))
                    write(out/'CHECKPOINT.json',dict(cases=len(rows),paid=b.run_used,latest_method=method))
        grouped=[]
        for p in range(3):
            for method in methods:
                r=[v for v in rows if v['parent']==p and v['method']==method]
                grouped.append(dict(parent=p,method=method,qualified=sum(v['qualified'] for v in r),cases=len(r),RHS=sum(v['RHS'] for v in r),
                    total_seconds=sum(v['total_seconds']+v['setup_seconds'] for v in r),solver_plus_validation_seconds=sum(v['solver_seconds']+v['validation_seconds']+v['setup_seconds'] for v in r)))
        result=dict(status='COMPLETE',rows=rows,grouped=grouped,paid=b.run_used,
            equal_trace_policy='All real physics and algebraic current residual trial arrays retained losslessly; per-endpoint batch flush; paid BEGIN remains durable before every solve',
            measurements='three deterministic reference materials x nine same free targets x three repeats; repeated timing units are not independent parents',
            stronger_control='ordinary frozen material Jacobian under identical full nonlinear residual and endpoint gate',
            shared_reference_RHS=27,scope='Separate optimized persistence/timing window; original exploratory window is not overwritten')
        write(out/'SUMMARY.json',result);b.end();print(json.dumps(dict(status=result['status'],paid=b.run_used,grouped=grouped),indent=2))
    except BaseException as e:
        if b.buffer:b.flush(out/'INTERRUPTED_TRACE.npz')
        write(out/'INTERRUPTED.json',dict(error=repr(e),rows=rows,paid=b.run_used));b.end('FAILED',error=repr(e));raise
if __name__=='__main__':main()

from contextlib import contextmanager
from pathlib import Path
import hashlib,json,os,resource,time
import numpy as np

P=Path(__file__).resolve().parents[1]
def encode(x):
    if isinstance(x,(np.integer,np.floating)):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,Path):return str(x)
    raise TypeError(type(x).__name__)
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(value,ensure_ascii=False,indent=2,default=encode)+'\n');os.replace(temp,path)
def identity(path):return hashlib.blake2b(Path(path).read_bytes()).hexdigest()

class BudgetStop(RuntimeError):pass
class Book:
    def __init__(self,run):
        self.run=run;self.path=P/'COST_LEDGER.jsonl';self.actual=0;self.formula=0;self.run_used=0;self.events=0
        self.t0=time.perf_counter();self.c0=time.process_time();self.previous_CPU=0.
        if (P/'CLOSURE.json').exists():raise BudgetStop('Campaign closed; no reset or continuation allowed')
        repeated_run=False
        if self.path.exists():
            for line in self.path.read_text().splitlines():
                row=json.loads(line);self.events+=1
                if row['event']=='RUN_START' and row['run']==run:repeated_run=True
                if row['event']=='BEGIN':
                    self.actual+=row.get('actual_complex_RHS',0);self.formula+=row.get('closed_form_source_RHS_equivalent',0)
                    if row['run']==run:self.run_used+=row.get('actual_complex_RHS',0)+row.get('closed_form_source_RHS_equivalent',0)
                if row['event']=='RUN_END' and row['run']==run:
                    self.previous_CPU+=row.get('controlled_CPU_seconds',0.)
        if repeated_run:raise BudgetStop('Run names are single-use; preserve interrupted/finished runs and register any separate repair explicitly')
        self.append(dict(event='RUN_START',new_physics_RHS=0,BLAS_threads=1))
    def append(self,row):
        self.events+=1;row=dict(row,seq=self.events,run=self.run,utc=time.time())
        with self.path.open('a') as f:f.write(json.dumps(row,default=encode)+'\n')
    @contextmanager
    def op(self,kind,rhs=0,formula=0,**cost):
        if self.actual+self.formula+rhs+formula>30000 or self.run_used+rhs+formula>10000:raise BudgetStop('Hard actual+formula-equivalent cap before action')
        if self.previous_CPU+time.process_time()-self.c0>7200 or resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>8*1024**3:raise BudgetStop('CPU or memory cap')
        self.actual+=rhs;self.formula+=formula;self.run_used+=rhs+formula
        action=self.events+1;self.append(dict(event='BEGIN',kind=kind,action=action,actual_complex_RHS=rhs,closed_form_source_RHS_equivalent=formula,**cost))
        t=time.perf_counter();c=time.process_time()
        try:
            yield action
        except BaseException as exc:
            self.append(dict(event='END',kind=kind,action=action,status='FAILED_PAID',error=repr(exc),wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-c));raise
        else:self.append(dict(event='END',kind=kind,action=action,status='COMPLETE',wall_seconds=time.perf_counter()-t,CPU_seconds=time.process_time()-c))
    def save(self,action,**arrays):
        path=P/'raw'/self.run/(str(action).zfill(6)+'.npz');path.parent.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(path,**arrays)
        self.append(dict(event='SAVE',action=action,path=path.relative_to(P).as_posix(),bytes=path.stat().st_size,blake2b=identity(path),array_dtypes={k:str(np.asarray(v).dtype) for k,v in arrays.items()}))
    def end(self,status='COMPLETE',**extra):
        self.append(dict(event='RUN_END',status=status,run_combined_budget_used=self.run_used,controlled_wall_seconds=time.perf_counter()-self.t0,controlled_CPU_seconds=time.process_time()-self.c0,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,**extra))

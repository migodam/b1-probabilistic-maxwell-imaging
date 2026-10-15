"""Prior-constrained pilot bookkeeping and material descriptors; no solver import."""
from pathlib import Path
from contextlib import contextmanager
from collections import Counter
import json, time, os, fcntl
import numpy as np
from . import prior

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'research/prior_constrained_20261015'
def plain(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    if isinstance(x, dict): return {str(k): plain(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [plain(v) for v in x]
    return x
def dump(path, obj):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plain(obj), ensure_ascii=False, indent=2, allow_nan=False)+'\n')
def pack(z): return np.concatenate((z.real.T,z.imag.T),axis=1).reshape(-1)
def ell(pred,y,sigma): return float(np.vdot(pred-y,pred-y).real / sigma**2)

def descriptor(chi, bg, amplitude=1.2):
    d=prior.structure_descriptor(prior.rho_field(chi,bg),amplitude=amplitude)
    d['centered']=d['smoothed']-d['smoothed'].mean(axis=(-2,-1),keepdims=True)
    return d
def separation(a,b):
    delta=np.mean(a['masks']!=b['masks'],axis=(-2,-1))
    cal=json.loads((OUT/'delegated/prior/STRUCTURE_CALIBRATION.json').read_text())['recommended_minima_max_of_scales']
    r=dict(shape_multiscale_rms=float(np.sqrt(np.mean((a['centered']-b['centered'])**2))),
        occupancy_multithreshold_mean=float(delta.mean()),
        boundary_multithreshold_mean=float(np.mean(a['boundaries']!=b['boundaries'])),
        persistent_levels=int(np.sum(delta.mean(axis=0)>=.05)),
        topology_differing_levels=int(np.sum((a['components']!=b['components'])|(a['holes']!=b['holes']))),
        changed_cells=int(np.sum(np.abs(a['normalized_field']-b['normalized_field'])>=.1)),
        dc_difference=abs(a['dc']-b['dc']))
    r['structure_qualified']=bool(all(r[k]>=v for k,v in cal.items()) and r['persistent_levels']>=3 and r['changed_cells']>=6)
    return r

class Limit(RuntimeError): pass
class Book:
    """One append-only authorization, stage limits cannot be traded or reset."""
    CAPS={'E0':(0,0),'E1':(64,2048),'E2':(192,8192)}
    def __init__(self,stage):
        self.stage=stage; self.path=OUT/'COST_LEDGER.jsonl'; self.counts=Counter(); self.context={}
        OUT.mkdir(parents=True,exist_ok=True)
    def rows(self): return [json.loads(x) for x in self.path.read_text().splitlines()] if self.path.exists() else []
    def append(self,row):
        with self.path.open('a') as f:
            f.write(json.dumps(plain(dict(utc=time.time(),stage=self.stage,**row)),allow_nan=False)+'\n'); f.flush(); os.fsync(f.fileno())
    def totals(self,stage=None):
        rows=[r for r in self.rows() if stage is None or r['stage']==stage]
        return {k:sum(r.get('amount',0) for r in rows if r['event']==event) for k,event in [('states','STATE'),('rhs','RHS')]}
    def check(self):
        if (OUT/'RUN_CLOSED.json').exists(): raise Limit('pilot closed')
        if self.stage!='E0' and not (OUT/'E0_COMPLETE.json').exists(): raise Limit('E0 must complete first')
        if self.stage=='E2':
            decision=json.loads((OUT/'E1_DECISION.json').read_text())
            if not decision.get('E2_authorized_by_evidence'): raise Limit('two qualified independent pairs required')
    def reserve(self,event,amount,**metadata):
        self.check(); col='states' if event=='STATE' else 'rhs'; cap=self.CAPS[self.stage][col=='rhs']
        with (OUT/'budget.lock').open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            totals=self.totals(self.stage); total=self.totals()
            if totals[col]+amount>cap or total[col]+amount>{'states':256,'rhs':10240}[col]: raise Limit(f'{self.stage} {col} hard cap')
            self.append(dict(event=event,amount=amount,cumulative_stage=totals[col]+amount,context=self.context,**metadata))
    def reserve_state(self,material,*,resolution,purpose): self.reserve('STATE',1,resolution=resolution,purpose=purpose)
    def reserve_rhs(self,count,role): self.reserve('RHS',int(count),role=role)
    @contextmanager
    def span(self,label,**counts):
        self.counts.update(counts); t=time.perf_counter(); cpu=time.process_time(); status='PASS'; err=None
        try: yield
        except BaseException as e: status='FAIL'; err=repr(e); raise
        finally: self.append(dict(event='ACTION',label=label,counts=counts,wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-cpu,status=status,error=err,context=self.context))

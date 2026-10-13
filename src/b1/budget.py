"""Append-only reservations: failures and interrupted calls are charged."""
from __future__ import annotations
import json, time, fcntl, os
from pathlib import Path
from contextlib import contextmanager

class PhysicsLimit(RuntimeError): pass

class Ledger:
    def __init__(self,root,cap=384):
        self.root=Path(root); self.root.mkdir(parents=True,exist_ok=True)
        self.path=self.root/'PHYSICS_LEDGER.jsonl'; self.cap=cap
        self.lock=self.root/'physics.lock'
        self.contract=self.root/'BUDGET_CONTRACT.json'
        with self.locked():
            if self.contract.exists():
                d=json.loads(self.contract.read_text())
                if d['cap']!=cap:raise PhysicsLimit('cumulative cap is immutable without user authorization')
                if self.count()<d['high_water_reservations']:raise PhysicsLimit('ledger missing/truncated; cannot reset charged physics')
            else:
                # Adoption of already append-only reservations never erases history.
                rows=self.rows()
                if any(r.get('cap',cap)!=cap for r in rows if r['event']=='BEGIN'):raise PhysicsLimit('recorded budget mismatch')
                self.contract.write_text(json.dumps(dict(schema='B1-budget-v1',cap=cap,high_water_reservations=self.count(),historical_2d_calls=74,permission='human instruction required for cap extension'),indent=2)+'\n')
    @contextmanager
    def locked(self):
        with self.lock.open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            try: yield
            finally: fcntl.flock(f,fcntl.LOCK_UN)
    def rows(self):
        return [json.loads(l) for l in self.path.read_text().splitlines()] if self.path.exists() else []
    def count(self):
        return sum(r['event']=='BEGIN' for r in self.rows())
    def append(self,row):
        with self.path.open('a') as f:
            f.write(json.dumps(row,ensure_ascii=False,allow_nan=False)+'\n'); f.flush(); os.fsync(f.fileno())
    def evaluate(self,kind,parameters,fn):
        with self.locked():
            n=self.count()
            if n>=self.cap:
                self.append(dict(event='REJECT',kind=kind,reason='PHYSICS_CAP',cap=self.cap))
                raise PhysicsLimit(f'B1 cumulative budget {n}/{self.cap}')
            identity=n+1
            self.append(dict(event='BEGIN',id=identity,kind=kind,parameters=parameters,cap=self.cap))
            contract=json.loads(self.contract.read_text());contract['high_water_reservations']=identity
            self.contract.write_text(json.dumps(contract,indent=2)+'\n')
        t=time.perf_counter()
        try:
            out=fn()
        except BaseException as e:
            with self.locked(): self.append(dict(event='END',id=identity,status='FAIL',error=type(e).__name__+': '+str(e),wall_seconds=time.perf_counter()-t))
            raise
        with self.locked(): self.append(dict(event='END',id=identity,status='PASS',wall_seconds=time.perf_counter()-t))
        return out
    def cache_hit(self,key):
        with self.locked(): self.append(dict(event='CACHE_HIT',key=key,new_physics=0))
    def summary(self):
        rows=self.rows(); begins=[r for r in rows if r['event']=='BEGIN']; ends=[r for r in rows if r['event']=='END']
        return dict(new_physics_evaluations=len(begins),cap=self.cap,remaining=self.cap-len(begins),
                    by_kind={k:sum(r['kind']==k for r in begins) for k in sorted({r['kind'] for r in begins})},
                    failures=sum(r['status']=='FAIL' for r in ends),unfinished=len(begins)-len(ends),
                    full_forward_wall_seconds=sum(r['wall_seconds'] for r in ends),
                    cache_hits=sum(r['event']=='CACHE_HIT' for r in rows), historical_package_2d_calls=74,
                    combined_package_and_b1_analytic_evaluations=74+len(begins),gpu_calls=0,training_updates=0,locked_arrays_read=0)

"""Unified E0/E1 state and RHS reservation with durable all-resolution counting."""
from __future__ import annotations
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
import json, math, os, time


class BudgetStop(RuntimeError):
    pass


class Book:
    def __init__(self, root, *, stage, state_limit=192, rhs_limit=16384):
        self.root=Path(root).resolve(); self.stage=str(stage)
        if stage not in {'e0','tests','e1','report'}:raise ValueError('REGISTERED_STAGE_REQUIRED')
        if not 0<=state_limit<=192 or not 0<=rhs_limit<=16384:raise ValueError('CAP_INCREASE_FORBIDDEN')
        self.state_limit=state_limit;self.rhs_limit=rhs_limit;self.counts=Counter()
        self.ledger=self.root/'results/a30/PHYSICAL_STATE_LEDGER.jsonl'
        self.ledger.parent.mkdir(parents=True,exist_ok=True)
        self.lock=self.ledger.with_suffix('.lock')

    def _reserve(self, event, amount, **metadata):
        if isinstance(amount,bool) or not isinstance(amount,int) or amount<1:raise ValueError('POSITIVE_INTEGER_RESERVATION')
        if self.stage!='e1':raise BudgetStop('REAL_PHYSICS_FORBIDDEN_IN_'+self.stage.upper())
        gate=self.root/'results/a30/E0_GATE.json'
        if not gate.exists() or json.loads(gate.read_text()).get('status')!='PASS':raise BudgetStop('E0_PASS_REQUIRED')
        try: fd=os.open(self.lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY,0o600)
        except FileExistsError:raise BudgetStop('ACCOUNTING_LOCK_BUSY') from None
        try:
            os.write(fd,str(os.getpid()).encode());os.close(fd)
            rows=[json.loads(x) for x in self.ledger.read_text().splitlines()] if self.ledger.exists() else []
            used=sum(r['amount'] for r in rows if r.get('event')==event)
            limit=self.state_limit if event=='STATE_RESERVED' else self.rhs_limit
            if used+amount>limit:raise BudgetStop('UNIFIED_'+event+'_LIMIT')
            record=dict(event=event,amount=amount,cumulative=used+amount,limit=limit,stage=self.stage,
                        timestamp_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**metadata)
            with self.ledger.open('a') as f:
                f.write(json.dumps(record,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
            return record
        finally:self.lock.unlink(missing_ok=True)

    def reserve_state(self, material=None, *, resolution=24, purpose='unspecified', count=1):
        # Counting is per evaluation, not per unique material nor resolution.
        shape=list(getattr(material,'shape',()))
        return self._reserve('STATE_RESERVED',count,resolution=int(resolution),purpose=purpose,material_shape=shape)

    def reserve_rhs(self,count,role='unspecified'):
        return self._reserve('RHS_RESERVED',int(count),role=role)

    @contextmanager
    def span(self,label,**counts):
        start=time.perf_counter();cpu=time.process_time();status='FAILED'
        self.counts.update(counts)
        try:yield;status='COMPLETE'
        finally:
            row=dict(event='ACTION',label=label,status=status,counts=counts,stage=self.stage,
                     wall_seconds=time.perf_counter()-start,process_cpu_seconds=time.process_time()-cpu)
            with (self.root/'results/a30/ACTIONS.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')


def physical_totals(root):
    p=Path(root)/'results/a30/PHYSICAL_STATE_LEDGER.jsonl'
    rows=[json.loads(x) for x in p.read_text().splitlines()] if p.exists() else []
    return dict(material_state_evaluations=sum(r['amount'] for r in rows if r['event']=='STATE_RESERVED'),
                complex_RHS=sum(r['amount'] for r in rows if r['event']=='RHS_RESERVED'),
                reservations=len(rows),counting='Every real invocation at all resolutions; attempts/repeats stay paid')

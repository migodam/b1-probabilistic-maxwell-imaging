"""Atomic project-local accounting. Reservations are charged even on failure."""
import contextlib
import fcntl
import json
import math
import os
from pathlib import Path
import tempfile
import time
import numpy as np


def json_default(value):
    """Scientific diagnostic metadata only; accounting counts stay native."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, complex):
        return [value.real, value.imag]
    raise TypeError(f"Unsupported JSON metadata type: {type(value).__name__}")

CAPS = {
    'E0E1': dict(rhs=4000, factorizations=1200, process_cpu_seconds=7200),
    'A0': dict(rhs=12000, factorizations=3000, process_cpu_seconds=14400),
    'E2E3': dict(rhs=12000, factorizations=3000, process_cpu_seconds=14400),
    'E4': dict(rhs=12000, factorizations=3000, process_cpu_seconds=14400),
    'E5': dict(rhs=8000, factorizations=2000, process_cpu_seconds=10800),
}


def stage_key(phase):
    key = str(phase).upper().replace('-', '').replace('_', '').replace('+', '').replace('/', '')
    key = {'E0':'E0E1','E1':'E0E1','E01':'E0E1','E2':'E2E3','E3':'E2E3','E23':'E2E3'}.get(key,key)
    if key not in CAPS:
        raise ValueError('unknown accounting stage: '+str(phase))
    return key


class BudgetExceeded(RuntimeError):
    pass


class BudgetBook:
    def __init__(self, path=None):
        self.path = Path(path or os.environ.get('WSLAM_LEDGER', Path(__file__).resolve().parents[2]/'results'/'budget_ledger.json')).resolve()
        self.path.parent.mkdir(parents=True,exist_ok=True)
        self.lock_path = self.path.with_suffix(self.path.suffix+'.lock')

    @contextlib.contextmanager
    def _locked(self):
        with self.lock_path.open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            try:
                if self.path.exists():
                    with self.path.open() as f:
                        data = json.load(f)
                    if data.get('schema') != 'wslam-budget-v1' or data.get('caps') != CAPS:
                        raise ValueError('incompatible ledger; old ledgers cannot be imported')
                else:
                    data = dict(schema='wslam-budget-v1',caps=CAPS,events=[])
                yield data
                fd, temporary = tempfile.mkstemp(prefix=self.path.name+'.',dir=self.path.parent)
                try:
                    with os.fdopen(fd,'w') as f:
                        # Lossless compact encoding avoids repeatedly writing
                        # large indentation overhead in a long physical bank.
                        # Schema, atomic replacement and charges are unchanged.
                        json.dump(data,f,separators=(',',':'),allow_nan=False,default=json_default)
                        f.flush()
                        os.fsync(f.fileno())
                    os.replace(temporary,self.path)
                finally:
                    if os.path.exists(temporary):
                        os.unlink(temporary)
            finally:
                fcntl.flock(lock,fcntl.LOCK_UN)

    @staticmethod
    def _usage(data,stage):
        used = dict(rhs=0,factorizations=0,process_cpu_seconds=0.0)
        for e in data['events']:
            if e.get('stage') == stage:
                for key in used:
                    used[key] += e.get(key,0)
        return used

    def snapshot(self,phase=None):
        with self._locked() as data:
            stages = [stage_key(phase)] if phase is not None else list(CAPS)
            return {s:dict(cap=dict(CAPS[s]),used=self._usage(data,s),remaining={k:max(0,CAPS[s][k]-self._usage(data,s)[k]) for k in CAPS[s]}) for s in stages}

    def totals(self,stage=None):
        """Charged counts plus measured process CPU; no refunds."""
        snapshots = self.snapshot(stage)
        total = dict(rhs=0, factorizations=0, cpu_seconds=0.0)
        for entry in snapshots.values():
            used = entry['used']
            total['rhs'] += used['rhs']
            total['factorizations'] += used['factorizations']
            total['cpu_seconds'] += used['process_cpu_seconds']
        return total

    def reserve(self,phase,label,rhs=0,factorizations=0,**metadata):
        stage = stage_key(phase)
        if any(not isinstance(n,int) or isinstance(n,bool) or n<0 for n in (rhs,factorizations)):
            raise ValueError('RHS and factorization charges must be nonnegative integers')
        if {'stage','kind','time','rhs','factorizations','process_cpu_seconds','id','label'} & metadata.keys():
            raise ValueError('metadata overrides accounting fields')
        with self._locked() as data:
            used = self._usage(data,stage)
            if used['process_cpu_seconds'] >= CAPS[stage]['process_cpu_seconds'] or used['rhs']+rhs > CAPS[stage]['rhs'] or used['factorizations']+factorizations > CAPS[stage]['factorizations']:
                raise BudgetExceeded(f'{stage}: pre-action reservation denied; used={used}, request rhs={rhs}, factorizations={factorizations}')
            record = dict(id=len(data['events'])+1,stage=stage,kind='reservation',time=time.time(),label=str(label),rhs=rhs,factorizations=factorizations,**metadata)
            data['events'].append(record)
            return dict(record)

    def event(self,phase=None,label='',**metadata):
        phase = phase if phase is not None else metadata.pop('stage',None)
        if phase is None:
            phase = os.environ.get('WSLAM_STAGE')
        stage = stage_key(phase)
        cpu = metadata.pop('process_cpu_seconds',0.0)
        if not isinstance(cpu,(int,float)) or not math.isfinite(cpu) or cpu<0:
            raise ValueError('process CPU must be finite and nonnegative')
        for key in ('rhs', 'factorizations'):
            if key in metadata:
                observed = 'observed_' + key
                if observed in metadata:
                    raise ValueError('conflicting informational accounting metadata')
                metadata[observed] = metadata.pop(key)
        if {'stage','id','time'} & metadata.keys():
            raise ValueError('event cannot change reserved accounting fields')
        with self._locked() as data:
            record = dict(id=len(data['events'])+1,stage=stage,time=time.time(),label=str(label),process_cpu_seconds=cpu,**metadata)
            record.setdefault('kind','event')
            data['events'].append(record)
            return dict(record)


Book = BudgetBook

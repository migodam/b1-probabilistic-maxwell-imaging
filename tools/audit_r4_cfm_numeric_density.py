#!/usr/bin/env python3
"""Authorized bounded offline numerical CFM-density audit; no physical calls."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import json
from pathlib import Path
import resource
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import torch
from scipy.special import ndtri
from p3r4.core import Book,OUT
from p3r4.density import load_density
from p3r4.evaluation import checkpoint_identity,write_new_json,read_online_measurement

TOLERANCE=.05
NFES=(16,32,64)
PARENTS=tuple(f'R4DEV_{i:05d}' for i in range(4))
MEMORY_CAP=2*1024**3


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)


def main():
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    checkpoint=OUT/'models/density/density_cfm_32768_s2026102501/BEST.pt'
    identity=checkpoint_identity(checkpoint)
    job='cfm_density_trace_audit_'+str(time.time_ns())
    output=OUT/'AUDIT'/job
    output.mkdir(parents=True,exist_ok=False)
    # This file is persisted before GT or density computation and never tuned.
    frozen=dict(status='FROZEN_BEFORE_NUMERICAL_DENSITY_COMPUTE',checkpoint=identity,
        parents=list(PARENTS),NFEs=list(NFES),absolute_agreement_tolerance_nats=TOLERANCE,
        primary_numeric_agreement='absolute logq32-logq64 <= .05 nat on all4 fixed parents',
        coordinate='raw original Gaussianized dx66',device='cpu',threads=1,
        trace='complete deterministic66-coordinate divergence, no Hutchinson estimator',
        RSS_cap_bytes=MEMORY_CAP,controlled_wall_cap_seconds=600,physics_calls_authorized=0,optimizer_updates_authorized=0,
        evidence_role='isolated offline audit of first4 already exposed DEV; not fresh holdout or model selection',
        density_scope='numerical continuous-ODE density, not exact finite-RK4-map density or exact-MH proposal')
    write_new_json(output/'FROZEN_AUDIT.json',frozen)
    book=Book('I',job,analysis_only=True,max_wall=600)
    wall,cpu=time.monotonic(),time.process_time()
    records=[];status='COMPLETE';error=None
    try:
        if rss()>MEMORY_CAP:raise MemoryError('RSS cap before model load')
        predictor=load_density(checkpoint,'cpu')
        if predictor.model.kind!='cfm' or predictor.model.dimension!=66 or predictor.model.context_dim!=512:
            raise ValueError('fixed raw512/full66 CFM checkpoint required')
        if predictor.payload.get('contract',{}).get('cfm_base')!='standard':
            raise ValueError('audit contract requires existing standard-Gaussian base')
        if checkpoint_identity(checkpoint)!=identity:raise ValueError('checkpoint changed after freeze')
        for parent in PARENTS:
            # No training, normalizer fitting, proposal selection, physical
            # forward, or live-thread GT writeback occurs in this process.
            y=read_online_measurement(OUT/'inputs/online'/f'{parent}.npz',parent)
            with np.load(OUT/'inputs/offline_only'/f'{parent}.npz',allow_pickle=False) as z:
                q=z['q']
            if q.shape!=(66,) or not np.isfinite(q).all() or np.any(q[64:]<=0) or np.any(q[64:]>=1):
                raise ValueError('finite original native GT and interior uniform properties required; no clipping')
            x=np.r_[q[:64],ndtri(q[64:])][None]
            for nfe in NFES:
                if rss()>MEMORY_CAP:raise MemoryError('RSS cap exceeded before density call')
                cw,cc=time.monotonic(),time.process_time()
                logq=float(predictor.log_prob(x,y,batch_size=1,allow_approximate=True,nfe=nfe)[0])
                finite=bool(np.isfinite(logq))
                row=dict(parent=parent,NFE=nfe,raw_logq_x66=logq if finite else None,
                         raw_NLL_x66=-logq if finite else None,status='FINITE_NUMERIC_DENSITY' if finite else 'NONFINITE_NUMERIC_DENSITY',
                         complete_divergence=True,velocity_evaluations=nfe,trace_component_VJPs=66*nfe,
                         wall_seconds=time.monotonic()-cw,cpu_seconds=time.process_time()-cc,peak_RSS_bytes=rss())
                records.append(row)
                with (output/'TRACE_RESULTS.jsonl').open('a') as file:file.write(json.dumps(row,allow_nan=False)+'\n')
                print(json.dumps(row,allow_nan=False),flush=True)
                if rss()>MEMORY_CAP:raise MemoryError('RSS cap exceeded after density call')
        if checkpoint_identity(checkpoint)!=identity:raise ValueError('checkpoint changed during audit')
    except BaseException as exc:
        status='FAILED';error=repr(exc)
        raise
    finally:
        convergence=[]
        for parent in PARENTS:
            entries={r['NFE']:r['raw_logq_x66'] for r in records if r['parent']==parent}
            fine=entries.get(32) is not None and entries.get(64) is not None
            all_three=all(entries.get(n) is not None for n in NFES)
            convergence.append(dict(parent=parent,
                absolute_delta_32_64=abs(entries[32]-entries[64]) if fine else None,
                absolute_delta_16_64=abs(entries[16]-entries[64]) if all_three else None,
                within_prefrozen_tolerance=bool(abs(entries[32]-entries[64])<=TOLERANCE) if fine else None))
        passed=(len(records)==12 and all(c['within_prefrozen_tolerance'] is True for c in convergence))
        result=dict(status=status,error=error,job=job,job_receipt=str(book.jobdir/'CLOSED.json'),
            frozen_contract=frozen,computations_completed=len(records),records=records,convergence=convergence,
            numeric_agreement_status='WITHIN_PREFROZEN_TOLERANCE_ON_FIXED4_EXPOSED_DEV' if passed else 'NOT_ESTABLISHED_AT_FROZEN_TOLERANCE',
            exact_discrete_sampling_map_density=False,exact_MH_proposal_density=False,posterior_qualification='NOT_ESTABLISHED',
            checkpoint_unchanged=checkpoint_identity(checkpoint)==identity,physics_calls=0,optimizer_updates=0,
            wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,peak_RSS_bytes=rss())
        write_new_json(output/'RESULT.json',result)
        book.close(status,error)
        print(json.dumps(dict(status=status,output=str(output),numeric_agreement=result['numeric_agreement_status'],physics_calls=0,optimizer_updates=0)),flush=True)


if __name__=='__main__':main()

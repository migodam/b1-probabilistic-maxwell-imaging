"""Authorized2-object exposed-DEV adaptive CFM density diagnostic;0physics."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,resource,sys,time,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import torch
from scipy.special import ndtri
from p3r4.core import Book,OUT,BudgetStop,write,save
from p3r4.density import load_density
from p3r4.anchor_density import load_anchor_density
from p3r4.probability_eval import read_online_observation,read_online_anchor_context
from p3r4.adaptive_cfm_audit import FullTraceODE,fixed32_double
from p3r4.train_physics_transport import identity,verify_identity

PARENTS=['R4DEV_00000','R4DEV_00001']
TOLS=[dict(atol=1e-5,rtol=1e-4),dict(atol=1e-6,rtol=1e-5)]
MODELS=[('raw_matched','raw',OUT/'models/density/extended_cfm_matched_s2026102501/BEST.pt'),
 ('anchor_none','anchor',OUT/'models/anchor_density/anchor_cfm_none_s2026102501/BEST.pt'),
 ('anchor_geometry','anchor',OUT/'models/geometry_weighted_cfm_v1/geometry_cfm_geometry_s2026102501/BEST.pt')]


def rss():return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)


def prepare_point(predictor,loader,parent):
    online=OUT/'inputs/online'/(parent+'.npz')
    y,_=read_online_observation(online,parent)
    condition={}
    context=OUT/'inputs/anchor_conditions/dev'/(parent+'.npz')
    if loader=='anchor':condition=read_online_anchor_context(context,parent)
    # Only these two explicitly authorized exposed DEV GT arrays are read.
    truth=OUT/'inputs/offline_only'/(parent+'.npz')
    with np.load(truth,allow_pickle=False) as z:q=z['q']
    if q.shape!=(66,) or not np.isfinite(q).all() or np.any(q[64:]<=0) or np.any(q[64:]>=1):raise ValueError('finite interior exposed GT; no clipping')
    raw=np.r_[q[:64],ndtri(q[64:])][None]
    public=float(predictor.log_prob(raw,y,batch_size=1,allow_approximate=True,nfe=32,**condition)[0])
    if loader=='anchor':
        features,center,base,coordinate=predictor._prepare(y,condition['context_x'],condition.get('base_mean'),condition.get('base_chol'),condition.get('context_J'))
        if coordinate.mode!='none':raise ValueError('fixed models require original none coordinate map')
        residual=coordinate.forward(raw-center);extra=coordinate.logabsdet[0];inner=predictor.inner
    else:features=y;center=np.zeros((1,66));base={};residual=raw;extra=0.;inner=predictor
    stats=inner.stats
    standardized=((residual@stats['frame']-stats['x_mean'])/stats['x_scale'])[0]
    cy=inner._y(features) # Original API-normalized float32 condition, then promoted.
    mean,chol=inner._base(base.get('base_mean'),base.get('base_chol'),1)
    outer=float(extra-np.log(stats['x_scale']).sum()+stats.get('frame_logabsdet',0.))
    return dict(raw=raw[0],standardized=standardized,public_standardized=standardized.astype(np.float32).astype(np.float64),
        condition=cy.astype(np.float64),base_mean=mean.astype(np.float64) if mean is not None else None,
        base_chol=chol.astype(np.float64) if chol is not None else None,outer_logdet=outer,
        center=center[0],x_scale=stats['x_scale'],x_mean=stats['x_mean'],frame_inverse=stats['frame_inverse'],
        public_fixed32_logq=public,provenance=dict(observation=identity(online),exposed_truth=identity(truth),
            context=identity(context) if loader=='anchor' else None))


def raw_inverse(s,p):return (np.asarray(s)*p['x_scale']+p['x_mean'])@p['frame_inverse']+p['center']


def main():
    torch.set_num_threads(1);torch.set_num_interop_threads(1)
    job='adaptive_cfm_density_'+str(time.time_ns());out=OUT/'AUDIT'/job;out.mkdir(parents=True)
    sources=[identity(ROOT/p) for p in ['src/p3r4/adaptive_cfm_audit.py','tools/audit_r4_adaptive_cfm.py',
        'src/p3r4/density.py','src/p3r4/anchor_density.py','src/p3r4/conditional_coordinates.py','src/p3r4/probability_eval.py']]
    frozen=dict(status='FROZEN_BEFORE_GT_OR_NUMERIC_DENSITY',parents=PARENTS,
        models=[dict(id=i,loader=l,checkpoint=identity(p)) for i,l,p in MODELS],source_identities=sources,
        solver='DOP853',tolerances=TOLS,agreement_nats=.05,cycle_error='raw maxabs/RMS/relative report only; no new acceptance threshold',
        backward_ODE='ds/dt=v(s,t,c), dl/dt=-trace(dv/ds), integrate1->0 from(s1,0)',
        density='logq_raw=logp_base(s0)-l_back+logabsdet(raw_to_network_map)',
        sign='l_back=+integral_0^1 div(v)dt; subtract it from base log density',
        precision='stored Float32 parameter/buffer values promoted to CPU Float64; same original affine maps and API Float32 normalized condition/base values promoted',
        target_point='true raw Gaussianized GT mapped by frozen affine in Float64; fixed32 also records source Float32-rounded input separately',
        divergence='all66 exact autograd diagonal VJPs; no stochastic trace',
        forward_cycle='velocity-only DOP8530->1 from recovered s0, invert frozen affine to rawdx66; no resampling',
        fixed_comparisons='source publicFloat32 RK4 NFE32; Float64 RK4 NFE32 at rounded s1; Float64 RK4 NFE32 at true s1',
        old_fixed_audit=str(OUT/'AUDIT/cfm_density_trace_audit_1792871479806924000/RESULT.json'),
        old_fixed_scope='older unmatched small rawCFM, different checkpoint; failed3/4 retained, not numerically interchanged with matched model',
        maximum_wall=900,CPU_threads=1,RSS_cap_bytes=2*1024**3,physics_RHS=0,optimizer_updates=0,
        heldout_GT='FORBIDDEN; exposed DEV0/1 only',scope='integration diagnostic only; no frozenNFE32/samplelaw/lockedNLL/model selection change',
        posterior_qualification='NOT_ESTABLISHED',no_fallback_or_expansion=True)
    write(out/'FROZEN_AUDIT.json',frozen)
    book=Book('I',job,analysis_only=True,max_wall=900);start,cpu=time.monotonic(),time.process_time()
    records=[];comparisons=[];status,error='COMPLETE',None
    def guard():
        if rss()>2*1024**3:raise BudgetStop('adaptive CFM2GiB cap')
        if time.monotonic()-start>895:raise BudgetStop('adaptive CFM preserve-close wall margin')
    try:
        for m in frozen['models']:
            guard();predictor=load_density(m['checkpoint']['path'],'cpu') if m['loader']=='raw' else load_anchor_density(m['checkpoint']['path'],'cpu')
            if predictor.model.kind!='cfm' or predictor.model.dimension!=66:raise ValueError('unchanged full66 CFM required')
            points={}
            # Public Float32 evaluations precede promotion of the same in-memory
            # network. No checkpoint/source parameter array is rewritten.
            for parent in PARENTS:
                guard();p=prepare_point(predictor,m['loader'],parent);points[parent]=p
                save(out/'inputs'/f"{m['id']}_{parent}.npz",**{k:v for k,v in p.items() if isinstance(v,np.ndarray)},
                    role=np.array('EXPOSED_DEV_NUMERICAL_DIAGNOSTIC_ONLY'),outer_logdet=np.array(p['outer_logdet']))
                write(out/'inputs'/f"{m['id']}_{parent}.json",dict(provenance=p['provenance'],public_fixed32_logq=p['public_fixed32_logq'],
                    public_input_rounding_maxabs=float(np.max(abs(p['standardized']-p['public_standardized'])))))
            predictor.model.double().eval().requires_grad_(False)
            for parent,p in points.items():
                guard();ode=FullTraceODE(predictor.model,p['condition'],guard)
                fixed_true=fixed32_double(ode.model,p['condition'],p['standardized'],p['base_mean'],p['base_chol'],p['outer_logdet'])
                fixed_rounded=fixed32_double(ode.model,p['condition'],p['public_standardized'],p['base_mean'],p['base_chol'],p['outer_logdet'])
                fixed=dict(model=m['id'],parent=parent,publicFloat32_NFE32_logq=p['public_fixed32_logq'],
                    Float64_NFE32_true_input_logq=fixed_true,Float64_NFE32_rounded_input_logq=fixed_rounded,
                    arithmetic_precision_delta=p['public_fixed32_logq']-fixed_rounded,input_rounding_delta=fixed_rounded-fixed_true)
                comparisons.append(fixed);write(out/'fixed32'/f"{m['id']}_{parent}.json",fixed)
                for index,tol in enumerate(TOLS):
                    row=dict(model=m['id'],parent=parent,tolerance_index=index,**tol,status='FAILED')
                    calls_before=ode.velocity_calls;traces_before=ode.trace_component_VJPs;tick=time.monotonic()
                    try:
                        backward=ode.backward_density(p['standardized'],p['base_mean'],p['base_chol'],p['outer_logdet'],**tol)
                        x0=backward.pop('x0_standardized')
                        save(out/'solutions'/f"{m['id']}_{parent}_tol{index}_backward.npz",x0_standardized=x0,
                            backward_logJac=np.array(backward['backward_accumulated_negative_divergence']))
                        row['backward']=backward
                        if not backward['success']:raise RuntimeError('DOP853 backward failed: '+backward['message'])
                        cycle=ode.forward_cycle(x0,**tol);s1=cycle.pop('x1_standardized')
                        reconstructed=raw_inverse(s1,p);difference=reconstructed-p['raw']
                        save(out/'solutions'/f"{m['id']}_{parent}_tol{index}_cycle.npz",x1_standardized=s1,raw_cycle=reconstructed,raw_target=p['raw'])
                        row.update(status='COMPLETE' if cycle['success'] else 'FORWARD_CYCLE_FAILED',forward_cycle=cycle,
                            raw_cycle_maxabs=float(abs(difference).max()),raw_cycle_RMS=float(np.sqrt(np.mean(difference**2))),
                            raw_cycle_relative_L2=float(np.linalg.norm(difference)/max(np.linalg.norm(p['raw']),1e-30)))
                    except BudgetStop:raise
                    except Exception as exc:row.update(status='FAILED',error=repr(exc))
                    finally:
                        row.update(total_velocity_calls=ode.velocity_calls-calls_before,
                            total_trace_component_VJPs=ode.trace_component_VJPs-traces_before,
                            wall_seconds=time.monotonic()-tick,RSS_peak_bytes=rss())
                        records.append(row)
                        with (out/'RECORDS.jsonl').open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
                        print(json.dumps(row,allow_nan=False),flush=True)
            del predictor
        for m in frozen['models']:verify_identity(m['checkpoint'])
        for s in sources:verify_identity(s)
    except BaseException as exc:status,error='PARTIAL_BUDGET_STOP' if isinstance(exc,BudgetStop) else 'FAILED',repr(exc)
    finally:
        agreement=[]
        for m in frozen['models']:
            for parent in PARENTS:
                pair=[next((r for r in records if r['model']==m['id'] and r['parent']==parent and r['tolerance_index']==i),None) for i in (0,1)]
                complete=all(r is not None and r['status']=='COMPLETE' for r in pair)
                delta=abs(pair[0]['backward']['raw_NLL_x66']-pair[1]['backward']['raw_NLL_x66']) if complete else None
                fixed=next((r for r in comparisons if r['model']==m['id'] and r['parent']==parent),None)
                agreement.append(dict(model=m['id'],parent=parent,all_tolerances_complete=complete,
                    absolute_tolerance_NLL_delta=delta,within_prefrozen_point05_nat=delta<=.05 if delta is not None else None,
                    tight_NLL=pair[1]['backward']['raw_NLL_x66'] if complete else None,
                    tight_raw_cycle_maxabs=pair[1]['raw_cycle_maxabs'] if complete else None,
                    fixed32_double_minus_adaptive_tight_NLL=(-fixed['Float64_NFE32_true_input_logq']-pair[1]['backward']['raw_NLL_x66']) if complete else None))
        write(out/'RESULT.json',dict(status=status,error=error,frozen=frozen,records=records,agreement=agreement,fixed32=comparisons,
            all6_within_prefrozen_tolerance=len(agreement)==6 and all(a['within_prefrozen_point05_nat'] is True for a in agreement),
            planned_solutions=12,recorded_solutions=len(records),physics_RHS=0,optimizer_updates=0,heldout_GT_reads=0,
            exact_discrete_map_density=False,exact_MH_density=False,frozen_algorithm_changed=False,posterior_qualification='NOT_ESTABLISHED',
            controlled_work_wall_seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,RSS_peak_bytes=rss(),receipt=str(book.jobdir/'CLOSED.json')))
        book.close(status,error);print(out,flush=True)
    if status!='COMPLETE':raise SystemExit(1)


if __name__=='__main__':main()

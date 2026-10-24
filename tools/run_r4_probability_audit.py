#!/usr/bin/env python3
"""Freeze generated samples -> bounded paid likelihoods -> isolated offline SBC.

No subcommand generates samples or contexts, trains, selects small residuals,
changes a sampler, drops duplicates, or silently replaces failed sample0.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import argparse,json,sys,time,resource,traceback
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
import numpy as np
from scipy.special import ndtr
from p3r4.evaluation import checkpoint_identity,write_new_json
from p3r4.probability_eval import safe_id
from p3r4.probability_physics import stable_importance,likelihood_rank,gamma256_summary
STAGE=ROOT/'research/p3l_r4_conditional_transport'


def load(path):return json.loads(Path(path).read_text())


def plain(value):
    if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [plain(v) for v in value]
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    return value


def freeze(args):
    spec=load(args.spec);cfg=load(STAGE/'configs/TARGET.json');records=[]
    if spec.get('likelihood_samples',0) not in (0,64):raise ValueError('fixed first64 or sample0 only')
    if not spec.get('generated'):raise ValueError('frozen generated directories required')
    group_ids=[g['id'] for g in spec['generated']]
    if len(group_ids)!=len(set(group_ids)) or not all(safe_id(i) for i in group_ids):raise ValueError('unique safe generated IDs')
    for group in spec['generated']:
        if group.get('density') not in ('explicit','none'):raise ValueError('density explicit/none must be declared')
        parents=group['parents']
        if not parents or len(parents)!=len(set(parents)) or not all(safe_id(p) for p in parents):raise ValueError('unique safe frozen parents')
        directory=Path(group['directory']).resolve();gatepath=Path(group['online_gate']).resolve();gate=load(gatepath)
        if gate.get('status')!='ONLINE_COMPLETE' or gate.get('truth_read') is not False:raise ValueError('complete observation-only generator gate required')
        if gate.get('freeze',{}).get('decoder',{}).get('config')!=cfg:raise ValueError('generated samples use a different declared statistical target')
        sealed={str(Path(v['identity']['path']).resolve()):v['identity'] for v in gate['artifacts']}
        extended=group.get('extended_parents',parents if spec.get('likelihood_samples')==64 else [])
        if not set(extended).issubset(set(parents)) or len(extended)!=len(set(extended)):raise ValueError('extended parent subset must be predeclared')
        group.update(directory=str(directory),online_gate=str(gatepath),gate_identity=checkpoint_identity(gatepath),extended_parents=extended,
                     observations=str(Path(group['observations']).resolve()))
        for parent in parents:
            source=directory/f'{parent}.npz';obs=Path(group['observations'])/f'{parent}.npz'
            if 'offline_only' in source.parts or 'offline_only' in obs.parts:raise ValueError('offline data forbidden in paid driver')
            index=list(range(64)) if parent in extended and spec.get('likelihood_samples')==64 else [0]
            record=dict(group=group['id'],parent=parent,source=str(source),observation=str(obs),sample_indices=index,density=group['density'],status='FROZEN')
            if source.exists():
                identity=checkpoint_identity(source)
                if sealed.get(str(source))!=identity:raise ValueError('sample is not sealed by declared online gate')
                record['source_identity']=identity
            else:record.update(status='GENERATOR_FAILED_MISSING_ARTIFACT',generator_failures=[v for v in gate.get('failures',[]) if v.get('parent')==parent])
            if obs.exists():record['observation_identity']=checkpoint_identity(obs)
            else:record['observation_failure']='MISSING_ORIGINAL_ONLINE_OBSERVATION'
            records.append(record)
    spec.update(status='FROZEN_PAID_PROBABILITY_AUDIT',freeze_pid=os.getpid(),freeze_time_ns=time.time_ns(),target_id=cfg['target_id'],sigma=cfg['noise']['sigma'],target_config=cfg,records=records,
        fixed_sample_index=0,planned_sample_evaluations=sum(len(r['sample_indices']) for r in records),
        planned_max_RHS=8*sum(len(r['sample_indices']) for r in records),truth_read=False,no_residual_selection=True,
        duplicate_policy='retain every fixed original index, including duplicate/rejected states; no dedup/resample',posterior_qualification='NOT_ESTABLISHED')
    write_new_json(args.out,spec);print(json.dumps(dict(status='FROZEN',records=len(records),planned_max_RHS=spec['planned_max_RHS'])) )


def online_inputs(record,frozen):
    if record['status']!='FROZEN':raise ValueError(record['status'])
    if checkpoint_identity(record['source'])!=record['source_identity'] or checkpoint_identity(record['observation'])!=record.get('observation_identity'):raise ValueError('sealed online inputs changed/missing')
    with np.load(record['source'],allow_pickle=False) as z:
        if str(z['parent'].item())!=record['parent'] or int(z['fixed_sample_index'])!=0 or not np.array_equal(z['sample_index'],np.arange(128)):raise ValueError('frozen sample identity/order')
        q=z['q_samples'];x=z['x_samples'];source_y=z['observation_features'];logq=z['log_q_x'] if record['density']=='explicit' else None
    if q.shape!=(128,66) or x.shape!=(128,66) or not np.isfinite(q).all() or not np.isfinite(x).all():raise ValueError('complete finite original128 samples required')
    if np.max(np.abs(q[:,:64]-x[:,:64]))>1e-10 or np.max(np.abs(q[:,64:]-ndtr(x[:,64:])))>2e-7:raise ValueError('q samples/probit x mismatch; no clipping permitted')
    if np.any((q[:,64:]<0)|(q[:,64:]>1)):raise ValueError('illegal native attributes')
    if logq is not None and (logq.shape!=(128,) or not np.isfinite(logq).all()):raise ValueError('explicit raw dx66 density required')
    with np.load(record['observation'],allow_pickle=False) as z:
        y=z['observation']
        if 'target_id' in z and str(z['target_id'].item())!=frozen['target_id']:raise ValueError('online target mismatch')
        if 'sigma' in z and float(z['sigma'])!=frozen['sigma']:raise ValueError('online sigma mismatch')
    if y.shape!=(32,8) or not np.iscomplexobj(y) or not np.isfinite(y).all():raise ValueError('original complex online observation required')
    from p3r4.probability_eval import source_major_pack
    if not np.array_equal(source_y.reshape(1,512),source_major_pack(y)):raise ValueError('generator condition differs from original online observation')
    return q,x,logq,y


def run(args,frozen):
    if args.stage is None or not args.job or args.max_rhs is None or args.max_wall is None:raise ValueError('explicit stage/job/hard max-rhs/max-wall required')
    if args.max_rhs<1 or args.max_wall<=15 or args.max_rhs<frozen['planned_max_RHS']:raise ValueError('hard caps must cover predeclared sample roster; reduce roster before freeze rather than run selective partial set')
    if args.out.exists():raise FileExistsError('fresh paid output only; no replay')
    # No GT argument exists in the physical execution function.
    if args.truth is not None:raise ValueError('truth path prohibited during paid online audit')
    from p3r4.core import Book,Target,BudgetStop,save,write
    book=Book(args.stage,args.job,max_wall=args.max_wall);original=book._update;costs=dict(RHS=0,states=0,LU=0);status='COMPLETE';error=None;all_records=[];artifacts=[];start=time.monotonic()
    def limited(event,**kw):
        if costs['RHS']+kw.get('rhs',0)>args.max_rhs:raise BudgetStop('probability audit hardRHS cap')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)>3*1024**3:raise BudgetStop('probability audit3GiB cap')
        answer=original(event,**kw)
        for k,a in [('RHS','rhs'),('states','states'),('LU','lu')]:costs[k]+=kw.get(a,0)
        return answer
    book._update=limited
    try:
        args.out.mkdir(parents=True);target=Target(book,keep_fields=True)
        if target.cfg!=frozen['target_config']:raise ValueError('statistical target changed since audit freeze')
        write(args.out/'FROZEN_EXECUTION.json',frozen)
        for record in frozen['records']:
            group,parent=record['group'],record['parent'];folder=args.out/group/parent;folder.mkdir(parents=True);indices=record['sample_indices']
            phi=np.full(len(indices),np.nan);failure=[]
            try:q,x,logq,y=online_inputs(record,frozen);input_error=None
            except Exception as exc:q=x=logq=y=None;input_error=repr(exc)
            for row,index in enumerate(indices):
                tick=time.monotonic();before=costs['RHS'];entry=dict(group=group,parent=parent,sample_index=index,status='FAILED',no_selection=True)
                try:
                    if input_error:raise ValueError(input_error)
                    book.context=dict(parent=parent,group=group,sample_index=index,purpose='frozen_probability_likelihood_audit')
                    point=target.evaluate(q[index],purpose=f'{group}_{parent}_fixed_sample_{index:03d}')
                    if costs['RHS']-before!=8:raise AssertionError('every original8-source Full evaluation must charge8RHS')
                    phi[row]=target.phi(point,y)
                    if not np.isfinite(phi[row]):raise FloatingPointError('nonfinite physical likelihood')
                    dest=folder/f'sample_{index:03d}.npz'
                    save(dest,index=np.array(index),x=x[index],q=q[index],chi=point.state.chi,current=point.state.current,field=point.state.field,prediction=point.state.data,
                         phi=np.array(phi[row]),log_q_x=np.array(logq[index] if logq is not None else np.nan),observation=y,endpoint=np.array(point.key))
                    artifacts.append(dict(relative_path=str(dest.relative_to(args.out)),identity=checkpoint_identity(dest)))
                    entry.update(status='COMPLETE',Phi=float(phi[row]),endpoint=point.key,paid_RHS=costs['RHS']-before)
                    del point;target.backend()._transfers.clear()
                except BudgetStop:raise
                except Exception as exc:phi[row]=np.nan;entry.update(error=repr(exc),paid_RHS=costs['RHS']-before);failure.append(entry.copy())
                entry['wall_seconds']=time.monotonic()-tick
                with (folder/'SAMPLE_RECORDS.jsonl').open('a') as f:f.write(json.dumps(entry,allow_nan=False)+'\n');f.flush();os.fsync(f.fileno())
            # Missing samples remain NaN; no importance normalization over only
            # surviving samples and no replacement of failed fixed index0.
            importance=dict(status='NOT_RUN_NOT_REQUESTED')
            if frozen.get('importance_weights') and record['density']=='explicit' and len(indices)==64 and x is not None:
                importance=stable_importance(phi,x[indices],logq[indices])
            elif frozen.get('importance_weights'):importance=dict(status='NOT_EVALUATED_EXPLICIT_DENSITY_OR64_FIXED_SAMPLES_UNAVAILABLE')
            save(folder/'LIKELIHOOD_VECTOR.npz',sample_indices=np.array(indices),phi=phi,x=np.empty((0,66)) if x is None else x[indices],
                 log_q_x=np.full(len(indices),np.nan) if logq is None else logq[indices],
                 log_weights=importance.get('log_weights',np.full(len(indices),np.nan)),normalized_weights=importance.get('normalized_weights',np.full(len(indices),np.nan)))
            record_result=dict(group=group,parent=parent,sample_indices=indices,successful=int(np.isfinite(phi).sum()),failed=len(failure),failures=failure,
                               fixed_index0_Phi=float(phi[0]) if np.isfinite(phi[0]) else None,
                               importance={k:v for k,v in importance.items() if k not in ('log_weights','normalized_weights')})
            write(folder/'SUMMARY.json',record_result);all_records.append(record_result)
            for name in ('LIKELIHOOD_VECTOR.npz','SUMMARY.json','SAMPLE_RECORDS.jsonl'):
                dest=folder/name;artifacts.append(dict(relative_path=str(dest.relative_to(args.out)),identity=checkpoint_identity(dest)))
            print(json.dumps({k:v for k,v in record_result.items() if k not in ('failures',)}),flush=True)
        write(args.out/'ONLINE_LIKELIHOOD_RESULTS.json',dict(status='COMPLETE',records=all_records,costs=costs,truth_read=False,no_residual_selection=True,posterior_qualification='NOT_ESTABLISHED'))
    except BaseException as exc:
        status='STOPPED' if isinstance(exc,BudgetStop) else 'ERROR';error=repr(exc)
        write(args.out/'ERROR.json',dict(status=status,error=error,traceback=traceback.format_exc(),costs=costs));raise
    finally:
        book.close(status,error)
        write(args.out/'PHYSICS_COMPLETE.json',dict(status=status,pid=os.getpid(),freeze=frozen,artifacts=artifacts,records=all_records,
            completed_time_ns=time.time_ns(),truth_read=False,no_residual_selection=True,costs=costs,receipt=str(book.jobdir/'CLOSED.json'),wall_seconds=time.monotonic()-start,
            failure_policy='all fixed indices retained; failures explicit, no residual selection/resample',posterior_qualification='NOT_ESTABLISHED'))


def score(args,frozen):
    gate=load(args.online/'PHYSICS_COMPLETE.json')
    if gate['status']!='COMPLETE' or gate['truth_read'] is not False or gate['freeze']!=frozen or gate['pid']==os.getpid():raise ValueError('closed physical gate and separate scoring process required')
    if load(gate['receipt'])['status']!='COMPLETE':raise ValueError('physical Book must close before offline score')
    for artifact in gate['artifacts']:
        if checkpoint_identity(args.online/artifact['relative_path'])!=artifact['identity']:raise ValueError('sealed physical output changed')
    if args.out.exists():raise FileExistsError('fresh offline output required')
    if args.truth is not None and 'offline_only' not in args.truth.resolve().parts:raise ValueError('truth must come from isolated offline_only directory')
    # First isolated truth access occurs only below, after all gates.
    groups={};rows=[]
    for i,record in enumerate(frozen['records']):
        folder=args.online/record['group']/record['parent'];summary=load(folder/'SUMMARY.json')
        with np.load(folder/'LIKELIHOOD_VECTOR.npz',allow_pickle=False) as z:phi=z['phi'];indices=z['sample_indices']
        groups.setdefault(record['group'],[]).append(phi[0])
        rank=dict(status='NOT_RUN_ISOLATED_TRUTH_PREDICTION_UNAVAILABLE')
        if args.truth is not None and len(indices)==64 and np.isfinite(phi).all():
            # Offline only reuses an already paid truth prediction/noise. No
            # truth q/chi is forwarded to any physical routine or constructor.
            try:
                with np.load(args.truth/f'{record["parent"]}.npz',allow_pickle=False) as z:truth_prediction=z['prediction']
                with np.load(folder/'sample_000.npz',allow_pickle=False) as z:observation=z['observation']
                if truth_prediction.shape!=(32,8) or not np.iscomplexobj(truth_prediction) or not np.isfinite(truth_prediction).all():raise ValueError('invalid isolated already-paid truth prediction')
                truth_phi=float(np.linalg.norm(truth_prediction-observation)**2/frozen['sigma']**2)
                rank=likelihood_rank(phi,truth_phi,np.random.SeedSequence([20261025099,i]))
            except Exception as exc:rank=dict(status='FAILED_ISOLATED_TRUTH_READ',error=repr(exc))
        rows.append(dict(group=record['group'],parent=record['parent'],physical_summary=summary,data_dependent_likelihood_rank=rank))
    result=dict(status='OFFLINE_LIKELIHOOD_SCORE_COMPLETE',scoring_pid=os.getpid(),physical_pid=gate['pid'],Gamma256={k:gamma256_summary(v) for k,v in groups.items()},
        records=rows,costs=gate['costs'],truth_access='isolated offline evaluator; pre-persisted prediction only',
        posterior_qualification='NOT_ESTABLISHED',no_residual_selection=True,physics_calls=0,optimizer_updates=0,
        ranks='fixed first64 only; duplicates/reject states kept; chain dependence limits rank/IS-ESS interpretation')
    args.out.mkdir(parents=True);write_new_json(args.out/'RESULTS.json',plain(result));print(json.dumps(dict(status=result['status'],groups=list(groups))))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=('freeze','run','score'))
    for key in ('spec','freeze','out','online','truth'):p.add_argument('--'+key,type=Path,required=key=='out')
    p.add_argument('--stage',choices=('I','II','III','IV'));p.add_argument('--job');p.add_argument('--max-rhs',type=int);p.add_argument('--max-wall',type=float)
    a=p.parse_args()
    if a.command=='freeze':freeze(a);return
    frozen=load(a.freeze)
    if frozen.get('status')!='FROZEN_PAID_PROBABILITY_AUDIT' or frozen.get('fixed_sample_index')!=0:raise ValueError('invalid frozen audit')
    if a.command=='run':run(a,frozen)
    else:score(a,frozen)

if __name__=='__main__':main()

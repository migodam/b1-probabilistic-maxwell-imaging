"""Online context generation; cannot load truth or training labels."""
from pathlib import Path
import sys,json,argparse,gc,os,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save
from p3l.inversion import reconstruct

def freeze_frame():
    dest=OUT/'D/FIXED_FRAME.npz'
    if dest.exists():return
    b=Book('D','fixed_train_frame',max_wall=600);status='COMPLETE';error=None
    try:
        t=Target(b);rows=json.loads((OUT/'configs/PARENTS.json').read_text())
        ids=[v['parent'] for v in rows if v['group']=='development' and v['role']=='train']
        risk=np.zeros((64,64))
        for parent in ids:
            p=t.load_cache(OUT/f'A/anchors/{parent}.npz',owner='paid P3L development train anchor; no GT')
            J=t.jacobian(p)[:,:64];D=t.prior.derivative(p.q)[:,:64]
            C=np.linalg.inv(J.T@J+np.eye(64));M=np.real(D.conj().T@D)/576
            risk+=C@M@C
        vals,V=np.linalg.eigh((risk+risk.T)/2/len(ids));selected=V[:,-8:][:,::-1]
        for k in range(8):
            if selected[np.argmax(np.abs(selected[:,k])),k]<0:selected[:,k]*=-1
        Q=np.linalg.qr(selected,mode='complete')[0]
        save(dest,Q=Q,Vw=Q[:,:8],Vs=Q[:,8:],risk_eigenvalues=vals[::-1])
        write(dest.with_suffix('.json'),dict(kind='fixed orthogonal Gaussian64 risk frame, not a posterior support restriction',parents=ids,training_only=True,GT_used=False,weak_dimension=8,chart_determinant_absolute=1,rule='leading eigenvectors of mean(C M C), C=(I+JzT Jz)^-1, M=Re(Dg*Dg)/576; complete Euclidean QR',bounded_properties='remain explicit; network labels use exact probit transport'))
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

def main():
    a=argparse.ArgumentParser();a.add_argument('group',choices=['train','validation','locked','development']);a.add_argument('--start',type=int,default=0);a.add_argument('--count',type=int,default=64);a.add_argument('--evaluations',type=int,default=16);a.add_argument('--resume-token',default='');args=a.parse_args()
    if not (OUT/f'AUDIT/{args.group}_DATA_GENERATOR_EXIT_REQUIRED.json').exists():raise RuntimeError('separate data generator must finish before online contexts')
    if args.group=='locked' and not (OUT/'METHODS_FROZEN_FOR_LOCKED.json').exists():raise RuntimeError('locked method freeze missing')
    freeze_frame();rows=json.loads((OUT/'configs/PARENTS.json').read_text());rows=[v for v in rows if v['group']==args.group][args.start:args.start+args.count]
    suffix=('_'+args.resume_token) if args.resume_token else ''
    if suffix and args.resume_token!='memory_resume1':raise RuntimeError('unregistered recovery token')
    b=Book('D',f'contexts_{args.group}_{args.start:04d}_{args.evaluations}'+suffix,max_wall=7200);status='COMPLETE';error=None
    try:
        t=Target(b,keep_fields=True)
        def memory_check(parent):
            t.backend()._transfers.clear();gc.collect()
            rss_kib=int(subprocess.check_output(['ps','-o','rss=','-p',str(os.getpid())],text=True).strip())
            write(b.jobdir/'MEMORY_CURRENT.json',dict(parent=parent,current_RSS_MiB=rss_kib/1024,native_transfer_entries=len(t.backend()._transfers),paid_workers=1))
            if rss_kib>2500*1024:raise RuntimeError('memory policy stop above2500MiB; preserve paid endpoints')
        for row in rows:
            parent=row['parent'];y=t.observation(parent)
            dest=OUT/f'D/contexts/{parent}.npz';cache_dest=OUT/f'D/contexts/cache/{parent}.npz'
            if dest.exists():
                if not args.resume_token:raise RuntimeError('context already present; no overwrite')
                if dest.with_suffix('.json').exists():old=json.loads(dest.with_suffix('.json').read_text())
                else:
                    receipts=json.loads((OUT/'AUDIT/MEMORY_INTERRUPTION_RECOVERY.json').read_text())['receipts']
                    origins=[r['job'] for r in receipts if r['stage']=='D' and (OUT/f'D/jobs/{r["job"]}/inversion/{parent}/coarse_TRF/best_point.json').exists()]
                    if len(origins)!=1:raise RuntimeError('partial context metadata provenance ambiguous')
                    with np.load(dest,allow_pickle=False) as z:old=dict(parent=parent,source_job=origins[0],GT_used=False,objective=float(z['objective']),context_evaluation_cap=args.evaluations,metadata_recovered_after_atomic_context_save=True)
                    write(dest.with_suffix('.json'),old)
                with np.load(dest,allow_pickle=False) as z:savedq=z['q'].copy();savedf=z['prediction'].copy()
                cache_was_missing=not cache_dest.exists()
                if cache_dest.exists():p=t.load_cache(cache_dest,owner='completed per-parent context from interrupted job; new equation audit')
                else:
                    p=t.evaluate(savedq,'restore_missing_context_factor_after_memory_stop');t.jacobian(p);t.persist_cache(p,cache_dest)
                if not np.array_equal(p.q,savedq) or not np.allclose(p.state.data,savedf,rtol=1e-12,atol=1e-14):raise RuntimeError('recovery cache/context disagreement')
                write(dest.with_suffix('.recovery.json'),dict(parent=parent,original_source_job=old['source_job'],validated_by_job=b.job,status='PAID_COMPLETED_CONTEXT_REUSED_AFTER_MEMORY_STOP',new_physics_for_validation=cache_was_missing,GT_used=False))
                del p;memory_check(parent);continue
            initial=None;previous_attempts=0;origin=None
            if args.resume_token:
                receipts=json.loads((OUT/'AUDIT/MEMORY_INTERRUPTION_RECOVERY.json').read_text())['receipts']
                for receipt in receipts:
                    if receipt['stage']!='D' or not receipt['job'].startswith('contexts_'+args.group+'_'):continue
                    path=OUT/f'D/jobs/{receipt["job"]}/inversion/{parent}/coarse_TRF/best_point.json'
                    if path.exists():
                        if initial is not None:raise RuntimeError('ambiguous interrupted optimizer provenance')
                        oldbest=json.loads(path.read_text());origin=receipt['job']
                        with path.with_name('trace.jsonl').open() as stream:previous_attempts=sum(json.loads(line).get('event')=='EVALUATION_ATTEMPT' for line in stream)
                        if previous_attempts+1>=args.evaluations:raise RuntimeError('interrupted cap leaves no restore-plus-optimizer allowance; stop before any new restore solve')
                        initial=t.evaluate(np.asarray(oldbest['q']),'restore_paid_best_after_memory_stop')
            remaining=args.evaluations-previous_attempts-int(initial is not None)
            if remaining<1:raise RuntimeError('interrupted optimizer cap exhausted; preserve endpoint for separately specified recovery, no implicit extra step')
            p,report=reconstruct(t,y,parent,method='coarse_TRF',max_evaluations=remaining,initial_point=initial,polish=False)
            J=t.jacobian(p);sv=np.linalg.svd(J[:,:64],compute_uv=False)
            scalar=np.r_[t.white(y),t.white(p.state.data-y),p.q,np.log1p(sv*sv)]
            D=t.prior.derivative(p.q)
            image=np.stack((p.state.chi.real,p.state.chi.imag,t.prior.A@p.q[:64],np.sqrt(np.sum(np.abs(D[:,:64])**2,axis=1)))).reshape(4,24,24)
            if dest.exists():raise RuntimeError('context already present; no overwrite')
            save(dest,scalar=scalar,image=image,q=p.q,prediction=p.state.data,J=J,objective=report['objective'],Phi=report['phi'])
            t.persist_cache(p,cache_dest)
            write(dest.with_suffix('.json'),dict(parent=parent,context_evaluation_cap=args.evaluations,GT_used=False,noise_ball_qualified=report['noise_ball_qualified'],objective=report['objective'],KKT=report['kkt']['scaled_KKT'],context_status='fixed cheap online estimate; optimized96 baseline evaluated separately',source_job=b.job,interruption_origin=origin,previous_paid_optimizer_attempts=previous_attempts,remaining_optimizer_cap=remaining,recovery_restore_is_paid=initial is not None))
            print(json.dumps(dict(parent=parent,context=True,Phi=report['phi'],fresh_states=report['fresh_forward_states'])),flush=True)
            del p,initial;memory_check(parent)
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()

"""Root-authorized zero-solve reuse of paid Original Native96 endpoints.

Only q/chi/current/field/prediction/receiver_T/J are loaded. LU/pivots are never
loaded or persisted. Exact Original receiver-J chain rule and cached Maxwell
arithmetic are checked; the original cold phase cost is retained in every report.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_name]='1'
import argparse,json,sys,time
from pathlib import Path
import numpy as np
TRACK=Path(__file__).resolve().parent;OUT=TRACK.parent;REPO=OUT.parents[1]
sys.path.insert(0,str(OUT))
from track_c.evaluate_dev import make_prior


def build():
    if not (TRACK/'REVISION2.json').exists():raise RuntimeError('freeze Revision2 first')
    inputs=TRACK/'matched_inputs_TRAIN512';prior=make_prior()
    r4=REPO.parent/'b1-p3l-r4-conditional-transport/research/p3l_r4_conditional_transport'
    with np.load(r4/'inputs/background/weak.npz',allow_pickle=False) as z:
        G,S,incident=(z[key].copy() for key in ('G','S','incident'))
    sigma=float(prior.config['noise']['sigma'])
    trainparents=[row['parent'] for row in json.loads((inputs/'PARENTS.json').read_text())['parents']]
    devparents=[row['parent'] for row in json.loads((OUT/'configs/parents_CDEV.json').read_text())['parents']]
    jobmap={}
    for cache in (r4/'I/jobs').glob('contexts_train_*/inversion/*/full_TRF/best_cache.npz'):
        parent=cache.parent.parent.name
        if parent in trainparents:
            if parent in jobmap:raise ValueError('duplicate historical native cache')
            jobmap[parent]=cache
    start,cpu=time.perf_counter(),time.process_time();records=[]
    for group,parents in [('TRAIN512',trainparents),('CDEV',devparents)]:
        root=TRACK/'native96_only_contexts'/group
        for i,parent in enumerate(parents):
            destination=root/parent/'native96_only'
            if (destination/'REPORT.json').exists():raise FileExistsError('native reuse never overwrites')
            if group=='TRAIN512':
                cache=jobmap[parent];cost=json.loads((cache.parent/'result.json').read_text())
                online_path=inputs/'online'/f'{parent}.npz'
            else:
                old=TRACK/'contexts/CDEV'/parent/'native96_probit64'
                cache=old/f'original_optimizer/inversion/{parent}/full_TRF/best_cache.npz'
                cost=json.loads((old/'REPORT.json').read_text())['original_native_report']
                online_path=OUT/'raw/online'/f'{parent}.npz'
            tick=time.perf_counter()
            with np.load(cache,allow_pickle=False) as z:
                a={key:z[key].copy() for key in ('q','chi','current','field','prediction','receiver_T')}
                saved_J=z['J'].copy() if 'J' in z.files and z['J'].size else None
            chi=prior.decode(a['q'])[0]
            if not np.array_equal(chi,a['chi']):raise ValueError('native cache decoder mismatch')
            if a['receiver_T'].shape!=(32,576):raise ValueError('native cached receiver transfer unavailable; no unbilled solve')
            field=incident+G@a['current'];rhs=chi[:,None]*incident
            errors=dict(current=float(np.linalg.norm(a['current']-chi[:,None]*(G@a['current'])-rhs)/max(1e-30,np.linalg.norm(rhs))),
                field=float(np.linalg.norm(field-a['field'])/max(1e-30,np.linalg.norm(a['field']))),
                prediction=float(np.linalg.norm(S@a['current']-a['prediction'])/max(1e-30,np.linalg.norm(a['prediction']))),
                receiver_transfer=float(np.linalg.norm(a['receiver_T']-(a['receiver_T']*chi[None,:])@G-S)/max(1e-30,np.linalg.norm(S))))
            if max(errors.values())>1e-10:raise ValueError(f'native paid cache arithmetic failed: {errors}')
            C=np.einsum('qp,ps->qsp',a['receiver_T'],a['field']).reshape(256,576)
            derivatives=C@prior.derivative(a['q'])
            # Physical derivative whitening retains Original float64 arithmetic.
            J=np.column_stack([np.concatenate((derivatives[:,k].reshape(32,8).real.T,
                                               derivatives[:,k].reshape(32,8).imag.T),axis=1).ravel()*np.sqrt(2)/sigma
                               for k in range(66)])
            if saved_J is not None:
                errors['saved_J_relative']=float(np.linalg.norm(J-saved_J)/max(1e-30,np.linalg.norm(saved_J)))
                if errors['saved_J_relative']>1e-12:raise ValueError('Original native savedJ differs from exact chain rule')
            with np.load(online_path,allow_pickle=False) as z:y=z['observation'].copy()
            destination.mkdir(parents=True)
            np.savez_compressed(destination/'CONTEXT.npz',parent=np.array(parent),observation=y,
                                residual=y-a['prediction'],prediction=a['prediction'],q0=a['q'],chi0=chi,J=J)
            report=dict(parent=parent,group=group,config='native96_only',new_anchor_definition=True,
                original_native_objective=cost['target'],original_native_endpoint=cost['endpoint'],
                cold_native_states=cost['fresh_forward_states'],cold_native_J32_calls=cost['receiver_J32_calls'],
                cold_native_RHS=8*cost['fresh_forward_states']+32*cost['receiver_J32_calls'],
                cold_native_wall_seconds=cost['wall_seconds'],new_reuse_RHS=0,new_LU=0,
                build_wall_seconds=time.perf_counter()-tick,validation_errors=errors,
                source_cache=str(cache),LU_read=False,GT_read=False,provider_or_head_updates=0,
                wall_claim='cold Native96 phase cost retained; arithmetic/cache reuse does not imply free inference')
            (destination/'REPORT.json').write_text(json.dumps(report,indent=2)+'\n');records.append(report)
        (root/'ANCHORS_COMPLETE.json').write_text(json.dumps(dict(status='COMPLETE',n=len(parents),parents=parents,config='native96_only',new_RHS=0),indent=2)+'\n')
    receipt=dict(status='NATIVE_CONTEXTS_BUILT',TRAIN512=512,CDEV=64,new_RHS=0,new_LU=0,LU_read=False,
        wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu,
        max_cached_arithmetic_error=max(max(row['validation_errors'].values()) for row in records),
        mean_cold_TRAIN_native_RHS=float(np.mean([r['cold_native_RHS'] for r in records if r['group']=='TRAIN512'])),
        records=records)
    path=TRACK/'native96_only_contexts/BUILD_RECEIPT.json';path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='records'}),flush=True)


if __name__=='__main__':build()

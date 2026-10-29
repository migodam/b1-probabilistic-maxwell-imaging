"""Zero-Maxwell audit replay of saved E1/B current corrector traces.

Uses saved G/S/E0 and original qualified flags solely to reproduce warm starts.
No physics module import or forward call. Preview unless --execute is supplied.
"""
import argparse
import json
import os
from pathlib import Path
from types import SimpleNamespace
for _v in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_v]='1'
import numpy as np
from scipy.linalg import null_space
from accounting import Book,P,write,identity
from charts import CurrentChart


class TraceChart(CurrentChart):
    def residual_jac(self,z,free_values,need_jac=True):
        r,D,decoded=super().residual_jac(z,free_values,need_jac)
        m,j,E=decoded
        self.trace.append(dict(z=z.copy(),free=np.asarray(free_values).copy(),m=m.copy(),j=j.copy(),E=E.copy(),r=r.copy(),
                               D=None if D is None else D.copy(),need_jac=bool(need_jac)))
        return r,D,decoded


def load(path):
    with np.load(path) as z:return {k:z[k].copy() for k in z.files}


def forbidden(*a,**kw):raise RuntimeError('No physical model or forward calls permitted in algebra trace replay')


def run(book,parents):
    destination=P/'raw/E1_CURRENT_TRACE_REPLAY'
    destination.mkdir(parents=True,exist_ok=False)
    sources={}
    def record(path):
        sources[str(path.relative_to(P))]=identity(path)
    for file in [Path(__file__),P/'src/charts.py',P/'src/accounting.py',P/'EXPERIMENT_CONTRACT.json']:
        record(file)
    receipts=[]
    for parent in parents:
        folder=P/'results/E1'/f'P{parent}'/'B'
        reference_path=folder/'REFERENCE.npz';summary_path=folder/'SUMMARY.json'
        record(reference_path);record(summary_path)
        ref=load(reference_path);summary=json.loads(summary_path.read_text())
        model=SimpleNamespace(book=book,G=ref['G'],S=ref['S'],E0=ref['E0'],
                              n=len(ref['m']),ns=ref['E0'].shape[1],nr=ref['S'].shape[0],
                              forward=forbidden,jacobian=forbidden)
        with book.op('AUDIT_REPLAY_RECEIVER_NULLSPACE',algebraic_complex_SVD=1,zero_new_Maxwell_RHS=True):
            model.N=null_space(model.S)
        chart=TraceChart(model,ref['y'],{'j':ref['j']},ref['dep'],ref['free'])
        chart.trace=[]
        zref=chart.encode(ref['j'],ref['m'])
        rootseed={'free':ref['m'][ref['free']].copy(),'z':zref,'key':'REFERENCE'}
        accepted=[rootseed];radial={}
        for row in summary['rows']:
            key=row['key'];values=np.asarray(row['free_values'])
            endpoint=folder/key/'ENDPOINT.npz'
            if row['kind']=='fixed_grid':
                seed=min(accepted,key=lambda a:np.linalg.norm(a['free']-values))
            elif row['kind']=='radial_continuation':
                direction=key.split('_')[1]
                seed=radial.get(direction,rootseed)
            else:raise RuntimeError('Unexpected original trajectory kind')
            chart.trace=[]
            with book.op('AUDIT_CURRENT_CORRECTOR_REPLAY',target=key,parent=parent,zero_new_Maxwell_RHS=True):
                result=chart.correct(seed['z'],values,'AUDIT_REPLAY_'+key)
            comparison={'original_endpoint_present':endpoint.exists(),'original_qualified_used_only_for_seed_policy':bool(row['qualified']),
                        'replayed_corrector_success':bool(result['success']),'original_corrector_success':bool(row['corrector_success'])}
            original=None
            if endpoint.exists():
                record(endpoint);original=load(endpoint)
                dm=float(np.max(abs(result['m']-original['m'])))
                dz=float(np.max(abs(result['z']-original['z'])))
                comparison.update(m_max_absolute_difference=dm,z_max_absolute_difference=dz,
                                  matches_original_at_1e_12=bool(dm<=1e-12 and dz<=1e-12))
            else:comparison['matches_original_at_1e_12']=None
            # Original independently qualified endpoints determine original warm starts.
            # Saved original z avoids turning reconstruction drift into a new path.
            if row['qualified']:
                if original is None:raise RuntimeError('Original qualified endpoint missing')
                entry={'free':values.copy(),'z':original['z'].copy(),'key':key}
                accepted.append(entry)
                if row['kind']=='radial_continuation':radial[direction]=entry
            calls=chart.trace
            arrays=dict(initial_z=seed['z'],free_values=values,final_z=result['z'],
                        final_m=result.get('m',np.full(model.n,np.nan)),reference_z=zref)
            if calls:
                for item in ['z','free','m','j','E','r']:
                    arrays[item]=np.stack([entry[item] for entry in calls])
                jacshape=(2*model.n*model.ns,chart.nc+len(chart.dep))
                arrays['D']=np.stack([entry['D'] if entry['D'] is not None else np.zeros(jacshape) for entry in calls])
                arrays['D_present']=np.asarray([entry['D'] is not None for entry in calls])
                arrays['need_jac']=np.asarray([entry['need_jac'] for entry in calls])
            output=destination/f'P{parent}'/(key+'.npz');output.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(output,**arrays)
            receipt=dict(parent=parent,key=key,warm_seed=seed['key'],trace_calls=len(calls),
                         output=str(output.relative_to(P)),output_blake2b=identity(output),
                         history=result['history'],comparison=comparison)
            write(output.with_suffix('.json'),receipt);receipts.append(receipt)
    checked=[r['comparison']['matches_original_at_1e_12'] for r in receipts if r['comparison']['original_endpoint_present']]
    status=dict(scope='Algebra audit replay only; not new science, physical validation, or cost advantage',
                new_actual_complex_RHS=0,new_closed_form_source_RHS_equivalent=0,parents=parents,
                targets=len(receipts),all_existing_endpoints_match=bool(all(checked)),
                mismatches=[{'parent':r['parent'],'key':r['key'],'comparison':r['comparison']} for r in receipts
                            if r['comparison']['matches_original_at_1e_12'] is False],
                sources_blake2b=sources,output=str(destination.relative_to(P)),
                seed_policy='Original summary qualification and original endpoint z reproduce original grid/radial warm starts; no new success selection',
                lossless_trace_fields=['z','free','m','j','E','r','D','D_present','need_jac'],receipts=receipts)
    write(destination/'REPRODUCIBILITY_RECEIPT.json',status)
    return status


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--execute',action='store_true')
    parser.add_argument('--parents',type=int,nargs='+',default=[0,1,2],choices=[0,1,2]);args=parser.parse_args()
    if not args.execute:
        print('PREVIEW: saved-matrix current traces only; zero new Maxwell RHS; --execute to replay serially');return
    book=Book('E1_CURRENT_TRACE_REPLAY');start=(book.actual,book.formula)
    try:
        status=run(book,args.parents)
        if start!=(book.actual,book.formula):raise RuntimeError('Unexpected physical accounting increment')
    except BaseException as exc:
        book.end('FAILED',error=repr(exc));raise
    else:
        book.end('COMPLETE',new_RHS=0,targets=status['targets'],all_existing_endpoints_match=status['all_existing_endpoints_match'])
        print(json.dumps({k:status[k] for k in ['targets','all_existing_endpoints_match','mismatches','output']},indent=2))

if __name__=='__main__':main()

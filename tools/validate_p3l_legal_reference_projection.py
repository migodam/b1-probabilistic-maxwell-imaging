"""Pay Full forwards for all32 fixed DEV Gaussian-reference projections.

No truth is imported. This exposed-only diagnostic is excluded from the
frozen128 locked estimators. It neither changes nor qualifies the reference.
"""
from pathlib import Path
import json
import shutil
import sys
import numpy as np
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,save


def main():
    folder=OUT/'B/legal_reference_projection'
    source=json.loads((folder/'SUMMARY.json').read_text())
    parents=[f'DEV_{i:03d}' for i in range(32)]
    if source['parents']!=parents:raise RuntimeError('exact all32 projection cohort required')
    closure=json.loads((OUT/'schedules/surrogate_force_hmc/CLOSED.json').read_text())
    if closure.get('status')!='COMPLETE':raise RuntimeError('existing physical sequence must finish first')
    book=Book('B','validate_legal_reference_projection_all32',max_wall=600)
    status,error='ERROR',None
    try:
        write(book.jobdir/'PROTOCOL.json',dict(parents=parents,physics='Full forward only, no Jacobian',
            maximum_new_RHS=256,maximum_states=32,Gaussian_reference='UNQUALIFIED',
            study='EXPOSED_DEV_ONLY_NOT_FROZEN_LOCKED_METHOD',GT_read=False,
            selection='all32 regardless of fit/objective, never according to GT',
            subsequent_scoring='separate evaluator after COMPLETE closure, retain all32 and failures'))
        shutil.copyfile(Path(__file__),book.jobdir/'VALIDATION_SOURCE.py')
        target=Target(book,keep_fields=True);rows=[]
        for parent in parents:
            row=dict(parent=parent,status='FAILURE_RETAINED')
            try:
                with np.load(folder/parent/'RESULT.npz',allow_pickle=False) as data:
                    q=data['q'].copy();start_material=data['material_start'].copy()
                if not target.prior.legal(q):raise ValueError('projection outside original prior')
                y=target.observation(parent);p=target.evaluate(q,'ALL32_EXPOSED_GAUSSIAN_REFERENCE_PROJECTION_VALIDATION')
                Phi=target.phi(p,y)
                target.persist_cache(p,book.jobdir/'cache'/f'{parent}.npz')
                row.update(status='FULL_EVALUATED_NOT_POSTERIOR_QUALIFIED',endpoint=p.key,
                    Phi=Phi,noise_ball=bool(Phi<=.5*chi2.ppf(.95,512)),
                    Gaussian_prior_cost=.5*float(q[:64]@q[:64]),
                    native_objective=Phi+.5*float(q[:64]@q[:64]),
                    properties_u=q[64:].tolist(),
                    material_RMS_from_anchor=float(np.sqrt(np.mean(np.abs(p.state.chi-start_material)**2))),
                    backward_residual=float(p.state.residual))
            except Exception as exc:row['error']=repr(exc)
            rows.append(row)
            write(folder/'FULL_VALIDATION.json',dict(rows=rows,GT_read=False,
                all32_preserved=len(rows)==32,posterior_qualification=False))
        status='COMPLETE'
    except BaseException as exc:error=repr(exc);raise
    finally:book.close(status,error)


if __name__=='__main__':main()

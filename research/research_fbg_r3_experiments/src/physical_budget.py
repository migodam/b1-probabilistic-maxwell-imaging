"""Experimental prepayment of full Maxwell LU/RHS, using the existing ledger.

This is a finite-run budget, not a service or alternative physics backend.
The public-data incumbent is retained when a budget ends between operations.
"""
from contextlib import contextmanager
import numpy as np


class PhysicalBudgetStop(RuntimeError):
    pass


def install_budget(model,Y,max_lu=None,max_rhs=None):
    before=dict(model.cost)
    original_operation=model._physics.operation
    original_state=model.state
    tracker=dict(incumbent=None,candidates=[],max_lu=max_lu,max_rhs=max_rhs)
    scale=float(np.vdot(Y,Y).real)

    @contextmanager
    def operation(kind,rhs=0,**cost):
        lu=model.cost['LU']-before.get('LU',0)
        actual_rhs=model.cost['actual_complex_RHS']-before.get('actual_complex_RHS',0)
        if max_lu is not None and lu+cost.get('LU',0)>max_lu:
            raise PhysicalBudgetStop('FULL_LU_CAP')
        if max_rhs is not None and actual_rhs+rhs>max_rhs:
            raise PhysicalBudgetStop('FULL_COMPLEX_RHS_CAP')
        with original_operation(kind,rhs=rhs,**cost):
            yield

    def state(m):
        z=original_state(m)
        r=z['Y']-Y
        score=float(.5*np.vdot(r,r).real/scale+model.reg.evaluate(m)[0])
        if tracker['incumbent'] is None or score<tracker['incumbent']['physical_score']:
            item=dict(m=np.asarray(m).copy(),prediction=z['Y'].copy(),physical_score=score,
                      relative_residual=float(np.linalg.norm(r)/np.sqrt(scale)),
                      model_cost={k:v-before.get(k,0) for k,v in model.cost.items()})
            tracker['incumbent']=item
            tracker['candidates'].append(item)
        return z

    model._physics.operation=operation
    model.state=state
    return tracker

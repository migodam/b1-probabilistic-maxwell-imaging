"""Actual original-Mie new-action evaluator; never grants a physics budget.

The consumer binds its separately authorized, enforcing ledger and a world's
original shared calibration gain. Imports and factory construction solve no
Maxwell problem. Calls are NOT executed by package export.
"""
import json
from pathlib import Path
import numpy as np
from b1.physics import Forward
from .target import ROOT

def make_mie_action_evaluator(radial_edges, external_budget, gain, *, target_identity,
                              cache_directory, allowed_frequencies=(1.6,2.4)):
    edges=np.asarray(radial_edges,float)
    if edges[0]!=0 or edges[-1]!=1 or np.any(np.diff(edges)<=0):raise ValueError('exact radial full-field edges required')
    if not .7<=gain<=1.3:raise ValueError('one original world gain in original support required')
    if not callable(getattr(external_budget,'evaluate',None)):raise TypeError('separately authorized enforcing physics ledger required')
    config=json.loads((ROOT/'configs/protocol.json').read_text());config['nmax']=14
    backend=Forward(config,external_budget,Path(cache_directory))
    def evaluate(material, action, authorization):
        if not authorization.get('authorized') or authorization.get('target_identity')!=target_identity:
            raise PermissionError('new action needs separate target-matched authorization')
        if authorization.get('remaining_attempts',0)<1 or not authorization.get('ledger_id'):
            raise PermissionError('missing action budget/ledger')
        k=float(action['wavenumber'])
        if k not in allowed_frequencies:raise ValueError('frequency outside frozen action contract')
        eps=np.asarray(material,dtype=complex).reshape(-1)
        if len(eps)!=len(edges)-1 or not np.isfinite(eps).all():raise ValueError('full radial array shape')
        nonvacuum=eps[np.abs(eps-1)>1e-12]
        if not len(nonvacuum) or np.max(abs(nonvacuum-nonvacuum[0]))>1e-12:
            raise ValueError('one shared complex material plus exact vacuum required')
        er,ei=nonvacuum[0].real,nonvacuum[0].imag
        if not (2.2-1e-12<=er<=2.8+1e-12 and .04-1e-12<=ei<=.16+1e-12):
            raise ValueError('material outside original support')
        # Merge identical neighbouring material intervals, preserving exact
        # full material everywhere; union grids need not add fictitious layers.
        last=np.r_[np.flatnonzero(eps[1:]!=eps[:-1]),len(eps)-1]
        radii=edges[last+1];layers=eps[last]
        lo=0.
        for hi,layer in zip(radii,layers):
            if abs(layer-1)<=1e-12:
                from b1_joint.material import RMAX
                central=lo==0 and hi<=RMAX[0]+1e-12
                annular=abs((lo+hi)/2-.65)<=1e-12 and (hi-lo)/2<=RMAX[1]+1e-12
                if not (central or annular):raise ValueError('void geometry outside frozen radial decoder support')
            lo=hi
        result=backend.run(-1,0,0,radii=radii,eps_layers=layers,k=k,nmax=14,
                           kind='RESEARCH_B_SEPARATELY_AUTHORIZED_NEW_ACTION',cache=False)
        f,_,_,meta=result
        return dict(target_identity=target_identity,action_id=action['id'],
                    realification_order='all-real-then-all-imaginary',prediction=(gain*f).tolist(),
                    shared_gain=gain,cost_receipt=dict(new_full_state_frequency_attempts=1,
                    vector_rhs=meta.get('dense_boundary_solves',0),ledger_id=authorization['ledger_id']),
                    validation_status=dict(boundary_backward_residual=meta.get('boundary_backward_residual'),
                    material_scope='same exact radial layers; new original vector-Mie frequency evaluation',
                    probability_semantics='candidate_only',new_action_physics_executed=True))
    return evaluate

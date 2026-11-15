"""R15 spectral material feedback with explicit, paid cyclic state interventions.

The current PR program and its complete-state codec are inherited from R14.
Projected-increment reflection is a new hybrid heuristic, not an exact DR
projection onto a frequency-separated Maxwell material set.
"""
from __future__ import annotations
from collections import Counter
from copy import deepcopy
from pathlib import Path
from types import FunctionType, MethodType
import json, os, time
import numpy as np
from scipy.linalg import svd
import dualmetric_core as inherited
from field_current import DenseFieldCurrentProjectors, _Port
from joint_spectral import conditional_four_spaces
from twofold_denominator import TwofoldDenominator, twofold_material_readout

METHODS=('g19','g19_s19','twofold','twofold_hop','twofold_dr',
         'twofold_ap','twofold_full_dr','twofold_reset')


def write(path,value):
    def encode(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value,default=encode,indent=2,allow_nan=False)+'\n')


class TwofoldG19(inherited.StatefulDualMetric):
    def __init__(self,model,Y,m0,config,method,out,audit_callback=None):
        if method not in METHODS:raise ValueError(method)
        self.method=method;self.out=Path(out);self.audit_callback=audit_callback
        self.feedback_operator=None
        numerical=dict(config,initialize_physical_current=True)
        # Same scalar current transform, original full-derivative material
        # helper substituted below. The internal '01' is dispatch, not A+.
        super().__init__(model,Y,m0,numerical,pair='01')
        self._r15_configure()
        self._check()
        geometry=conditional_four_spaces(self.G,self.S,diagnostics=False)
        self.geometry_metadata=geometry['metadata']
        self.basis_G19=geometry['G_right_frame'][:,:19].copy()
        self.basis_S19=geometry['S_right_frame'][:,:19].copy()
        self.U11=geometry['U11'];self.H_basis=geometry['U01']
        for k,v in self.geometry_metadata['operation_counts'].items():
            self.cost['joint_setup_'+k]+=v
        U=self.basis_G19 if method=='g19' else self.basis_S19 if method=='g19_s19' else self.U11
        self.feedback_operator=TwofoldDenominator(self.G,U,self.shape,
            one_hop=method not in ('g19','g19_s19','twofold'),cost=self.cost)
        self._bind_native_step()
        self.memory_X=np.zeros_like(self.J);self.memory_Q=np.zeros_like(self.Q)
        self.cycle_index=0;self.pending_cycle=None;self.cycle_history=[]
        self.projection_cache={};self.projector=None
        self.cache_path=self.out/'FIXED_OPERATOR_CACHE.npz'
        self._save_fixed_cache()
        write(self.out/'TWOFOLD_SETUP.json',dict(method=method,
            geometry=self.geometry_metadata,denominator=self.feedback_operator.setup_receipt,
            global_G19_is_original_strong_rank=int(self.geometry_metadata['ranks']['global_G_strong']),
            selected_rank=self.feedback_operator.rank,material_real_DOF=self.n,
            exclusion='Excluded from extra rational/routing feedback, not from Maxwell data/material DOF; local residual hop retained',
            reflection='P_H on BOTH scaled current and field increments; full Maxell projectors; not frequency-decoupled feasibility',
            r15_protocol=self.protocol,truth_access='NONE',setup_complete_process_CPU=time.process_time()))
        self._check()

    def _r15_configure(self):
        self.protocol=dict(self.config.get('r15_protocol',{}))
        self.period=int(self.protocol.get('period',16))
        self.scale=float(self.protocol.get('DRscale',6.))
        self.omega=float(self.protocol.get('omega',1.))
        if (self.period,self.scale,self.omega)!=(16,6.,1.):raise ValueError('Frozen period/scale/omega')

    def _bind_native_step(self):
        def readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,config,cost):
            local=dict(config,matrix_denominator=self.feedback_operator,matrix_inner_maxiter=8)
            return twofold_material_readout(blocks,J,Q,m,beta,Os,Ls,rho,bounds,reg,local,cost)
        function=inherited.StatefulDualMetric.step
        scope=dict(function.__globals__,mobius_material_readout=readout)
        cloned=FunctionType(function.__code__,scope,function.__name__,function.__defaults__,function.__closure__)
        self._native_step=MethodType(cloned,self)

    def _P(self,m,rhs,family):
        if family=='0':return super()._P(m,rhs,family)
        return self.feedback_operator.transform(m,self.beta1,rhs)

    def _PH(self,m,rhs,family):
        if family=='0':return super()._PH(m,rhs,family)
        return self.feedback_operator.adjoint(m,self.beta1,rhs)

    def value_material(self,m=None,J=None,Q=None):
        m=self.m if m is None else m
        J,Q=(self.J,self.Q) if J is None else (J,Q)
        return self.feedback_operator.value_gradient(J,Q,m,self.beta1,self.Os,self.rho,self.reg)[0]

    def _save_fixed_cache(self):
        operator=self.feedback_operator
        arrays=dict(U=operator.U,B=operator.B,U11=self.U11,H_basis=self.H_basis,
            basis_G19=self.basis_G19,basis_S19=self.basis_S19,
            retained_neighbor=getattr(operator,'retained_neighbor',np.empty((0,0),complex)))
        np.savez_compressed(self.cache_path,**arrays)

    def _ensure_projector(self):
        if self.projector is not None:return
        self._check()
        self.projector=DenseFieldCurrentProjectors(self.G,self.E,self.S,self.Y,
            float(self.config['r15_data_delta']),scale=self.scale,
            lower=self.lo,upper=self.hi,components=1,grouped_ports=True)
        p=self.projector;arrays=dict(H=p.H[0],H_lower=np.asarray(p.H[1]),
            T=p.ports[0].T,V=p.ports[0].V,eigenvalues=p.ports[0].eigenvalues,
            U=p.ports[0].U,positive=p.ports[0].positive,
            targets=np.column_stack(p.targets),energy_scale=np.asarray(p._energy_scale),
            ball_tolerance=np.asarray(p._ball_tolerance),floor=np.asarray(p.minimum_reachable_residual))
        np.savez_compressed(self.out/'BACKGROUND_PROJECTOR_CACHE.npz',**arrays)
        self._check()

    def _Hproject(self,value):
        if self.method=='twofold_full_dr':return value.copy()
        self.cost['routing_spectral_projection_actions']+=1
        self.cost['routing_spectral_projection_complex_RHS']+=value.shape[1]
        return self.H_basis@(self.H_basis.conj().T@value)

    def _inject_shadow(self,cx,cm):
        """Original affine DPCS/Fourier state, with explicit PR reset and fee."""
        key=(self.stage,self.U.k)
        if key not in self.projection_cache:
            matrix=self.U.apply(np.eye(self.U.k,dtype=complex))
            left,singular,right=svd(matrix,full_matrices=False,check_finite=False)
            rank=int(np.count_nonzero(singular>1e-12*singular[0]))
            self.projection_cache[key]=(left[:,:rank],singular[:rank],right[:rank])
            self.cost['injection_Fourier_projection_SVD']+=1
            self.cost['injection_Fourier_projection_columns']+=self.U.k
        left,singular,right=self.projection_cache[key]
        target=self.scale*cx-self.JD
        self.A=right.conj().T@((left.conj().T@target)/singular[:,None])
        self.cost['injection_projection_SVD_apply_RHS']+=self.ns
        fluctuation=self.U.apply(self.A)
        self.J=self.JD+fluctuation
        self.Q=self.E+self.GJD+self.green.apply(fluctuation)
        self.m=np.asarray(cm,float).copy()
        self.old_grad=self.old_direction=None
        if self.stage_pending or self.done:
            self.old_coefficients={k:a.copy() for k,a in zip(self.frequencies,self.A)}
        gap=self.J-self.scale*cx
        self.cost['injection_receiver_actions']+=1
        self.cost['injection_receiver_complex_RHS']+=self.ns
        return dict(current_projection_relative_loss=float(np.linalg.norm(gap)/max(np.linalg.norm(self.scale*cx),1e-300)),
            projected_port_change_relative=float(np.linalg.norm(self.S@gap)/max(np.linalg.norm(self.Y),1e-300)),
            post_projection_constitutive_defect=float(np.linalg.norm(self.J-self.m[:,None]*self.Q)/max(np.linalg.norm(self.J),1e-300)),
            source_current_real_DOF=int(2*self.ns*len(singular)),PR_history_reset=True)

    def step(self):
        if self.pending_cycle is not None:return self._complete_cycle()
        row=self._native_step()
        row.update(pair_internal_dispatch='01',method=self.method,
            material_transform='twofold_residual_first_hop' if self.feedback_operator.one_hop else 'rank_projected_G19')
        self.history[-1].update(pair_internal_dispatch='01',method=self.method,
            material_transform=row['material_transform'])
        due=self.iterations%self.period==0 or row.get('stage_end',False)
        row['cycle_event']=bool(due)
        if not due:return row
        self.cycle_index+=1
        self.pending_cycle=dict(index=self.cycle_index,row=deepcopy(row),status='BEFORE_CYCLE')
        if self.method in ('twofold_dr','twofold_ap','twofold_full_dr'):
            self._ensure_projector()
        if self.cycle_index<=4:self.checkpoint(self.out/f'CYCLE_{self.cycle_index:03d}_BEFORE.npz')
        return self._complete_cycle()

    def _complete_cycle(self):
        self._check();pending=self.pending_cycle;row=pending['row'];idx=pending['index']
        cpu,wall=time.process_time(),time.perf_counter()
        event=dict(cycle=idx,iteration=self.iterations,stage=self.stage,method=self.method,
            status='RUNNING',memory_policy='retained weak-band normal memory with moving G19 complement anchor',
            routing_dimension=self.n if self.method=='twofold_full_dr' else self.H_basis.shape[1],truth_access='NONE')
        try:
            if self.method in ('twofold_dr','twofold_ap','twofold_full_dr'):
                self._ensure_projector();p=self.projector
                gx,gq=self.J/self.scale,self.Q.copy()
                x,q=gx+self.memory_X,gq+self.memory_Q
                cx,cq,cm=p.material(x,q)
                if self.method=='twofold_ap':
                    ax,aq=p.affine_ball(cx,cq)
                    self.memory_X=self._Hproject(ax-gx)
                    self.memory_Q=self._Hproject(aq-gq)
                else:
                    ax,aq=p.affine_ball(2*cx-x,2*cq-q)
                    self.memory_X+=self.omega*self._Hproject(ax-cx)
                    self.memory_Q+=self.omega*self._Hproject(aq-cq)
                self._check()
                full_gap=float(np.sqrt(np.linalg.norm(ax-cx)**2+np.linalg.norm(aq-cq)**2))
                hx,hq=self._Hproject(ax-cx),self._Hproject(aq-cq)
                block_gap=float(np.sqrt(np.linalg.norm(hx)**2+np.linalg.norm(hq)**2))
                nx,nq=gx+self.memory_X,gq+self.memory_Q
                sx,sq,sm=p.material(nx,nq)
                event.update(full_shadow_gap=full_gap,routed_shadow_gap=block_gap,
                    reflected_memory_norm=float(np.sqrt(np.linalg.norm(self.memory_X)**2+np.linalg.norm(self.memory_Q)**2)),
                    shadow_material_change=float(np.linalg.norm(sm-self.m)),
                    shadow_constitutive_defect=float(np.linalg.norm(self.scale*sx-sm[:,None]*sq)),
                    projector_cost=dict(p.costs),
                    A_multiplier=float(p.last_A_info['multiplier']) if np.isfinite(p.last_A_info.get('multiplier',np.nan)) else None,
                    A_multiplier_pseudoinverse_limit=bool(p.last_A_info.get('pseudoinverse_limit',False)))
                # No GT or immediate-Phi accept/reject. Global paid candidate
                # retention remains the original physical selection rule.
                if self.audit_callback is not None:self.audit_callback(sm,dict(event,source='cyclic post-reflection material shadow'))
                event.update(self._inject_shadow(sx,sm))
            else:
                if self.method=='twofold_reset':
                    self.old_grad=self.old_direction=None;self.cost['scheduled_PR_only_resets']+=1
                    event['PR_history_reset']=True
                if self.audit_callback is not None:self.audit_callback(self.m,dict(event,source='matched cycle state audit'))
            event['status']='COMPLETE'
            self.pending_cycle=None
            row.update(m=self.m.copy(),J=self.J.copy(),Q=self.Q.copy(),coefficients=self.A.copy(),cycle_event=True,cycle_metadata=event)
            return row
        finally:
            event.update(CPU_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall)
            if event['status']=='RUNNING':event['status']='INTERRUPTED_RETAINED_NOT_AN_ELIGIBLE_CANDIDATE'
            self.cycle_history.append(event)
            with (self.out/'CYCLES.jsonl').open('a') as f:f.write(json.dumps(event,allow_nan=False)+'\n')
            if idx<=4 and event['status']=='COMPLETE':self.checkpoint(self.out/f'CYCLE_{idx:03d}_AFTER.npz')

    def state_dict(self):
        state=super().state_dict()
        state['r15']=dict(method=self.method,geometry_metadata=self.geometry_metadata,
            cycle_index=self.cycle_index,pending_cycle=deepcopy(self.pending_cycle),
            memory_X=self.memory_X.copy(),memory_Q=self.memory_Q.copy(),
            cycle_history=deepcopy(self.cycle_history),projection_cache=deepcopy(self.projection_cache),
            projector_cost=None if self.projector is None else dict(self.projector.costs),
            one_hop=self.feedback_operator.one_hop,selected_rank=self.feedback_operator.rank)
        return state

    def checkpoint(self,path):
        path=Path(path);state=self.state_dict()
        state['r15']['fixed_cache_relative']=os.path.relpath(self.cache_path,path.parent)
        background=self.out/'BACKGROUND_PROJECTOR_CACHE.npz'
        state['r15']['background_cache_relative']=os.path.relpath(background,path.parent) if background.is_file() else None
        inherited._save_state_npz(path,state);return path

    @classmethod
    def from_checkpoint(cls,model,Y,config,path,*,out=None,audit_callback=None,
                        method_override=None,reset_reflection_memory=False):
        path=Path(path);state=inherited._load_state_npz(path);details=state.pop('r15')
        obj=super().from_state_dict(model,Y,config,state)
        obj.method=details['method'] if method_override is None else method_override
        if obj.method not in METHODS:raise ValueError(obj.method)
        if obj.method!=details['method']:
            shared={'twofold_hop','twofold_dr','twofold_ap','twofold_full_dr','twofold_reset'}
            if not {obj.method,details['method']}<=shared:
                raise ValueError('Causal override requires the same frozen denominator and one-hop objective')
        obj.out=path.parent if out is None else Path(out);obj.out.mkdir(exist_ok=True,parents=True)
        obj.audit_callback=audit_callback;obj._r15_configure()
        obj.cache_path=(path.parent/details['fixed_cache_relative']).resolve()
        with np.load(obj.cache_path,allow_pickle=False) as z:cache={k:z[k].copy() for k in z.files}
        obj.U11=cache['U11'];obj.H_basis=cache['H_basis'];obj.basis_G19=cache['basis_G19'];obj.basis_S19=cache['basis_S19']
        operator=TwofoldDenominator.__new__(TwofoldDenominator)
        operator.blocks=None;operator.cost=obj.cost;operator.U=cache['U'];operator.B=cache['B']
        operator.n,operator.rank=operator.U.shape;operator.projection='galerkin';operator.one_hop=details['one_hop']
        operator.setup_receipt=dict(restored_from_checkpoint=str(path),selected_rank=operator.rank,one_hop=operator.one_hop)
        if operator.one_hop:operator.blocks=obj.blocks;operator.retained_neighbor=cache['retained_neighbor']
        obj.feedback_operator=operator;obj._bind_native_step();obj.geometry_metadata=details['geometry_metadata']
        obj.memory_X=details['memory_X'];obj.memory_Q=details['memory_Q']
        if reset_reflection_memory:obj.memory_X[:]=0;obj.memory_Q[:]=0
        obj.cycle_index=details['cycle_index'];obj.pending_cycle=details['pending_cycle']
        obj.cycle_history=details['cycle_history'];obj.projection_cache=details['projection_cache'];obj.projector=None
        relative=details['background_cache_relative']
        if relative is not None:
            with np.load(path.parent/relative,allow_pickle=False) as z:a={k:z[k].copy() for k in z.files}
            p=DenseFieldCurrentProjectors.__new__(DenseFieldCurrentProjectors)
            p.G=obj.G;p.E=obj.E;p.B=obj.scale*obj.G;p.scale=obj.scale;p.delta=float(config['r15_data_delta'])
            p.components=1;p.cells=obj.n;p.sources=obj.ns;p.lower=obj.lo;p.upper=obj.hi;p.grouped_ports=True
            p.H=(a['H'],bool(a['H_lower']));port=_Port(a['T'],a['V'],a['eigenvalues'],a['U'],a['positive'])
            p.ports=[port]*obj.ns;p.port_groups=[(port,np.arange(obj.ns))];p.targets=[a['targets'][:,i] for i in range(obj.ns)]
            p._energy_scale=float(a['energy_scale']);p._ball_tolerance=float(a['ball_tolerance']);p.minimum_reachable_residual=float(a['floor'])
            p.costs=details['projector_cost'];p.last_A_info={};p.public_metadata={};obj.projector=p
        return obj

"""Protocol adapter for the fixed-operator research prototype."""
from dataclasses import dataclass, replace
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from time import perf_counter

import numpy as np
from flint import acb_mat, arb, ctx

PAPER = Path(__file__).resolve().parents[1]
A, B, C = (PAPER/'experiments'/name for name in 'ABC')
for directory in (C, B, A, PAPER/'src'):
    sys.path.insert(0, str(directory))

import bounds as a_bounds
import compressed as a_chart
import certify as b_checks
import design as c_design
import fixed_trace as trace
import model as c_model

if __package__:
    from .rectangular_chart import AffineEqualities, evaluate_chart, pullback_border, rectangular_equality_charts
else:
    from rectangular_chart import AffineEqualities, evaluate_chart, pullback_border, rectangular_equality_charts

OPERATORS = ('EXACT_STORED_FLOAT','INTERVAL_DISCRETIZATION','SAMPLE_MODEL_ERROR')
OUTCOMES = ('EXCLUDED','WITNESS_FOUND','UNRESOLVED')


def _module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _immutable(a, dtype=None):
    data = np.array(a, dtype=dtype, copy=True)
    if not np.all(np.isfinite(data)):
        raise ValueError('All array inputs must be finite.')
    data.flags.writeable = False
    return data


def _nonnegative(x, name):
    x = float(x)
    if not np.isfinite(x) or x < 0:
        raise ValueError(f'{name} must be finite and nonnegative.')
    return x


def _coords(grid, n):
    if not grid or int(np.prod(grid)) != n:
        raise ValueError('Grid must account for every modeled cell.')
    return np.indices(grid).reshape(len(grid),n).T


@dataclass(frozen=True)
class ModelContract:
    G: np.ndarray
    S: np.ndarray
    Bdict: np.ndarray
    grid: tuple
    units: str
    operator: str = OPERATORS[0]
    time_convention: str = 'exp(-i omega t)'

    def __post_init__(self):
        for name in ('G','S','Bdict'):
            object.__setattr__(self,name,_immutable(getattr(self,name),complex))
        n = self.G.shape[0]
        if self.G.shape != (n,n) or self.S.ndim != 2 or self.Bdict.ndim != 2:
            raise ValueError('Incompatible fixed operator dimensions.')
        if self.S.shape[1] != n or self.Bdict.shape[0] != n:
            raise ValueError('Receiver and source operators must use the same cells.')
        if self.operator not in OPERATORS or self.time_convention != 'exp(-i omega t)' or not self.units:
            raise ValueError('Unrecognized operator contract, time convention, or units.')
        _coords(self.grid,n)

    @property
    def model_hash(self):
        meta = [*map(str,self.grid),self.units,self.operator,self.time_convention]
        return c_model.hash_arrays(self.G,self.S,self.Bdict,np.asarray(meta))


@dataclass(frozen=True)
class MaterialNode:
    node_id: str
    contract_id: str
    ell: np.ndarray
    u: np.ndarray
    parent_id: str | None = None
    candidate_centers: tuple = ()
    children: tuple = ()
    region_label: str = 'all modeled cells'
    project_version: str = '1'
    status: str = 'UNRESOLVED'

    def __post_init__(self):
        for name in ('ell','u'):
            a = np.asarray(getattr(self,name))
            if np.iscomplexobj(a):
                raise ValueError('Only real contrast boxes are implemented.')
            object.__setattr__(self,name,_immutable(a,float))
        if self.ell.ndim != 1 or self.u.shape != self.ell.shape or np.any(self.ell>self.u):
            raise ValueError('Material boxes require ordered endpoints on each cell.')
        if self.status not in OUTCOMES:
            raise ValueError('Unrecognized material node outcome.')

    @property
    def box_hash(self):
        return c_model.hash_arrays(self.ell,self.u)


@dataclass(frozen=True)
class NoiseContract:
    r: float
    source: str
    per_column_r: float | None = None

    def __post_init__(self):
        object.__setattr__(self,'r',_nonnegative(self.r,'Frobenius radius'))
        if not self.source:
            raise ValueError('The Frobenius radius requires a proof source.')
        if self.per_column_r is not None:
            object.__setattr__(self,'per_column_r',_nonnegative(self.per_column_r,'Column radius'))

    def omega(self, columns):
        # E*E ≼ ||E||_F² I, including all cross terms.
        return trace.diag([arb(self.r)**2]*columns)


@dataclass(frozen=True)
class Acquisition:
    model: ModelContract
    V: np.ndarray
    Y: np.ndarray
    noise: NoiseContract
    receiving_basis_id: str
    energies: tuple
    power_gram: np.ndarray | None = None
    block_hashes: tuple = ()

    def __post_init__(self):
        for name in ('V','Y'):
            object.__setattr__(self,name,_immutable(getattr(self,name),complex))
        if (not self.receiving_basis_id or self.V.ndim != 2 or self.Y.ndim != 2
                or self.V.shape[0] != self.model.Bdict.shape[1]
                or self.Y.shape != (self.model.S.shape[0],self.V.shape[1])):
            raise ValueError('Acquisition source and receiving dimensions do not agree.')
        if len(self.energies)!=self.V.shape[1]:
            raise ValueError('Acquisition requires one declared energy per measured column.')
        object.__setattr__(self,'energies',tuple(_nonnegative(x,'Beam energy') for x in self.energies))
        if self.power_gram is not None:
            power = _immutable(self.power_gram,complex)
            if power.shape != (self.V.shape[0],)*2 or not np.array_equal(power,power.conj().T):
                raise ValueError('A declared power Gram must be Hermitian and index the source dictionary.')
            object.__setattr__(self,'power_gram',power)

    @property
    def data_hash(self):
        gram = np.empty((0,0)) if self.power_gram is None else self.power_gram
        return c_model.hash_arrays(self.V,self.Y,gram,np.asarray(self.energies),
                                   np.asarray([self.receiving_basis_id,self.noise.source,
                                               repr(self.noise.r),repr(self.noise.per_column_r)]))


@dataclass(frozen=True)
class BorderProposal:
    K: acb_mat
    L: acb_mat
    C: acb_mat
    full: acb_mat
    errors: dict
    hashes: dict
    ell: np.ndarray
    u: np.ndarray
    w: np.ndarray
    optical: np.ndarray

    @property
    def n(self):
        return self.K.nrows()


@dataclass(frozen=True)
class Kernel:
    base: BorderProposal
    trace_rows: acb_mat
    affine: acb_mat
    reduced: trace.FixedKernel
    extra: acb_mat | None = None
    pair_indices: tuple | None = None


@dataclass(frozen=True)
class PairRows:
    T: acb_mat
    E: acb_mat
    indices: tuple
    data_hash: str
    model_hash: str

    def __iter__(self):
        return iter((self.T,self.E))


@dataclass(frozen=True)
class AChartKernel:
    bank: a_chart.TraceBank
    geometry_hash: str


@dataclass(frozen=True)
class ReducedBorder:
    full: acb_mat
    n: int
    eliminated: int
    errors: dict
    hashes: dict


def _stored(model):
    if model.operator!=OPERATORS[0]:
        raise ValueError('This adapter requires exact stored arrays; use the analytic C adapter for its interval model.')


def _weights(weights,n):
    raw, optical=weights
    if np.iscomplexobj(raw) or np.iscomplexobj(optical):
        raise ValueError('Sector and optical weights must be real.')
    w, v=[np.broadcast_to(np.asarray(k,float),(n,)) for k in (raw,optical)]
    if np.any(w<0) or not (np.all(np.isfinite(w)) and np.all(np.isfinite(v))):
        raise ValueError('The real sector requires finite nonnegative weights and real optical multipliers.')
    return w,v


def assemble_border(model,box,weights,Z):
    """26. Assemble the exact stored real-sector border. Z is an Acquisition."""
    _stored(model)
    n=model.G.shape[0]
    if Z.model.model_hash!=model.model_hash or box.ell.shape!=(n,):
        raise ValueError('Box, acquired data, and model must share their declared grid and operator.')
    w, optical=_weights(weights,n)
    G,S,B,V,Y=map(trace.mat,(model.G,model.S,model.Bdict,Z.V,Z.Y))
    full=trace.assemble_border(G,S,B,V,Y,box.ell,box.u,w,optical)
    state=range(n)
    source=range(n,full.nrows())
    return BorderProposal(trace.block(full,state,state),-trace.block(full,state,source),
                          trace.block(full,source,source),full,
                          dict(method='Arb enclosure of every entry',precision_bits=ctx.prec),
                          dict(model_hash=model.model_hash,data_hash=Z.data_hash,box_hash=box.box_hash,
                               weights_hash=c_model.hash_arrays(w,optical)),box.ell,box.u,
                          _immutable(w),_immutable(optical))


def local_pair_rows(model,indices,acquisition):
    """27. Affine [j_i,e_i] rows in caller's fixed pixel order."""
    _stored(model)
    ids=tuple(int(i) for i in indices)
    if len(set(ids))!=len(ids) or any(i<0 or i>=model.G.shape[0] for i in ids):
        raise ValueError('Trace indices must be unique in-grid cells.')
    if acquisition.model.model_hash!=model.model_hash:
        raise ValueError('Source dictionary and acquisition differ.')
    if not ids:
        return PairRows(acb_mat(0,model.G.shape[0]),acb_mat(0,acquisition.V.shape[1]),ids,
                        acquisition.data_hash,model.model_hash)
    rows=trace.local_pair_rows(trace.mat(model.G),trace.mat(model.Bdict),trace.mat(acquisition.V),ids)
    return PairRows(*rows,ids,acquisition.data_hash,model.model_hash)


def build_fixed_kernel(base,trace_rows,verification=None):
    """28. Split off only a verified positive exact trace-null complement."""
    T,E=trace_rows
    indices=None
    if isinstance(trace_rows,PairRows):
        if trace_rows.data_hash!=base.hashes['data_hash'] or trace_rows.model_hash!=base.hashes['model_hash']:
            raise ValueError('Trace rows came from another acquisition.')
        indices=trace_rows.indices
    extra=None if verification is None else verification.get('extra_exact_covectors')
    if not T.nrows():
        return dict(status='UNRESOLVED',reason='The adapter requires retained physical trace rows.')
    try:
        reduced=trace.build_fixed_kernel(base.full,T,E,base.n,extra=extra)
    except (ArithmeticError, ValueError, ZeroDivisionError, np.linalg.LinAlgError) as exc:
        return dict(status='UNRESOLVED',reason=str(exc))
    return Kernel(base,T,E,reduced,extra,indices)


def update_box(kernel,local_D,Z):
    """29. Keep the entire affine trace quadratic including the source corner."""
    if Z.data_hash!=kernel.base.hashes['data_hash'] or Z.model.model_hash!=kernel.base.hashes['model_hash']:
        raise ValueError('A changed acquisition requires a new base border.')
    delta=local_D if isinstance(local_D,acb_mat) else trace.mat(local_D)
    if delta.nrows()!=kernel.trace_rows.nrows() or delta.ncols()!=delta.nrows():
        raise ValueError('Local update must index the retained physical trace pair rows.')
    return ReducedBorder(trace.update_box(kernel.reduced,delta),kernel.reduced.dim,
                         kernel.reduced.eliminated,kernel.base.errors,kernel.base.hashes)


def certify_inertia(border,threshold,frame_proposal=None):
    """30. Verify both signs independently on supplied or proposed frames."""
    t,j,r=threshold
    t,r=_nonnegative(t,'Threshold'),_nonnegative(r,'Radius')
    if int(j)!=j or not 1<=j<=border.full.nrows()-border.n:
        raise ValueError('j must index a retained source dimension.')
    shifted=border.full-trace.diag([0]*border.n+[arb(t)]*(border.full.nrows()-border.n))
    start=perf_counter()
    values,vectors=np.linalg.eigh(b_checks.herm(b_checks.midpoint(shifted)))
    proposals={'positive':vectors[:,values>0],'negative':vectors[:,values<0]}
    frames = proposals if frame_proposal is None else frame_proposal
    rec=dict(t=t,j=int(j),radius=r,state_dimension=border.n,precision_bits=ctx.prec,
             proposal_positive=int(sum(values>0)),proposal_negative=int(sum(values<0)),
             eliminated_positive=getattr(border,'eliminated',0),status='UNRESOLVED',**border.hashes)
    used={}
    for sign,value in (('positive',1),('negative',-1)):
        f=np.asarray(frames.get(sign,np.empty((len(values),0))),complex)
        if f.ndim!=2 or f.shape[0]!=len(values) or not np.all(np.isfinite(f)):
            raise ValueError('Sign frame must have finite columns in the threshold form.')
        gap, gram=b_checks.sign_floor(shifted,f,value)
        rec[f'verified_{sign}_at_least']=f.shape[1] if gap is not None else 0
        rec[f'{sign}_gap_lower']=gap
        rec[f'{sign}_frame_gram_lower']=gram
        if gap is not None:
            used[sign]=f
    sq=arb(int(j))*arb(t)
    if rec['verified_positive_at_least']>=border.n+j and sq>=arb(r)**2:
        rec.update(status='EXCLUDED',residual_lower=max(0.,b_checks.lower(sq.sqrt())))
    rec['unresolved_signs']=len(values)-rec['verified_positive_at_least']-rec['verified_negative_at_least']
    rec['interval_check_seconds']=perf_counter()-start
    return rec, used


def certify_box(node,acquisition,kernel_bank,budget):
    """31. Finite predeclared weight/threshold attempts, or independently check a witness."""
    records=[]
    if node.ell.shape!=(acquisition.model.G.shape[0],):
        raise ValueError('Material node must specify every modeled cell.')
    trials=budget.get('trials',())
    for index,trial in enumerate(trials):
        if index>=budget.get('max_trials',len(trials)):
            break
        if 'kernel' in trial:
            kernel=kernel_bank[trial['kernel']]
            border=_node_update(node,acquisition,kernel,trial['weights'])
            if border is None:
                border=assemble_border(acquisition.model,node,trial['weights'],acquisition)
        else:
            border=assemble_border(acquisition.model,node,trial['weights'],acquisition)
        result,_=certify_inertia(border,(trial['t'],trial['j'],acquisition.noise.r))
        records.append(dict(node_id=node.node_id,**result))
        if result['status']=='EXCLUDED':
            return dict(status='EXCLUDED',records=records)
    if 'material' in budget:
        witness=verify_witness(node,acquisition,budget['material'])
        if witness['status']=='WITNESS_FOUND':
            return dict(status='WITNESS_FOUND',records=records,witness=witness)
    return dict(status='UNRESOLVED',records=records)


def _node_update(node,acquisition,kernel,weights):
    if not isinstance(kernel,Kernel) or not kernel.pair_indices:
        return None
    base=kernel.base
    w,v=_weights(weights,base.n)
    ids=np.asarray(kernel.pair_indices,dtype=int)
    changed=np.flatnonzero((node.ell!=base.ell)|(node.u!=base.u)|(w!=base.w)|(v!=base.optical))
    if not np.all(np.isin(changed,ids)):
        return None
    q=len(ids)
    identity=np.eye(2*q)
    jj,ee=map(trace.mat,(identity[::2],identity[1::2]))
    def local(ell,upper,weight,opt):
        return a_bounds.sector_ball(jj,ee,ell[ids],upper[ids],weight[ids],opt[ids])
    delta=local(node.ell,node.u,w,v)-local(base.ell,base.u,base.w,base.optical)
    border=update_box(kernel,delta,acquisition)
    return replace(border,hashes=dict(border.hashes,box_hash=node.box_hash,
                                      weights_hash=c_model.hash_arrays(w,v)))


def extend_traces(kernel,new_rows):
    """32. Rebuild one joint exact kernel rather than combining separate facts."""
    add_t,add_e=new_rows
    start=perf_counter()
    T=trace.mat(kernel.trace_rows.tolist()+add_t.tolist())
    E=trace.mat(kernel.affine.tolist()+add_e.tolist())
    joint=(T,E)
    if kernel.pair_indices is not None and isinstance(new_rows,PairRows):
        ids=kernel.pair_indices+new_rows.indices
        if len(set(ids))!=len(ids):
            return dict(kernel=kernel,elapsed_seconds=perf_counter()-start,
                        previous_rows=kernel.trace_rows.nrows(),joint_rows=T.nrows(),
                        status='UNRESOLVED',reason='Overlapping named rows require explicit deduplication.')
        joint=PairRows(T,E,ids,new_rows.data_hash,new_rows.model_hash)
    updated=build_fixed_kernel(kernel.base,joint,dict(extra_exact_covectors=kernel.extra))
    return dict(kernel=updated,elapsed_seconds=perf_counter()-start,
                previous_rows=kernel.trace_rows.nrows(),joint_rows=T.nrows())


def _measured_incident(acquisition):
    _stored(acquisition.model)
    return trace.mat(acquisition.model.Bdict)*trace.mat(acquisition.V)


def prepare_a_chart_kernel(acquisition,background_cap,weight):
    """Bind A's existing certified bank to the actual fixed geometry."""
    return AChartKernel(a_chart.TraceBank(acquisition,background_cap,weight),
                        c_model.hash_arrays(acquisition.G,acquisition.S))


def exact_row_oracle(node,selected_rows,acquisition,kernel,budget):
    """33. Reuse A's verified identity-pivot chart on its frozen joint trace region."""
    a = acquisition
    pos={int(cell):i for i,cell in enumerate(a_chart.REGION)}
    if not isinstance(a,a_bounds.Acquisition) or not isinstance(kernel,AChartKernel):
        return dict(status='UNRESOLVED',reason='This implementation uses the A acquisition and TraceBank.')
    if kernel.geometry_hash!=c_model.hash_arrays(a.G,a.S):
        raise ValueError('The A chart bank and acquisition geometry differ.')
    bank=kernel.bank
    if not bank.available or not selected_rows or len(set(selected_rows))!=len(selected_rows) or not all(i in pos for i in selected_rows):
        return dict(status='UNRESOLVED',reason='No verified joint A trace bank for the requested rows.')
    outside=np.ones(a.n,dtype=bool)
    outside[a_chart.REGION]=False
    if (node.ell.shape!=(a.n,) or not np.array_equal(node.ell[outside],bank.ell[outside])
            or not np.array_equal(node.u[outside],bank.upp[outside])):
        return dict(status='UNRESOLVED',reason='Outside-trace material bounds differ from the shared A bank.')
    stored={}
    def save(frame):
        key=f'frame_{len(stored)}'
        stored[key]=frame
        return key
    params=tuple(budget['weights'])
    _weights((params[0]*params[1],params[0]*params[2]),1)
    if float(params[0])!=float(bank.mu):
        raise ValueError('The shared trace bank and chart must use the same base weight.')
    args=(a,bank.base(a,False),bank.base(a,True),node.ell,node.u,
          tuple(pos[i] for i in selected_rows),params,save)
    proof=a_chart.interval_search(*args,perf_counter()+budget.get('seconds',1.0),
                                 max_depth=budget.get('max_depth',4),max_leaves=budget.get('max_leaves',16))
    return dict(status='EXCLUDED' if proof['excluded'] else 'UNRESOLVED',proof=proof,frames=stored)


def build_task_cover(task_contract):
    """34. B's finite lattice-center cover and one-step closed refinement."""
    b_run=_module('_zalg_b_run',B/'run.py')
    spec=dict(task_contract)
    contract_id=spec.pop('contract_id')
    contract=b_run.Contract(**spec)
    x0,x1,y0,y1=contract.window
    if not (0<=x0<=x1<contract.side and 0<=y0<=y1<contract.side):
        raise ValueError('Task window must enumerate valid grid centers.')
    def wrap(win, parent=None):
        l,u=contract.cover_box(win)
        candidates=tuple(b_run.members(win))
        child_ids=tuple(','.join(map(str,w)) for w in b_run.quadrants(win)) if len(candidates)>1 else ()
        return MaterialNode(','.join(map(str,win)),contract_id,l,u,parent,
                            candidates,child_ids,contract.mode)
    def refine(node):
        c=node.candidate_centers
        win=(min(x for x,y in c),max(x for x,y in c),min(y for x,y in c),max(y for x,y in c))
        children=tuple(wrap(w,node.node_id) for w in b_run.quadrants(win)) if len(c)>1 else ()
        if children and set(c)!=set().union(*(set(x.candidate_centers) for x in children)):
            raise ArithmeticError('Center refinement lost a registered center.')
        return children
    return dict(root=wrap(tuple(contract.window)),refine=refine,
                coverage='Declared finite lattice centers; endpoints are coordinatewise outer bounds.')


def verify_witness(node,acquisition,material):
    """35b. Independent exact original forward solve, including all known rows."""
    a=acquisition
    if np.iscomplexobj(material):
        raise ValueError('The physical witness requires real material at every modeled cell.')
    material=_immutable(material,float)
    if material.shape!=node.ell.shape or material.shape!=(a.model.G.shape[0],):
        raise ValueError('Witness must specify every modeled cell.')
    status=dict(status='UNRESOLVED',node_id=node.node_id,model_hash=a.model.model_hash,
                data_hash=a.data_hash,box_hash=node.box_hash,material=material,
                precision_bits=ctx.prec,radius=a.noise.r,material_in_box=bool(np.all(
                    (material>=node.ell)&(material<=node.u))))
    if not status['material_in_box']:
        return status
    start=perf_counter()
    G,S,Y=map(trace.mat,(a.model.G,a.model.S,a.Y))
    m=[arb(x) for x in material]
    try:
        state=(trace.diag([1]*len(m))-trace.diag(m)*G).solve(trace.diag(m)*_measured_incident(a))
        residual=c_model.norm(S*state-Y)
        status.update(status='WITNESS_FOUND' if residual<=arb(a.noise.r) else 'UNRESOLVED',
                      verified_original_square_solve=True,residual_upper=b_checks.upper(residual))
    except (ArithmeticError,ValueError,ZeroDivisionError) as exc:
        status['reason']=str(exc)
    return dict(status,interval_seconds=perf_counter()-start)


def search_witness(node,acquisition,fixed_budget):
    """35a. Proposal via B's bounded least-squares search, then independent verification."""
    a=acquisition
    _stored(a.model)
    if len(a.model.grid)!=2 or a.model.grid[0]!=a.model.grid[1] or np.any(node.ell<0):
        return dict(status='UNRESOLVED',reason='B search adapter uses square grids and nonnegative boxes.')
    candidate=(node.ell+node.u)/2
    search=dict(method='box midpoint',forward_evaluations=0)
    started=perf_counter()
    if a.noise.r>0 and np.linalg.norm(a.Y)>0 and fixed_budget.get('max_evaluations',0)>0:
        try:
            g=b_checks.Geometry(a.model.grid[0],a.model.S,a.model.Bdict@a.V,green=a.model.G)
            candidate,_,search=b_checks.witness_search(g,a.Y,node.ell,node.u,a.noise.r,
                                      max_evaluations=fixed_budget['max_evaluations'])
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as exc:
            search['reason']=str(exc)
    audit=verify_witness(node,a,candidate)
    return dict(candidate=candidate,search=search,witness=audit,status=audit['status'],
                seconds=perf_counter()-started)


def build_response_envelope(node,certificate,acquisition):
    """36. C's analytic model; exact old physical columns and proved Frobenius bound."""
    a=acquisition
    g=c_model.geometry(a['layout'])
    y=np.asarray(a['Y'],complex)
    if (y.shape!=(g.s.shape[0],len(c_model.OLD))
            or tuple(a['source_indices'])!=c_model.OLD or node.ell.shape!=(g.g.shape[0],)):
        return dict(status='UNRESOLVED',reason='Outside the C analytic geometry or old acquisition definition.')
    if certificate.get('status')=='EXCLUDED':
        return dict(status='UNRESOLVED',reason='An excluded class has no prediction.')
    energy=_nonnegative(a['per_beam_energy'],'Old beam energy')
    if energy==0:
        return dict(status='UNRESOLVED',reason='Old physical columns must have positive fluence.')
    z=c_design.source_pair(c_model.OLD,arb(energy).sqrt(),y)
    # C's function takes per-column scale and reconstructs the old Frobenius bound exactly.
    sigma=arb(a['noise'].r)/arb(len(c_model.OLD)).sqrt()
    weight=_nonnegative(certificate['weight'],'Sector weight')
    border, state=c_design.form_border(g,node.ell,node.u,weight)
    if border is None:
        return dict(status='UNRESOLVED',proof=state)
    record, proof=c_design.envelope(border,g,z,sigma)
    return dict(status='BOUNDED_ENVELOPE' if record else 'UNRESOLVED',envelope=record,proof=proof,
                node_id=node.node_id,box_hash=node.box_hash,layout=g.layout,
                data_hash=c_model.hash_arrays(y,np.asarray([energy,a['noise'].r])))


def score_action(action,envelopes,noise_contract):
    """37. Verify one native two-class direction, preserving the shared source Gram."""
    valid=len(envelopes)==2 and all(e.get('envelope') is not None for e in envelopes)
    if not valid or noise_contract.per_column_r is None:
        return dict(status='UNRESOLVED',worst_survivors_upper=len(envelopes))
    if not (envelopes[0]['layout']==envelopes[1]['layout']==action['layout'] and
            envelopes[0]['data_hash']==envelopes[1]['data_hash']):
        raise ValueError('The two envelopes and action must use one old acquisition and geometry.')
    ids=tuple(action['source_indices'])
    if len(ids) not in c_model.BATCHES or len(set(ids))!=len(ids) or not all(
            i in range(len(c_model.ANGLES)) and i not in c_model.OLD for i in ids):
        raise ValueError('Action is outside C registered source menu.')
    sigma=max(arb(noise_contract.per_column_r),
              (arb(noise_contract.r)/arb(len(ids)).sqrt()).upper())
    g=c_model.geometry(action['layout'])
    V=c_model.minor(g.ex_b,range(g.g.shape[0]),(*c_model.OLD,*ids))
    rank_floor=c_model.low(c_model.spectrum(c_model.star(V)*V)[0][0])
    if rank_floor<=0:
        return dict(status='UNRESOLVED',reason='Joint incident source rank is unresolved.',worst_survivors_upper=2)
    amp=arb(_nonnegative(action['per_beam_energy'],'New beam energy')).sqrt()
    if not amp>0:
        return dict(status='UNRESOLVED',reason='New physical columns must have positive fluence.',worst_survivors_upper=2)
    pair=[e['envelope'] for e in envelopes]
    _,w=c_design.split_query(pair,ids,amp,sigma)
    if w is None:
        return dict(status='UNRESOLVED',worst_survivors_upper=2)
    proof=c_design.certify_pair(pair,w,ids,amp,sigma)
    return dict(status='STRICT_EDGE' if proof['strict_edge'] else 'UNRESOLVED',
                worst_survivors_upper=proof['worst_survivors_bound'],proof=proof,
                direction=w,direction_hash=c_model.hash_arrays(w),source_gram_floor=rank_floor)


def clique_cover_bound(labels,pair_checks,areas=None):
    """Greedy disjoint clique cover using only independently verified strict edges."""
    labels=tuple(labels)
    if len(set(labels))!=len(labels):
        raise ValueError('Duplicate class identifiers.')
    adjacency={label:set() for label in labels}
    for (a,b),check in pair_checks.items():
        if a==b or a not in adjacency or b not in adjacency:
            raise ValueError('A pair must reference two distinct declared classes.')
        candidate=check.get('proof',check)
        margin=candidate.get('noisy_separation_gap_lower',0.)
        if candidate.get('strict_edge') is True and np.isfinite(margin) and margin>0.:
            adjacency[a].add(b)
            adjacency[b].add(a)
    cliques=[]
    for label in sorted(labels,key=lambda v:-len(adjacency[v])):
        group=next((g for g in cliques if all(other in adjacency[label] for other in g)),None)
        if group is None:
            cliques.append([label])
        else:
            group.append(label)
    receipt=dict(cliques=cliques,verified_edges=sum(map(len,adjacency.values()))//2,
                 worst_survivors_upper=len(cliques))
    if areas is not None:
        if set(areas)!=set(labels):
            raise ValueError('Provide one nonnegative area per class.')
        bounds=[max(_nonnegative(areas[v],'Class area') for v in group) for group in cliques]
        receipt['worst_disjoint_area_upper']=c_model.high(sum(map(arb,bounds),arb(0)))
    return receipt


def select_and_freeze_action(experiment_path=C,diagnostic=False):
    """38. Call the registered C planner in its own process with no future acquisition."""
    p=Path(experiment_path).resolve()
    base=p/'diagnostic' if diagnostic else p
    if (base/'acquired').exists() and not (base/'proposal.json').exists():
        raise RuntimeError('Future data appeared before the proposal.')
    if not (base/'proposal.json').exists():
        env=dict(os.environ,ZALG_EXPERIMENT_C_DIAGNOSTIC='1' if diagnostic else '0',
                 OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
        subprocess.run([sys.executable,'-B',str(p/'design.py')],env=env,check=True,capture_output=True,text=True)
    filename=base/'proposal.json'
    record=c_model.read_json(filename)
    if record['status']!='SEALED_READ_ONLY_BEFORE_ANY_NEW_DATA':
        raise ValueError('Proposal is not in the experiment C sealed format.')
    if (base/'acquired/acquisition.json').exists():
        receipt=c_model.read_json(base/'acquired/acquisition.json')
        if receipt['frozen_proposal_hash']!=c_model.hash_file(filename):
            raise ValueError('Archived acquisition references a different proposal.')
    return dict(path=str(filename),proposal=record,proposal_hash=c_model.hash_file(filename))


def append_observation(old,new,expected_old_hash):
    """39. Append original physical measurements as recorded and combine disjoint error budgets."""
    if old.data_hash!=expected_old_hash:
        raise ValueError('Previously frozen acquisition changed.')
    if old.model.model_hash!=new.model.model_hash or old.receiving_basis_id!=new.receiving_basis_id:
        raise ValueError('Old and appended data must share the original receiving basis and model.')
    if (old.power_gram is None)!=(new.power_gram is None) or old.power_gram is not None and not np.array_equal(
            old.power_gram,new.power_gram):
        raise ValueError('Incident power Gram changed.')
    radius=c_model.high((arb(old.noise.r)**2+arb(new.noise.r)**2).sqrt())
    noise=NoiseContract(radius,'disjoint old/new source columns and their existing Frobenius bounds')
    result=Acquisition(old.model,np.column_stack((old.V,new.V)),np.column_stack((old.Y,new.Y)),
                       noise,old.receiving_basis_id,old.energies+new.energies,old.power_gram,
                       old.block_hashes+(old.data_hash,new.data_hash))
    m=old.V.shape[1]
    if not np.array_equal(old.Y,result.Y[:,:m]) or not np.array_equal(old.V,result.V[:,:m]):
        raise ArithmeticError('Old observation block changed on append.')
    return result


def replay_certificate(record,original_inputs):
    """40. Reassemble the original full stored-model border and check saved frames."""
    model,node,weights,acquisition=(original_inputs[name] for name in ('model','node','weights','acquisition'))
    _stored(model)
    n=model.G.shape[0]
    w,v=_weights(weights,n)
    hashes=dict(model_hash=model.model_hash,data_hash=acquisition.data_hash,box_hash=node.box_hash,
                weights_hash=c_model.hash_arrays(w,v))
    if record.get('eliminated_positive',0)!=0:
        return dict(status='UNRESOLVED',reason='A compressed frame requires its dedicated A/B replay.')
    if acquisition.model.model_hash!=model.model_hash or any(record.get(key)!=value for key,value in hashes.items()) or record.get('radius')!=acquisition.noise.r:
        raise ValueError('Stored proof input or threshold contract changed.')
    G,S,Y=map(trace.mat,(model.G,model.S,acquisition.Y))
    m=Y.ncols()
    current=c_model.hjoin(trace.diag([1]*n),acb_mat(n,m))
    field=c_model.hjoin(G,_measured_incident(acquisition))
    receiver=c_model.hjoin(S,-Y)
    full=trace.adj(receiver)*receiver+a_bounds.sector_ball(current,field,node.ell,node.u,w,v)
    border=ReducedBorder(full,n,0,dict(precision_bits=ctx.prec),hashes)
    frames=original_inputs.get('frames',{})
    receipt,_=certify_inertia(border,(record['t'],record['j'],acquisition.noise.r),frames)
    receipt['independent_reassembly']=True
    return receipt

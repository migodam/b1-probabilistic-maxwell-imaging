"""Paid P1B GL65-offset x GL5-contrast full-Maxwell material node bank.

Import does no work. main uses the SAME append-only P1B Book, plus a persistent
15000-RHS scope wrapper. Exact P1A float-hex nodes may be reused only after
literal public-family/walls/physics-source equality; no P1A Book is constructed.
No observations or truth are read. This is a finite quadrature, not an exact
continuous posterior or a new material prior.
"""
from __future__ import annotations
import fcntl
import json
import os
from pathlib import Path
import numpy as np
from .accounting import ROOT,Book,BudgetStop,write
from .probability_family import Family

OUT=ROOT/'results/P1B/full_material_rule'
SCOPE='P1B_FULL_MATERIAL_GL65X5_V1'
SCOPE_RHS_CAP=15000
NU,NC=65,5


def material_rule():
    zu,wu=np.polynomial.legendre.leggauss(NU);zc,wc=np.polynomial.legendre.leggauss(NC)
    u=.01*zu;c=3+.06*zc;b,uu,cc=np.meshgrid(np.arange(2),u,c,indexing='ij')
    weights=.5*np.broadcast_to((wu/2)[None,:,None]*(wc/2)[None,None,:],(2,NU,NC))
    centers=np.array([-.9,-.8])
    return dict(branches=b.ravel(),u=uu.ravel(),contrast=cc.ravel(),x=centers[b.ravel()]+uu.ravel(),
        weights=weights.ravel(),u_axis=u,contrast_axis=c,u_weights=wu/2,contrast_weights=wc/2,shape=np.array([2,NU,NC]))


def _key(b,u,c):return f'{int(b)}|{float(u).hex()}|{float(c).hex()}'


def physics_sources():
    files=[ROOT/'src/wr'/n for n in ('physics.py','native_world.py','physical_witness.py','probability_family.py')]
    files += [ROOT.parent/'research_w_continuation/src/wq'/n for n in ('fem.py','meshes.py')]
    files += [ROOT.parent/'research_w_continuation/vendor/p0/geometry.py']
    return files


def _contract(directory,family,*,create=False):
    directory=Path(directory);snapshot=directory/'PUBLIC_FAMILY_SNAPSHOT.json';walls=directory/'PUBLIC_WALLS_SNAPSHOT.npz'
    source_dir=directory/'SOURCE_CONTRACT'
    if create:
        directory.mkdir(parents=True,exist_ok=True);source_dir.mkdir(exist_ok=True)
    if snapshot.exists():
        if json.loads(snapshot.read_text())!=family.config:raise RuntimeError('Bank public-family/acquisition mismatch')
    elif create:write(snapshot,family.config)
    else:raise RuntimeError('Bank public-family snapshot absent')
    if walls.exists():
        with np.load(walls,allow_pickle=False) as z:
            if any(not np.array_equal(z[f'q{i}'],q) for i,q in enumerate(family.walls)):
                raise RuntimeError('Bank public wall vectors differ')
    elif create:np.savez_compressed(walls,q0=family.walls[0],q1=family.walls[1])
    else:raise RuntimeError('Bank public wall snapshot absent')
    for i,path in enumerate(physics_sources()):
        saved=source_dir/f'{i}_{path.name}'
        if saved.exists():
            if saved.read_bytes()!=path.read_bytes():raise RuntimeError('Literal physics source changed:'+str(path))
        elif create:saved.write_bytes(path.read_bytes())
        else:raise RuntimeError('Literal physics contract incomplete:'+str(saved))


def _rows(path):
    if not path.exists():return []
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def _append(path,row):
    with Path(path).open('a') as stream:
        fcntl.flock(stream,fcntl.LOCK_EX);stream.write(json.dumps(row,allow_nan=False)+'\n')
        stream.flush();os.fsync(stream.fileno())


def _response(path,node):
    path=Path(path).resolve()
    if not (path.is_relative_to(OUT.resolve()) or path.is_relative_to((ROOT/'results/P1A/bank').resolve())):
        raise RuntimeError('Response path must be inside a public paid-node bank')
    with np.load(path,allow_pickle=False) as z:
        F=z['F'].copy()
        if (int(z['branch'])!=node['branch'] or float(z['u'])!=node['u'] or float(z['contrast'])!=node['contrast']
            or F.shape!=(5,4,12) or not np.iscomplexobj(F) or not np.isfinite(F).all()):
            raise RuntimeError('Public response coordinates/shape/finite contract mismatch')
    return F


class ScopeBook:
    """Extra scope quota only; stage closure/budgets remain in existing Book.

    A separate scope lock serializes quota checks and parent reservations.
    Prior attempts are counted from the SAME persistent append-only ledger.
    No stage guard/cap is copied or bypassed here.
    """
    def __init__(self,parent,directory):
        if parent.stage!='P1B':raise ValueError('Full material bank must use Book(P1B)')
        self.parent=parent;self.lock=Path(directory)/'SCOPE.lock';self.node_key=None;self.attempt=None
        self.scope=SCOPE;self.scope_cap_rhs=SCOPE_RHS_CAP
    def __getattr__(self,name):return getattr(self.parent,name)
    def used_rhs(self):
        return sum(r.get('rhs',0) for r in _rows(self.parent.path)
            if r.get('event')=='reserve' and r.get('stage')=='P1B' and r.get('scope')==SCOPE)
    def reserve(self,phase,label,rhs=0,factorizations=0,**metadata):
        with self.lock.open('a') as stream:
            fcntl.flock(stream,fcntl.LOCK_EX);used=self.used_rhs()
            if used+rhs>SCOPE_RHS_CAP:
                self.parent.event(event='scope_budget_reject',scope=SCOPE,scope_cap_rhs=SCOPE_RHS_CAP,
                    scope_used_rhs=used,requested_rhs=rhs,node_key=self.node_key,attempt=self.attempt)
                raise BudgetStop('Persistent full-material scope15000 RHS exceeded')
            return self.parent.reserve(phase,label,rhs=rhs,factorizations=factorizations,
                scope=SCOPE,scope_cap_rhs=SCOPE_RHS_CAP,scope_node_key=self.node_key,scope_attempt=self.attempt,**metadata)
    def event(self,**row):
        self.parent.event(scope=SCOPE,scope_node_key=self.node_key,scope_attempt=self.attempt,**row)


class FullBank:
    def __init__(self):
        self.family=Family.load();OUT.mkdir(parents=True,exist_ok=True)
        paid=bool(_rows(OUT/'INDEX.jsonl')) or any(OUT.glob('physical_*'))
        _contract(OUT,self.family,create=not paid)
        self.book=ScopeBook(Book('P1B'),OUT);self.index=OUT/'INDEX.jsonl'
        self.entries=_rows(self.index);self.cache={r['key']:r for r in self.entries}
        # Read P1A PUBLIC cache without constructing Bank() or any P1A Book.
        self.p1a_root=ROOT/'results/P1A/bank';p1a_rows=_rows(self.p1a_root/'INDEX.jsonl')
        if p1a_rows:_contract(self.p1a_root,self.family,create=False)
        self.p1a={r['key']:r for r in p1a_rows}
        self._recover_completed_attempts()
    def _recover_completed_attempts(self):
        for folder in sorted(OUT.glob('physical_*')):
            if not (folder/'POINT_COMPLETE.json').exists() or not (folder/'FULL_RESPONSE.npz').exists():continue
            public=json.loads((folder/'PUBLIC_PRIOR_NODE.json').read_text())['coordinates']
            node=self.family.coordinates(public['branch'],public['u'],public['contrast'])
            key=_key(node['branch'],node['u'],node['contrast'])
            if key in self.cache:continue
            _response(folder/'FULL_RESPONSE.npz',node)
            row=dict(key=key,branch=node['branch'],u=node['u'],contrast=node['contrast'],
                path=str((folder/'FULL_RESPONSE.npz').relative_to(ROOT)),source_stage='P1B',
                original_point_rhs=20,recovered_completed_attempt=True,response_kind='FULL_CACHED_NODE')
            _append(self.index,row);self.entries.append(row);self.cache[key]=row
    def point(self,b,u,c):
        node=self.family.coordinates(b,u,c);key=_key(node['branch'],node['u'],node['contrast'])
        if key in self.cache:
            row=self.cache[key];F=_response(ROOT/row['path'],node)
            return F,{**row,'access_incremental_rhs':0,'cache_hit':True}
        if key in self.p1a:
            old=self.p1a[key];source=self.p1a_root/old['path'];F=_response(source,node)
            row=dict(key=key,branch=node['branch'],u=node['u'],contrast=node['contrast'],path=str(source.relative_to(ROOT)),
                source_stage='P1A',original_point_rhs=old.get('original_point_rhs',old.get('incremental_rhs',20)),
                response_kind='FULL_CACHED_NODE',exact_floathex_reuse=True)
            _append(self.index,row);self.entries.append(row);self.cache[key]=row
            return F,{**row,'access_incremental_rhs':0,'cache_hit':True}
        existing=[int(p.name.split('_')[-1]) for p in OUT.glob('physical_*') if p.name.split('_')[-1].isdigit()]
        folder=OUT/f'physical_{max(existing,default=-1)+1:05d}'
        self.book.node_key=key;self.book.attempt=folder.name
        try:
            F=self.family.point(node['branch'],node['u'],node['contrast'],self.book,str(folder.relative_to(ROOT/'results')))
        finally:self.book.node_key=None;self.book.attempt=None
        # Family.point already persisted complete per-frequency fields and F.
        F=_response(folder/'FULL_RESPONSE.npz',node)
        row=dict(key=key,branch=node['branch'],u=node['u'],contrast=node['contrast'],
            path=str((folder/'FULL_RESPONSE.npz').relative_to(ROOT)),source_stage='P1B',original_point_rhs=20,
            response_kind='FULL_CACHED_NODE',recovered_completed_attempt=False)
        _append(self.index,row);self.entries.append(row);self.cache[key]=row
        return F,{**row,'access_incremental_rhs':20,'cache_hit':False}


def load_full_rule():
    """Read-only completed public full bank, without constructing any Book."""
    if not (OUT/'COMPLETE.json').exists():raise RuntimeError('P1B full material bank incomplete')
    family=Family.load();_contract(OUT,family,create=False)
    with np.load(OUT/'FULL_RULE.npz',allow_pickle=False) as z:rule={k:z[k].copy() for k in z.files}
    expected=material_rule()
    for key,value in expected.items():
        if not np.array_equal(rule[key],value):raise RuntimeError('Full material quadrature changed:'+key)
    if rule['F'].shape!=(2*NU*NC,5,4,12) or not np.iscomplexobj(rule['F']) or not np.isfinite(rule['F']).all():
        raise RuntimeError('Invalid complete full-Maxwell material node arrays')
    metadata=json.loads((OUT/'FULL_RULE.json').read_text())
    if metadata.get('response_kind')!='FULL_CACHED_NODE' or metadata.get('material_rule')!=[2,NU,NC]:
        raise RuntimeError('Full rule metadata contract mismatch')
    return rule,metadata


def main():
    gate=ROOT/'results/P1A/reference/GATE.json'
    if not gate.exists() or json.loads(gate.read_text()).get('probability_reference_qualified') is not True:
        raise RuntimeError('Paid P1B full bank blocked: P1A probability_reference_qualified must beTrue')
    protocol=json.loads((ROOT/'configs/P1_PROBABILITY_PROTOCOL.json').read_text())
    spec=protocol['P1B']['finite_material_rule']
    if spec['offset_GL_order']!=NU or spec['contrast_GL_order']!=NC or spec['branches']!=2:
        raise RuntimeError('Pre-observation protocol must specify full GL65xGL5x2')
    bank=FullBank();before=bank.book.used_rhs();rule=material_rule();data=[];nodes=[]
    if (OUT/'FULL_RULE.npz').exists():
        # Validate an existing final aggregate before reusing or finishing marker.
        with np.load(OUT/'FULL_RULE.npz',allow_pickle=False) as z:
            for key,value in rule.items():
                if not np.array_equal(z[key],value):raise RuntimeError('Existing full material rule differs:'+key)
            if z['F'].shape!=(650,5,4,12) or not np.isfinite(z['F']).all():raise RuntimeError('Invalid existing aggregate F')
    for index,(b,u,c) in enumerate(zip(rule['branches'],rule['u'],rule['contrast'])):
        F,row=bank.point(int(b),float(u),float(c));data.append(F);nodes.append(row)
        write(OUT/'PROGRESS.json',dict(completed=index+1,total=len(rule['weights']),
            scope_used_rhs=bank.book.used_rhs(),scope_cap_rhs=SCOPE_RHS_CAP,
            current_run_new_rhs=bank.book.used_rhs()-before,
            reused_nodes_so_far=sum(r['access_incremental_rhs']==0 for r in nodes),response_kind='FULL_CACHED_NODE'))
    F=np.asarray(data)
    if (OUT/'FULL_RULE.npz').exists():
        with np.load(OUT/'FULL_RULE.npz',allow_pickle=False) as z:
            if not np.array_equal(z['F'],F):raise RuntimeError('Existing aggregate differs from bound node cache')
    else:
        temporary=OUT/'FULL_RULE.npz.tmp'
        with temporary.open('wb') as stream:
            np.savez_compressed(stream,F=F,**rule);stream.flush();os.fsync(stream.fileno())
        os.link(temporary,OUT/'FULL_RULE.npz');temporary.unlink()
    metadata=dict(family_id=bank.family.config['family_id'],response_kind='FULL_CACHED_NODE',
        source_rule=str((OUT/'FULL_RULE.npz').relative_to(ROOT)),material_rule=[2,NU,NC],nodes=nodes,
        physical_prior_changed=False,interpolation_used=False,
        interpretation='paid complete Maxwell node vectors; finite quadrature approximates unchanged continuous controlled prior',
        scope=SCOPE,scope_cap_rhs=SCOPE_RHS_CAP,scope_used_rhs=bank.book.used_rhs(),
        current_run_new_rhs=bank.book.used_rhs()-before,
        original_node_cost_rhs=sum(r['original_point_rhs'] for r in nodes),
        current_run_reuse_count=sum(r['access_incremental_rhs']==0 for r in nodes),
        source_contract='SOURCE_CONTRACT literal physics bytes plus public family/wall snapshots',
        original_minimum_feature='UNRESOLVED',truth_access=False)
    write(OUT/'FULL_RULE.json',metadata)
    write(OUT/'COMPLETE.json',dict(status='FULL_MATERIAL_RULE_COMPLETE',nodes=len(F),material_rule=[2,NU,NC],
        scope=SCOPE,scope_used_rhs=bank.book.used_rhs(),current_run_new_rhs=bank.book.used_rhs()-before,
        global_P1B_totals=bank.book.totals()['P1B'],response_kind='FULL_CACHED_NODE',
        continuous_posterior_exact_claim=False,truth_access=False))
    return metadata


if __name__=='__main__':main()

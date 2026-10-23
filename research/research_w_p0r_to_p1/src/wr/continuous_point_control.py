"""Additional observation-only continuous point baseline; import does no work.

Array interpolation proposes ONE point per branch. Only actual Full responses
select the winner. Neither proxy stationarity nor a paid point certifies a
continuous global MAP or Full-objective KKT. No truth/label files are accessed.
"""
from __future__ import annotations
import fcntl
import json
import shutil
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from .accounting import ROOT, Book, BudgetStop, write
from .prob_reference import lagrange
from .progressive_inference import source_statistics, profile_shared_gamma

SCOPE='P1B_CONTINUOUS_POINT_CONTROL_V1'
CAP=1280
OUT=ROOT/'results/P1B/continuous_point_control'
SCALE=np.array([.01,.06])
OFFSET=np.array([0.,3.])


def key(b,u,c):return f'{int(b)}|{float(u).hex()}|{float(c).hex()}'
def rows(path):
    return [json.loads(s) for s in Path(path).read_text().splitlines() if s.strip()] if Path(path).exists() else []


def profile(F,y):
    """Same unique shared gain box across all frequencies, sources, receivers."""
    p=profile_shared_gamma(source_statistics(np.asarray(F)[None],np.asarray(y)))
    a=float(p['amplitude'][0]);phase=float(p['phase'][0]);g=a*np.exp(1j*phase)
    # Direct residual avoids subtracting nearly equal large squared norms.
    phi=float(.5*np.sum(abs(g*F-y)**2)/.01**2)
    return phi,a,phase


def basis(nodes,z):
    """Existing barycentric Lagrange basis and its exact polynomial derivative."""
    nodes=np.asarray(nodes,float);L=lagrange(nodes,np.array([z]))[0]
    delta=nodes[:,None]-nodes[None,:];np.fill_diagonal(delta,1.)
    w=1/np.prod(delta,axis=1);hits=np.flatnonzero(abs(z-nodes)<1e-14)
    if len(hits):
        i=int(hits[0]);d=np.zeros(len(nodes));use=np.arange(len(nodes))!=i
        d[use]=w[use]/(w[i]*(nodes[i]-nodes[use]));d[i]=-d.sum()
    else:
        inv=1/(z-nodes);d=L*(np.dot(L,inv)-inv)
    return L,d


class ArrayProxy:
    def __init__(self,rule,branch,y):
        self.y=np.asarray(y);self.branch=int(branch);use=np.asarray(rule['branches'])==branch
        us=np.unique(rule['u'][use]);cs=np.unique(rule['contrast'][use])
        if len(us)!=9 or len(cs)!=9 or np.sum(use)!=81:raise ValueError('Expected paid full9 tensor rule')
        ids=np.flatnonzero(use)
        order=[next(i for i in ids if rule['u'][i]==u and rule['contrast'][i]==c) for u in us for c in cs]
        self.cube=np.asarray(rule['F'])[order].reshape(9,9,5,4,12)
        if not np.isfinite(self.cube).all():raise ValueError('Nonfinite paid responses')
        self.un=us/.01;self.cn=(cs-3)/.06;self.evaluations=0
    def response(self,z):
        lu,du=basis(self.un,float(z[0]));lc,dc=basis(self.cn,float(z[1]))
        F=np.einsum('i,j,ijfrs->frs',lu,lc,self.cube)
        dF=np.stack([np.einsum('i,j,ijfrs->frs',du,lc,self.cube),np.einsum('i,j,ijfrs->frs',lu,dc,self.cube)])
        return F,dF
    def fun(self,z):
        self.evaluations+=1;F,dF=self.response(z);phi,a,p=profile(F,self.y)
        gamma=a*np.exp(1j*p);res=gamma*F-self.y
        grad=np.array([np.real(np.vdot(res,gamma*d))/.01**2 for d in dF])
        return phi,grad
    def optimize(self):
        start_cpu,start_wall=time.process_time(),time.perf_counter()
        grid=[]
        for u in self.un:
            for c in self.cn:
                z=np.array([u,c]);grid.append((self.fun(z)[0],z))
        zg=min(grid,key=lambda t:t[0])[1]
        # One property-boundary start, chosen deterministically from two endpoint
        # response evaluations at the best grid offset; no labels involved.
        boundary=[np.array([zg[0],s]) for s in (-1.,1.)]
        zb=min(boundary,key=lambda z:self.fun(z)[0])
        starts=[zg,np.zeros(2),zb];candidates=[];records=[]
        for label,z in zip(('full9_best','center','property_boundary'),starts):
            cpu,wall=time.process_time(),time.perf_counter();initial=self.fun(z)[0]
            opt=minimize(self.fun,z,jac=True,method='L-BFGS-B',bounds=[(-1,1)]*2,
                         options=dict(maxiter=150,maxfun=1000,ftol=1e-12,gtol=1e-7,maxls=30))
            if not np.isfinite(opt.fun) or not np.isfinite(opt.x).all():raise RuntimeError('Nonfinite array optimizer result')
            candidates.extend([(initial,z.copy()),(float(opt.fun),opt.x.copy())])
            records.append(dict(start=label,z=z.tolist(),initial_Phi=initial,final_Phi=float(opt.fun),
                success=bool(opt.success),message=str(opt.message),nit=int(opt.nit),nfev=int(opt.nfev),
                cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall))
        _,z=min(candidates,key=lambda t:t[0]);phi,g=self.fun(z);F,_=self.response(z)
        pg=float(np.linalg.norm(z-np.clip(z-g,-1,1),ord=np.inf))
        uc=np.clip(OFFSET+SCALE*z,[-.01,2.94],[.01,3.06])
        return dict(u=float(uc[0]),contrast=float(uc[1]),Phi=phi,z=z.tolist(),
            proxy_projected_gradient_scaled=pg,proxy_gradient_scaled=g.tolist(),
            proxy_KKT_only=True,actual_Full_KKT_certified=False,starts=records,
            array_evaluations=self.evaluations,cpu_seconds=time.process_time()-start_cpu,
            wall_seconds=time.perf_counter()-start_wall),F


class ControlBook:
    """Persistent root-ledger sublimit, including failed attempts; never resets."""
    def __init__(self,parent):
        self.parent=parent;self.node_key=None
        self.lock=ROOT/'P1B_CONTINUOUS_POINT_CONTROL_SCOPE.lock'
    def __getattr__(self,name):return getattr(self.parent,name)
    def local(self):return [r for r in rows(self.parent.path) if r.get('stage')=='P1B' and r.get('scope')==SCOPE]
    def used(self):return sum(r.get('rhs',0) for r in self.local() if r.get('event')=='reserve')
    def begin(self,k):
        with self.lock.open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX);rr=self.local()
            attempts=[r for r in rr if r.get('event')=='continuous_full_vector_attempt']
            if any(r['point_key']==k for r in attempts):raise RuntimeError('Previously attempted point: automatic retry forbidden')
            if len(attempts)>=64 or self.used()+20>CAP:raise BudgetStop('Persistent continuous control allowance insufficient for whole vector')
            self.node_key=k
            self.parent.event(event='continuous_full_vector_attempt',scope=SCOPE,point_key=k,vector_max=64,scope_rhs_cap=CAP)
    def reserve(self,phase,label,rhs=0,factorizations=0,**metadata):
        with self.lock.open('a') as f:
            fcntl.flock(f,fcntl.LOCK_EX)
            if self.node_key is None:raise RuntimeError('Missing full-vector attempt marker')
            if self.used()+rhs>CAP:raise BudgetStop('Persistent continuous control RHS cap')
            return self.parent.reserve(phase,label,rhs,factorizations,scope=SCOPE,scope_node_key=self.node_key,**metadata)
    def event(self,**metadata):self.parent.event(scope=SCOPE,scope_node_key=self.node_key,**metadata)


def _point_record(b,u,c,F,y,source,response_error=0.,potential_error=0.):
    phi,a,p=profile(F,y)
    return dict(branch=int(b),u=float(u),contrast=float(c),amplitude=a,phase=p,
                x=float((-.9 if b==0 else -.8)+u),Phi=phi,physical_source=str(source),
                response_error_white=float(response_error),potential_error=float(potential_error),
                actual_full_response=True,actual_Full_KKT_certified=False)


def main():
    """Explicit serial-queue entry only. Existing/partial output is never retried."""
    from .probability_family import Family
    from .progressive_full_bank import _contract,physics_sources
    start_cpu,start_wall=time.process_time(),time.perf_counter()
    addpath=ROOT/'configs/P1B_STRONG_POINT_CONTROL_ADDENDUM.json';add=json.loads(addpath.read_text())
    if add['persistent_scope']!=SCOPE or add['scope_rhs_cap']!=CAP or add['new_full_vector_max']!=64:raise RuntimeError('Frozen allowance mismatch')
    gate=json.loads((ROOT/'results/P1A/reference/GATE.json').read_text())
    if not gate.get('probability_reference_qualified'):raise RuntimeError('P1A probability gate closed')
    inf=ROOT/'results/P1B/inference';pv=ROOT/'results/P1B/point_validation'
    if not ((inf/'COMPLETE.json').exists() or (inf/'GATE.json').exists()):raise RuntimeError('Inference not frozen')
    complete=json.loads((pv/'COMPLETE.json').read_text())
    if complete.get('status')!='FULL_POINT_VALIDATION_COMPLETE' or complete.get('observations')!=32:raise RuntimeError('Full original point validation incomplete')
    inference=json.loads((inf/'RESULTS.json').read_text());validated=json.loads((pv/'RESULTS.json').read_text())
    with np.load(ROOT/'results/P1B/data/OBSERVATIONS.npz',allow_pickle=False) as obs:
        ys=obs['y'].copy();parents=obs['parent_ids'].copy()
    if ys.shape!=(32,5,4,12) or not np.isfinite(ys).all() or len(inference)!=32 or len(validated)!=32:raise RuntimeError('Expected 32 frozen observations')
    _,counts=np.unique(parents,return_counts=True)
    if len(counts)!=16 or not np.all(counts==2):raise RuntimeError('Expected 16 parents with two observations each')
    for records in (inference,validated):
        if sorted(r['observation'] for r in records)!=list(range(32)):raise RuntimeError('Observation identity mismatch')
        if any(int(parents[r['observation']])!=r['parent'] for r in records):raise RuntimeError('Parent identity mismatch')
    family=Family.load();bank=ROOT/'results/P1A/bank';_contract(bank,family,create=False)
    with np.load(ROOT/'results/P1A/rules/full_9.npz',allow_pickle=False) as z:rule={k:z[k].copy() for k in z.files}
    # Index only public paid Full endpoints; no Bank object or P1A Book.
    sources={}
    for directory in (bank,ROOT/'results/P1B/full_material_rule'):
        entries=rows(directory/'INDEX.jsonl')
        if entries:_contract(directory,family,create=False)
        for r in entries:
            path=Path(r['path']);path=path if path.is_absolute() else ((ROOT/path) if str(path).startswith('results/') else directory/path)
            sources[r['key']]=path
    # Existing successful point-validation attempts are also eligible exact caches.
    # Read only tiny endpoint responses/metadata, never field arrays or labels.
    for source in sorted(pv.glob('observation_*/*/FULL_RESPONSE.npz')):
        if not (source.parent/'POINT_COMPLETE.json').exists():continue
        with np.load(source,allow_pickle=False) as z:
            sources[key(int(z['branch']),float(z['u']),float(z['contrast']))]=source
    for source in sorted(pv.glob('observation_*/*_FINAL.npz')):
        with np.load(source,allow_pickle=False) as z:
            pp=z['params'];k=key(*pp[:3])
        sources.setdefault(k,source)
    baseline={}
    for r in validated:
        i=int(r['observation']);path=pv/f'observation_{i:03d}/direct_joint_point_FINAL.npz'
        with np.load(path,allow_pickle=False) as z:F=z['F'].copy();p=z['params'].copy()
        m=r['methods']['direct joint point']
        if key(*p[:3])!=key(m['branch'],m['u'],m['contrast']):raise RuntimeError('Validated baseline coordinate mismatch')
        if F.shape!=(5,4,12) or not np.isfinite(F).all():raise RuntimeError('Invalid actual baseline response')
        baseline[i]=(p,F,path);sources.setdefault(key(*p[:3]),path)
    OUT.mkdir(parents=True,exist_ok=False) # partial/failure requires parent review; no automatic retry
    _contract(OUT,family,create=True)
    sourceout=OUT/'CONTROL_SOURCE';sourceout.mkdir()
    for n,path in enumerate([Path(__file__),ROOT/'src/wr/prob_reference.py',ROOT/'src/wr/progressive_inference.py',*physics_sources(),addpath,ROOT/'configs/P1_PROBABILITY_PROTOCOL.json']):
        shutil.copyfile(path,sourceout/f'{n}_{path.name}')
    write(OUT/'INPUT_PROVENANCE.json',dict(inference=str(inf/'RESULTS.json'),point_validation=str(pv/'RESULTS.json'),
        observation_file='results/P1B/data/OBSERVATIONS.npz',setup_rule='results/P1A/rules/full_9.npz',
        setup_paid_stage='P1A',setup_is_not_free=True,truth_access=False,inference_frozen=True,
        source_equality='literal bytes, public family JSON and public wall arrays; no hashes'))
    book=ControlBook(Book('P1B'));answer=[];active=None
    try:
        for i,y in enumerate(ys):
            directory=OUT/f'observation_{i:03d}';directory.mkdir();p,F,path=baseline[i]
            original=_point_record(int(p[0]),float(p[1]),float(p[2]),F,y,path)
            candidates=[];active=dict(observation=i,parent=int(parents[i]),candidate_branch=None)
            for b in (0,1):
                active['candidate_branch']=b;proposal,Fproxy=ArrayProxy(rule,b,y).optimize()
                u,c=proposal['u'],proposal['contrast'];k=key(b,u,c)
                dest=directory/f'branch_{b}';dest.mkdir();write(dest/'PROPOSAL.json',proposal)
                np.savez_compressed(dest/'PROXY.npz',F=Fproxy,q=family.walls[b],u=u,contrast=c)
                cpu,wall=time.process_time(),time.perf_counter();before=book.used();cache=k in sources
                if cache:
                    source=sources[k]
                    with np.load(source,allow_pickle=False) as z:
                        actual=z['F'].copy()
                        if 'params' in z:
                            pp=z['params'];savedkey=key(*pp[:3])
                        else:savedkey=key(int(z['branch']),float(z['u']),float(z['contrast']))
                    if savedkey!=k:raise RuntimeError('Floathex cache mismatch')
                else:
                    # Contract is rechecked directly before each newly paid call.
                    _contract(bank,family,create=False);book.begin(k)
                    physical=dest/'physical'
                    actual=family.point(b,u,c,book,str(physical.relative_to(ROOT/'results')))
                    source=physical/'FULL_RESPONSE.npz';sources[k]=source
                if actual.shape!=(5,4,12) or not np.iscomplexobj(actual) or not np.isfinite(actual).all():raise RuntimeError('Invalid Full response')
                err=float(np.linalg.norm(actual-Fproxy)/.01);actualphi=profile(actual,y)[0]
                pe=abs(actualphi-proposal['Phi']);qualified=err<=.05 and pe<=.05
                point=_point_record(b,u,c,actual,y,source,err,pe)
                point.update(proxy=proposal,exact_floathex_cache_hit=cache,new_rhs=book.used()-before,
                    cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,
                    proxy_geometry_qualified=qualified,
                    proxy_geometry_status='QUALIFIED_AT_THIS_POINT_ONLY' if qualified else 'UNQUALIFIED_PROXY_GEOMETRY',
                    fields_source=str(source.parent),full_candidate_eligible_despite_proxy_failure=True)
                np.savez_compressed(dest/'FULL_EVALUATED_POINT.npz',F=actual,F_proxy=Fproxy,q=family.walls[b],
                    params=[b,u,c,point['amplitude'],point['phase']],Phi=actualphi)
                write(dest/'RESULT.json',point);candidates.append(point)
                write(directory/'PROGRESS.json',dict(**active,completed_candidates=candidates))
            best=min([original,*candidates],key=lambda r:r['Phi']).copy()
            best['selection_source']='original_validated_direct_grid' if best['physical_source']==original['physical_source'] else 'continuous_branch_candidate'
            row=dict(observation=i,parent=int(parents[i]),best=best,candidates=candidates,
                baseline_direct_grid=original,actual_Phi_improvement=original['Phi']-best['Phi'],
                truth_access=False,inference_frozen=True,global_continuous_MAP_certified=False,
                actual_Full_KKT_certified=False,additional_control_not_replacement=True,
                cost=dict(new_rhs=sum(p['new_rhs'] for p in candidates),
                    array_cpu_seconds=sum(p['proxy']['cpu_seconds'] for p in candidates),
                    array_wall_seconds=sum(p['proxy']['wall_seconds'] for p in candidates),
                    full_or_cache_cpu_seconds=sum(p['cpu_seconds'] for p in candidates),
                    setup_provenance='Paid P1A full9 and prior Full-validated direct points; shared setup not free'))
            answer.append(row);write(directory/'RESULT.json',row);write(OUT/'RESULTS.json',answer)
        write(OUT/'COMPLETE.json',dict(status='FULL_CONTINUOUS_POINT_CONTROL_COMPLETE',truth_access=False,
            inference_frozen=True,observations=32,new_rhs=book.used(),scope=SCOPE,scope_rhs_cap=CAP,
            cpu_seconds=time.process_time()-start_cpu,wall_seconds=time.perf_counter()-start_wall,
            new_full_vectors=sum(r.get('event')=='continuous_full_vector_attempt' for r in book.local()),
            physical_budget=book.totals()['P1B'],automatic_GO=False))
    except BaseException as exc:
        write(OUT/'FAILURE.json',dict(status='INCOMPLETE_NO_AUTOMATIC_RETRY',error_type=type(exc).__name__,
            error=str(exc),active=active,completed_observations=len(answer),scope_used_rhs=book.used(),
            cpu_seconds=time.process_time()-start_cpu,wall_seconds=time.perf_counter()-start_wall,
            truth_access=False,inference_frozen=True,partial_endpoints_preserved=True))
        raise


if __name__=='__main__':main()

"""Pure norm/cache equivalence, saved-state replay, and lossless I/O tests."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,time,resource,sys,zipfile,importlib.util,tempfile
import numpy as np
from p3r4.decision_refresh_tight import DecisionRefreshTight
from p3r4.decision_refresh_fast import DecisionRefreshFast
from p3r4.interval_chain import IntervalState
from p3r4.core import save,array_io_profile
root=Path.cwd();stage=root/'research/p3l_r4_conditional_transport';out=stage/'AUDIT/fast_chain';start=time.perf_counter();cpu=time.process_time()
with np.load(stage/'inputs/background/weak.npz') as z:G=z['G'];S=z['S'];incident=z['incident']
with np.load(stage/'results/contexts/cache/R4DEV_00032.npz') as z:base={k:z[k] for k in ('chi','current','field')}
with np.load(stage/'results/II/interval_chain_tight_r1/interval/R4DEV_00032/PAID_ORIGINAL_R0.npz') as z:R0=z['R0']
sigma=json.loads((stage/'configs/TARGET.json').read_text())['noise']['sigma']
args=(base['chi'],base['current'],base['field'],G,S,incident,R0,sigma,576,True)
tight=DecisionRefreshTight(*args);fast=DecisionRefreshFast(*args)
rng=np.random.default_rng(20261025132);maxdelta=0.;maxeta=0.
for _ in range(30):
 chi=base['chi']+.001*(rng.normal(size=576)+1j*rng.normal(size=576));a=tight.bound_constants(chi);b=fast.bound_constants(chi)
 maxdelta=max(maxdelta,abs(a['perturbation_operator_F']-b['perturbation_operator_F']));maxeta=max(maxeta,abs(a['eta_upper']-b['eta_upper']));assert np.isclose(a['eta_upper'],b['eta_upper'],atol=1e-14,rtol=1e-13)
cachevalue=fast.bound_constants(chi);assert fast.bound_constants(chi) is cachevalue
try:cachevalue['eta_upper']=0
except TypeError:pass
else:raise AssertionError('cachedconstantmutated')
try:fast.G[0,0]=0
except ValueError:pass
else:raise AssertionError('cachedGmutable')
assert len(fast._bound_cache)==2
with np.load(stage/'results/II/interval_chain_tight_r1/full/R4DEV_00032/CHAIN.npz') as z:truthphi=z['candidate_phi']
with np.load(stage/'results/contexts/dev/R4DEV_00032.npz') as z:y=z['observation']
replay=[]
for i in range(4):
 with np.load(stage/f'results/II/interval_chain_tight_r1/interval/R4DEV_00032/candidates/s{i:04d}.npz') as z:point={k:z[k] for k in ('x','q','chi','current','field','prediction','residual')}
 a=tight.inspect(point['chi'],point['current'],y,tight.bound_constants(point['chi']));b=fast.inspect(point['chi'],point['current'],y,fast.bound_constants(point['chi']))
 assert np.allclose([a['phi_lower'],a['phi_upper']],[b['phi_lower'],b['phi_upper']],atol=1e-10,rtol=1e-13)
 assert b['phi_lower']<=truthphi[i]<=b['phi_upper']
 replay.append(dict(candidate=i,max_interval_difference=max(abs(a['phi_lower']-b['phi_lower']),abs(a['phi_upper']-b['phi_upper'])),contains_saved_independentFull=True))
oldenv=os.environ.get('R4_ARRAY_CODEC');io=Path(tempfile.mkdtemp(prefix='SCHEMA_ONLY_',dir=out))
try:
 os.environ.pop('R4_ARRAY_CODEC',None);save(io/'compressed.npz',**point)
 os.environ['R4_ARRAY_CODEC']='stored';save(io/'stored.npz',**point)
 with zipfile.ZipFile(io/'compressed.npz') as z:assert all(i.compress_type==zipfile.ZIP_DEFLATED for i in z.infolist())
 with zipfile.ZipFile(io/'stored.npz') as z:assert all(i.compress_type==zipfile.ZIP_STORED for i in z.infolist())
 with np.load(io/'compressed.npz') as za,np.load(io/'stored.npz') as zb:assert all(np.array_equal(za[k],zb[k]) for k in point)
 os.environ['R4_ARRAY_CODEC']='invalid'
 try:save(io/'invalid.npz',x=np.zeros(1))
 except ValueError:pass
 else:raise AssertionError('unrecognizedcodecaccepted')
 os.environ['R4_ARRAY_CODEC']='stored'
 spec=importlib.util.spec_from_file_location('fast_runner',root/'tools/run_r4_interval_chain_fast.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 class BookStub:jobdir=io/'fakeBook'
 writer=module.EndpointWriter(io/'writer',BookStub(),save,incident)
 state=IntervalState(point['x'],point['q'],point['chi'],dict(b),False,0,'');state.endpoint_id='c0000'
 writer.final(state);first=np.load(io/'writer/endpoints/c0000.npz')['current'];writer.final(state)
 assert writer.stats['array_writes']==1
 state.point=fast.inspect(point['chi'],point['current']+fast.R0@b['residual'],y,fast.bound_constants(point['chi']));state.refinements=1
 writer.proof(state,'current','refined',1);writer.final(state)
 assert writer.stats['array_overwrites']==1 and writer.versions['c0000']==2
 raw=io/'fakeBook/endpoints/paid001.npz';save(raw,q=state.q,chi=state.chi,current=state.point['current'],field=state.point['field'],prediction=state.point['prediction'],backward_residual=np.array(0.))
 save(io/'fakeBook/attempts/paid001.npz',q=state.q)
 exact=IntervalState(state.x,state.q,state.chi,state.point,True,0,'paid001');exact.endpoint_id='c0001';writer.final(exact)
 assert writer.stats['array_writes']==2 and writer.stats['Target_endpoint_references']==1 and not (io/'writer/endpoints/c0001.npz').exists()
finally:
 if oldenv is None:os.environ.pop('R4_ARRAY_CODEC',None)
 else:os.environ['R4_ARRAY_CODEC']=oldenv
result=dict(status='PASS',perturbation_norm_max_difference=maxdelta,eta_max_difference=maxeta,saved_replay=replay,cache_immutable_twoentries=True,
 stored_lossless=True,default_compressed_preserved=True,invalid_codec_rejected=True,canonical_samepoint_rewrite_versioned=True,Full_fields_referenced_not_duplicated=True,
 fast_setup=fast.fast_setup_profile,fast_counts=dict(fast.counts),array_io_profile=array_io_profile(),
 physics_RHS=0,physics_states=0,training_updates=0,paidlaunch=False,wall_seconds=time.perf_counter()-start,cpu_seconds=time.process_time()-cpu,
 peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
(out/'FAST_EQUIVALENCE_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('fast_setup','fast_counts','array_io_profile','saved_replay')}))

"""Synthetic schema-only test. No context-generation or training claim."""
import importlib.util,json,time,resource,sys
from pathlib import Path
import numpy as np
root=Path.cwd();spec=importlib.util.spec_from_file_location('packer',root/'tools/pack_r4_context_training.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
base=root/'research/p3l_r4_conditional_transport/AUDIT/risk/CONTEXT_PACK_SCHEMA_ONLY';base.mkdir(exist_ok=True);(base/'contexts/train').mkdir(parents=True,exist_ok=True);(base/'outputs').mkdir(exist_ok=True)
rng=np.random.default_rng(20261025042);parents=m.expected_parents('train',2);identities=[];observations=[]
wall,cpu=time.monotonic(),time.process_time()
for i,parent in enumerate(parents):
 obs=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8));observations.append(obs)
 q=np.zeros(66);q[64:]=.5
 row=dict(observation=obs,q0=q,chi0=np.full(576,.1+.04j),J=np.zeros((512,66)),prediction=obs*.9,residual=obs*.1,phi=np.array(1.),x=np.zeros(66),H=np.eye(66),parent=np.array(parent))
 path=base/'contexts/train'/f'{parent}.npz';np.savez(path,**row);s=path.stat();identities.append(dict(parent=parent,path=str(path),bytes=s.st_size,mtime_ns=s.st_mtime_ns))
x=np.zeros((2,66),dtype='float32');y=m.packed_y(np.array(observations));material=np.full((2,576),.1+.04j,dtype='complex64')
np.savez(base/'synthetic_joint_train.npz',x=x,y=y,material=material,parent=np.array(parents))
row=m.pack_split(base/'synthetic_joint_train.npz',base/'contexts','train',parents,identities,base/'outputs',.25)
from p3r2.residual_risk import read_anchor
from p3r4.train_density import load_data
anchor=read_anchor(base/'outputs/anchor_train.npz',supervised=True);density=load_data(base/'outputs/density_train.npz')
assert anchor['q0'].shape==(2,66) and np.array_equal(anchor['chi_GT'],material)
assert np.array_equal(density['x'],x) and np.array_equal(density['y'],y)
assert np.array_equal(density['base_mean'],np.zeros((2,66)))
assert np.array_equal(np.diagonal(density['base_chol'],axis1=1,axis2=2),np.full((2,66),.25))
with np.load(base/'outputs/online_train.npz') as z:
 assert all(key not in z for key in ('chi_GT','q_GT','x_GT','material'))
result=dict(status='PASS',scope='synthetic schema only; zero claimed generated strong contexts',physics_calls=0,training_updates=0,
            read_anchor_supervised='PASS',density_full66_loader='PASS',online_truth_keys_absent='PASS',
            base_fullrank_diagonal=.25,wall_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
            peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)*(1 if sys.platform=='darwin' else 1024))
(base/'SCHEMA_ONLY_RESULT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

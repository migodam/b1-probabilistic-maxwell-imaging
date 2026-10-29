"""Separate GT-generating process. Never import this module from a solver."""
import os,sys,json,time,secrets
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(P/'src'))
from physics import Model

def main():
 c=json.loads((P/'experiments/contract.json').read_text()); private=P/'private'; public=P/'experiments/public'
 if (private/'sealed_targets.npz').exists(): raise RuntimeError('Do not replace hidden target set')
 seed=secrets.randbits(128);rng=np.random.default_rng(seed);targets={};cases=[];t=time.perf_counter();cpu=time.process_time()
 for preset,g in c['presets'].items():
  model=Model(None,**g)
  for i in range(3):
   case=f'{preset}_{i:02d}';m=rng.uniform(.5,4,model.n);state=model.forward(m,independent=True)
   targets[case]=m;np.savez(public/f'{case}.npz',y=state['y']);cases.append({'case':case,'preset':preset,'observation_file':f'{case}.npz'})
 np.savez(private/'sealed_targets.npz',**targets)
 (private/'generator_secret.json').write_text(json.dumps({'seed':str(seed),'created_unix':time.time(),'wall_seconds':time.perf_counter()-t,'CPU_seconds':time.process_time()-cpu,'generation_source_RHS':sum(len(c['presets'][x['preset']]['sources']) for x in cases)},indent=2)+'\n')
 os.chmod(private/'generator_secret.json',0o600);os.chmod(private/'sealed_targets.npz',0o600)
 (public/'manifest.json').write_text(json.dumps({'schema':c['schema'],'cases':cases,'generator':'separate process; hidden targets and seed private'},indent=2)+'\n')
 print('Generated 6 observation-only public cases; private truth sealed, values suppressed.')
if __name__=='__main__':main()

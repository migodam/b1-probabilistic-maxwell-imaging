"""Track D closed-by-default contract helpers; no physics/model/truth imports."""
from pathlib import Path
import hashlib,json,time
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent
GROUPS={'DLOW':.5,'DHIGH':2.};SEEDS=(2026102601,2026102602,2026102603)
def identity(path):
 p=Path(path);h=hashlib.blake2b()
 with p.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return dict(path=str(p.resolve()),bytes=p.stat().st_size,blake2b=h.hexdigest())
def write(path,obj):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def freeze():
 p=ROOT/'D_FREEZE.json'
 if not p.exists():raise RuntimeError('NOT_RUN: root D_FREEZE absent')
 d=json.loads(p.read_text())
 if d.get('status')!='ROOT_FROZEN' or d.get('execution_authorized') is not True or d.get('C_finalist_qualified') is not True:raise RuntimeError('NOT_RUN: root qualification/execution freeze required')
 if d.get('groups')!=list(GROUPS) or d.get('seeds')!=list(SEEDS):raise ValueError('two fixed noise cohorts and three frozen seeds required')
 if d.get('statistics')!={'bootstrap':20000,'rng_seed':202610270699,'familywise_alpha_per_condition':.05,'comparisons_per_condition':3,'NI_fraction':.01,'wall_fraction':.8}:raise ValueError('freeze statistical contract exactly; any amendment requires versioned code')
 for group in GROUPS:
  p=OUT/f'configs/parents_{group}.json'
  if identity(p)!=d['manifest_identities'][group]:raise ValueError('manifest differs from root freeze')
 arms=d['matched_arms']
 if not {'matched_full','matched_selected'}<=set(arms) or arms['matched_full']['context']!='full':raise ValueError('matched primary pair required')
 if arms['matched_selected']['context'] not in ('native','low'):raise ValueError('selected candidate must be a reduced anchor')
 if set(arms)&{'strongCNN','distilled_w128_n32768','deterministic_matched'}:raise ValueError('matched names cannot shadow original baselines')
 for name,arm in arms.items():
  if arm['context'] not in ('full','native','low') or set(arm['checkpoints'])!=set(map(str,SEEDS)):raise ValueError('three checkpoint seeds and allowed context required')
  for value in arm['checkpoints'].values():
   if identity(value['path'])!=value:raise ValueError('checkpoint differs from root freeze')
 return d

def parents(group):
 if group not in GROUPS:raise ValueError('source_subset NOT_RUN; DLOW/DHIGH only')
 rows=json.loads((OUT/f'configs/parents_{group}.json').read_text())['parents']
 ids=[r['parent'] for r in rows]
 if len(ids)!=256 or len(set(ids))!=256:raise ValueError('exactly256 unique parents required')
 return ids

def context_path(group,parent,kind,d):
 config={'full':'native96_probit64','native':'native96_only','low':d['selected_low_config']}[kind]
 return ROOT/'contexts'/group/parent/config/'CONTEXT.npz'

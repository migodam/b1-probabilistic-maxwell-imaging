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
def _json_value(value):
 # Original optimizer reports contain NumPy arrays/scalars. This changes only
 # report serialization, never a forward, optimizer, prediction or score.
 if hasattr(value,'tolist'):return value.tolist()
 if isinstance(value,Path):return str(value)
 if isinstance(value,complex):return [value.real,value.imag]
 raise TypeError(f'Unsupported report type {type(value).__name__}')
def write(path,obj):
 p=Path(path);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False,default=_json_value)+'\n');tmp.replace(p)
def freeze():
 p=ROOT/'D_FREEZE.json'
 if not p.exists():raise RuntimeError('NOT_RUN: root D_FREEZE absent')
 d=json.loads(p.read_text())
 if d.get('status')!='ROOT_FROZEN' or d.get('execution_authorized') is not True or d.get('C_finalist_qualified') is not True:raise RuntimeError('NOT_RUN: root qualification/execution freeze required')
 if d.get('groups')!=list(GROUPS) or d.get('seeds')!=list(SEEDS):raise ValueError('two fixed noise cohorts and three frozen seeds required')
 if d.get('statistics')!={'bootstrap':20000,'rng_seed':202610270699,'familywise_alpha_global':.05,'conditions':2,'familywise_alpha_per_condition':.025,'comparisons_per_condition':4,'comparisons_global':8,'per_comparison_alpha':.05/8,'NI_fraction':.01,'wall_fraction':.8}:raise ValueError('freeze statistical contract exactly; any amendment requires versioned code')
 for group in GROUPS:
  p=OUT/f'configs/parents_{group}.json'
  if identity(p)!=d['manifest_identities'][group]:raise ValueError('manifest differs from root freeze')
 arms=d['estimator_arms']
 full,selected,distill=d['primary_full'],d['primary_selected'],d['distill_reference'];actual=d['actualgradient_reference']
 if full not in arms or selected not in arms or distill not in arms or actual not in arms:raise ValueError('declared reference arms missing')
 if arms[full]!={'original':'deterministic_matched','context':'full'} or arms[selected]!={'original':'deterministic_matched','context':'low','feature_mode':'kkt_q'}:raise ValueError('frozen original deterministic Full/Low finalist required')
 if arms[distill]!={'original':'distilled_w128_n32768','context':'full'}:raise ValueError('full historical Distill reference required')
 if not any(v=={'original':'strongCNN','context':'full'} for v in arms.values()):raise ValueError('full Strong reference required')
 if arms[actual]!={'original':'deterministic_matched','context':'low'}:raise ValueError('ordinary low actual-gradient control required')
 if arms.get('full_kkt_deterministic')!={'original':'deterministic_matched','context':'full','feature_mode':'kkt_q'}:raise ValueError('full q-closure control required')
 if d['selected_low_config']!='native24_probit64':raise ValueError('approved Revision1 low anchor required')
 for arm in arms.values():
  if arm.get('context') not in ('full','low') or arm.get('original') not in ('strongCNN','distilled_w128_n32768','deterministic_matched') or arm.get('feature_mode','actual') not in ('actual','kkt_q'):raise ValueError('original frozen estimators only; Native/matched training candidate excluded')
 if not d.get('frozen_checkpoint_identities'):raise ValueError('root must freeze actual original heads and Gaussian/provider dependencies')
 for value in d['frozen_checkpoint_identities']:
  if identity(value['path'])!=value:raise ValueError('original checkpoint differs from root freeze')
 for value in d.get('frozen_source_identities',[]):
  if identity(value['path'])!=value:raise ValueError('frozen source changed; record real bug before amendment')
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

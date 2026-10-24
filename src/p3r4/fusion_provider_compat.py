"""Additive recovery of six closed provider headers; no tensor/feature changes.

Known trainer count-shadowing wrote parameter_count6418 into moment_K. Recovery
requires a checkpoint-identity seal and independently closed TRAIN/validation
provider gates. Ordinary and legacy checkpoints keep the original strict loader.
"""
from pathlib import Path
import json,torch
from .risk_fusion import FusionPredictor,FusionRisk,matching_moment_gates
from .evaluation import checkpoint_identity
AUDIT=Path(__file__).resolve().parents[2]/'research/p3l_r4_conditional_transport/AUDIT/fusion/providers/CLOSED_K_METADATA_CLARIFICATION.json'

def load_fusion_compatible(path,device='cpu'):
 p=torch.load(path,map_location='cpu',weights_only=False)
 if p.get('moment_K')!=6418:return FusionPredictor(path,device)
 audit=json.loads(AUDIT.read_text());identity=checkpoint_identity(path)
 rows=[r for r in audit['records'] if r['checkpoint_identity']==identity]
 if len(rows)!=1:raise ValueError('unknown malformed provider checkpoint; no inference recovery')
 r=rows[0];gates=[]
 for role in ('train','validation'):
  seal=r[role+'_gate_identity']
  if checkpoint_identity(seal['path'])!=seal:raise ValueError('closed source gate changed')
  gates.append(json.loads(Path(seal['path']).read_text()))
 provider,count,source=matching_moment_gates(*gates)
 if source!=p['moment_checkpoint_identity'] or provider!=p['moment_provider'] or count!=r['true_K']:raise ValueError('provider source recovery mismatch')
 if sum(t.numel() for k,t in p['state_dict'].items() if not k=='roi')!=6418:raise ValueError('known architecture tensor count required')
 model=FusionPredictor.__new__(FusionPredictor);model.p=p.copy();model.p['moment_K']=count
 model.geometry=p['geometry'];model.device=device;model.moment_provider=provider;model.moment_K=count
 model.model=FusionRisk(model.geometry).to(device);model.model.load_state_dict(p['state_dict']);model.model.eval()
 model.metadata_recovery=r
 return model

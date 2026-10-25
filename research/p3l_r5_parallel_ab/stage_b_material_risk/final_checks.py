import json,sys,datetime
from pathlib import Path
OUT=Path(__file__).resolve().parent;sys.path.insert(0,str(OUT))
import numpy as np,torch
from inference import Action,OLD
from p3r2.residual_risk import read_anchor
from p3r4.risk import full_material
from p3r2.fusion import passive_box
torch.set_num_threads(1);data=dict(np.load(OUT/'prepared_validation.npz'));online=read_anchor(OLD/'research/p3l_r4_conditional_transport/inputs/strong_context_1024_256/anchor_validation.npz');g=torch.load(OUT/'geometry.pt',map_location='cpu',weights_only=False)['geometry'];rows=[]
for method in ['deterministic_matched','distilled_w128_n32768']:
 for seed in [2026102601,2026102602,2026102603]:
  first=Action(method,seed).predict(data['y'],online,batch_size=64);second=Action(method,seed).predict(data['y'],online,batch_size=64);rows.append(dict(method=method,seed=seed,exact_array_equal=bool(np.array_equal(first,second)),max_abs_difference=float(np.max(abs(first-second)))))
obj=dict(status='PASS' if all(r['exact_array_equal'] for r in rows) else 'FAIL',scope='CPU1 twofreshActionloads all256VALparents; distinctfromcrossdevice3e-6gate',records=rows,physics_RHS=0,original_CPU_MPS_gate='FAIL_RETAINED_NO_TOLERANCE_CHANGE');(OUT/'CPU_EXACT_REPRODUCIBILITY.json').write_text(json.dumps(obj,indent=2)+'\n')
preds={}
for key,name in [('strong','strong_CNN'),('K32','original_K32_frozen')]:preds[name]=np.array([full_material(data[f'{key}_s{s}'],g) for s in [2026102501,2026102502,2026102503]])
for name,pattern in [('deterministic_matched','experiments/deterministic_matched_s{}/validation.npz'),('distilled_w128_n32768','learned_context/distilled_w128_n32768_s{}/validation.npz')]:preds[name]=np.array([np.load(OUT/pattern.format(s))['prediction'] for s in [2026102601,2026102602,2026102603]])
rng=np.random.default_rng(2026102643);ids=rng.integers(256,size=(20000,256));records=[]
for clipping in [False,True]:
 risk={name:np.mean(abs(np.array([passive_box(v,g['config'],g['roi']) for v in pp]) if clipping else pp-data['truth'][None])**2,axis=2).mean(0) for name,pp in []}
 risk={}
 for name,pp in preds.items():
  if clipping:pp=np.array([passive_box(v,g['config'],g['roi']) for v in pp])
  risk[name]=np.mean(abs(pp-data['truth'][None])**2,axis=2).mean(0)
 B=risk['strong_CNN'];O=risk['original_K32_frozen'];den=(B-O)[ids].mean(1)
 for name in ['deterministic_matched','distilled_w128_n32768']:
  C=risk[name];gain=(O-C)[ids].mean(1);ret=(B-C)[ids].mean(1)/den;positive=den>0
  records.append(dict(action='passive_clipped' if clipping else 'raw',candidate=name,mean_reference_minus_candidate=float(np.mean(O-C)),paired_gain_CI_98_333333=np.quantile(gain,[1/120,119/120]).tolist(),retained_increment_point=float(np.mean(B-C)/np.mean(B-O)),retained_increment_CI_98_333333_condition_positive_denominator=np.quantile(ret[positive],[1/120,119/120]).tolist(),nonpositive_denominator_bootstrap_count=int((~positive).sum())))
boot=dict(status='EXPLORATORY_SELECTION_DATA_ONLY',seed=2026102643,replicates=20000,parent_count=256,aggregation='parentbootstrap on3seedaverageloss; seedsnotresampled; ensembleimagelossnotused',CI_level=.9833333333333333,records=records,not_independent_confirmation=True);(OUT/'BOOTSTRAP_EXPLORATORY.json').write_text(json.dumps(boot,indent=2)+'\n');print(obj['status']);print(json.dumps(records,indent=2))

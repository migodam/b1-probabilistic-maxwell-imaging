import sys,json,time,csv
from pathlib import Path
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2];OLD=ROOT.parent/'b1-p3l-r4-conditional-transport';sys.path.insert(0,str(OLD/'src'))
import numpy as np,torch
from p3r4.risk import full_material
from p3r2.fusion import passive_box
from p3r4.risk_fusion import FusionRisk
S=OLD/'research/p3l_r4_conditional_transport';torch.set_num_threads(1);d=dict(np.load(OUT/'prepared_validation.npz'));g=torch.load(OUT/'geometry.pt',weights_only=False,map_location='cpu')['geometry'];inventory=json.loads((OUT/'ASSET_INVENTORY.json').read_text());A=np.asarray(g['A']);cfg=g['config'];truth=d['truth'];pred={};newseeds=[2026102601,2026102602,2026102603];oldseeds=[2026102501,2026102502,2026102503];replay=[]
for family,key in [('strong_CNN','strong'),('original_K32_frozen','K32'),('deterministic_frozen','det_frozen')]:pred[family]=np.array([full_material(d[f'{key}_s{s}'],g) for s in oldseeds])
for arm in ['deterministic_matched','direct_moments_matched']:
 preds=[]
 for seed in newseeds:
  path=OUT/'experiments'/f'{arm}_s{seed}';p=dict(np.load(path/'validation.npz'));preds.append(p['prediction'])
  if arm=='deterministic_matched':
   payload=torch.load(path/'BEST.pt',weights_only=False,map_location='cpu');model=FusionRisk(g);model.load_state_dict(payload['state_dict']);model.eval();mom=(d['det']-payload['moment_mean'])/payload['moment_scale'];rows=[]
   with torch.inference_mode():
    for b in range(0,len(mom),64):rows.append(model(torch.tensor(d['image'][b:b+64]),torch.tensor(mom[b:b+64]),torch.tensor(d['anchor'][b:b+64])).numpy())
   cp=full_material(np.concatenate(rows),g);replay.append(dict(seed=seed,CPU_MPS_max_abs_prediction=float(np.max(abs(cp-p['prediction']))),CPU_MPS_mean_risk_difference=float(np.mean(abs(cp-truth)**2)-np.mean(p['risk']))))
 pred[arm]=np.asarray(preds)
for arm in ['raw_y','decoder_only','physical_residual','same_spectrum_random']:pred['RRR_'+arm]=np.asarray([np.load(OUT/'low_rank_correction'/f'{arm}_s{s}'/'validation.npz')['prediction'] for s in newseeds])
risks={k:np.mean(abs(v-truth[None])**2,2).mean(0) for k,v in pred.items()};summary=[]
for k,v in pred.items():
 cr=np.array([passive_box(x,cfg,g['roi']) for x in v]);summary.append(dict(family=k,seeds=3,seed_mean_parent_MSE=float(risks[k].mean()),passive_clipped_seed_mean_MSE=float(np.mean(abs(cr-truth[None])**2)),ensemble_image_MSE=float(np.mean(abs(v.mean(0)-truth)**2)),seed_MSE=[float(np.mean(abs(x-truth)**2)) for x in v]))
rng=np.random.default_rng(2026102641);idx=rng.integers(0,256,size=(20000,256));comparisons=[]
for candidate,reference in [('deterministic_frozen','original_K32_frozen'),('deterministic_matched','original_K32_frozen'),('direct_moments_matched','original_K32_frozen')]+[(k,'strong_CNN') for k in pred if k.startswith('RRR_')]+[('RRR_physical_residual',k) for k in ['RRR_raw_y','RRR_decoder_only','RRR_same_spectrum_random']]:
 difference=risks[reference]-risks[candidate];boot=difference[idx].mean(1);ci=np.quantile(boot,[.0083333333333333,.9916666666666667]);comparisons.append(dict(candidate=candidate,reference=reference,mean_reference_minus_candidate=float(difference.mean()),relative_gain=float(difference.mean()/risks[reference].mean()),CI_98_333333=ci.tolist(),scope='exploratory VAL selection data; parent bootstrap averages losses acrossseed; not locked confirmatory'))
field=(cfg['field_mean']+d['q0'][:,:64]@A.T)/cfg['phase_transition_halfwidth'];mask=np.zeros(576,bool);mask[g['roi']]=True;regions=dict(pure_host=(field<=-1)&mask,transition=(abs(field)<1)&mask,pure_inclusion=(field>=1)&mask);phase=[]
for region,sel in regions.items():
 rec=dict(region=region,classification='frozen original y-only anchor field; not true posterior phase',ROI_cell_fraction=float(sel.sum()/(256*len(g['roi']))),full_grid_cell_fraction=float(sel.mean()),cell_count=int(sel.sum()),risks={})
 for family,v in pred.items():
  err=np.mean(abs(v-truth[None])**2,0);rec['risks'][family]=float(err[sel].mean())
 rec['incremental_gain_per_cell']={k:rec['risks']['strong_CNN']-v for k,v in rec['risks'].items()};phase.append(rec)
base=float(risks['strong_CNN'].mean());original=float(risks['original_K32_frozen'].mean());gates=[];bench=json.loads((OUT/'COST_BENCHMARK.json').read_text());cost={r['provider']+'_'+str(r['batch']):r['warm_per_parent_seconds'] for r in bench['rows']}
for family,provider in [('deterministic_matched','gaussian_deterministic'),('direct_moments_matched','no_density_direct')]:
 keep=(base-float(risks[family].mean()))/(base-original);ratios={str(b):cost[provider+'_'+str(b)]/cost['gaussian_K32_'+str(b)] for b in [1,64]};gates.append(dict(family=family,VAL_retained_incremental_gain_fraction=keep,risk_retention_VAL_gate=bool(keep>=.8),provider_cost_ratios=ratios,provider_cost_single_gate=bool(ratios['1']<=1/3),provider_cost_batch64_gate=bool(ratios['64']<=1/3),final_compression_gate='NOT_ESTABLISHED: no independent confirmation; fullpipeline cost reported separately'))
res=dict(status='DEVELOPMENT_COMPLETE_INDEPENDENT_CONFIRMATION_RESOURCE_LIMITED',data_regime='head/RRR1024TRAIN256VAL; restoredproviders32768TRAIN2048VAL; samejointproviderdata',representation_arms=8,shared_reference='3frozen strongCNN seeds notcountednewarm',summary=summary,comparisons=comparisons,phase=phase,compression_gates=gates,CPU_MPS_replay=replay,RRR_caveats=['OOF correction conditional on unsupervised fullTRAIN PCA and fixed baseline trained on fullTRAIN; not end-to-end crossfit','rawy decoder physical RRR deterministic conditional on baseline; seed variations pair3baseline seeds; onlyrandomfeatures/folds random','RRR failure limited to frozen128diminputPCA ridge rankgrid; not proof no nonlinear compressiblecontext'],new_confirmation=dict(status='RESOURCE_LIMITED_NOT_RUN',required_parents=512,RHS_cap=32768,expected_R4_full_context_RHS_512=1934840,selection_used_no_exposed_R4_TEST=True,new_truth_generated=False),physics_RHS=0)
(OUT/'RESULTS.json').write_text(json.dumps(res,indent=2)+'\n')
with (OUT/'METRICS.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=['family','seed_mean_parent_MSE','passive_clipped_seed_mean_MSE','ensemble_image_MSE']);wr.writeheader();wr.writerows([{k:r[k] for k in wr.fieldnames} for r in summary])
print(json.dumps(summary,indent=2));print(json.dumps(gates,indent=2))

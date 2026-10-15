from pathlib import Path
import json,numpy as np
from b1_joint.accounting import dump
from research_a.budget import Budget
from research_a.target import OUT
from research_a.search import structure_count
with Budget(OUT).job('confirmation-analysis'):
 rows=json.loads((OUT/'results/CONFIRM_METRICS.json').read_text());campaign=json.loads((OUT/'results/CONFIRM_CAMPAIGN.json').read_text())
 dump(OUT/'results/METRICS.json',rows)
 arms=['finite_static','multistart','population','mdn','adaptive'];summary={}
 for a in arms:
  rr=[r for r in rows if r['arm']==a];cc=[m for o in campaign for m in o['methods'] if m['arm']==a]
  summary[a]=dict(Best4=float(np.mean([r['by_k']['4']['best']['normalized_mae'] for r in rr])),Best4_complex=float(np.mean([r['by_k']['4']['best']['complex_mae'] for r in rr])),Best4_Re=float(np.mean([r['by_k']['4']['best']['real_mae'] for r in rr])),Best4_Im=float(np.mean([r['by_k']['4']['best']['imag_mae'] for r in rr])),Selected=float(np.mean([r['selected']['normalized_mae'] for r in rr if r['selected']])),base=float(np.mean([r['base']['normalized_mae'] for r in rr])),Coverage4=float(np.mean([r['by_k']['4'].get('coverage',False) for r in rr])),distinct=float(np.mean([r['distinct_archive'] for r in rr])),structural_distinct=float(np.mean([structure_count(m['archive']) for m in cc])),wall_mean=float(np.mean([r['wall'] for r in rr])),wall_median=float(np.median([r['wall'] for r in rr])),proxy_calls_mean=float(np.mean([r['proxy_evaluations'] for r in rr])),new_full_calls=384,rejected=sum(not c['valid'] for m in cc for c in m['attempts']),duplicate=sum(c['discovery_status']=='DUPLICATE_MATERIAL' for m in cc for c in m['attempts']))
 rng=np.random.default_rng(261015104);boot=rng.integers(0,32,(2000,32));comparisons={}
 for a,b in [('finite_static','multistart'),('adaptive','finite_static'),('mdn','multistart')]:
  x=np.array([r['by_k']['4']['best']['normalized_mae'] for r in rows if r['arm']==a]);y=np.array([r['by_k']['4']['best']['normalized_mae'] for r in rows if r['arm']==b]);d=x-y
  comparisons[a+'-'+b]=dict(mean=float(d.mean()),paired_parent95=np.quantile(d[boot].mean(1),[.025,.975]).tolist(),relative_improvement=float(1-x.mean()/y.mean()))
 evidence=dict(nparents=32,noise_per_parent=1,methods=summary,comparisons=comparisons,threshold=json.loads((OUT/'configs/CONFIRM_FREEZE.json').read_text())['coverage_normalized_mae'],scoring='offline only; exact radial-volume field MAE; no truth selection online',fresh_after_freeze=True)
 dump(OUT/'results/SUMMARY_E1.json',evidence)
 improvement=comparisons['finite_static-multistart']['mean']<0
 entry=dict(E2_allowed=bool(improvement and summary['finite_static']['structural_distinct']>1),basis='Frozen E2-entry consistency +full-validstructuralfamilies +Best4improvement signal; not fullalgorithm/TAP acceptance',evidence=evidence,full_superiority_claim=False)
 dump(OUT/'results/E1_ENTRY_DECISION.json',entry)
 print(json.dumps(evidence,indent=2))

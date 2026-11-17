from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv
book=Ledger('accounting_and_figures');plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
scores=json.loads((ROOT/'results/scores.json').read_text());qual=json.loads((ROOT/'results/qualification.json').read_text());totals=Counter();cert_counts=Counter();rows=[]
for case in scores:
 receipt=json.loads((ROOT/'results'/case/'receipt.json').read_text());totals.update(receipt['counts'])
 for arm in ('clean','noisy'):
  for method in ('Q','TRF64','reciprocal64'):
   d=json.loads((ROOT/'results'/case/arm/(method+'.json')).read_text())
   for cr in d['certificates']:
    cert=cr.get('certificate',{});cert_counts.update(cert.get('counter',{}))
   rows.append(dict(case=case,arm=arm,method=method,discovery_CPU_seconds=d['discovery_CPU_total'],discovery_wall_seconds=d['wall_seconds'],root_count=len(d['accepted_roots'])))
# Every top-level numerical receipt is represented, but inclusive CPU is not added twice.
ex=json.loads((ROOT/'results/executor_receipt.json').read_text());stage_cpu={'six_case_processes_including_children':sum(r['CPU_seconds_including_reaped_children'] for r in ex['receipts'])}
for name,key in [('debug','historical_debug'),('formula_tests','formula_tests'),('generation_receipt','generation'),('qualification_before_serialization_audit_receipt','qualification_v1'),('qualification_receipt','qualification_final'),('scoring_receipt','scoring')]:
 d=json.loads((ROOT/'results'/(name+'.json')).read_text());r=d.get('receipt',d);stage_cpu[key]=r['CPU_seconds'];totals.update(r.get('counts',{}))
# Worker cumulative costs include saved failed engineering tests in tests history where available.
stage_cpu['worker_certification_tests_conservative_allowance']=10.
stage_cpu['qualification_failed_attempt_upper_bound']=2.
stage_cpu['scoring_failed_attempt_upper_bound']=2.
stage_cpu['misc_numerical_process_startup_reserve']=10.
account=dict(stage_CPU_seconds=stage_cpu,charged_CPU_seconds_before_plots=sum(stage_cpu.values()),CPU_budget_seconds=7200,peak_case_rss_bytes=max(r['peak_rss_bytes'] for r in ex['receipts']),physical_and_reduced_operations=dict(totals),certification_operations_separate=dict(cert_counts),note='Case process CPU includes discovery, validation, certificates, design, correction and broker child startup. Inclusive nested event times MUST NOT be summed. Worker/missing startup timing charged conservative explicit reserve, not falsely exact measurement. Audit/literature/LLM service CPU not measured locally as scientific experiment CPU.',unmetered_components='Tiny metadata commands and LLM remote inference; numerical failed attempts covered by conservative allowances.')
save(ROOT/'results/cost_accounting.json',account)
with (ROOT/'results/cost_ledger.csv').open('x') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
# 1 historical and new branches, GT distinct from reconstructions
fig,ax=plt.subplots(1,2,figsize=(11,4))
for col,arm in enumerate(('clean','noisy')):
 a=ax[col];case='T0';gt=scores[case]['truth'];q=json.loads((ROOT/'results'/case/arm/'Q.json').read_text())
 a.plot(range(1,7),gt,'k*-',lw=2,label='Synthetic GT (scoring only)')
 for i,m in enumerate(q['accepted_roots']):a.plot(range(1,7),m,'o--',label=f'Q20 old-data root {i}')
 a.set(title='T0 '+arm,xlabel='Independent material cell',ylabel='Material contrast');a.legend(fontsize=8)
fig.tight_layout();fig.savefig(ROOT/'figures/01_material_branches.png');plt.close(fig)
#2 source separation certified local roots not noisy feasible tubes
fig,ax=plt.subplots(figsize=(11,4));d=json.loads((ROOT/'results/T0/clean/design_freeze.json').read_text());gap=d['min_pair_gaps'];lower=qual['T0']['clean']['pairs'][0]['source_gap_lower_bounds'];xx=np.arange(len(gap));ax.plot(xx,gap,'o',label='Center separation');ax.plot(xx,lower,'s',label='Certified local-box lower bound');ax.axhline(2*(EPS+MODEL_ALLOWANCE),color='red',ls='--',label='2(noise + assumed model allowance)');ax.set_xticks(xx);ax.set_xticklabels([v['label'] for v in d['menu']],rotation=65,ha='right',fontsize=7);ax.set(ylabel='Complex response separation',title='T0 clean: source menu; old-source combinations cannot separate');ax.legend();fig.tight_layout();fig.savefig(ROOT/'figures/02_source_separability.png');plt.close(fig)
#3 noisy estimates versus strong baselines
fig,ax=plt.subplots(figsize=(10,4));methods=['fixed','maxmin','fisher','reciprocal64_active'];x=np.arange(6)
for j,m in enumerate(methods):ax.plot(x,[scores[f'T{i}']['arms']['noisy']['strategies'][m]['corrected_RMS'] for i in x],'o-',label=m)
ax.plot(x,[scores[f'T{i}']['arms']['noisy']['direct']['joint_RMS'] for i in x],'s-',lw=2,label='basis completion + joint TRF64')
ax.axhline(.1,color='k',ls='--',label='Frozen RMS success threshold');ax.set_yscale('log');ax.set_xticks(x,[f'T{i}' for i in x]);ax.set(ylabel='Material RMS error (not data residual)',title='Fixed noise 1e-4 + radius mismatch 0.1%');ax.legend(ncols=2,fontsize=8);fig.tight_layout();fig.savefig(ROOT/'figures/03_noisy_material_errors.png');plt.close(fig)
#4 discovery cost vs measurement tradeoff
fig,ax=plt.subplots(1,2,figsize=(10,4));names=['Q','TRF64','reciprocal64'];vals=[sum(r['discovery_CPU_seconds'] for r in rows if r['method']==m) for m in names];ax[0].bar(names,vals);ax[0].set(ylabel='CPU seconds, all six targets / both arms',title='Root discovery only; certificates extra')
ax[1].bar(x-.17,np.ones(6),.34,label='One-source methods');ax[1].bar(x+.17,[scores[f'T{i}']['arms']['noisy']['direct']['new_measurements'] for i in x],.34,label='Direct completion');ax[1].set_xticks(x,[f'T{i}' for i in x]);ax[1].set(ylabel='Added unit-energy / unit-duration measurements',title='Fewer measurements do not imply same accuracy');ax[1].legend(fontsize=8);fig.tight_layout();fig.savefig(ROOT/'figures/04_cost_tradeoffs.png');plt.close(fig)
#5 hard filtering and mismatch failure
fig,ax=plt.subplots(1,2,figsize=(11,4));hard=[len(scores[f'T{i}']['arms']['noisy']['strategies']['maxmin']['hard_retained']) for i in x];tube=[len(scores[f'T{i}']['arms']['noisy']['strategies']['maxmin']['tube_retained']) for i in x];ax[0].bar(x-.17,hard,.34,label='Point hard filter');ax[0].bar(x+.17,tube,.34,label='Local numerical tube fit');ax[0].set_xticks(x,[f'T{i}' for i in x]);ax[0].set(ylabel='Retained candidates',title='Retention is not global uniqueness');ax[0].legend(fontsize=8)
mm=[max(scores[f'T{i}']['arms']['noisy']['old_model_error']) for i in x];ax[1].bar(x,mm);ax[1].axhline(MODEL_ALLOWANCE,color='r',ls='--',label='Predeclared allowance');ax[1].set_xticks(x,[f'T{i}' for i in x]);ax[1].set(ylabel='Actual old-data discrepancy at GT',title='Declared model allowance fails in 4 / 6 targets');ax[1].legend(fontsize=8);fig.tight_layout();fig.savefig(ROOT/'figures/05_failures.png');plt.close(fig)
save(ROOT/'results/plot_receipt.json',book.receipt());print(json.dumps(account,indent=2))

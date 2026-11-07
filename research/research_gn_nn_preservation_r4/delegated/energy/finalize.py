"""Complete postprocessing of the finished E0 run; no physical solves."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[k]='1'
import csv,json,sys,time
from pathlib import Path
from collections import defaultdict
import numpy as np
R=Path(__file__).resolve().parents[4];sys.path.insert(0,str(R))
from research.research_gn_nn_preservation_r4.experiments.e0_energy import residual_relation,writecsv

HERE=Path(__file__).resolve().parents[2];REPO=HERE.parents[1]
old=REPO/'research/research_gn_nn_completion_r3';out=HERE/'results/e0'
start=time.perf_counter()
with np.load(old/'data/operators.npz')as z:S=np.array(z['S'])
relation=residual_relation(old,out,S)
keys=('run','arm','seed','passes','noise','level','stage','quantity')
groups=defaultdict(lambda:defaultdict(list))
with (out/'per_case_energy.csv').open()as f:
    for r in csv.DictReader(f):
        for band in ('S_strong9','S_weak39'):
            groups[tuple(r[k]for k in keys)+(band,)][r['parent_id']].append(float(r[band+'__delta']))
rng=np.random.default_rng(2026110705);extra=[]
for key,mothers in sorted(groups.items()):
    v=np.array([np.mean(x)for _,x in sorted(mothers.items())]);n=len(v)
    means=v[rng.integers(0,n,(1000,n))].mean(1)
    extra.append(dict(zip(keys+('projection',),key),mothers=n,delta_mean=float(v.mean()),delta_ci_low=float(np.quantile(means,.025)),delta_ci_high=float(np.quantile(means,.975)),mother_improvement_fraction=float(np.mean(v<0)),bootstrap_repetitions=1000,bootstrap_unit='mother parent, paired after-before absolute energy; seeds not pooled'))
with (out/'mother_paired_bootstrap.csv').open()as f:boot=list(csv.DictReader(f))
boot=[r for r in boot if r['projection']not in ('S_strong9','S_weak39')]+extra
writecsv(out/'mother_paired_bootstrap.csv',boot)
selected=[r for r in boot if str(r['seed'])=='7' and str(r['passes'])=='1' and r['level']=='low' and r['quantity']=='Maxwell_current' and r['projection'] in ('P_in_strong19','S_strong9','S_visible48','D509','total','G_field','S_field') and r['stage']in ('anchor_to_raw_unclipped','anchor_to_final_top1')]
writecsv(out/'DECISIVE_TRUE_CURRENT_SEED7.csv',selected)
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
arms=('v_image','v_physics','j_full','j_null','j_null_current')
for ni,noise in enumerate((0.,.001)):
    for pi,projection in enumerate(('P_in_strong19','S_strong9')):
        ax=axs[ni,pi]
        for stage,offset in [('anchor_to_raw_unclipped',-.13),('anchor_to_final_top1',.13)]:
            rs=[next(r for r in selected if r['arm']==arm and float(r['noise'])==noise and r['projection']==projection and r['stage']==stage)for arm in arms]
            y=np.array([float(r['delta_mean'])for r in rs]);lo=np.array([float(r['delta_ci_low'])for r in rs]);hi=np.array([float(r['delta_ci_high'])for r in rs])
            ax.errorbar(np.arange(5)+offset,y,yerr=[y-lo,hi-y],fmt='o',label='raw unclipped'if'raw_unclipped'in stage else'final top1',capsize=3)
        ax.axhline(0,color='k',linewidth=.7);ax.set_xticks(range(5),arms,rotation=25)
        ax.set_yscale('symlog',linthresh=1e-7 if noise==0 else 1e-6)
        ax.set_title(('clean'if noise==0 else'noisy')+' / '+projection)
        ax.set_ylabel('absolute true-current energy Δ vs GN');ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Actual Maxwell-current strong errors: seed 7, 64 low-amplitude mothers\n95% mother-paired bootstrap intervals; symlog scale')
fig.savefig(HERE/'figures/e0/true_current_strong_paired_delta.png',dpi=160);plt.close(fig)
(HERE/'figures/e0/signed_cross_term.png').unlink(missing_ok=True)
with (out/'grouped_absolute_energy.csv').open()as f:means=list(csv.DictReader(f))
with (out/'real_J_diagnostic.csv').open()as f:j=list(csv.DictReader(f))
jgroup=defaultdict(dict)
for r in j:jgroup[r['run'],r['obs_id']][r['stage']]=r
jagg=defaultdict(list)
for (run,obs),states in jgroup.items():
    a=states['anchor']
    for stage in ('raw_unclipped','final_top1'):
        if stage in states:
            b=states[stage]
            jagg[run,a['noise'],stage].append(float(b['J_strong_error_energy'])-float(a['J_strong_error_energy']))
writecsv(out/'real_J_diagnostic_delta_summary.csv',[dict(run=k[0],noise=k[1],stage=k[2],observations=len(v),mean_real_J_strong_material_error_delta=float(np.mean(v)),improvement_fraction=float(np.mean(np.array(v)<0)),scope='deterministic diagnostic subset only; not all mothers')for k,v in sorted(jagg.items())])
receipt=json.loads((out/'RECEIPT.json').read_text());receipt['data_residual_relation']=relation
receipt['postprocessing_seconds']=time.perf_counter()-start
receipt['total_recorded_wall_seconds']=receipt['wall_seconds']+receipt['postprocessing_seconds']
(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2))
def get(arm,noise,projection,stage):
    return next(r for r in selected if r['arm']==arm and float(r['noise'])==noise and r['projection']==projection and r['stage']==stage)
lines=['# E0 absolute energy audit','',f"**COMPLETE:** {receipt['counts']['frozen_model_observations']}/3264 model-observations: all 64 low-amplitude mothers (clean/noisy), all 32 existing A2.5 mothers (clean/noisy), all 17 confirmation models (three seeds plus two cascades). No training, no R4 test truth access.",'',
    '- Main rows below are **actual nonlinear Maxwell-current errors**, obtained by solving the unchanged frozen Model at the historical raw material and final top1. Auxiliary-head discrepancies and material-coordinate Green energies remain separate quantities.',
    '- For clean low-amplitude seed-7 results, final internal-strong19 and receiver-strong9 error increases in all five arms; D509 decreases. For noisy results, final internal-strong19 decreases in all five arms, while receiver-strong9 still increases. Visible48 and D509 decrease. Thus the absolute-current results do not support a blanket strong-error improvement statement.',
    '- The raw→final pipeline is recorded directly, including update energy and signed cross term. It includes finite-candidate correction and selection and does not isolate a TRF effect. All 1278 unclipped/clipped differences are float32 lower-bound roundoff only (max Linf 4.47035e-10, max RMS 3.56345e-10; zero above 1e-7), so no material clipping confound is evidenced. Every band records ΔE = update energy + signed cross term; a large update alone does not establish improvement.',
    '- Receiver true-error and observation-misfit have different references. data_residual_relation.csv records E_true = E_misfit + E_noise + signed cross, using frozen physical residual receipts and stored historical GT current, with zero additional LU. Clean GN can fit far below the .001 NN correction gate; stopping accuracy can partly account for increased true receiver-field error. This audit does not establish causality.',
    '- Complex Green projections of real material error are labelled material_complex_Green_embedding; they are not material observability or real-Jacobian null spaces. The deterministic real-J diagnostic is separate (32 observations, 16 mothers), never promoted to the full test population.','',
    '## Actual-current final minus GN absolute squared error (seed 7, one pass)','',
    '|noise|arm|internal strong19|receiver strong9|visible48|D509|',
    '|---:|---|---:|---:|---:|---:|']
for noise in (0.,.001):
    for arm in ('v_image','v_physics','j_full','j_null','j_null_current'):
        vals=[float(get(arm,noise,p,'anchor_to_final_top1')['delta_mean'])for p in ('P_in_strong19','S_strong9','S_visible48','D509')]
        lines.append('|'+f'{noise}|{arm}|'+('|'.join(f'{v:.6g}'for v in vals))+'|')
lines+=['','Negative ΔE is improvement. All entries are absolute squared norms, never shares. DECISIVE_TRUE_CURRENT_SEED7.csv includes raw/final and mother-paired 95% bootstrap intervals. The full bootstrap file keeps seeds, noise and low/stress distinct; grouped_absolute_energy.csv also separates amplitude.','',
    'Fixed geometry: P_in=19/557, S=9/39/528, visible48; B67 is the union of internal strong19 and receiver visible48, D509 its orthogonal complement. Sequential I→O and O→I, true numerical intersections (weak/weak rank20 and weak/null rank509), and the numerical-intersection remainder are separate columns.','',
    f"Cost: {receipt['wall_seconds']:.1f} s main audit; {receipt['physical_costcounts']['LU']} LU, {receipt['physical_costcounts']['actual_complex_RHS']} complex RHS; one CPU worker with BLAS=1. Full-VIE max equation defect {receipt['max_full_VIE_equation_defect']:.3g}. All {receipt['counts']['final_top1_available']} frozen final top1 outputs were available. No missing E0 absolute-energy cases. Exact unclipped-raw observed-misfit relations are formally NOT_RUN for {relation['NOT_RUN']} float32 lower-bound-roundoff cases; the actual verified raw_1 residual is retained in a separately labelled column. No material clipping effect or missing absolute-current energy follows from these formal mismatches.",'',
    'Files: per_case_energy.csv, grouped_absolute_energy.csv, mother_paired_bootstrap.csv, DECISIVE_TRUE_CURRENT_SEED7.csv, data_residual_relation.csv, real_J_diagnostic.csv, real_J_diagnostic_delta_summary.csv, GEOMETRY.json, RECEIPT.json. Figures: absolute_raw_final.png, noisy_D509_paired_delta.png, true_current_strong_paired_delta.png. Final scientific interpretation remains with the parent review.']
text='\n'.join(lines)+'\n'
(out/'E0_FINDINGS.md').write_text(text);(HERE/'delegated/energy/E0_FINDINGS.md').write_text(text)
print(json.dumps(dict(residual_relation=relation,bootstrap_rows=len(boot),main_receipt=receipt),indent=2))

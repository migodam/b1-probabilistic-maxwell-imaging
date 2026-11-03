"""Reproduce offline descriptive R5 figures from saved outputs; no physics imports."""
from pathlib import Path
from datetime import datetime, timezone
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]
OUT=ROOT/'FIGURES'
OUT.mkdir(exist_ok=True)
RESULTS=ROOT/'RESULTS'
DEV='obs_b20d5a55677d4eca'
CASES=[DEV,'obs_42cb233e4f224bb3','obs_63abafa750fc483c','obs_100bd65498644322']
NAMES=['Old curved bar','R4 failure 42cb','R4 failure 63ab','Ring 100b']
COLORS={'global':'#3269ad','onehop':'#d55e00'}
INPUTS=set()
OUTPUTS=[]
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180,'pdf.fonttype':42})

def read(p):
    p=Path(p);INPUTS.add(p);return json.loads(p.read_text())
def array(p,key):
    p=Path(p);INPUTS.add(p)
    with np.load(p,allow_pickle=False) as z:return z[key].copy()
def title(fig,s,sub):
    fig.text(.5,.995,s,fontsize=15,ha='center',va='top')
    fig.text(.5,.995-.34/fig.get_figheight(),sub,ha='center',fontsize=10,color='#444444')
def save(fig,name):
    fig.get_layout_engine().set(rect=(0,0,1,.995-.7/fig.get_figheight()))
    for ext in ('png','pdf'):
        p=OUT/(name+'.'+ext);fig.savefig(p,bbox_inches='tight');OUTPUTS.append(str(p.relative_to(ROOT)))
    plt.close(fig)
def diag(case,arm,stage=1,iteration=500):
    return RESULTS/'trajectory_diagnostics'/f'{case}__{arm}'/f'checkpoint_s{stage}_i{iteration:03d}.json'
def image(ax,a,label,vmax=20):
    im=ax.imshow(np.asarray(a).reshape(24,24),origin='lower',vmin=0,vmax=vmax,cmap='viridis',interpolation='nearest')
    ax.set_title(label);ax.set_xticks([]);ax.set_yticks([]);return im
def truth(case):
    p=REPO/'research/research_max_struct_r2/sealed_truth'/f'{case}.npz' if case==DEV else REPO/'research/research_max_struct_r4/sealed_truth/confirmation'/f'{case}.npz'
    return array(p,'m')

traj=read(RESULTS/'OFFLINE_TRAJECTORY_SCORES.json')
fact=read(RESULTS/'OFFLINE_FACTORIAL_SCORES.json')
ext=read(RESULTS/'OFFLINE_FACTORIAL_EXTENSION_SCORES.json')
rows=[r for r in traj['rows'] if r['source']=='sealed_replay_checkpoint']

# Only sampled after-step states; lines are visual guides, never a dense trace.
fig,axes=plt.subplots(4,2,figsize=(12,11),sharex=True,layout='constrained')
for i,(case,name) in enumerate(zip(CASES,NAMES)):
    for arm in ('global','onehop'):
        rr=sorted([r for r in rows if r['case']==case and r['arm']==arm],key=lambda r:(r['stage'],r['iteration']))
        x=[r['stage']*500+r['iteration'] for r in rr]
        phi=[read(r['diagnostic_path'].replace('.npz','.json'))['physical']['Phi_after'] for r in rr]
        axes[i,0].semilogy(x,phi,'o--',ms=3,lw=1,color=COLORS[arm],label=arm)
        axes[i,1].plot(x,[r['after']['amplitude_normalized_nrmse'] for r in rr],'o--',ms=3,lw=1,color=COLORS[arm])
    axes[i,0].set_ylabel(name+'\nOriginal $\\Phi$')
    axes[i,1].set_ylabel('Material RMSE / 12')
    for ax in axes[i]:
        for v in (500,1000):ax.axvline(v,color='.65',ls=':',lw=.8)
        ax.grid(alpha=.18)
axes[0,0].legend(frameon=False,ncol=2)
for ax in axes[-1]:ax.set_xlabel('Cumulative iteration (stage boundaries at 500 and 1000)')
title(fig,'R5 auxiliary trajectories: physical fit and material error','Scalar TM · public development · 21 sampled after-step states per arm; dashed segments are guides only')
fig.get_layout_engine().set(rect=(0,0,1,.95))
save(fig,'01_sampled_trajectories')

# Common amplitude and cross limits across all four objects, no individual rescaling.
zs=[]
for case in CASES:
    p=diag(case,'onehop').with_suffix('.npz')
    zs.append({k:array(p,k) for k in ('m_before','u_cell_amplitude','Lu_cell_amplitude','cross_cell')})
ampmax=max(max(z['u_cell_amplitude'].max(),z['Lu_cell_amplitude'].max()) for z in zs)
crossmax=max(abs(z['cross_cell']).max() for z in zs)
fig,axes=plt.subplots(4,4,figsize=(12,11),layout='constrained')
for i,(name,z) in enumerate(zip(NAMES,zs)):
    image(axes[i,0],z['m_before'],name+'\n$m_{before}$')
    for j,k in enumerate(('u_cell_amplitude','Lu_cell_amplitude'),1):
        imamp=axes[i,j].imshow(z[k].reshape(24,24),origin='lower',vmin=0,vmax=ampmax,cmap='magma')
        axes[i,j].set_title(('$\\|u_i\\|_{\\Omega_s}$' if j==1 else '$\\|(Lu)_i\\|_{\\Omega_s}$') if i==0 else '')
    imcross=axes[i,3].imshow(z['cross_cell'].reshape(24,24),origin='lower',cmap='RdBu_r',norm=SymLogNorm(linthresh=crossmax*.01,vmin=-crossmax,vmax=crossmax))
    axes[i,3].set_title('Signed cross (symlog)' if i==0 else '')
    for ax in axes[i]:ax.set_xticks([]);ax.set_yticks([])
fig.colorbar(imamp,ax=axes[:,1:3],shrink=.65,label='Source-weighted amplitude (shared linear scale)')
fig.colorbar(imcross,ax=axes[:,3],shrink=.65,label='$2\\,\\mathrm{Re}\\sum_s\\bar u_{is}(Lu)_{is}\\Omega_{s}$')
title(fig,'Fixed onehop state: spatial amplitudes and interference','Scalar TM · public development · zero-based stage 1 / iteration 500, evaluated at m_before · material scale 0–20')
fig.get_layout_engine().set(rect=(0,0,1,.95))
save(fig,'02_fixed_state_spatial')

# Factorial positions retain source labels and actual historical gate result.
fr={r['arm']:r for r in fact['rows']}
fig=plt.figure(figsize=(11,7),layout='constrained');gs=fig.add_gridspec(2,3,width_ratios=[.85,1,1])
ax=fig.add_subplot(gs[:,0]);im=image(ax,truth(DEV),'Offline GT\nold curved bar')
for i,(arm,lab) in enumerate([('global','Current scalar / Material scalar'),('material_only','Current scalar / Material onehop'),('current_only','Current onehop / Material scalar'),('onehop','Current onehop / Material onehop')]):
    ax=fig.add_subplot(gs[i//2,i%2+1]);r=fr[arm]
    image(ax,array(r['prediction_path'],'m'),lab)
    ax.set_xlabel(f"{arm} ({'new R5' if r['source']=='new_R5_run' else 'historical R4'})\nNRMSE {r['amplitude_normalized_nrmse']:.3f} · IoU {r['otsu_iou']:.3f} · {r['historical_gate_status']}\nPhi {r['physical_objective']:.3g}",fontsize=9)
fig.colorbar(im,ax=fig.axes,shrink=.65,label='Material (shared 0–20)')
title(fig,'Old curved bar: current / material position factorial','Scalar TM · public development · saved final outputs · historical joint contract; no new blind cases')
fig.get_layout_engine().set(rect=(0,0,1,.95))
save(fig,'03_old_bar_factorial')

er={(r['case'],r['arm']):r for r in ext['rows']}
fig,axes=plt.subplots(3,4,figsize=(12,9),layout='constrained')
for i,(case,name) in enumerate(zip(CASES[1:],NAMES[1:])):
    im=image(axes[i,0],truth(case),name+'\nOffline GT')
    for j,arm in enumerate(('global','onehop','material_only'),1):
        r=er[case,arm];image(axes[i,j],array(r['prediction_path'],'m'),arm+(' (new)' if j==3 else ' (history)'))
        axes[i,j].set_xlabel(f"NRMSE {r['amplitude_normalized_nrmse']:.3f} · IoU {r['otsu_iou']:.3f}\nPhi {r['physical_objective']:.3g} · {r['historical_gate_status']}",fontsize=9)
fig.colorbar(im,ax=axes,shrink=.6,label='Material (shared 0–20)')
title(fig,'Material-only extension: two failures and a ring','Scalar TM · public development · historical endpoints used a different platform; not same-run or equal-compute controls')
fig.get_layout_engine().set(rect=(0,0,1,.95))
save(fig,'04_extension_materials')

# Same saved states for all readouts; directional quantities do not use GT.
fig,axes=plt.subplots(4,3,figsize=(13,10),layout='constrained')
kinds=['ordinary','self','plus','minus'];kc=['#777777','#3269ad','#d55e00','#8c6bb1']
for i,(case,name) in enumerate(zip(CASES,NAMES)):
    rr=next(r for r in rows if r['case']==case and r['arm']=='onehop' and r['stage']==1 and r['iteration']==500)
    d=read(diag(case,'onehop'))
    vals=[[d['counterfactual_readouts'][k]['deltaPhi'] for k in kinds],
          [rr['counterfactuals'][k]['delta_nrmse'] for k in kinds],
          [-d['gradients'][k]['exact_state']['cosine_negative_physical_gradient'] for k in kinds]]
    for j in range(3):
        ax=axes[i,j];ax.bar(kinds,vals[j],color=kc);ax.axhline(0,color='.4',lw=.7);ax.grid(axis='y',alpha=.15)
        if j<2:ax.ticklabel_format(axis='y',style='sci',scilimits=(-2,2))
        else:ax.set_ylim(-1,1)
    axes[i,0].set_ylabel(name)
for ax,t in zip(axes[0],['Readout $\\Delta\\Phi$ (lower is better)','Readout $\\Delta$ NRMSE (offline GT)','cos($-g_{state}$, $-g_{physical}$)']):ax.set_title(t)
title(fig,'Same-state readout counterfactuals and gradient directions','Onehop m_before · zero-based stage 1 / iteration 500 · fixed J/Q, prior and original inner allowances · descriptive only')
fig.get_layout_engine().set(rect=(0,0,1,.95))
save(fig,'05_same_state_counterfactuals')

handoff_path=RESULTS/'OFFLINE_HANDOFF_SCORES.json'
if handoff_path.exists():
    handoff=read(handoff_path)
    finals=[r for r in handoff['rows'] if r['state_kind']=='final_prediction']
    fig,axes=plt.subplots(1,4,figsize=(12,4.6),layout='constrained')
    im=image(axes[0],truth(DEV),'Offline GT')
    r=fr['material_only'];image(axes[1],array(r['prediction_path'],'m'),'Full material-only')
    axes[1].set_xlabel(f"NRMSE {r['amplitude_normalized_nrmse']:.3f} · {r['historical_gate_status']}\nPhi {r['physical_objective']:.3g}")
    for ax,r in zip(axes[2:],finals):
        label='Late only\nglobal → material-only' if r['handoff_id']=='global_stage0_to_material_only' else 'Early only\nmaterial-only → global'
        image(ax,array(r['prediction_path'],'m'),label)
        ax.set_xlabel(f"NRMSE {r['amplitude_normalized_nrmse']:.3f} · {r['historical_gate_status']}\nPhi {r['physical_objective']:.3g}")
    fig.colorbar(im,ax=axes,shrink=.65,label='Material (0–20)')
    title(fig,'Natural stage handoffs on the old curved bar','Public development · split after mf3; continue mf6/mf10 · final correctors hit 1000 iterations; no convergence claim')
    fig.get_layout_engine().set(rect=(0,0,1,.90))
    save(fig,'06_natural_stage_handoffs')

captions='''# R5 figure captions and scientific scope

All figures are **OFFLINE_DIAGNOSTIC**, scalar-TM, public-development descriptions, not fresh blind evidence. Lower original physical Phi does not establish material recovery. GT was accessed only after the corresponding prediction/diagnostic seals, and is used here only for offline display and scoring. No online selection, optimization, physical evaluations, algorithm changes or new gates occur in this script.

1. **01_sampled_trajectories**: Four preselected cases, global and onehop auxiliary replay, before the common physical corrector. Each arm contains exactly 21 fixed callback checkpoints (1,2,10,50,150,300,500 in each of three stages); plotted values refer to m_after. Dots are observations; dashed segments are interpolation guides, not continuously observed trajectories. Vertical dotted lines separate stages. These auxiliary checkpoints are distinct from final corrected predictions shown in Figures 3/4/6. NRMSE means RMSE divided by amplitude 12.
2. **02_fixed_state_spatial**: All four onehop replays at zero-based stage 1, iteration 500; every array is evaluated at m_before, before the last material readout. The u and Lu maps use the saved source-weighted cell amplitudes with one shared linear range, without per-object rescaling. Signed cross is 2 Re sum_s(conj(u) Lu Omega_s), with one symmetric symlog range across objects; the linear central band is 1% of the largest absolute cross. Amplitudes and cross have different dimensions and separate colorbars. The material column has shared 0–20 scale. Cross sign alone is not a material-gradient or GT certificate.
3. **03_old_bar_factorial**: Offline GT beside a 2×2 factorial: rows select scalar/onehop current transformation, columns select scalar/onehop material readout. The mixed-position arms are new R5 runs; global/onehop are historical R4 endpoints. PASS/FAIL retains the full historical joint contract (material, prediction-only Otsu IoU, data and box), not a new threshold. One disclosed example does not give a population success rate.
4. **04_extension_materials**: GT, historical global/onehop and new material-only final predictions for two R4 failures and a normal ring. Every material map shares 0–20; annotations report NRMSE, IoU, original Phi and historical joint gate. Historical references ran on a different platform and are not same-execution or equal-compute controls. The two new failure outcomes are retained.
5. **05_same_state_counterfactuals**: At the same four stage-1/iteration-500 onehop m_before states, ordinary/self/plus/minus material readouts share saved J/Q, starting material, prior and inherited inner allowances. Delta means counterfactual minus starting-state value; negative Phi/error deltas mean improvement. The third column uses exact negative state-gradient alignment with the negative physical gradient and no GT. These are finite readouts and local directions, not alternate trajectories or convergence guarantees; prior gradient is not folded into the state-gradient column. The material-error deltas are offline only.
6. **06_natural_stage_handoffs** (when scores exist): Offline GT, full material-only endpoint, and two frozen natural-stage handoffs. Splitting after the original mf3 stage preserves the natural PR reset and coefficient zero-padding before mf6/mf10. Both final handoff outcomes fail the historical joint contract and do not reproduce full material-only success on this single old development case. Both correctors reached 1000 iterations; no convergence or mathematical basin claim is supported. Candidate states are not counted as independent runs.

Reproduce from the repository root with the project Python environment: `python research/research_max_struct_r5/experiments/make_figures.py`. The script imports only NumPy/Matplotlib and standard-library modules. It reads scored JSON and saved arrays; no Maxwell solver is imported. Sealed source artifacts are never written. Figure outputs are regenerated in FIGURES/.
'''
(OUT/'FIGURE_CAPTIONS.md').write_text(captions)
receipt=dict(status='COMPLETE',classification='OFFLINE_DIAGNOSTIC_PUBLIC_DEVELOPMENT_NOT_BLIND',generated_utc=datetime.now(timezone.utc).isoformat(),script=str(Path(__file__).relative_to(ROOT)),input_paths=[str(p) for p in sorted(INPUTS)],outputs=OUTPUTS+['FIGURES/FIGURE_CAPTIONS.md'],new_physical_LU=0,new_complex_RHS=0,physics_imported=False,algorithm_changes=False,GT_use='Offline display/scoring only, after existing seals',checkpoint_count_per_arm=21,checkpoint_sampling='fixed sparse observations; line segments are visual guides',handoff_scores_available=(RESULTS/'OFFLINE_HANDOFF_SCORES.json').exists(),limitations=['Public development, not fresh blind evidence','Historical endpoints are not same-platform/equal-compute controls','Lower Phi is not material recovery'])
(ROOT/'audit/FIGURE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'COMPLETE','figures':len(OUTPUTS)//2,'handoff_scores_available':receipt['handoff_scores_available']}))

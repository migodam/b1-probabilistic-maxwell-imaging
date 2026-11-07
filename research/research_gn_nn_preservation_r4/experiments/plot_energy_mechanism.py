"""Plot saved E0 CSVs only. No physical solver, SVD, training or GT file read."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_k]='1'
import csv
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

HERE=Path(__file__).resolve().parents[1]
RESULTS=HERE/'results/e0'
FIGURES=HERE/'figures/e0'
ARMS=('j_full','j_null','v_physics')
NAMES={'j_full':'J-full','j_null':'J-null','v_physics':'voxel physics'}
PHYSICAL='Maxwell_current'
AUXILIARY='auxiliary_current_head_OFFLINE_not_Maxwell_solution'
RAW='anchor_to_raw_unclipped'
FINAL='anchor_to_final_top1'
AUX='anchor_to_auxiliary_head'

def fmt(x,signed=False):
    if x==0:return '0'
    return f'{x:+.2e}'if signed else f'{x:.2e}'

def saved_tables():
    with (RESULTS/'grouped_absolute_energy.csv').open()as f:
        rows=[r for r in csv.DictReader(f)if r['arm']in ARMS and r['seed']=='7'and r['passes']=='1'and r['level']=='low']
    groups=defaultdict(list)
    for r in rows:groups[r['arm'],float(r['noise']),r['quantity'],r['stage'],r['projection']].append(r)
    def value(arm,noise,quantity,stage,projection,metric):
        rs=groups.get((arm,noise,quantity,stage,projection),[])
        if not rs:return None
        n=sum(int(r['cases'])for r in rs)
        if n!=64:raise ValueError('Saved grouped cohort is not 64 observations: '+str((arm,noise,quantity,stage,projection,n)))
        return sum(float(r[metric+'_sum'])for r in rs)/n
    # Compare identities, not just counts, before connecting unlike quantities.
    sets=defaultdict(set)
    with (RESULTS/'per_case_energy.csv').open()as f:
        for r in csv.DictReader(f):
            if r['arm']in ARMS and r['seed']=='7'and r['passes']=='1'and r['level']=='low':
                sets[r['arm'],float(r['noise']),r['quantity'],r['stage']].add((r['obs_id'],r['parent_id']))
    for arm in ARMS:
        for noise in (0.,.001):
            raw=sets[arm,noise,PHYSICAL,RAW]
            final=sets[arm,noise,PHYSICAL,FINAL]
            if len(raw)!=64 or raw!=final:raise ValueError('Raw/final saved cohort mismatch')
            aux=sets.get((arm,noise,AUXILIARY,AUX),set())
            if arm in ('j_full','j_null')and aux!=raw:raise ValueError('Auxiliary/physical saved cohort mismatch')
            if arm=='v_physics'and aux:raise ValueError('Unexpected voxel auxiliary-head rows')
    for noise in (0.,.001):
        reference=sets['j_full',noise,PHYSICAL,RAW]
        if any(sets[arm,noise,PHYSICAL,RAW]!=reference for arm in ARMS):
            raise ValueError('Selected arms do not share the identical saved mother/observation cohort')
    return value,sets

def twofold(value):
    definitions=('I_to_O','O_to_I','intersection')
    labels=('Sequential I → O: ‖P_O P_I e‖²','Sequential O → I: ‖P_I P_O e‖²','Numerical intersections + remainder')
    ibands=('strong','weak');obands=('strong','weak','null')
    layouts=[];largest=0.
    for noise in (0.,.001):
        for arm in ARMS:
            row=[]
            for definition in definitions:
                before=np.empty((2,3));raw=np.empty((2,3));final=np.empty((2,3))
                for i,ib in enumerate(ibands):
                    for j,ob in enumerate(obands):
                        projection=f'{definition}__{ib}__{ob}'
                        before[i,j]=value(arm,noise,PHYSICAL,RAW,projection,'before')
                        raw[i,j]=value(arm,noise,PHYSICAL,RAW,projection,'delta')
                        final[i,j]=value(arm,noise,PHYSICAL,FINAL,projection,'delta')
                largest=max(largest,float(np.max(abs(final))))
                row.append((before,raw,final))
            layouts.append((noise,arm,row))
    fig,axes=plt.subplots(6,3,figsize=(18,21))
    fig.subplots_adjust(left=.12,right=.97,bottom=.075,top=.91,hspace=.83,wspace=.30)
    norm=SymLogNorm(linthresh=1e-8,linscale=.65,vmin=-largest,vmax=largest,base=10)
    image=None
    for r,(noise,arm,panels)in enumerate(layouts):
        for c,(before,raw,final)in enumerate(panels):
            ax=axes[r,c];image=ax.imshow(final,cmap='RdBu_r',norm=norm,aspect='auto')
            ax.set_xticks(range(3),('S strong9','S weak39','S null528'),fontsize=9)
            ax.xaxis.tick_top();ax.tick_params(axis='x',length=0,pad=5)
            ax.set_yticks(range(2),('I strong19','I weak557'),fontsize=9)
            ax.tick_params(axis='y',length=0)
            ax.set_xticks(np.arange(-.5,3,1),minor=True);ax.set_yticks(np.arange(-.5,2,1),minor=True)
            ax.grid(which='minor',color='white',linewidth=1.5);ax.tick_params(which='minor',bottom=False,left=False)
            if r==0:ax.set_title(labels[c],fontsize=12,pad=32)
            for i in range(2):
                for j in range(3):
                    text=f'GN {fmt(before[i,j])}\nΔraw {fmt(raw[i,j],True)}\nΔfinal {fmt(final[i,j],True)}'
                    rgba=image.cmap(image.norm(final[i,j]));luminance=.2126*rgba[0]+.7152*rgba[1]+.0722*rgba[2]
                    ax.text(j,i,text,ha='center',va='center',fontsize=8.3,color='white'if luminance<.48 else'black',linespacing=1.45)
            if c==0:
                ax.text(-.34,.5,('clean'if noise==0 else'noise 0.001')+'\n'+NAMES[arm],transform=ax.transAxes,ha='center',va='center',fontsize=11,fontweight='bold')
            if c==2:
                before_rem=value(arm,noise,PHYSICAL,RAW,'intersection_remainder','before')
                raw_rem=value(arm,noise,PHYSICAL,RAW,'intersection_remainder','delta')
                final_rem=value(arm,noise,PHYSICAL,FINAL,'intersection_remainder','delta')
                ax.text(.5,-.10,f'Remainder rank47: GN {fmt(before_rem)}; Δraw {fmt(raw_rem,True)}; Δfinal {fmt(final_rem,True)}',transform=ax.transAxes,ha='center',va='top',fontsize=8)
                ax.text(.5,-.23,'Intersection ranks: [0,0,0] / [0,20,509]',transform=ax.transAxes,ha='center',va='top',fontsize=8,color='.25')
            else:
                ax.text(.5,-.10,'Sequential operator energies; not orthogonal joint cells',transform=ax.transAxes,ha='center',va='top',fontsize=8,color='.25')
    fig.suptitle('Historical R3 disclosed diagnostic — NOT new TEST\nActual Maxwell-current error: seed7, one pass, low amplitude, 64 mothers per noise',fontsize=16,y=.985)
    fig.text(.5,.943,'Each cell: original GN absolute energy, signed Δraw and signed Δfinal versus GN. Units: discrete squared norm Σᵢ,ₛ |eᵢₛ|².',ha='center',fontsize=11)
    cax=fig.add_axes((.22,.039,.56,.012));bar=fig.colorbar(image,cax=cax,orientation='horizontal')
    bar.set_label('Cell color: signed absolute Δfinal (symmetric log color scale); blue = decrease, red = increase',fontsize=10)
    fig.text(.5,.009,'Two sequential orders and numerical intersections are different definitions. Negative values are retained; no error-energy shares or physical-power interpretation.',ha='center',fontsize=10)
    fig.savefig(FIGURES/'twofold_absolute_delta.png',dpi=180);plt.close(fig)

def auxiliary_chain(value):
    bands=('total','P_in_strong19','S_strong9','D509')
    titles=('Total current discrepancy','Internal strong19','Receiver strong9','D509')
    colors={'j_full':'#1261a0','j_null':'#d2691e','v_physics':'#238443'}
    fig,axes=plt.subplots(2,4,figsize=(19,9))
    fig.subplots_adjust(left=.065,right=.985,bottom=.16,top=.79,hspace=.40,wspace=.28)
    handles=[]
    for ni,noise in enumerate((0.,.001)):
        for bi,band in enumerate(bands):
            ax=axes[ni,bi]
            for arm in ARMS:
                vals=[value(arm,noise,PHYSICAL,RAW,band,'before'),
                      value(arm,noise,AUXILIARY,AUX,band,'after'),
                      value(arm,noise,PHYSICAL,RAW,band,'after'),
                      value(arm,noise,PHYSICAL,FINAL,band,'after')]
                if arm=='v_physics':
                    # No inferred or interpolated auxiliary stage for voxel arm.
                    ax.plot([0],[vals[0]],'o',color=colors[arm],markersize=5)
                    h,=ax.plot([2,3],vals[2:],'o-',color=colors[arm],label=NAMES[arm],linewidth=1.7)
                else:
                    h,=ax.plot(range(4),vals,'o-',color=colors[arm],label=NAMES[arm],linewidth=1.7)
                if ni==0 and bi==0:handles.append(h)
            ax.set_xticks(range(4),('GN\nphysical','auxiliary\nhead','raw-material\nphysical','final top1\nphysical'),fontsize=8.5)
            ax.set_yscale('log');ax.grid(alpha=.23,which='major')
            ax.set_title(('clean'if noise==0 else'noise 0.001')+' / '+titles[bi],fontsize=11)
            if bi==0:ax.set_ylabel('mean absolute discrepancy energy\n(discrete squared norm)',fontsize=10)
    fig.suptitle('Historical R3 disclosed diagnostic — NOT new TEST\nAuxiliary head → real Maxwell current at raw material → real Maxwell current at final top1',fontsize=15,y=.98)
    fig.text(.5,.875,'Seed7, one pass; identical 64 low-amplitude mothers per noise and identical GT-current reference at every available point.',ha='center',fontsize=11)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.86),ncol=3,frameon=False)
    fig.text(.5,.085,'The auxiliary column is a head discrepancy, not a physical Maxwell-solution current error. Voxel has no auxiliary head (N/A, not zero).',ha='center',fontsize=11)
    fig.text(.5,.045,'Lines connect associated saved stage energies only: they do not isolate a readout/corrector effect or supply missing head→physical vector update/cross terms.',ha='center',fontsize=10)
    fig.savefig(FIGURES/'auxiliary_to_physical_chain.png',dpi=180);plt.close(fig)

def documentation(value,sets):
    lines=['# E0 mechanism figures','',
        '**Scope:** read-only plotting of completed historical R3 E0 CSVs. No new TEST truth, raw GT file, LU, Jacobian, SVD, optimizer or training call is made. The plot script imports no project physics or network module.','',
        'Run `python research/research_gn_nn_preservation_r4/experiments/plot_energy_mechanism.py` from the repository. It writes only the two requested PNGs and this note. Inputs are `results/e0/grouped_absolute_energy.csv` and `results/e0/per_case_energy.csv`. The former supplies saved energy sums; the latter is read only to confirm matching observation/mother identities before connecting stages.','',
        'Both figures use seed7, one pass, low amplitude, clean/noisy separately, with all 64 mothers per noise. The selected arms are J-full, J-null and voxel physics (`v_physics`), the voxel arm with physical training loss. `v_image`, other seeds, cascades and A2.5 are outside these figures, not removed from E0. Amplitude groups are pooled with their saved case counts; means are not unweighted means of amplitude means.','',
        'The unit is the discrete squared current norm Σ over 576 cells and 12 sources of |error|², averaged across the 64 observations. It is not physical power, flux, a continuum integral or an error-energy share. All signs and numeric magnitudes come from saved CSV fields. No new confidence intervals or significance claim are produced.','',
        '## Twofold absolute delta','',
        '`figures/e0/twofold_absolute_delta.png` has six rows (noise × arm) and three strictly separate definition columns. Each small matrix is 2×3: internal strong19/weak557 by receiver strong9/weak39/null528. Cell text gives original GN energy, Δraw versus GN and Δfinal versus GN. Cell color shows signed Δfinal on a shared symmetric-log color scale; negative decreases remain negative and blue. Color scaling changes display only, not the saved values.','',
        '- I→O means ‖P_O P_I e‖². O→I means ‖P_I P_O e‖². They are sequential fixed-operator energies. The resulting joint cells are not labelled mutually orthogonal subspaces or a true intersection partition; the two orders must not be interchanged.',
        '- The numerical-intersection column is a separate definition. Nonzero intersection ranks are weak/weak20 and weak/null509 (D509); the other four cells have rank0. The rank47 orthogonal remainder is printed explicitly below every matrix. Its GN value and both signed deltas are retained. This column does not substitute TSOM projected spans or sequential products for intersections.',
        '- Errors are the actual nonlinear Maxwell-current errors already computed in E0. Raw is the frozen unclipped proposal convention used by that historical audit. All differences between that proposal and the verified raw_1 material were previously documented as float32 lower-box roundoff only (max Linf4.47035e-10), not substantial clipping. No new solve is made here.','',
        '## Auxiliary-to-physical chain','',
        '`figures/e0/auxiliary_to_physical_chain.png` shows total, internal strong19, receiver strong9 and D509 energies, separately for each noise level. The four positions are GN physical current, auxiliary head, actual Maxwell current at the raw material, and actual Maxwell current at the selected final top1.','',
        'The scalar energies are comparable as discrepancies to the same stored GT current under the same fixed projections, same discrete norm and exactly matched 64 observation identities. J-full/J-null have all four saved positions. Voxel physics has no auxiliary-head rows: its missing point is N/A and is neither zero-filled nor joined by an invented intermediate value. The auxiliary point is always labelled a head discrepancy; only the GN/raw/final points are true nonlinear Maxwell-solution current errors.','',
        '**Limits/gaps:** the saved E0 CSVs do not contain the full head→raw-physical vector displacement or its signed projected cross term, a constitutive-defect table, or an isolated readout-only intervention. Final top1 also includes finite-candidate physical correction and selection. Consequently these connecting lines do not establish a causal readout contraction, explain the source of an increase, or isolate a TRF effect. Only the available absolute stage energies are plotted; no missing vector or defect value is reconstructed.','',
        '## Saved scalar anchors for interpreting the chain','',
        '|noise|arm|band|GN physical|auxiliary head|raw physical|final physical|',
        '|---:|---|---|---:|---:|---:|---:|']
    for noise in (0.,.001):
        for arm in ARMS:
            for band in ('S_strong9','D509'):
                vals=[value(arm,noise,PHYSICAL,RAW,band,'before'),value(arm,noise,AUXILIARY,AUX,band,'after'),value(arm,noise,PHYSICAL,RAW,band,'after'),value(arm,noise,PHYSICAL,FINAL,band,'after')]
                lines.append('|'+f'{noise}|{NAMES[arm]}|{band}|'+('|'.join('N/A'if v is None else f'{v:.7g}'for v in vals))+'|')
    lines+=['',
        'For example, clean J-null receiver-strong9 head discrepancy stays near the GN baseline, while the corresponding raw-material real Maxwell-current error is much larger; physical correction lowers that raw error but leaves it above the extremely small GN baseline. This is a descriptive historical stage association, not a new TEST result or a theorem about arbitrary readouts. D509 can improve at the physical endpoints while these strong-direction comparisons differ. The parent retains mechanism interpretation and final scientific judgment.','',
        '## Coverage and writes','',
        'Every plotted physical raw/final group and available auxiliary group has 64 identical observation/mother pairs. Voxel auxiliary data are intentionally absent. No E0 source table or pre-existing source file is edited. The only new outputs are this note, `twofold_absolute_delta.png`, and `auxiliary_to_physical_chain.png`; the new plotting source is `experiments/plot_energy_mechanism.py`.']
    (HERE/'delegated/energy/MECHANISM_FIGURES.md').write_text('\n'.join(lines)+'\n')

def main():
    value,sets=saved_tables()
    twofold(value);auxiliary_chain(value);documentation(value,sets)
    print('COMPLETE: two saved-table figures and MECHANISM_FIGURES.md; no physics calls or new TEST truth reads')

if __name__=='__main__':main()

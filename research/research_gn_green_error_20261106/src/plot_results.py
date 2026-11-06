"""Plot closed offline diagnostic tables/arrays; performs no physical solves."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Patch
P=Path(__file__).resolve().parents[1];R=P/'results';F=P/'figures';F.mkdir(exist_ok=True)
summary=json.loads((R/'SUMMARY.json').read_text())
if not summary['status'].startswith('COMPLETED'):raise RuntimeError('Analysis must close first')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.dpi':180})
BANDS=['strong','weak','complement']; COLORS=['#d57241','#3a89ad','#8a80b3'];AMPS=[.4,2.5,12.];Q=['material','current'];METHODS=['native_trf','reciprocal_trf'];ML=['Native TRF','Reciprocal TRF']
def read(name):
 with (R/name).open() as f:rows=list(csv.DictReader(f))
 for r in rows:
  for k in r:
   if k in ['amplitude','tau','pooled_fraction','allocated_fraction','unassigned_fraction','component_energy','total_energy'] and r[k]!='':r[k]=float(r[k])
 return rows
def select(rows,**kw):return [r for r in rows if all(r.get(k)==v for k,v in kw.items())]
def save(fig,name):
 fig.savefig(F/name,bbox_inches='tight',facecolor='white');plt.close(fig)
single=read('aggregate_single.csv');two=read('aggregate_twofold.csv');rem=read('aggregate_remainder.csv');cases=read('case_metrics.csv')
# 1. Single Green input-energy partitions.
fig,axs=plt.subplots(2,2,figsize=(12,8),sharey=True)
for qi,q in enumerate(Q):
 for oi,op in enumerate(['G_in','G_ob']):
  ax=axs[qi,oi]
  for mi,m in enumerate(METHODS):
   x=np.arange(3)+(mi-.5)*.34;bottom=np.zeros(3)
   for b,c in zip(BANDS,COLORS):
    y=np.array([select(single,quantity=q,operator=op,method=m,noise='clean',tau=.1,amplitude=a,band=b)[0]['pooled_fraction']*100 for a in AMPS])
    ax.bar(x,y,.32,bottom=bottom,color=c,edgecolor='white',linewidth=.6,hatch=None if mi==0 else '//');bottom+=y
  ax.set(title=f'{q.capitalize()} error | {op}',xticks=np.arange(3),xticklabels=['0.4','2.5','12'],ylim=(0,100),xlabel='Contrast amplitude');ax.yaxis.set_major_formatter(PercentFormatter());ax.grid(axis='y',alpha=.2);ax.set_axisbelow(True)
  if oi==0:ax.set_ylabel('Pooled input error energy')
fig.suptitle('Bounded GN-family TRF endpoints: single Green partitions',fontsize=16,y=1.01)
fig.legend(handles=[Patch(facecolor=c,label=b.capitalize()) for b,c in zip(BANDS,COLORS)]+[Patch(facecolor='white',edgecolor='#555',label=ML[0]),Patch(facecolor='white',edgecolor='#555',hatch='//',label=ML[1])],loc='lower center',ncol=5,bbox_to_anchor=(.5,-.02),frameon=False)
fig.text(.5,.025,'Clean observations; 12 shapes per contrast and method; complex Green right-domain projections; primary threshold = 0.1',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.06,1,.98));save(fig,'single_green_error_energy.png')
# 2/3. Joint cells: never normalize joint definitions to 100.
def matrix(q,a,definition):
 rows=select(two,quantity=q,method='native_trf',noise='clean',tau=.1,amplitude=a,definition=definition)
 return np.array([[next(r['pooled_fraction']*100 for r in rows if r['in_band']==i and r['ob_band']==j) for j in BANDS] for i in BANDS])
def heat(ax,z,title):
 im=ax.imshow(z,vmin=0,vmax=100,cmap='YlGnBu');ax.set(title=title,xticks=range(3),yticks=range(3),xticklabels=['Strong','Weak','Complement'],yticklabels=['Strong','Weak','Complement']);ax.set_xlabel('G_ob band');ax.set_ylabel('G_in band')
 for i in range(3):
  for j in range(3):ax.text(j,i,f'{z[i,j]:.2f}%',ha='center',va='center',fontsize=9,color='white' if z[i,j]>48 else '#172b36')
 return im
fig,axs=plt.subplots(2,3,figsize=(13,8),layout='constrained')
for qi,q in enumerate(Q):
 for ai,a in enumerate(AMPS):im=heat(axs[qi,ai],matrix(q,a,'in_then_ob'),f'{q.capitalize()} | amplitude {a:g}')
fig.colorbar(im,ax=axs,shrink=.75,label='Pooled input error energy (%)')
fig.suptitle('Sequential joint energies: G_in then G_ob\nNative TRF, clean observations, threshold = 0.1; each panel sums to 100%',fontsize=15)
save(fig,'twofold_energy_9cells.png')
fig,axs=plt.subplots(2,3,figsize=(14,9),layout='constrained');defs=['true_intersection','in_then_ob','tsom_projected_span'];names=['Numerical intersections','G_in then G_ob','TSOM projected spans']
for qi,q in enumerate(Q):
 for di,(d,name) in enumerate(zip(defs,names)):
  z=matrix(q,12.,d);r=select(rem,quantity=q,method='native_trf',noise='clean',tau=.1,amplitude=12.,definition=d)[0]
  extra=f"Unassigned: {100*r['unassigned_fraction']:.2f}%" if d=='true_intersection' else ('Partition: 100%' if d=='in_then_ob' else f"Overlapping span sum: {100*r['allocated_fraction']:.2f}%")
  im=heat(axs[qi,di],z,f'{q.capitalize()} | {name}\n{extra}')
fig.colorbar(im,ax=axs,shrink=.75,label='Original input error energy (%)')
fig.suptitle('Joint definitions differ: no renormalization across cells\nNative TRF, amplitude 12, clean observations; threshold = 0.1',fontsize=15)
save(fig,'twofold_definitions_comparison.png')
# 4. Fixed spectra and cumulative energy, denominator includes complement.
with np.load(R/'fixed_spectra.npz') as z:spec={k:z[k].copy() for k in z.files}
with np.load(R/'cumulative_mode_energies.npz') as z:cum={k:z[k].copy() for k in z.files}
metadata={r['task_id']:r for r in cases}
fig,axs=plt.subplots(2,3,figsize=(14,8),layout='constrained');ac=['#31947e','#d19527','#9d4375']
for oi,op in enumerate(['G_in','G_ob']):
 ax=axs[oi,0];sv=spec['in_s' if oi==0 else 'ob_s'];ax.semilogy(np.arange(1,577),np.maximum(sv/sv[0],1e-16),color='#244f6b');ax.axhline(.1,color='#c65d38',ls='--',label='Strong cutoff 0.1');ax.axhline(1e-12,color='#777',ls=':',label='Rank cutoff 1e-12');ax.set(title=f'{op}: normalized singular values',xlabel='Right-domain mode index',ylabel='Sigma / sigma_max',ylim=(1e-16,2));ax.legend(fontsize=8);ax.grid(alpha=.2)
 for qi,q in enumerate(Q):
  ax=axs[oi,qi+1]
  if op=='G_ob':
   rank=int(np.count_nonzero(sv>1e-12*sv[0]));ax.axvspan(rank,576,color='#8a80b3',alpha=.09);ax.axvline(rank,color='#777',ls=':',lw=1);ax.text(.51,.05,'Complement mode order is arbitrary',transform=ax.transAxes,ha='center',fontsize=8,color='#555')
  for a,col in zip(AMPS,ac):
   mask=np.array([qq==q and oo==op and metadata[str(t)]['method']=='native_trf' and metadata[str(t)]['noise']=='clean' and metadata[str(t)]['amplitude']==a for t,qq,oo in zip(cum['task_id'],cum['quantity'],cum['operator'])])
   energies=cum['mode_energy'][mask].sum(axis=0);curve=np.cumsum(energies)/energies.sum()*100
   ax.plot(np.arange(1,577),curve,color=col,label=f'Amplitude {a:g}')
  ax.set(title=f'{op}: {q} error',xlabel='Right-domain mode index',ylabel='Cumulative input error energy (%)',ylim=(0,101));ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Fixed complex Green spectra and pooled error mode energies\nNative TRF; clean observations; cumulative denominator includes every input mode',fontsize=15)
save(fig,'operator_spectra_and_error_cumulative.png')
# 5. Predefined digit7 only, magnitudes of complex projected material-error components.
with np.load(R/'error_arrays.npz') as z:err={k:z[k].copy() for k in z.files}
index={str(t):i for i,t in enumerate(err['task_id'])};V=spec['ob_V'];sv=spec['ob_s'];groups=[np.flatnonzero(sv>=.1*sv[0]),np.flatnonzero((sv<.1*sv[0])&(sv>1e-12*sv[0])),np.flatnonzero(sv<=1e-12*sv[0])]
chosen=[select(cases,shape='digit7',method='native_trf',noise='clean',amplitude=a)[0] for a in AMPS]
ims=[]
for r in chosen:
 i=index[r['task_id']];d=err['delta_m'][i];components=[V[:,g]@(V[:,g].conj().T@d) for g in groups];ims.append([err['truth_m'][i],err['estimate_m'][i],np.abs(d),*[np.abs(c) for c in components]])
materialmax=max(float(np.max(x)) for row in ims for x in row[:2]);errmax=max(float(np.max(x)) for row in ims for x in row[2:]);fig,axs=plt.subplots(3,6,figsize=(17,9),layout='constrained');titles=['Ground truth m','Native TRF m','Total |delta m|','G_ob strong |component|','G_ob weak |component|','G_ob complement |component|']
for i,row in enumerate(ims):
 for j,x in enumerate(row):
  im=axs[i,j].imshow(x.reshape(24,24),vmin=0,vmax=materialmax if j<2 else errmax,cmap='viridis' if j<2 else 'magma',interpolation='nearest');axs[i,j].set_xticks([]);axs[i,j].set_yticks([])
  if i==0:axs[i,j].set_title(titles[j],fontsize=10)
  if j==0:axs[i,j].set_ylabel(f'Amplitude {AMPS[i]:g}')
matim=axs[0,0].images[0];eim=axs[0,2].images[0];fig.colorbar(matim,ax=axs[:,:2],shrink=.6,label='Real material contrast m');fig.colorbar(eim,ax=axs[:,2:],shrink=.6,label='Material error / complex component magnitude')
fig.suptitle('Predefined digit7: clean Native TRF endpoints and receiver-Green error components\nShared scales across contrasts; projected components are complex magnitudes, not material candidates',fontsize=15)
save(fig,'digit7_material_error_maps.png')
print(json.dumps({'status':'PLOTTED_FROM_CLOSED_RESULTS','cases':len(cases),'figures':[str(f) for f in sorted(F.glob('*.png'))]}))

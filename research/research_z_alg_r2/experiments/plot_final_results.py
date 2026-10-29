"""Paper figures from already-paid, frozen results. No forward solves."""
import json,pathlib,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle
P=pathlib.Path(__file__).resolve().parents[1]
def read(n):return json.loads((P/n).read_text())
def decode(x):
 if isinstance(x,dict) and set(x)=={'real','imag'}:return complex(x['real'],x['imag'])
 if isinstance(x,list):return [decode(v) for v in x]
 if isinstance(x,dict):return {k:decode(v) for k,v in x.items()}
 return x
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
s=read('results/quotient_holdout_summary.json');methods=s['methods'];colors=['#1e7d71','#947442','#546c96','#8e5a83']
fig,ax=plt.subplots(1,3,figsize=(14.8,4.4),gridspec_kw={'width_ratios':[1.3,1.1,1.5]})
labels=['Quotient Q6','Six-path','Reduced TRF ×8','Reduced TRF ×64'];times=[d['total_wall_including_target_logging'] for d in methods]
ax[0].barh(labels,times,color=colors);ax[0].invert_yaxis();ax[0].set_xscale('log');ax[0].set_xlabel('Total seconds, 4 fresh targets (log scale)');ax[0].set_title('A  Same returned physical root sets')
for i,t in enumerate(times):ax[0].text(t*1.08,i,f'{t:.4f} s',va='center',fontsize=9)
ax[0].set_xlim(.01,3)
casecounts=np.array([[c['qualified_counts'][m['method']] for c in s['cases']] for m in methods])
for i in range(4):ax[1].bar(np.arange(4)+(i-1.5)*.18,casecounts[i],width=.17,label=labels[i],color=colors[i])
ax[1].set_xticks(range(4),[c['case'].split('_')[-1] for c in s['cases']]);ax[1].set_ylim(0,2.65);ax[1].set_yticks([0,1,2]);ax[1].set_xlabel('Fresh target');ax[1].set_ylabel('Qualified material roots');ax[1].set_title('B  Truth contained in 4/4 sets')
ax[1].legend(fontsize=7,ncol=2,loc='upper left',frameon=False)
for method in methods[:1]:
 file=P/'results/quotient_holdout_QH4_03_quotient.json';d=decode(json.loads(file.read_text()))
 result=d.get('algebraic_result',d.get('homotopy',d.get('quotient')))
 if result is None:
  result=next(v for v in d.values() if isinstance(v,dict) and 'paths' in v)
 x=np.asarray([p['x'] for p in result['paths']]);ax[2].scatter(x[:,0].real,x[:,0].imag,s=65,c='#287e8e',edgecolor='white')
 for i,z in enumerate(x[:,0]):ax[2].annotate(str(i+1),(z.real,z.imag),xytext=(5,5),textcoords='offset points')
 ax[2].set_xlabel('Re(c₁)');ax[2].set_ylabel('Im(c₁)');ax[2].set_title('C  Six complex candidates before filtering')
fig.suptitle('Fresh blind holdout: identical material coverage, cheaper finite algebraic enumeration',fontsize=13)
fig.text(.5,.01,'All methods: 16 full RHS each (6 setup + 10 validation). Small-system operations differ. One timed batch; numerical, not certified completeness.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.055,1,.93])
for ext in ['png','pdf','svg']:fig.savefig(P/f'figures/quotient_fresh_holdout.{ext}',dpi=190)
plt.close(fig)
fig,ax=plt.subplots(figsize=(14,3.6));ax.axis('off')
steps=[('Measured data','s complex amplitudes'),('Exact compatibility','2s paired bilinear equations'),('Eliminate conjugate coordinates','maximal (s+1)-minors'),('Fixed quotient basis','monomials of degree ≤ s'),('Enumerate candidates','D × D eigensolve; D = binom(2s,s)'),('Physical materials','reality + box + full Maxwell')]
for i,(title,sub) in enumerate(steps):
 x=.01+i*.166
 ax.add_patch(FancyBboxPatch((x,.35),.145,.4,boxstyle='round,pad=.008',facecolor='#edf5f3' if i in [3,4] else '#f2f4f7',edgecolor='#63928a' if i in [3,4] else '#99a5b3',transform=ax.transAxes))
 words=title.split();line='';parts=[]
 for w in words:
  if len(line+w)>17:parts.append(line);line=w+' '
  else:line+=w+' '
 parts.append(line)
 ax.text(x+.0725,.61,'\n'.join(p.strip() for p in parts),ha='center',va='center',transform=ax.transAxes,fontsize=10,fontweight='bold')
 ax.text(x+.0725,.4,sub,ha='center',va='center',transform=ax.transAxes,fontsize=8)
 if i<5:ax.annotate('',xy=(x+.158,.55),xytext=(x+.147,.55),xycoords='axes fraction',arrowprops={'arrowstyle':'->','color':'#425466'})
ax.text(.5,.88,'Why the representation changes the computation',ha='center',fontsize=16,fontweight='bold',transform=ax.transAxes)
ax.text(.5,.18,'Conditions: fixed Green · 2s independent real materials · one receiver / s sources · full source rank · nonsingular leading block T',ha='center',fontsize=10,transform=ax.transAxes)
ax.text(.5,.07,'Implemented: four materials / D=6 and six materials / D=20. Known algebraic tools; physical specialization and measured results are the claim candidate.',ha='center',fontsize=10,color='#52606d',transform=ax.transAxes)
fig.tight_layout()
for ext in ['png','pdf','svg']:fig.savefig(P/f'figures/algebraic_mechanism.{ext}',dpi=190)
print('Saved two figure sets from frozen arrays, 0 RHS.')

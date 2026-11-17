from complete_quotient import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
begin=time.process_time();h=json.loads((OUT/'results/H12_probe.json').read_text());ct=json.loads((OUT/'results/H12_contract.json').read_text());fig,ax=plt.subplots(1,2,figsize=(11,4.7),constrained_layout=True);c=np.array(ct['centers']);ax[0].scatter(c[:6,0],c[:6,1],s=140,label='6 independent unknown cells');ax[0].scatter(c[6:,0],c[6:,1],s=140,marker='s',label='6 known host cells')
for i,(xx,yy) in enumerate(c):ax[0].text(xx,yy,str(i),ha='center',va='center',fontsize=8,color='white')
ax[0].axis('equal');ax[0].set(title='H12: n = 12, k = 6, s = 3',xlabel='x',ylabel='y');ax[0].legend(fontsize=8,loc='upper left')
for i,p in enumerate(h['predictions']):ax[1].plot(range(6),p['material'],'o-',label=f'Certified explanation {i+1}')
ax[1].set(title='Two equal-old-response explanations',xlabel='Active cell',ylabel='Material contrast');ax[1].legend(fontsize=8);fig.supxlabel(f"One realizable source gap: {h['source']['predicted_gap']:.5f}; reduced/full difference: {h['full_reduced_error']:.2e}.\nReduced-model certificate only; no noise or host-error guarantee in this probe.",fontsize=10);fig.savefig(OUT/'figures/07_known_host.png',dpi=180);save(OUT/'results/host_plot_revision_receipt.json',dict(CPU_seconds=time.process_time()-begin))

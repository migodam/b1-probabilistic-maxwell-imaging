from pathlib import Path
import json,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parents[1]; start=time.perf_counter()
fit=json.loads((R/'results/active/FIT.json').read_text());score=json.loads((R/'results/active/SCORE.json').read_text())
coverage=json.loads((R/'results/active/COVERAGE.json').read_text())
fig=plt.figure(figsize=(12,4.5),layout='constrained')
ax=fig.add_subplot(121)
rows=sorted(fit['trials'],key=lambda r:r['t'])
ax.plot([r['t'] for r in rows],[r['residual'] for r in rows],'o-',label='New source: evaluated candidates')
ax.axhline(6,color='r',ls='--',label='New-source noise radius')
ax.axhline(fit['old_family_residual'],color='gray',ls=':',label='Old-source residual (all t)')
ax.axvline(score['truth_t'],color='k',ls=':',label='Truth (post-seal only)')
ax.set(xlabel='Known-family parameter t',ylabel='Whitened residual norm',title='One additional illumination; one scalar unknown')
ax.legend(fontsize=8)
ax=fig.add_subplot(122,projection='3d')
a=-.12+.03/2+np.arange(8)*.03;p=np.stack(np.meshgrid(a,a,a,indexing='ij'),-1).reshape(-1,3)
q=np.sum(p*p,axis=1)/.1**2;b=np.zeros(len(p));ok=q<1;b[ok]=np.exp(1-1/(1-q[ok]))
s=ax.scatter(*p[ok].T,c=1+fit['t']*b[ok],vmin=1,vmax=4,cmap='viridis',s=16)
ax.set(xlabel='x / lambda',ylabel='y / lambda',zlabel='z / lambda',title='Fitted epsilon_zz: KNOWN spatial family')
fig.colorbar(s,ax=ax,shrink=.6)
fig.savefig(R/'figures/ACTIVE_FAMILY.png',dpi=180)
plt.close(fig)
(R/'audit/ACTIVE_PLOT.json').write_text(json.dumps(dict(wall_s=time.perf_counter()-start,LU=0,RHS=0,
    truth_access='POST_SEAL_SCORING',interpolation='lines connect evaluated points only, not confidence bounds'),indent=2))

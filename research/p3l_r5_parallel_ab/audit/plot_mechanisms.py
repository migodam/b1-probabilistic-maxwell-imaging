from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parents[1]
seam=json.loads((P/'audit/seam_maxwell/RESULTS.json').read_text())
bounds=json.loads((P/'audit/curvature_bounds/RESULTS.json').read_text())
f,ax=plt.subplots(1,2,figsize=(12,4.3),layout='constrained')
r=seam['rows'];x=np.array([v['h'] for v in r]);y=np.array([v['raw_response_distance'] for v in r])
ax[0].loglog(x,y,'o-',label='Full Maxwell secant')
ax[0].loglog(x,y[0]*(x/x[0])**3,'--',label='Cubic reference')
ax[0].set(xlabel='Latent displacement from host seam',ylabel='Raw real response distance',title=f"Finite response slope = {seam['loglog_slope']:.4f}")
ax[0].legend();ax[0].grid(True,which='both',alpha=.2)
for center in bounds['config']['centers']:
    a=[v for v in bounds['rows'] if v['center']==center]
    ax[1].loglog([v['noise_multiplier'] for v in a],[v['max_sufficient_response_R'] for v in a],'o-',label=center)
ax[1].axhline(3,color='#b33426',ls='--',label='Requested response radius R=3')
ax[1].set(xlabel='Noise / original sigma',ylabel='Maximum radius from sufficient bound',title='Norm bound is too conservative here')
ax[1].legend(fontsize=8);ax[1].grid(True,which='both',alpha=.2)
f.suptitle('Full TM24 numerical probes; fixed submodel/direction, not a posterior certificate',fontsize=12)
dest=P/'figures';dest.mkdir(exist_ok=True)
f.savefig(dest/'SEAM_AND_BOUND.png',dpi=180)
plt.close(f)

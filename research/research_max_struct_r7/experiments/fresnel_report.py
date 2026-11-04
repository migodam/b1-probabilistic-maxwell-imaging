from pathlib import Path
import os
os.environ['MPLBACKEND']='Agg'
import numpy as np,json,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
rows=json.loads((ROOT/'RESULTS/fresnel_summary.json').read_text())
fig,axs=plt.subplots(2,3,figsize=(13,8),constrained_layout=True)
for i,case in enumerate(['dielTM_dec8f','twodielTM_8f']):
 for j,arm in enumerate(['cie_dual','global_dual','a_material']):
  p=ROOT/'RESULTS'/f'fresnel_{case}_{arm}';d=np.load(p/'endpoint.npz');r=json.loads((p/'RECEIPT.json').read_text())
  im=axs[i,j].imshow(d['m'].reshape(24,24),origin='lower',extent=(-80,80,-80,80),vmin=0,vmax=20,cmap='viridis')
  axs[i,j].set(title=f'{case} / {arm}\ntrain {r["train_relative_complex_residual"]:.3f}; withheld {r["withheld_all_relative_complex_residual"]:.3f}',xlabel='x [mm]',ylabel='y [mm]')
fig.colorbar(im,ax=axs,label='Reconstructed disk-pulse contrast (no pixel GT)');fig.suptitle('Actual Fresnel 2001 measurements, 2 GHz; same physical endpoint in all three arms')
fig.savefig(ROOT/'FIGURES/fresnel_material_endpoints.png',dpi=160);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(12,4),constrained_layout=True)
for ax,case in zip(axs,['dielTM_dec8f','twodielTM_8f']):
 d=np.load(ROOT/'RESULTS'/f'fresnel_{case}_cie_dual'/'endpoint.npz');mask=np.isfinite(d['measured_Y']);meas=d['measured_Y'][mask];pred=d['prediction'][mask]
 ax.scatter(meas.real,pred.real,s=6,alpha=.4,label='real');ax.scatter(meas.imag,pred.imag,s=6,alpha=.4,label='imag')
 lim=max(abs(meas).max(),abs(pred).max());ax.plot([-lim,lim],[-lim,lim],'k--',lw=.7);ax.set(title=case,xlabel='Measured scattered field [unstated units]',ylabel='Predicted scattered field [same units]');ax.legend()
fig.suptitle('All measured channels; no object-scattered gain refit on held-out channels');fig.savefig(ROOT/'FIGURES/fresnel_complex_predictions.png',dpi=160)

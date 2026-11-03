"""Plot saved finite-vector-DDA mechanism only; no physics/GT imports or solves."""
from pathlib import Path
from datetime import datetime, timezone
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/'RESULTS/vector_pair'
OUT=ROOT/'FIGURES';OUT.mkdir(exist_ok=True)
cfg=json.loads((ROOT/'configs/vector_pair_frozen.json').read_text())
saved_cfg=json.loads((SRC/'FROZEN_CONFIG.json').read_text())
spec=json.loads((SRC/'SPECTRUM.json').read_text())
original_receipt=json.loads((SRC/'RECEIPT.json').read_text())
with np.load(SRC/'FIELDS.npz',allow_pickle=False) as z:
    centers=z['centers'].copy();directions=z['directions'].copy();polarizations=z['polarizations'].copy()
plt.rcParams.update({'font.size':11,'pdf.fonttype':42,'savefig.dpi':180,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(12,6));gs=fig.add_gridspec(1,2,width_ratios=[1,1.1],left=.02,right=.98,bottom=.24,top=.76,wspace=.16)
ax=fig.add_subplot(gs[0],projection='3d')
h=saved_cfg['h']
for center in centers:
    v=center+np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]])*h/2
    faces=[v[idx] for idx in ([0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7])]
    ax.add_collection3d(Poly3DCollection(faces,facecolor='#90b5d7',edgecolor='#496d8b',alpha=.20,lw=.8))
ax.plot(centers[:,0],centers[:,1],centers[:,2],'o--',color='#444444',lw=1,ms=4)
ax.text(-.025,.025,-.034,'pair axis: x',fontsize=10)
# Incident vector and two separate transverse polarizations are shown upstream.
origin=np.array([0.,-.075,0.])
ax.quiver(*origin,*(directions[0]*.045),color='#333333',arrow_length_ratio=.2,lw=2)
ax.text(.005,-.075,-.012,'k along +y',fontsize=10)
for pol,color,label in zip(polarizations,['#d55e00','#3269ad'],['E along x','E along z']):
    ax.quiver(*origin,*(pol*.035),color=color,arrow_length_ratio=.2,lw=2)
    endpoint=origin+pol*.041
    ax.text(*endpoint,label,color=color,fontsize=10)
ax.set(xlim=(-.055,.055),ylim=(-.085,.045),zlim=(-.04,.055),xlabel='$x / \\lambda_0$',ylabel='$y / \\lambda_0$',zlabel='$z / \\lambda_0$')
ax.set_xticks([-.04,0,.04]);ax.set_yticks([-.08,0,.04]);ax.set_zticks([-.04,0,.04]);ax.tick_params(labelsize=8,pad=0)
ax.set_box_aspect((1.1,1.3,1));ax.view_init(elev=23,azim=-58)
ax.set_title('Two known isotropic cells in 3D\n$h=0.04\\lambda_0$, $\\epsilon_r=3$, $\\mu_r=1$',pad=12)

ax=fig.add_subplot(gs[1]);pols=['x_longitudinal','y_transverse','z_transverse'];colors=['#d55e00','#3269ad']
mode_rows=[r for r in spec['spatial_polarization_modes'] if r['beta']==.2]
weights={}
for j,mode in enumerate(['coherent','antisymmetric']):
    rr=[next(r for r in mode_rows if r['polarization']==p and r['spatial_mode']==mode) for p in pols]
    vals=[r['residual_weight_abs_1_plus_a_lambda_squared'] for r in rr]
    weights[mode]=vals
    xx=np.arange(3)+(j-.5)*.33
    ax.bar(xx,vals,width=.30,color=colors[j],label=mode+' (+,+)' if j==0 else mode+' (+,−)',alpha=.9)
    for x,y,r in zip(xx,vals,rr):
        ax.text(x,y+.025,f'{y:.3f}'+(' *' if r['physically_excited_by_fixed_waves'] else ''),ha='center',fontsize=10)
ax.axhline(1,color='#444444',ls='--',lw=1,label='Unmodified metric = 1')
ax.set_xticks(range(3),['x\nparallel to pair','y\ntransverse to pair','z\ntransverse to pair'])
ax.set_ylim(0,1.9);ax.set_ylabel('Residual metric weight $|1+a\\lambda|^2$')
ax.set_title('Six fixed vector modes at $\\beta=0.2$\nCoherent x enhanced; coherent y/z attenuated')
ax.legend(loc='upper left',fontsize=9,frameon=False,ncol=2);ax.grid(axis='y',alpha=.16)
fig.text(.5,.97,'Vector DDA two-cell residual metric; not imaging',ha='center',va='top',fontsize=17)
fig.text(.5,.90,'Saved full-vector DDA: one 6×6 solve / two RHS; plot uses saved arrays',ha='center',fontsize=11,color='#444444')
fig.text(.5,.12,'“Longitudinal x” means parallel to the cell-pair axis, not parallel to k. Both incident E vectors are transverse to k.',ha='center',fontsize=10)
fig.text(.5,.075,'* Excited by the two saved plane waves: coherent x and z only. The y and antisymmetric bars are metric modes, not extra illuminations.',ha='center',fontsize=10)
fig.text(.5,.03,'No inverse recovery, 3D imaging, continuum-PDE qualification or convergence claim.',ha='center',fontsize=10,color='#444444')
paths=[]
for suffix in ('png','pdf'):
    path=OUT/('07_vector_pair_modes.'+suffix);fig.savefig(path,bbox_inches='tight');paths.append(str(path.relative_to(ROOT)))
plt.close(fig)
caption='''# Figure 7: finite vector-DDA pair mechanism

**Vector DDA two-cell residual metric; not imaging.** Left: true saved centers (−0.02,0,0) and (0.02,0,0), cubic cells of side 0.04 wavelength, known isotropic epsilon_r=3 and mu_r=1. Two plane waves propagate along +y and have transverse electric polarizations x and z. The arrows show the two illuminations separately, not their coherent superposition. The x mode is longitudinal only relative to the cell-pair separation; it is not an inadmissible longitudinal incident plane wave. Right: all six saved coherent/antisymmetric × Cartesian polarization residual-metric modes at beta=0.2, using the actual complex retarded Green and its original self term, not the electrostatic approximation. The horizontal reference is unmodified weight 1. Coherent x has weight 1.342837, while coherent y/z have weight 0.856193; antisymmetric weights reverse the directional trend (x 0.707605, y/z 1.154966). Asterisks identify coherent x/z, the only modes physically excited by these two fixed plane waves. Coherent y and antisymmetric modes are displayed as operator modes and do not represent additional physical illuminations.

This is a finite two-cell CM+RR DDA mechanism example, a counterexample to universal residual attenuation, not a demonstration of 3D imaging, unknown-material recovery, continuum-PDE validation or algorithm convergence. Plotting imports only saved data and plotting libraries, performs zero new physics and accesses no GT. Reproduce with the project Python: `python research/research_max_struct_r5/experiments/plot_vector_pair.py`. The original six figures are unchanged.
'''
(OUT/'07_VECTOR_PAIR_CAPTION.md').write_text(caption)
receipt=dict(status='COMPLETE',classification='FIXED_FINITE_VECTOR_DDA_MECHANISM_NOT_IMAGING',generated_utc=datetime.now(timezone.utc).isoformat(),source_paths=[str(p.relative_to(ROOT)) for p in [ROOT/'configs/vector_pair_frozen.json',SRC/'FROZEN_CONFIG.json',SRC/'FIELDS.npz',SRC/'SPECTRUM.json',SRC/'RECEIPT.json']],source_status=original_receipt['status'],beta=.2,mode_count=6,weights=weights,outputs=paths+['FIGURES/07_VECTOR_PAIR_CAPTION.md'],new_physical_LU=0,new_complex_RHS=0,physics_imported=False,GT_access='NONE',original_six_figures_modified=False)
(ROOT/'audit/VECTOR_PAIR_FIGURE_RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':'COMPLETE','outputs':paths}))

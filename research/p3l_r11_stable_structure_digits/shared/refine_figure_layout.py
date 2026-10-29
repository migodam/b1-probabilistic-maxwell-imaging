"""Presentation-only grid layout using sealed arrays; no physics or inference."""
from core import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
METHODS=['FULL256','AGENT_A','AGENT_B','AGENT_C'];NAMES=['GT','Full256','A exact support','B total field','C global bound']
z=np.load(P/'DIGIT_RECONSTRUCTION_ARRAYS.npz')
for ds in ['N','P']:
 prefix='' if ds=='N' else 'P_';data=[z[ds+'_GT_occupancy']]+[z[ds+'_'+m+'_occupancy'] for m in METHODS]
 for binary in [False,True]:
  fig,axs=plt.subplots(9,5,figsize=(11,18));fig.subplots_adjust(left=.055,right=.99,bottom=.01,top=.95,wspace=.13,hspace=.13)
  for i in range(9):
   for j in range(5):
    a=data[j][i].reshape(24,24);a=a>=.5 if binary else a;ax=axs[i,j];ax.imshow(a,vmin=0,vmax=1,cmap='gray_r',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
    if i==0:ax.set_title(NAMES[j],fontsize=10,pad=9)
    if j==0:ax.set_ylabel(str(i+1),rotation=0,labelpad=12,fontsize=14)
  fig.suptitle(f'{ds}: '+('binary occupancy masks; fixed threshold 0.5' if binary else 'same inverse objective and budget; physical 24 x 24 grid'),y=.985,fontsize=13)
  fig.savefig(P/'figures'/f'{prefix}DIGITS_{"BINARY_MASKS" if binary else "SIDE_BY_SIDE"}.png',dpi=140);plt.close(fig)
 gm=z[ds+'_GT_material'];pm=z[ds+'_FULL256_material'];dd=[gm.real,pm.real,gm.imag,pm.imag];names=['GT Re(chi)','Full Re(chi)','GT Im(chi)','Full Im(chi)']
 fig,axs=plt.subplots(9,4,figsize=(9,18));fig.subplots_adjust(left=.055,right=.99,bottom=.01,top=.95,wspace=.13,hspace=.13)
 for i in range(9):
  for j in range(4):
   ax=axs[i,j];ax.imshow(dd[j][i].reshape(24,24),vmin=.1 if j<2 else .04,vmax=1.9 if j<2 else .22,cmap='viridis',interpolation='nearest');ax.set_xticks([]);ax.set_yticks([])
   if i==0:ax.set_title(names[j],fontsize=11,pad=9)
   if j==0:ax.set_ylabel(str(i+1),rotation=0,labelpad=12,fontsize=14)
 fig.suptitle(f'{ds}: actual complex material, common real / imaginary ranges',y=.985,fontsize=12);fig.savefig(P/'figures'/f'{prefix}DIGITS_MATERIAL_REAL_IMAG.png',dpi=140);plt.close(fig)
write(P/'audit/PRESENTATION_REFINEMENT.json',dict(utc=time.time(),new_physics_RHS=0,changes='separatetitleandcolheaders; addbinarymaskgalleryfromsamearrays; nocolororvaluealteration',GT_or_prediction_arrays_changed=False))

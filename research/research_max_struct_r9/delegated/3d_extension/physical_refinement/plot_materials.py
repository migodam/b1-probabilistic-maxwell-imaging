"""Offline figure-only revision: fixed known paper-cube z=-0.65lambda plane."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import time
T=time.perf_counter();C=time.process_time()
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'RESULTS/original_vector3d/physical_refinement';FIG=ROOT/'FIGURES/original_vector3d/physical_refinement'
assert (OUT/'SEALED_ALL_THREE.json').exists() and (OUT/'SCORES.json').exists()
truth=np.load(ROOT/'RESULTS/original_vector3d/sealed_truth/truth.npz',allow_pickle=False)['inverse_material'].reshape(30,30,30)
rows=json.loads((OUT/'SCORES.json').read_text())['rows'];materials=[truth]+[np.load(OUT/r['start']/'prediction.npz',allow_pickle=False)['material'].reshape(30,30,30) for r in rows]
first=FIG/('materials_attempt2_colorbar_clip.png' if (FIG/'materials_attempt1_central_plane.png').exists() else 'materials_attempt1_central_plane.png');assert not first.exists();(FIG/'materials.png').rename(first)
fig,axs=plt.subplots(2,4,figsize=(12,6))
for col,(name,m) in enumerate(zip(['Known truth']+[r['start'] for r in rows],materials)):
    im=axs[0,col].imshow(m.max(2),vmin=0,vmax=4,cmap='viridis');axs[0,col].set_title(name+' max-z',fontsize=11)
    axs[1,col].imshow(m[:,:,8],vmin=0,vmax=4,cmap='viridis');axs[1,col].set_title('z=-0.65 lambda; '+('truth' if col==0 else 'NRMSE %.3f'%rows[col-1]['contrast_NRMSE']),fontsize=10)
for a in axs.ravel():a.set_xticks([]);a.set_yticks([])
fig.suptitle('Exposed-development CM30 finite-DDA materials; chi=[0,4]')
fig.tight_layout(rect=[0,0,.94,.96]);cb=fig.add_axes([.95,.13,.012,.72]);fig.colorbar(im,cax=cb,label='contrast chi')
fig.savefig(FIG/'materials.png',dpi=170,bbox_inches='tight');plt.close(fig)
p=ROOT/'delegated/3d_extension/physical_refinement/REPORT_ZH.md';s=p.read_text();s+='\n材料图下排使用固定 z=-0.65λ（index8）切片，对应已知论文立方体 z=-0.7λ 附近；中心 z=0.05λ 原本位于两层目标间隙，首版图另存 materials_attempt1_central_plane.png。图仅离线展示，不选迭代/材料/评分；颜色范围固定 χ∈[0,4]。\n';p.write_text(s) if '材料图下排使用固定' not in p.read_text() else None
(OUT/'LOCAL_FIGURE_REVISION_COST.json').write_text(json.dumps(dict(platform='Mac',CPU_seconds=time.process_time()-C,wall_seconds=time.perf_counter()-T,additive=True,included_in=None,operation='OFFLINE_FIXED_PAPER_CUBE_PLANE_FIGURE_REVISION',scope='Plot-only after allthree sealed/scored; no new metric or solve; attempt1 retained',counts=dict(new_physical_RHS=0)),indent=2)+'\n')

"""Plot stored phase-best checkpoints, never synthesize optimizer trajectories."""
import json,pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=pathlib.Path(__file__).resolve().parent
data=json.loads((P/'analysis/COARSE_PATH_ATTRIBUTION.json').read_text())
colors={'native_trf':'#346c9c','reciprocal_trf':'#d4823f','data_lbfgs':'#508c65','globaldual':'#b6427c'}
fig,axes=plt.subplots(2,2,figsize=(12,7.2))
for i,shape in enumerate(('ring','curved_bar')):
    for row in data:
        if row['shape']!=shape:continue
        points=[]
        for c in row['phase_candidates']:
            x={8:0,12:1}.get(c['stage'])
            if x is None:x=2 if c['origin'].endswith(('proposal','phase1')) else 3
            points.append((x,c))
        points.sort(key=lambda pair:pair[0])
        for j,key in enumerate(('relative_residual','NRMSE_OFFLINE')):
            axes[i,j].semilogy([p[0] for p in points],[p[1][key] for p in points],'-o',ms=4,
                color=colors[row['method']],label=row['method'])
    for j in range(2):
        axes[i,j].set_xticks(range(4),['8x8','12x12','24x24 phase1','24x24 corrected'],fontsize=8)
        axes[i,j].set_title(shape+' | '+('full data residual' if j==0 else 'offline NRMSE'))
        axes[i,j].grid(alpha=.2);axes[i,j].legend(fontsize=8)
fig.suptitle('Intermediate contrast A=2.5: large errors already present before final GlobalDual phase',fontsize=12)
fig.text(.01,.015,'Only frozen phase-best checkpoints are connected. Material scores are offline. Lines are not continuous trajectories, exact roots or independent branches.',fontsize=8)
fig.tight_layout(rect=(0,.04,1,.94));fig.savefig(P/'figures/coarse_path_attribution.png',dpi=180);plt.close(fig)

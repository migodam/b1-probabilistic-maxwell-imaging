"""Pure artifact aliases from already sealed estimates; no PDE or truth input."""
from pathlib import Path
import json,shutil
import numpy as np
def complete_output_aliases(out):
    out=Path(out)
    if not (out/'endpoint.npz').exists():return
    if (out/'target.png').exists():shutil.copyfile(out/'target.png',out/'reconstructed_target.png')
    with np.load(out/'endpoint.npz') as p:maps={k:p[k].copy() for k in ['pixel_x','pixel_y','image_eps','image_loss']}
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(9,3),constrained_layout=True)
    for ax,key,title in zip(axes,['image_eps','image_loss'],['Reconstructed wall real permittivity','Reconstructed wall loss']):
        im=ax.imshow(maps[key],origin='lower',extent=[maps['pixel_x'][0],maps['pixel_x'][-1],maps['pixel_y'][0],maps['pixel_y'][-1]],aspect='equal',cmap='magma');ax.set(xlim=(-1.5,1.5),ylim=(.5,1.25),xlabel='x / wavelength',ylabel='y / wavelength',title=title);fig.colorbar(im,ax=ax,shrink=.75)
    fig.savefig(out/'reconstructed_wall.png',dpi=180,bbox_inches='tight');plt.close(fig)
    (out/'ARTIFACT_ALIASES.json').write_text(json.dumps(dict(source='already sealed endpoint.npz and target.png; no new physics',target='reconstructed_target.png',wall='reconstructed_wall.png',wall_extent=[-1.5,1.5,.5,1.25]),indent=2)+'\n')
if __name__=='__main__':
    import sys
    for p in sys.argv[1:]:complete_output_aliases(p)

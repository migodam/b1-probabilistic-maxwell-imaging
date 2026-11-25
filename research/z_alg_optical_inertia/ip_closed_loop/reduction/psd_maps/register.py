"""Freeze selected comparisons from already public map arrays."""
from pathlib import Path
import json

import numpy as np

HERE = Path(__file__).resolve().parent
MAPS = HERE.parents[1]/'maps'


def register():
    nodes=[]
    for obj in ('M32-R01','M32-R02','M32-R03','M32-R04'):
        oid = f'ipcl-fov-f1__{obj}__weak__7tx9rx__n03'
        f=MAPS/'results'/f'{oid}__b0.03125__broad_sector.npz'
        with np.load(f,allow_pickle=False) as z:
            margin=z['margin_sector']; nu=z['negative_index']
            choose = np.flatnonzero((margin>0)&(nu>0))
            if not len(choose):
                raise ValueError(f'no archived indefinite positive at {oid}')
            ordered = choose[np.argsort(margin[choose],kind='stable')]
            if obj=='M32-R01':
                selected=[('weakest_archived_indefinite_exclusion',ordered[0]),
                          ('median_archived_indefinite_exclusion',ordered[len(ordered)//2])]
            else:
                selected=[('median_archived_indefinite_exclusion',ordered[len(ordered)//2])]
            for why,k in selected:
                nodes.append(dict(side=32,oid=oid,shape='ring',family='broad',index=int(k),reason=why,
                                  old_strict_margin=float(margin[k]), old_negative_index=int(nu[k])))
    for obj,shape,family in [('M96-R01','ring','broad'),('M96-R01','ring','primary'),('M96-P01','double','primary')]:
        oid=f'ipcl-fov-f1__{obj}__weak__7tx9rx__n03'
        with np.load(MAPS/'results'/f'{oid}__b0.03125__primary.npz',allow_pickle=False) as z:
            residual=z['reference_error'];order=np.argsort(residual,kind='stable')
            i=int(order[len(order)//2] if family=='broad' else order[-1])
            why='median_public_midpoint_fit_residual' if family=='broad' else 'largest_public_midpoint_fit_residual'
            nodes.append(dict(side=96,oid=oid,shape=shape,family=family,index=i,reason=why,
                              public_midpoint_fit_residual=float(residual[i]),old_strict_margin=None,old_negative_index=None))
    payload={'scope':'archive-conditioned node comparisons; no truth files accessed',
             'b':.03125,'broad_core':[3.,5.],'broad_shell':[0.,5.],'primary_core':[3.75,4.25],
             'ratio_bounds':[2.**-19,2.**5], 'scale_bounds':[2.**-20,2.**8], 'nodes':nodes}
    (HERE/'freeze.json').write_text(json.dumps(payload,indent=2)+'\n')


if __name__=='__main__':
    register()

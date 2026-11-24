"""Evaluate the complete frozen 16-tile menu for one observed matrix."""

import numpy as np

from ablation import one_sided
from branches import as_branch,expand
from native_box import clock,elapsed


def all_sites():
    for row,col in np.ndindex(4,4):
        mask=np.zeros((4,4))
        mask[row,col]=4.
        yield (row,col),as_branch(f'tile_{row}_{col}',mask)


class SiteMenu:
    def __init__(self,checks):
        self.checks=checks
        self.sites=list(all_sites())
        setup=clock()
        g=self.checks.geometry
        self.unit_born={}
        for coords,branch in self.sites:
            active=expand(branch,g.side)>0
            self.unit_born[coords]=g.s[:,active]@g.b[active]
        self.born_setup=elapsed(setup)

    def score(self,y,radius):
        y=np.asarray(y,complex)
        start=clock()
        born_start=clock()
        born=[]
        for tile,signature in self.unit_born.items():
            norm=np.vdot(signature,signature).real
            alpha=np.vdot(signature,y)/norm if norm>0 else 0j
            misfit=float(np.linalg.norm(y-alpha*signature))
            born.append(dict(tile=tile,fitted_scale=[float(alpha.real),float(alpha.imag)],
                             residual=misfit,relative_residual=misfit/max(np.linalg.norm(y),1e-30)))
        born.sort(key=lambda row:(row['residual'],row['tile']))
        born_cost=elapsed(born_start)
        scores=[]
        for coords,branch in self.sites:
            primary=self.checks.score(y,radius,branch)
            app=self.checks.tile_projection
            ablations=[one_sided(app,y,radius,branch,method)
                       for method in ('receive_only','illumination_only')]
            scores.append(dict(tile=coords,primary=primary,one_sided=ablations))
        count=lambda method:[s['tile'] for s in scores if next(c for c in s['one_sided'] if c['arm']==method)['decision']=='RETAIN']
        return dict(radius=float(radius),tile_scores=scores,
                    survivors=[s['tile'] for s in scores if s['primary']['decision']=='RETAIN'],
                    receive_only_survivors=count('receive_only'),
                    illumination_only_survivors=count('illumination_only'),
                    born_ranking=born,born_comparator_cost=born_cost,**elapsed(start))

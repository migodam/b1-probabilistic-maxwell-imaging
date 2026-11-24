"""Replay the frozen single-profile menu on development only; no truth inputs."""

import json

import numpy as np
from threadpoolctl import threadpool_limits

from branches import primary_menu
from development import SLIT, SHELL, inputs
from native_box import HERE, OPTICS, clock, elapsed
from screen import Checks


def main():
    start=clock()
    with np.load(OPTICS/'synthetic/public/fine_geometry.npz',allow_pickle=False) as p:
        checks=Checks(p['spectrum'],p['PW180_10x20__S'],p['PW180_10x20__B'])
    rows={}
    for case in (SLIT,SHELL):
        t=clock()
        y,r,_,radius=inputs(case)
        archive=HERE.parent/'diagnosis/trajectories'/case/'step_0064.npz'
        with np.load(archive,allow_pickle=False) as p:
            menu=primary_menu(p['m'],p['J'])
        rows[case]=dict(radius_scope=radius,preparation=elapsed(t),branches=[])
        for branch in menu['branches']:
            row=checks.score(y,r,branch)
            row['tile_caps']=branch.caps
            rows[case]['branches'].append(row)
            (HERE/'frozen_development.json').write_text(json.dumps(dict(geometry_setup=checks.setup_cost,
                 cases=rows,**elapsed(start)),indent=2,allow_nan=False)+'\n')
            print(case,branch.name,row['decision'],row['residual_lower'],row['check_time'],flush=True)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

"""Check the cheaper translated-tile proof on the two frozen development views."""

import json

import numpy as np
from threadpoolctl import threadpool_limits

from branches import primary_menu
from development import inputs,SLIT,SHELL
from fast_support import TileProjection
from native_box import HERE,OPTICS,NativeBox,clock,elapsed


def main():
    start=clock()
    with np.load(OPTICS/'synthetic/public/fine_geometry.npz',allow_pickle=False) as p:
        geometry=NativeBox(p['spectrum'],p['PW180_10x20__S'],p['PW180_10x20__B'])
    app=TileProjection(geometry)
    rows={}
    for case in (SLIT,SHELL):
        y,r,_,_=inputs(case)
        with np.load(HERE.parent/'diagnosis/trajectories'/case/'step_0064.npz',allow_pickle=False) as p:
            menu=primary_menu(p['m'],p['J'])
        rows[case]=[]
        for b in menu['branches'][1:3]:
            q=app.score(y,r,b)
            rows[case].append(q)
            print(case,b.name,q['residual_lower'],q['decision'],q['wall_seconds'],flush=True)
    (HERE/'fast_development.json').write_text(json.dumps(dict(cases=rows,**elapsed(start)),indent=2)+'\n')


if __name__=='__main__':
    with threadpool_limits(1):
        main()

"""Receiver-only and illumination-only versions of the translated-tile proof."""

import argparse
import json

import numpy as np
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from branches import expand,primary_menu
from development import SLIT,SHELL,inputs
from fast_support import TileProjection
from native_box import HERE,OPTICS,NativeBox,ball,clock,elapsed,frobenius_squared,gershgorin,herm,lo,star


def one_sided(app,y,radius,branch,arm):
    start=clock()
    app.prepare()
    g=app.geometry
    ids=np.flatnonzero(expand(branch,g.side))
    s,b=g.s[:,ids],g.b[ids]
    if arm=='receive_only':
        _,v=eigh(herm(s@s.conj().T),check_finite=False)
        w=np.eye(g.ns)
        ds=np.sqrt(np.cumsum(np.sum(abs(v.conj().T@s)**2,axis=1)))
        dy=np.sqrt(np.cumsum(np.sum(abs(v.conj().T@y)**2,axis=1)))
        screen=dy-float(app.factor.mid())*ds*np.linalg.norm(b)
        rank=int(np.argmax(screen))+1
        v=v[:,:rank]
    elif arm=='illumination_only':
        _,w=eigh(herm(b.conj().T@b),check_finite=False)
        v=np.eye(g.nr)
        db=np.sqrt(np.cumsum(np.sum(abs(b@w)**2,axis=0)))
        dy=np.sqrt(np.cumsum(np.sum(abs(y@w)**2,axis=0)))
        screen=dy-float(app.factor.mid())*db*np.linalg.norm(s)
        rank=int(np.argmax(screen))+1
        w=w[:,:rank]
    else:
        raise ValueError('Choose one-sided comparator.')
    vb,wb=ball(v),ball(w)
    norm=lambda a:frobenius_squared(a).sqrt()
    value=norm(star(vb)*ball(y)*wb).lower()
    upper=norm(star(vb)*ball(s)).upper()*app.factor*norm(ball(b)*wb).upper()
    vv=gershgorin(star(vb)*vb,upper=True).sqrt()
    ww=gershgorin(star(wb)*wb,upper=True).sqrt()
    lower=max(0.,lo((value-upper)/(vv*ww)))
    return dict(arm=arm,left_modes=v.shape[1],right_modes=w.shape[1],residual_lower=lower,
                radius=radius,decision='PRUNE' if lower>radius else 'RETAIN',**elapsed(start))


def compare(app,y,radius,branch):
    two=app.score(y,radius,branch)
    return dict(branch=branch.name,two_sided=two,
                one_sided=[one_sided(app,y,radius,branch,arm)
                           for arm in ('receive_only','illumination_only')])


def main(output):
    start=clock()
    with np.load(OPTICS/'synthetic/public/fine_geometry.npz',allow_pickle=False) as f:
        g=NativeBox(f['spectrum'],f['PW180_10x20__S'],f['PW180_10x20__B'])
    app=TileProjection(g)
    rows={}
    for obs in (SLIT,SHELL):
        y,r,_,_=inputs(obs)
        with np.load(HERE.parent/'diagnosis/trajectories'/obs/'step_0064.npz',allow_pickle=False) as f:
            menu=primary_menu(f['m'],f['J'])
        rows[obs]=[compare(app,y,r,branch) for branch in menu['branches'][1:3]]
        for q in rows[obs]:
            print(obs,q['branch'],[(r.get('arm','two_sided'),r['residual_lower'],r['decision'])
                  for r in [q['two_sided'],*q['one_sided']]],flush=True)
    output.write_text(json.dumps(dict(cases=rows,full_cost=elapsed(start)),indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=type(HERE),default=HERE/'ablation_development.json')
    with threadpool_limits(1):
        main(parser.parse_args().output)

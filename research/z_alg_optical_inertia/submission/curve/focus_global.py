"""Densify the numerical global control around object 04's 1% local crossing."""

from time import perf_counter, process_time

import numpy as np
from flint import ctx
from threadpoolctl import threadpool_limits

from curve_inputs import HERE, full_ring_inputs
from envelope import FourierEnvelope
from native import Geometry, Native, Weights
from search import Projected, MU, RATIOS, times, save_json

CAP = .75
MU_LOCAL = np.geomspace(0.000035,.0004,77)
RATIO_LOCAL = np.linspace(-.4,.65,106)


def main():
    ctx.prec,ctx.threads=112,2
    start,cpu=perf_counter(),process_time()
    spec,s,b,rows=full_ring_inputs()
    y=next(row for row in rows if row.object_id=='optfresh_04' and row.noise_pct==1)
    projected=Projected(spec,s,b,8)
    proposals=[]
    for ratios, mus in ((RATIOS,MU),(RATIO_LOCAL,MU_LOCAL)):
        for t in ratios:
            for mu,(nu,dd,front) in projected.scalar(CAP,float(t),[y],mus).items():
                if nu < len(y.y.T):
                    proposals.append(dict(t=float(t),mu=mu,nu=nu,d=float(dd[0])))
    proposals.sort(key=lambda r:r['d'],reverse=True)
    print('04 best dense global projected',proposals[:8],flush=True)
    search_cost=times(start,cpu)
    out=dict(scope='Finite numerical control, no universal global-joint claim',
             id=y.id, radius=y.padded_radius, cap=CAP,
             full_ring_search_state_dimension=projected.dimension,
             attempted_profiles=len(RATIOS)*len(MU)+len(MU_LOCAL)*len(RATIO_LOCAL),
             local_mu_grid=MU_LOCAL.tolist(),local_ratio_grid=RATIO_LOCAL.tolist(),
             best_projected_100=proposals[:100],search_cost=search_cost,controls=[])
    geo=Geometry(spec,s,b)
    env=FourierEnvelope(geo)
    choices=[('dense_global_winner',Weights.global_(CAP,r['mu'],r['t']))
             for r in proposals[:3]]
    choices.extend([
        ('uniform_mean_optical',Weights(CAP,.0001,(.0001,)*4,(.000013,)*4)),
        ('unperturbed_global',Weights.global_(CAP,.0001,.1)),
        ('local_optical',Weights(CAP,.0001,(.0001,)*4,(4.000000000000001e-06,1.6000000000000003e-05,
                                                      1.6000000000000003e-05,1.6000000000000003e-05)))])
    for name,w in choices:
        local_start,local_cpu=perf_counter(),process_time()
        actual=Native(geo,w)
        modes=actual.modes(4)
        _,negative_gap=actual.ensure_negative(actual.nu,return_gap=True)
        proof=env.verify(w,actual.nu)
        trial,v,vals,numeric=actual.numerical_distance(y.y,actual.nu)
        strict=actual.ensure_source(y.y,actual.nu,proof['positive_gap_lower'],trial,v,vals,
                                    negative_gap=negative_gap)
        out['controls'].append(dict(name=name,profile=w.as_dict(),nu=actual.nu,
                                    native_numeric=numeric,inertia_proof=proof,
                                    two_sided_source=strict,**times(local_start,local_cpu)))
        out.update(times(start,cpu))
        print(name,strict['distance_lower'],strict['distance_upper'],flush=True)
        save_json(HERE/'focus_global.json',out)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

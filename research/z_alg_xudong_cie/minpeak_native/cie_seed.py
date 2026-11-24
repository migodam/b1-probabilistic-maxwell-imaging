"""Run inherited Dual-CIE on the saved native operators to initialize box searches."""

import argparse
from collections import Counter
import json
import sys

import numpy as np
from threadpoolctl import threadpool_limits

from search import HERE, Model, full_ring_inputs, clocks, elapsed

RESEARCH = HERE.parents[1]
sys.path[:0] = [str(RESEARCH/'srap_cie_array_study_20261121'),
                str(RESEARCH/'delegated/max_struct_r1/r13_worktree')]
from engine import legacy, run_cie_arm
from research.research_z_alg_r4_diversity.src.physics import Regularizer


class Green:
    def __init__(self, native):
        self.native = native

    def apply(self, x):
        return self.native.green(x)

    def adjoint(self, x):
        return self.native.green(x,adjoint=True)


class Adapter:
    def __init__(self, spectrum, sensor, incident, measurement, limit):
        self.clock=legacy.Clock(limit)
        self.E=np.asarray(incident)
        self.n,self.ns=incident.shape
        self.shape=(96,96)
        self.ports=legacy.Ports(sensor)
        self.native=Model(spectrum,sensor,incident,measurement)
        self.green=Green(self.native)
        self.reg=Regularizer(self.shape)


def one(spectrum,sensor,incident,obs,args):
    start=clocks()
    model=Adapter(spectrum,sensor,incident,obs.y,args.wall)
    initial=np.zeros(model.n)
    best_residual=float('inf')
    best=None
    audits=[]
    case=dict(public_schedule=dict(modes=[5,8,12],betas=[8.,2.,.5],rank=args.rank))
    status='completed inherited CIE schedule'
    for step_i in [0]:
        residual, ports, _=model.native.full(initial,gradient=False)
        best_residual, best = residual,initial.copy()
        audits.append(dict(step=step_i,residual=float(residual),peak=0.,cap=args.clip))
    try:
        for mat,step,stage,at_event in run_cie_arm(model,obs.y,case,'dual',initial,obs.padded_radius):
            if step%25 and not at_event:
                continue
            material=np.minimum(mat,args.clip)
            residual,ports,_=model.native.full(material,gradient=False)
            if residual<best_residual:
                best_residual= residual
                best=material.copy()
            audits.append(dict(step=step,stage=stage,residual=float(residual),peak=float(material.max()),
                               raw_peak=float(mat.max()),cap=args.clip))
            print(obs.id,audits[-1],flush=True)
    except legacy.TimeLimit:
        status='inherited CIE wall limit'
    np.savez_compressed(args.out/(obs.id+'.npz'),m_seed=best)
    (args.out/(obs.id+'.json')).write_text(json.dumps({
        'observation':obs.id,'rank':args.rank,'schedule':case['public_schedule'],
        'status':status,'native_audits':audits,'seed_residual':float(best_residual),
        'native_counts':dict(model.native.count),'cie_counts':dict(model.clock.counts),
        **elapsed(start)},indent=2,allow_nan=False)+'\n')
    print('DONE',obs.id,best_residual,status,flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=type(HERE),default=HERE/'cie')
    parser.add_argument('--object',action='append')
    parser.add_argument('--noise',type=float,action='append')
    parser.add_argument('--wall',type=float,default=50.)
    parser.add_argument('--clip',type=float,default=4.)
    parser.add_argument('--rank',type=int,default=20)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    s,p,b,ob=full_ring_inputs()
    for obs in ob:
        if obs.noise_pct in (args.noise or [1.,2.]) and (not args.object or obs.object_id in args.object):
            one(s,p,b,obs,args)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

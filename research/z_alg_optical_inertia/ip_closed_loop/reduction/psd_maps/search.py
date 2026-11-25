"""Continuous, all-source positive-state control of the exact weak-background map boxes."""
import argparse
import json
import os
from pathlib import Path
import sys

for variable in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'

import numpy as np
from scipy.linalg import eigvalsh
from scipy.optimize import minimize_scalar
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
MAPS = HERE.parents[1]/'maps'
sys.path.insert(0,str(MAPS))
from common import Geometry, clock, elapsed, geom_path, pattern, place, save_json
from sector_map import SectorGeometry, hermitian, local_certificate, margin_distances, response_reduction


def read_spec():
    return json.loads((HERE/'freeze.json').read_text())


def box(node):
    side,shape,idx = node['side'],node['shape'],node['index']
    offset,prof = pattern(side,shape)
    y,x = divmod(idx,side)
    ids,profile = place(side,(x,y),offset,prof)
    if node['family']=='broad':
        lo = np.where(profile==1,3.,0.)/profile
        hi = 5./profile
    else:
        lo,hi = 3.75,4.25
    return ids,profile,lo,hi


def safe_number(x):
    return float(x) if np.isfinite(x) else None


class Objective:
    def __init__(self, sector,node,weights_cache,scale_bounds):
        self.sector = sector
        self.node = node
        self.a,self.profile,self.lo,self.hi = box(node)
        with np.load(MAPS/'public'/f'{node["oid"]}.npz',allow_pickle=False) as z:
            self.Y,self.r = z['Y'],float(z['r'])
        self.scale_bounds = scale_bounds
        self.weights_cache = weights_cache
        self.local = {}
        self.trials = []
        self.cache = {}
        # Fixed geometric clipping is enough for the electrostatic covariance cache;
        # the measured border is formed independently at each cut point.
        xy = sector.g.xy[self.a]
        self.cov_key = (node['shape'],node['family'],tuple(map(tuple,xy-xy[0])))

    def prepare(self,ratio):
        if ratio in self.local:
            return self.local[ratio]
        start = clock()
        key = (self.cov_key,ratio)
        if key not in self.weights_cache:
            self.weights_cache[key] = local_certificate(self.sector,self.a,self.profile,self.lo,self.hi,ratio)
        p = self.weights_cache[key]
        prep = None if p is None else response_reduction(self.sector,self.a[None,:],p)
        self.local[ratio] = prep
        self.trials.append(dict(stage='ratio_prepare',ratio=ratio,physical_negative_index=(None if p is None else p[-1]),
                                resolved=p is not None,cost=elapsed(start)))
        return prep

    def evaluate(self,ratio,scale,stage, positive_only=True):
        key = ratio,scale,positive_only
        if key in self.cache:
            return self.cache[key]
        started=clock()
        data = self.prepare(ratio)
        raw,bound,raw_nu,verified_nu = None,None,None,None
        if data is not None:
            P,yp,theta,*_,nu = data
            port=hermitian(np.eye(self.Y.shape[0])+P[0]/scale)
            vals,Q = np.linalg.eigh(port)
            raw_nu = nu-int(np.count_nonzero(vals<0))
            if 0<=raw_nu<self.Y.shape[1] and np.min(abs(vals))>1.e-13:
                if raw_nu==0 or not positive_only:
                    e=Q.conj().T@(self.Y-yp[0])
                    h=hermitian(e.conj().T@(e/vals[:,None])-scale*theta[0])
                    eig=eigvalsh(h,check_finite=False)[:self.Y.shape[1]-raw_nu]
                    raw=float(np.sqrt(np.maximum(eig,0).sum())-self.r)
                    expected_margin,verified,_ = margin_distances(data,self.Y,self.r,[scale])
                    if np.isfinite(expected_margin[0]) and (verified[0]==0 or not positive_only):
                        bound=float(expected_margin[0]); verified_nu=int(verified[0])
        rec=dict(stage=stage,ratio=ratio,scale=scale,raw_positive_part_margin=raw,
                 verified_margin=bound,raw_state_negative_index=raw_nu,
                 verified_state_negative_index=verified_nu,positive_state_required=positive_only,cost=elapsed(started))
        self.cache[key]=rec
        self.trials.append(rec)
        return rec

    def best_scale(self,ratio):
        data=self.prepare(ratio)
        if data is None:
            return None
        P,_,_,pe,_,_,nu=data
        poles=-eigvalsh(hermitian(P[0]),check_finite=False)
        rp=len(poles)
        if nu>rp:
            return None
        padding=8.*float(pe[0])+1.e-12
        lower=max(self.scale_bounds[0],poles[nu]+padding) if nu<rp else self.scale_bounds[0]
        upper=min(self.scale_bounds[1],poles[nu-1]-padding) if nu>0 else self.scale_bounds[1]
        if lower>=upper:
            self.trials.append(dict(stage='scale_interval_empty',ratio=ratio,numeric_scale_interval=[lower,upper]))
            return None
        logs=np.linspace(np.log2(lower),np.log2(upper),9)
        def sample(z):
            e=self.evaluate(ratio,float(2.**z),'continuous_scale')
            val=e['raw_positive_part_margin']
            return -val if val is not None else 1.e2
        values=np.array([sample(z) for z in logs])
        k=int(np.argmin(values)); left,right=logs[max(0,k-1)],logs[min(len(logs)-1,k+1)]
        ans=minimize_scalar(sample,bounds=(left,right),method='bounded',options={'maxiter':12,'xatol':.005})
        self.trials.append(dict(stage='scale_optimizer_stop',ratio=ratio,iterations=ans.nfev,
                                returned_log2_scale=float(ans.x),success=bool(ans.success)))
        entries=[z for (w,m,p),z in self.cache.items() if w==ratio and p and z['raw_positive_part_margin'] is not None]
        return max(entries,key=lambda z:z['raw_positive_part_margin']) if entries else None

    def optimize(self,ratio_bounds,old_weight_bank,old_scale_bank):
        started=clock()
        bank=[self.evaluate(r,s,'fixed_positive_bank') for r in old_weight_bank for s in old_scale_bank]
        comparable_unrestricted=[self.evaluate(r,s,'fixed_unrestricted_bank',positive_only=False)
                                 for r in old_weight_bank for s in old_scale_bank]
        bank_time=elapsed(started)
        st=clock()
        low,high=map(np.log2,ratio_bounds)
        points=np.linspace(low,high,13)
        def search(log):
            p=self.best_scale(float(2.**log))
            return -p['raw_positive_part_margin'] if p is not None else 1.e2
        values=np.array([search(point) for point in points]); k=int(np.argmin(values))
        left,right=points[max(0,k-1)],points[min(len(points)-1,k+1)]
        ans=minimize_scalar(search,bounds=(left,right),method='bounded',options={'maxiter':14,'xatol':.04})
        self.trials.append(dict(stage='ratio_optimizer_stop',evaluations=ans.nfev,
                                returned_log2_ratio=float(ans.x),success=bool(ans.success)))
        def best(values,tag):
            valid=[v for v in values if v[tag] is not None]
            return max(valid,key=lambda v:v[tag]) if valid else None
        all_candidates=[e for (_,_,p),e in self.cache.items() if p]
        return dict(positive_fixed_verified_best=best(bank,'verified_margin'),
                    unrestricted_fixed_verified_best=best(comparable_unrestricted,'verified_margin'),
                    positive_continuous_verified_best=best(all_candidates,'verified_margin'),
                    positive_continuous_raw_optimum=best(all_candidates,'raw_positive_part_margin'),
                    original_numeric_stop_is_global=False,
                    source_columns=self.Y.shape[1],noise_radius=self.r,
                    fixed_bank_cost=bank_time,additional_continuous_cost=elapsed(st),
                    trials=self.trials)


def run(side):
    spec=read_spec(); start=clock()
    g=Geometry(geom_path(side)); sector=SectorGeometry(g,spec['b'])
    preparation=elapsed(start)
    cache={}; result=[]
    weight_bank=(.125,.5,1.,2.,8.) if side==32 else (.125,.03125,.5)
    scale_bank=tuple(float(2.**x) for x in range(-12,5)) if side==32 else (.5,1.,.25,2.,.125,4.,.0625,.03125,.015625)
    for node in (n for n in spec['nodes'] if n['side']==side):
        st=clock(); o=Objective(sector,node,cache,spec['scale_bounds'])
        found=o.optimize(spec['ratio_bounds'],weight_bank,scale_bank)
        row={**node,**found,'whole_node_cost':elapsed(st),'physical_cells':o.a.tolist()}
        result.append(row)
        save_json(HERE/f'result_{side}.json',{'geometry_prep':preparation,'weights_bank':weight_bank,
                                            'scale_bank':scale_bank,'elapsed':elapsed(start),'nodes':result})
        print(side,node['oid'],node['index'],node['family'],
              [None if found[k] is None else (found[k]['verified_margin'],found[k]['ratio'],found[k]['scale']) for k in
               ('positive_fixed_verified_best','unrestricted_fixed_verified_best','positive_continuous_verified_best')],
              row['whole_node_cost'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('side',type=int,choices=(32,96)); args=p.parse_args()
    with threadpool_limits(limits=1):run(args.side)

"""Resolve narrow positive-state intervals missed by the first broad ratio grid."""
import numpy as np
from scipy.optimize import minimize_scalar

from search import (HERE,Geometry,Objective,SectorGeometry,clock,elapsed,geom_path,read_spec,
                    save_json,threadpool_limits)


def run():
    spec=read_spec()
    start=clock()
    g=Geometry(geom_path(32)); sector=SectorGeometry(g,spec['b'])
    prep=elapsed(start)
    cutouts={}; results=[]
    for node in (r for r in spec['nodes'] if r['side']==32 and r['family']=='broad'):
        started=clock()
        o=Objective(sector,node,cutouts,spec['scale_bounds'])
        log_ratios=np.linspace(-13.,0.,105)
        def trial(z):
            ret=o.best_scale(float(2.**z))
            return -ret['raw_positive_part_margin'] if ret is not None else 1.e2
        values=np.array([trial(x) for x in log_ratios])
        if np.min(values)<1.e2:
            i=int(np.argmin(values)); left,right=log_ratios[max(0,i-1)],log_ratios[min(len(log_ratios)-1,i+1)]
            solve=minimize_scalar(trial,bounds=(left,right),method='bounded',options={'maxiter':14,'xatol':.01})
            stop={'raw_stop_success':bool(solve.success),'reported_log2_ratio':float(solve.x),'extra_proposals':solve.nfev}
        else:
            stop={'reason':'no positive_state_scale_interval'}
        entries=list(o.cache.values())
        raw=[x for x in entries if x['raw_positive_part_margin'] is not None]
        exact=[x for x in entries if x['verified_margin'] is not None]
        found={**node,'positive_continuous_raw_optimum':max(raw,key=lambda q:q['raw_positive_part_margin']) if raw else None,
               'positive_continuous_verified_best':max(exact,key=lambda q:q['verified_margin']) if exact else None,
               'evaluated_ratio_count':len(o.local),'scale_optimizer_stop':stop,'whole_node_cost':elapsed(started),'trials':o.trials}
        results.append(found)
        print(node['index'],found['positive_continuous_raw_optimum'],found['positive_continuous_verified_best'],flush=True)
        save_json(HERE/'refined_broad_32.json',{'setup_cost':prep,'scope':'32 squared; exploratory finer ratio grid near PSD boundaries',
                                             'elapsed':elapsed(start),'nodes':results})


if __name__=='__main__':
    with threadpool_limits(limits=1):run()

"""Resolve the gap between earlier tested caps, using the unchanged rank bound.
Each source extrapolation is checked on the native complex data before reporting.
"""
from pathlib import Path
from time import perf_counter,process_time
import sys, argparse
from threadpoolctl import threadpool_limits
from flint import ctx
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[2]/'physical_revision/array'))
from screen import inputs,frozen_profiles,grid_profile,save_json
from verify_native import Geometry,FourierEnvelope,check_profile
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--transition',action='store_true',help='Resolve intermediate caps using two already published weights.')
    args=parser.parse_args()
    t,c=perf_counter(),process_time()
    spec,ss,bs,obs=inputs()
    g={2:Geometry(spec,ss[32],bs[32])}
    env={2:FourierEnvelope(g[2])}
    picked={32:[(x,y) for x,y in obs[32][32] if x['noise_pct']==2]}
    ctx.prec,ctx.threads=96,1
    base=next(frozen_profiles())
    results=[]
    checks=[(base,cap) for cap in [.78125,.8125,.875,.9375,3.]]
    if args.transition:
        checks=[(profile,cap) for cap in [.84375,.859375,.875,.90625] for profile in [base,grid_profile(601)]]
    for profile,cap in checks:
        t1,c1=perf_counter(),process_time()
        try:
            row=check_profile(profile,cap,g,env,ss[32],{32:bs[32]},picked,None, max(8,int(cap*25-8)))
        except (ArithmeticError,RuntimeError,ValueError) as ex:
            row=dict(cap=cap,error=str(ex),wall_seconds=perf_counter()-t1,cpu_seconds=process_time()-c1)
        results.append(row)
        output='cap_resolution_transition.json' if args.transition else 'cap_resolution.json'
        save_json(HERE/output,dict(profile=base,rows=results,wall_seconds=perf_counter()-t,cpu_seconds=process_time()-c,
            purpose='Cap sampling resolution for the same rank-only method, unchanged six native data balls.'))
        print(cap,[(x['object_id'],round(x.get('distance_lower',0)/x['radius'],3)) for x in row.get('observations',[])],flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()

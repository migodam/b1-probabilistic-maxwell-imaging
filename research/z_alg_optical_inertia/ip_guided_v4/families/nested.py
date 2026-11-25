"""Transport a separating sector from [0,6] to its three nonconcentric subintervals."""
import json
import numpy as np
from run import (HERE, SITES, WEIGHTS, BACKGROUNDS, clock, elapsed, inputs, ball,
                 ctx, kernel_from_symbol, star, arrays, patch, get_profile, trace_lower,
                 low, arb, threadpool_limits, save)

def main():
    ctx.prec=128
    original=json.loads((HERE/'results.json').read_text())
    tick=clock()
    spec,sensors,sources,data=inputs()
    s,b=sensors[32],sources[32]
    rows=[(z['id'],y,z['padded_radius']) for z,y in data[32][32] if z['noise_pct']==2]
    for z in original['controls']:
        a=np.load(HERE/(z['id']+'.npz'));rows.append((z['id'],a['Y'],z['radius']))
    sym=kernel_from_symbol(spec);v=ball(np.column_stack((s.conj().T,b)));gram=star(v)*v
    out={'method':'transported parent sector weight h_parent/h_child, matching unchanged exterior','runs':[],'errors':[]}
    for j,site in enumerate(SITES):
        g,_,bi=arrays(sym,patch(site),(s,b),gram,v)
        for bg in BACKGROUNDS:
            for mu in WEIGHTS:
                for half in (2,1,.5):
                    start=clock()
                    try:
                        form,nu,details=get_profile(g,bi,bg,4,half,mu,float(spec.real.max()),32,32,3/half)
                        endpoints={key:low(arb(trace_lower(form,nu,ball(y)))/arb(float(r))) for key,y,r in rows}
                        out['runs'].append(dict(site=j+1,b=bg,mu=mu,center=4,width=half,nu=nu,results=endpoints,facts=details,time=elapsed(start)))
                    except ArithmeticError as err:
                        out['errors'].append(dict(site=j+1,b=bg,mu=mu,width=half,error=str(err)))
        out['total']=elapsed(tick); save('nested.json',out)
if __name__=='__main__':
    with threadpool_limits(1):main()

import json
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from material_moments import bvn
records=[]
for rho in [-.9999999,-.99,-.5,0,.5,.99,.9999999]:
 for t,u in [(0,0),(.5,-.4),(4.,4.),(-4.,-3.),(.001,.001001)]:
  angle=np.arcsin(rho)
  fun=lambda theta:np.exp(-(t*t-2*np.sin(theta)*t*u+u*u)/(2*np.cos(theta)**2))/(2*np.pi)
  reference=ndtr(t)*ndtr(u)+quad(fun,0,angle,epsabs=3e-13,epsrel=3e-13,limit=300)[0];value=float(bvn(t,u,rho));records.append(dict(rho=rho,t=t,u=u,value=value,reference=reference,abs_error=abs(value-reference)))
maxerr=max(r['abs_error'] for r in records);out=dict(status='PASS' if maxerr<1e-9 else 'FAIL',max_abs_error=maxerr,records=records);Path(__file__).with_name('JOINT_CDF_VALIDATION.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],maxerr)

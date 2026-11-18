"""Cost of the conditional bandwidth choice, including the physical coefficient evaluations."""
from time import perf_counter,process_time
import json
import numpy as np
from scipy.optimize import brentq
from scipy.linalg import null_space
from freq_model import FrequencyTM,ri
from mirror_spectrum import Mirror,ROOT
from risk_design import conditional_risk

if __name__=='__main__':
    records=[]
    for case in json.loads((ROOT/'results/risk_design.json').read_text()):
        start=perf_counter(); cpu=process_time()
        model=Mirror([1,2]) if case['label']=='mirror' else FrequencyTM(case['geometry'])
        m=np.array(case['center']);v=np.array(case['v']);W=null_space(v[None,:])
        y,J,H=model.response(m,hessian=True)
        N=np.linalg.svd(ri(J)@W,full_matrices=True)[0][:,W.shape[1]:]
        c=N.T@ri(np.einsum('kij,i,j->k',H,v,v))/2
        d=ri(model.response(m,derivative=True)[3]@v)
        a=case['amplitude'];sigma=case['unit_power_sigma_real']/a**2
        theta=np.linalg.norm(c)/np.linalg.norm(d)
        def ev(t):
            return conditional_risk(c,d,N,0,(theta*t)**2/2,1,sigma)
        root=brentq(lambda t:ev(t)-.1,.5,4.,full_output=True)
        records.append(dict(label=case['label'],band=a*theta*root[0],LU=model.factorizations,
                            risk_integrals=root[1].function_calls,cpu_s=process_time()-cpu,wall_s=perf_counter()-start,
                            assumption='fold center, kernel and half separation already available'))
    (ROOT/'results/risk_cost.json').write_text(json.dumps(records,indent=2)+'\n')

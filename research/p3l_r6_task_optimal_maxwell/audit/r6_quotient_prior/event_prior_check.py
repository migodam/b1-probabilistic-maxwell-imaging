import importlib.util,json
from pathlib import Path
from scipy.integrate import quad
from scipy.stats import norm
import numpy as np
HERE=Path(__file__).resolve().parent;sp=importlib.util.spec_from_file_location('a',HERE/'audit.py');a=importlib.util.module_from_spec(sp);sp.loader.exec_module(a);ns=a.ns;A=ns['A'];Z0=ns['Z0'];threshold=.02
z1=float(ns['inverse_fraction'](threshold/ns['RHI']));z2=float(ns['inverse_fraction'](threshold/ns['RLO']))
def f(z):
 q=float(a.saturation(1.3*z));p=0 if q<threshold/ns['RHI'] else np.clip((ns['RHI']-threshold/q)/(ns['RHI']-ns['RLO']),0,1)
 return norm.pdf(z)/Z0*p
exact=ns['PH']+quad(f,-A,A,points=[z1,z2],epsabs=1e-12,epsrel=1e-12)[0];records=[]
for name,make in [('S1_12_12',lambda:ns['source_nodes'](2,12,12)),('S1_32_24',lambda:ns['source_nodes'](2,32,24)),('S3_10_8',lambda:ns['source_nodes'](3,10,8)),('S3_16_9',lambda:ns['source_nodes'](3,16,9)),('quotient12_6',lambda:ns['quotient_nodes'](12,6,48)),('quotient20_8',lambda:ns['quotient_nodes'](20,8,64)),('quotient64_32',lambda:ns['quotient_nodes'](64,32,128))]:
 n=make();val=float(np.sum(n['prior_weights']*(n['V'].real>threshold)));records.append(dict(method=name,event_prior_numerator=val,error=val-exact))
(HERE/'MATERIAL_EVENT_PRIOR_QA.json').write_text(json.dumps(dict(exact_event_prior_probability=exact,threshold=threshold,records=records,physics_RHS=0,MC_draws=0,GT_read=False),indent=2)+'\n')
print(exact,records)

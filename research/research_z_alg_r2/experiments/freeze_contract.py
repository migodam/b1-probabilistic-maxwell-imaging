"""Freeze public geometry and observation-only interface before any hidden draw."""
import json, math, datetime
from pathlib import Path
P=Path(__file__).resolve().parent
c={'schema':'Z-ALG-R2-BLIND-v1','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'materials':{'real_box':[.5,4.],'imaginary':0,'initialization':'fixed box midpoint 2.25, all cells','hidden_distribution':'independent uniform on [0.5,4.0], three targets per preset'},'physics':{'k':2*math.pi,'radius':.04,'noise':'none','basis':'constant disk pulse Green VIE'},'input_contract':{'allowed':['preset geometry','material box','complex observations y','case identifier','declared noise model'],'forbidden':['hidden materials','generator seed','hidden current or field','GT-derived initialization or support choices'],'packing':'source major, receiver major, real then imaginary'},'baseline':{'method':'scipy.optimize.least_squares TRF exact analytic full-J adjoint','bounds':[.5,4.],'max_nfev':120,'ftol':1e-12,'xtol':1e-12,'gtol':1e-12,'loss':'linear','x_scale':1.0,'initial':2.25,'frozen_J':'same TRF settings, one midpoint Jacobian'},'validation':{'independent':'fresh numpy.linalg.solve at every returned endpoint','response_relative_max':1e-8,'state_relative_max':1e-10,'material_RMS_success':1e-5,'legal_box_tolerance':1e-10,'resolvent_sigma_min_floor':1e-6,'nonunique_policy':'E1 response-fit and material error separate; never infer uniqueness'},'budget':{'max_nfev_per_method_case':120,'actual_complex_rhs_per_method_case':15000,'wall_seconds_per_method_case':120,'single_thread':True},'presets':{}}
c['presets']['E1_LOW6']={'centers':[[-.29,-.15],[.04,-.23],[.33,-.08],[-.22,.22],[.09,.13],[.38,.27]],'sources':[[-.61,-.72],[.63,-.57]],'receivers':[[.17,.86]]}
c['presets']['INFO12']={'centers':[[x,y] for y in [-.24,0.,.24] for x in [-.30,-.10,.10,.30]],'sources':[[-.65,-.70],[.72,-.48],[-.67,.42]],'receivers':[[.95*math.cos(2*math.pi*j/16+.071),.95*math.sin(2*math.pi*j/16+.071)] for j in range(16)]}
p=P/'contract.json'
if p.exists(): raise RuntimeError('Frozen contract already exists')
p.write_text(json.dumps(c,indent=2)+'\n')
print(p)

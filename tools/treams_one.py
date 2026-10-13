"""Reuse A15's plain helicity multilayer coefficient map in isolated ABI.
Caller MUST already reserve a B1 physical evaluation before launching this.
"""
import sys,json,os
from pathlib import Path
VENDOR=Path(os.environ.get("B1_TREAMS_VENDOR",str(Path(__file__).resolve().parents[3]/"A6/vendor")))
if 'numpy' in sys.modules or 'scipy' in sys.modules:raise RuntimeError('fresh subprocess required')
sys.path.insert(0,str(VENDOR))
import numpy as np
import scipy
import treams.coeffs
p=json.loads(sys.argv[1]);rr=np.array(p['radii']);ee=np.array([complex(*z) for z in p['eps_layers']]+[1+0j]);a=[];b=[]
# Unmodified coefficient-map lines from A15 run_reference_checks.py:
for l in range(1,p['nmax']+1):
 h=np.asarray(treams.coeffs.mie(l,p['k']*rr,ee,np.ones(len(ee)),np.zeros(len(ee))))
 if not np.isfinite(h).all() or np.linalg.norm(h-h.T)>1e-12 or abs(h[0,0]-h[1,1])>1e-12:raise ValueError('helicity symmetry')
 a.append(-(h[0,0]+h[0,1]));b.append(-(h[0,0]-h[0,1]))
print(json.dumps(dict(a=[[z.real,z.imag] for z in a],b=[[z.real,z.imag] for z in b],environment=dict(numpy=np.__version__,scipy=scipy.__version__,treams=str(treams.__file__)))))

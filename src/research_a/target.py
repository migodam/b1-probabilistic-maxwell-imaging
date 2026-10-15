import json,time,os
from pathlib import Path
import numpy as np
from scipy.stats import chi2
from b1.physics import Forward
from b1_joint.material import from_geometry,decode,geometry,RMAX,image,prior_sample
from b1_joint.response import Response

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'research/research_a_20261015'
HISTORY=Path(os.environ.get('B1_P2_RESULTS',str(ROOT/'research/joint_p2/results')))
if not (HISTORY/'COMMON_RESPONSE_BANK.npz').exists():HISTORY=ROOT.parent/'b1-joint-p2/research/joint_p2/results'
BOUNDS=np.array([[0,RMAX[0]],[0,RMAX[1]],[2.2,2.8],[.04,.16]])
TAU=float(np.sqrt(chi2.ppf(.95,32)))

def to_x(p):return from_geometry(*p)
def to_p(x):return np.r_[geometry(x),np.asarray(x)[2:]]
def profile(f,y,sigma):
    f=np.asarray(f);w=np.asarray(sigma)**-2
    g=np.clip(np.sum(f*y*w,axis=-1)/np.sum(f*f*w,axis=-1),.7,1.3)
    residual=(g[...,None]*f-y)/sigma
    return g,residual

class Target:
    def __init__(self,budget):
        self.budget=budget
        config=json.loads((ROOT/'configs/protocol.json').read_text());config['nmax']=14
        self.solver=Forward(config,budget,OUT/'results/new_forward_cache')
        self.proxy=Response(HISTORY/'COMMON_RESPONSE_BANK.npz')
        q=json.loads((HISTORY/'PHYSICS_QUALIFICATION.json').read_text())
        self.sigma=np.repeat(q['noise_sigma_per_frequency'],32)
    def forward(self,p,k=1.6,kind='ENDPOINT',cache=False,nmax=14,implementation='matching'):
        x=to_x(p);d=decode(x)
        return self.solver.run(-1,p[2],p[3],radii=d['radii'],eps_layers=d['eps_layers'],k=k,kind=kind,nmax=nmax,implementation=implementation,cache=cache)
    def predict_proxy(self,p):
        p=np.atleast_2d(p);return self.proxy.predict(np.array([to_x(z) for z in p]))
    def validate_material(self,p):
        d=decode(to_x(p));return dict(passive=bool(np.all(d['eps_layers'].imag>=0)),shared_across_sources=True,constitutive_relative_defect=None,reason='analytic concentric isotropic material is shared by both polarizations; current-field constitutive defect not exposed by angular Mie interface')
    def materialize(self,p,n=151):return image(to_x(p),n)
    def physics_identity(self):
        return dict(id='Research-A-original-P2-vector-Mie-n14-C',support='original four continuous coefficients',source_commit='7ed4fae',frequency=1.6,packing='Re16 then Im16',noise_sigma=self.sigma[:32].tolist(),gain='real bounded .7..1.3',tau=TAU)
    def jvp(self,p,direction,step=1e-5):
        """Paid central full finite difference, not an analytic/free adjoint."""
        p=np.asarray(p);d=np.asarray(direction)
        if np.any(p-step*d<BOUNDS[:,0]) or np.any(p+step*d>BOUNDS[:,1]):raise ValueError('central JVP outside physical support')
        a=self.forward(p+step*d,kind='JVP_PLUS')[0];b=self.forward(p-step*d,kind='JVP_MINUS')[0]
        return (a-b)/(2*step)
    def vjp(self,p,data_vector,step=1e-5):
        """Explicit eight-paid-solve transpose; no adjoint speed claim."""
        return np.array([np.dot(self.jvp(p,np.eye(4)[j],step),data_vector) for j in range(4)])

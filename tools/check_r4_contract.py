"""Necessary cache/algebra audit: no physical LU or solver RHS."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3r4.core import OUT,Book,write,pack
from p1_material.prior import MaterialPrior,saturation
import numpy as np
from scipy.special import ndtr,ndtri
from scipy.integrate import quad
from scipy.stats import kstest

def main():
 b=Book('0','probability_and_cached_operator_contract_r1',max_wall=300,analysis_only=True)
 try:
  cfg=json.loads((OUT/'configs/TARGET.json').read_text());pts=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];p=MaterialPrior(cfg,pts)
  old=ROOT.parent/'b1-p3l-r2-risk-posterior/research/p3l_r2_risk_posterior_exploration'
  with np.load(OUT/'inputs/background/weak.npz') as d:G=d['G'];S=d['S'];incident=d['incident']
  with np.load(old/'A/anchors/DEV_010.npz') as d:a={k:d[k] for k in d.files}
  chi=p.decode(a['q'])[0];L=np.eye(576)-chi[:,None]*G;rhs=chi[:,None]*incident
  err=dict(material=float(np.max(abs(chi-a['chi']))),state=float(np.linalg.norm(L@a['current']-rhs)/np.linalg.norm(rhs)),
     field=float(np.linalg.norm(incident+G@a['current']-a['field'])/np.linalg.norm(a['field'])),
     prediction=float(np.linalg.norm(S@a['current']-a['prediction'])/np.linalg.norm(a['prediction'])),
     receiver=float(np.linalg.norm(a['receiver_T']@L-S)/np.linalg.norm(S)))
  if max(err.values())>1e-10:raise ValueError('cached native operator consistency')
  rng=np.random.default_rng(2026102500);r=rng.normal(size=(32,8))+1j*rng.normal(size=(32,8));white=pack(r)*np.sqrt(2)/cfg['noise']['sigma']
  noise_error=float(abs(.5*white@white-np.vdot(r,r).real/cfg['noise']['sigma']**2)/(.5*white@white))
  marker=np.arange(256).reshape(32,8)+1j*(1000+np.arange(256).reshape(32,8))
  expected=np.concatenate([np.r_[marker[:,s].real,marker[:,s].imag] for s in range(8)])
  assert np.array_equal(pack(marker),expected)
  gauss=rng.normal(size=(100000,2));u=ndtr(gauss)
  jacobian_error=float(np.max(abs(-np.log(np.exp(-gauss*gauss/2)/np.sqrt(2*np.pi)).sum(1)-((gauss*gauss).sum(1)/2+np.log(2*np.pi)))))
  ks=[dict(statistic=float(kstest(u[:,i],'uniform').statistic),pvalue=float(kstest(u[:,i],'uniform').pvalue)) for i in range(2)]
  layer=[]
  for s in np.linspace(-1.2,1.2,49):
   val=quad(lambda t:15*(1-t*t)**2/16,-1,min(1,max(-1,s)),epsabs=1e-13)[0]
   layer.append(abs(float(saturation(s))-val))
  # Decoder-only interior derivative test, not a modified physical endpoint.
  # The exposed cache has a property at a native box boundary; an unconstrained
  # central difference there is correctly rejected by the prior.
  q=a['q'].copy();q[64:]=.5
  v=rng.normal(size=66);v[64:]*=.01;steps=[1e-3,1e-4,1e-5]
  fd=[dict(h=h,relative=float(np.linalg.norm((p.decode(q+h*v)[0]-p.decode(q-h*v)[0])/(2*h)-p.derivative(q)@v)/np.linalg.norm(p.derivative(q)@v))) for h in steps]
  assert noise_error<1e-14 and max(layer)<1e-13 and fd[-1]['relative']<1e-6
  write(OUT/'AUDIT/STAGE0_CONTRACT.json',dict(status='PASS',new_physics_RHS=0,new_physics_LU=0,
   cached_anchor='DEV_010 EXPOSED',cached_errors=err,source_major_pack=True,whitening_error=noise_error,
   probit_prior_Jacobian_error=jacobian_error,uniform_pushforward_KS=ks,layer_cake_max_error=max(layer),decoder_FD=fd,
   limited_scope='cached Maxwell equations and decoder finite differences only; new physical derivative checks charged at first paid job'))
  b.close();print('STAGE0_CONTRACT PASS; zero new physics')
 except BaseException as e:b.close('FAILED',repr(e));raise

if __name__=='__main__':main()

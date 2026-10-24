"""Decoder-only calculus and mock nonlinear-constraint control; zero Maxwell."""
import unittest,json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from p1_material.prior import MaterialPrior
from p3r4.core import to_probit,from_probit
from p3r4.representative import material_objective_gradient,solve_representative,LocalBudgetStop

ROOT=Path(__file__).resolve().parents[2]

class Budget:
 def __init__(self,states=16,rhs=256):self.maxstates=states;self.maxrhs=rhs;self.states=0;self.rhs=0
 def require(self,*,states,RHS):
  if self.states+states>self.maxstates or self.rhs+RHS>self.maxrhs:raise LocalBudgetStop('mock hard budget')
 def used(self):return dict(states=self.states,RHS=self.rhs,LU=self.states)

class MockPrior:
 def legal(self,q):return bool(np.shape(q)==(66,) and np.isfinite(q).all() and np.all(np.asarray(q)[64:]>=0) and np.all(np.asarray(q)[64:]<=1))
 def decode(self,q):return np.array([q[0]+1j*q[1]]),{}
 def derivative(self,q):
  D=np.zeros((1,66),complex);D[0,0]=1;D[0,1]=1j;return D

class MockTarget:
 def __init__(self,budget):self.prior=MockPrior();self.budget=budget;self.serial=0;self.trials=[]
 def point(self,q,key):return SimpleNamespace(q=np.asarray(q).copy(),key=key,state=SimpleNamespace(chi=self.prior.decode(q)[0],residual=0.),J=None,T=None,gradient=None,gradient_observation=None)
 def evaluate(self,q,purpose=''):
  self.budget.require(states=1,RHS=8);self.budget.states+=1;self.budget.rhs+=8;self.serial+=1;p=self.point(q,f'mock_{self.serial}');self.trials.append(p);return p
 def phi(self,p,y):return float(p.q[0]**2)
 def phi_gradient(self,p,y):
  if p.gradient is None:
   self.budget.require(states=0,RHS=8);self.budget.rhs+=8;p.gradient=np.zeros(66);p.gradient[0]=2*p.q[0];p.gradient_observation=np.asarray(y).copy()
  return p.gradient

class RepresentativeTests(unittest.TestCase):
 def test_true_decoder_material_gradient_including_probit(self):
  cfg=json.loads((ROOT/'research/p1_material_target_20261018/configs/ONLINE_TARGET.json').read_text());points=json.loads((ROOT/'research/p1_material_target_20261018/inputs/geometry.json').read_text())['points'];prior=MaterialPrior(cfg,points)
  rng=np.random.default_rng(2026102517);q=prior.draw(rng,1)[0];x=to_probit(q);mu=prior.decode(prior.draw(rng,1)[0])[0]
  value,g,_=material_objective_gradient(prior,x,mu)
  for k in (0,3,19,64,65):
   a=x.copy();b=x.copy();a[k]+=1e-6;b[k]-=1e-6
   numeric=(material_objective_gradient(prior,a,mu)[0]-material_objective_gradient(prior,b,mu)[0])/2e-6
   self.assertAlmostEqual(g[k],numeric,delta=2e-8)

 def test_mock_constraint_retains_only_feasible_final(self):
  budget=Budget();target=MockTarget(budget);x0=np.zeros(66);anchor=target.point(from_probit(x0),'anchor');events=[]
  result,point=solve_representative(target,anchor,np.zeros(1),np.array([2.+0j]),x0=x0,require=budget.require,used=budget.used,event=events.append,phi_threshold=1.)
  self.assertLessEqual(result['Phi'],1.);self.assertLess(result['material_objective'],4.);self.assertLessEqual(budget.states,16);self.assertLessEqual(budget.rhs,256)
  self.assertTrue(any(e.get('Phi',0)>1. for e in events if e['event']=='EVALUATED'))
  self.assertFalse(any(point.key==p.key and target.phi(p,None)>1. for p in target.trials))
  self.assertFalse(result['posterior_mean']);self.assertIn('stationarity_inf',result['kkt'])

 def test_mock_budget_fallback_no_rejected_replacement(self):
  budget=Budget(states=2,rhs=40);target=MockTarget(budget);x0=np.zeros(66);anchor=target.point(from_probit(x0),'anchor')
  result,point=solve_representative(target,anchor,np.zeros(1),np.array([2.+0j]),x0=x0,require=budget.require,used=budget.used,phi_threshold=1.)
  self.assertTrue(result['fallback_to_original_anchor']);self.assertEqual(point.key,'anchor');self.assertLessEqual(budget.rhs,40)

 def test_unqualified_anchor_rejected_before_query(self):
  budget=Budget();target=MockTarget(budget);x0=np.zeros(66);x0[0]=2;anchor=target.point(from_probit(x0),'bad')
  with self.assertRaises(ValueError):solve_representative(target,anchor,np.zeros(1),np.array([2.+0j]),x0=x0,require=budget.require,used=budget.used,phi_threshold=1.)
  self.assertEqual(budget.rhs,0)

if __name__=='__main__':unittest.main()

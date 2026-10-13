import json,tempfile,unittest
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from src.b1.budget import Ledger,PhysicsLimit
from src.b1.inference import infer,quadrature_nodes,log_likelihood,top_k,material_image,refinement
from src.b1.physics import topology_layers
C=json.loads(Path('configs/protocol.json').read_text())

class ProbabilityTests(unittest.TestCase):
 def bank(self,n=5):
  th,w=quadrature_nodes(n,C);return dict(H=np.repeat(np.arange(3),len(th)),theta=np.tile(th,(3,1)),weights=np.tile(w,3),F=np.zeros((3*len(th),4)))
 def test_proper_uniform_quadrature(self):
  for n in [3,5,7]:
   th,w=quadrature_nodes(n,C);self.assertAlmostEqual(w.sum(),1,14);np.testing.assert_allclose(w@th,[2.5,.1],atol=1e-14)
 def test_uninformative_restores_prior(self):
  b=self.bank();r=infer(b,np.ones(4),1,[.8,.15,.05]);np.testing.assert_allclose(r['topology_probabilities'],[.8,.15,.05],atol=1e-14)
 def test_topk_not_renormalized(self):
  b=self.bank();r=infer(b,np.zeros(4),1,[.52,.33,.15]);k=top_k(r,b,2);self.assertAlmostEqual(sum(q['posterior_mass'] for q in k['items']),.85);self.assertAlmostEqual(k['remaining_modeled_mass'],.15);self.assertIsNone(k['unmodeled_geometry_mass'])
 def test_gain_marginalized_matches_dense_covariance(self):
  F=np.array([[.1,.5,-1.],[.2,.3,.4]]);y=np.array([1.,-.2,.3]);sigma=.2;sd=.3
  L,_=log_likelihood(F,y,sigma,sd)
  for i,f in enumerate(F):
   c=sigma**2*np.eye(3)+sd**2*np.outer(f,f);r=y-f
   expected=-.5*(r@np.linalg.solve(c,r)+np.linalg.slogdet(c)[1]+3*np.log(2*np.pi));self.assertAlmostEqual(L[i],expected,12)
 def test_evidence_not_best_residual(self):
  b=self.bank(); b['F'][b['H']==1]=2.;b['F'][np.flatnonzero(b['H']==1)[0]]=0.
  r=infer(b,np.zeros(4),.1,[.5,.5,1e-100]);self.assertGreater(r['topology_probabilities'][0],.9)
 def test_permittivity_void_two_regions(self):
  _,im,v=material_image(2,[2.5,.1],C); self.assertTrue(v[75,75]);self.assertTrue(v[75,110]);self.assertFalse(v[75,97]);self.assertGreater(im.imag.max(),0)
 def test_real_packing_noise(self):
  from src.b1.physics import angular_data
  self.assertEqual(angular_data(np.zeros(10),np.zeros(10),C['angles_degrees'],C['wavenumber']).shape,(32,))
 def test_layer_definitions(self):
  rr,ee=topology_layers(2,2.5+.1j,C);np.testing.assert_allclose(rr,[.24,.46,.6,1]);np.testing.assert_allclose(ee,[1,2.5+.1j,1,2.5+.1j])
 def test_posterior_weights_sum(self):
  b=self.bank();r=infer(b,np.zeros(4),1,[.2,.3,.5]);self.assertAlmostEqual(r['weights'].sum(),1,13)
 def test_reference_requires_both_metrics(self):
  a=dict(topology_probabilities=[.99,.005,.005],log_bayes_factors_vs_solid=[0,-6,-6]); b=dict(topology_probabilities=[.99,.005,.005],log_bayes_factors_vs_solid=[0,-7,-7]);self.assertEqual(refinement(a,b,C)['status'],'REFERENCE_UNRESOLVED')
 def test_persistent_failure_charged(self):
  with tempfile.TemporaryDirectory() as tmp:
   l=Ledger(tmp,1)
   def bad():raise RuntimeError('injected failure, no physical solver')
   with self.assertRaises(RuntimeError):l.evaluate('SYNTHETIC',{},bad)
   self.assertEqual(Ledger(tmp,1).count(),1)
   with self.assertRaises(PhysicsLimit):l.evaluate('SYNTHETIC',{},lambda:0)
   self.assertEqual(l.summary()['failures'],1)
 def test_support_exclusion(self):
  p=np.array([.5,.3,.2]);q=np.array([.625,.375,0]);self.assertGreaterEqual(.5*abs(p-q).sum()+1e-15,p[2])
 def test_projective_probability_signature(self):
  rng=np.random.default_rng(20);J=rng.normal(size=(5,4));Y=rng.normal(size=(5,2));N=rng.normal(size=(2,4));G=np.array([[1.2,.1],[.2,1.3]]);s=rng.normal(size=(40,4));r=rng.normal(size=5)
  J1=J+Y@np.linalg.solve(G,N);c=np.linalg.solve(G,(s@N.T).T).T; h=np.sum((r@Y+s@J.T@Y)*c,axis=1)+.5*np.sum((c@Y.T)**2,axis=1)
  phi0=.5*np.sum((r+s@J.T)**2,axis=1);phi1=.5*np.sum((r+s@J1.T)**2,axis=1);np.testing.assert_allclose(phi1-phi0,h,atol=5e-13)
 def test_density_correction_identity(self):
  prior=np.array([.1,.2,.7]);q=np.array([.4,.3,.3]);lik=np.array([.2,.8,.4]);p=prior*lik;p/=p.sum();corr=q*(prior*lik/q);corr/=corr.sum();np.testing.assert_allclose(corr,p)
 def test_constant_potential_offset_cancels(self):
  x=np.array([-2.,.3,4.]);p=np.exp(x-logsumexp(x));q=np.exp(x+8-logsumexp(x+8));np.testing.assert_allclose(p,q)
 def test_field_offset_does_not_cancel(self):
  F=np.array([[-1.],[1.]]);l,_=log_likelihood(F,np.array([.4]),1);l2,_=log_likelihood(F+2,np.array([.4]),1);self.assertGreater(abs((l[1]-l[0])-(l2[1]-l2[0])),1)

if __name__=='__main__':unittest.main()

class RefinementTests(unittest.TestCase):
 def test_polynomial_reproduces_paid_training_fields(self):
  from src.b1.refine import fit,predict
  b=dict(np.load('results/banks/order7.npz'));coef=fit(b,7)
  for h in range(3):np.testing.assert_allclose(predict(coef,h,b['theta'][b['H']==h],b['bounds']),b['F'][b['H']==h],atol=5e-15)
 def test_constant_likelihood_continuous_quantiles(self):
  b=ProbabilityTests().bank();b['bounds']=np.array([[2.2,2.8],[.04,.16]])
  r=infer(b,np.zeros(4),1,[1/3]*3)
  np.testing.assert_allclose(r['conditional_material'][0]['quantiles_05_50_95'],[[2.23,2.5,2.77],[.046,.1,.154]],atol=1e-14)
 def test_budget_cap_immutable(self):
  with tempfile.TemporaryDirectory() as tmp:
   l=Ledger(tmp,1);l.evaluate('SYNTHETIC',{},lambda:0)
   with self.assertRaises(PhysicsLimit):Ledger(tmp,2)
 def test_budget_ledger_deletion_refused(self):
  with tempfile.TemporaryDirectory() as tmp:
   l=Ledger(tmp,1);l.evaluate('SYNTHETIC',{},lambda:0);l.path.unlink()
   with self.assertRaises(PhysicsLimit):Ledger(tmp,1)
 def test_real_experiment_complete_reservations(self):
  rows=[json.loads(x) for x in Path('results/PHYSICS_LEDGER.jsonl').read_text().splitlines()];begins=[r['id'] for r in rows if r['event']=='BEGIN'];ends=[r['id'] for r in rows if r['event']=='END']
  self.assertEqual(begins,list(range(1,385)));self.assertEqual(sorted(ends),begins);self.assertEqual(len(ends),len(set(ends)))

class TargetMetadataTests(unittest.TestCase):
 def test_gain_snapshot_contains_exact_likelihood(self):
  q=np.load('results/cases/refined_gain_marginalized/posterior.npz');self.assertEqual(float(q['gain_sd']),.15)
 def test_gain_cli_replay_matches_frozen_target(self):
  a=json.loads(Path('results/cli_gain_regression_corrected/probabilities.json').read_text());b=json.loads(Path('results/cases/refined_gain_marginalized/probabilities.json').read_text());np.testing.assert_allclose(a['result']['topology_probabilities'],b['result']['topology_probabilities'],atol=1e-14);self.assertEqual(a['measurement_target']['gain_sd'],.15)
 def test_stale_material_prior_is_rejected(self):
  from src.b1.cli_infer import validate_frozen
  import copy
  bad=copy.deepcopy(C);bad['material_prior']['eps_real']=[2.,3.]
  with self.assertRaises(ValueError):validate_frozen(Path.cwd(),bad)
  validate_frozen(Path.cwd(),C)

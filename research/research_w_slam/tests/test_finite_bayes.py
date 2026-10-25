"""Pure probability checks; no forward solver or physical simulation."""
import unittest
import numpy as np
from wslam.finite_bayes import (
    build_grid,build_coarse_grid,summarize,update_logweights,
    conditional_target_only,posterior_metrics,importance_sample,tempered_smc)


def tiny_grid():
    # Two walls, each with the same two target positions.
    return dict(label='TEST_DISCRETE_PRIOR',prior_weight=np.array([.1,.1,.4,.4]),
        support_id=np.array([0,0,1,1]),support=np.array([[0,0,0,0],[0,0,0,0],[1,0,0,0],[1,0,0,0]]),
        epsinf=np.array([3.75,3.75,2.,2.]),ell=np.array([.125,.125,.1,.1]),
        pose=np.array([[-1.,2.],[1.,2.],[-1.,2.],[1.,2.]]),wall_group_id=np.array([0,0,1,1]),
        legal_supports=np.array([[[0,0],[0,0]],[[1,0],[0,0]]]))


class FiniteBayesTests(unittest.TestCase):
    def test_prior_physical_count_mass_and_pose_outer_order(self):
        g=build_grid()
        self.assertEqual(len(g['prior_weight']),2950)
        self.assertEqual(g['wall_group_count'],118)
        self.assertEqual(len(g['legal_supports']),14)
        self.assertAlmostEqual(g['prior_weight'].sum(),1.)
        self.assertTrue(np.all(g['prior_weight']>0))
        self.assertFalse(g['continuous_uniform_integration_certificate'])
        self.assertTrue(np.all(g['pose'][:118]==g['pose'][0]))
        self.assertFalse(np.array_equal(g['pose'][117],g['pose'][118]))
        empty=~g['support'].any(axis=1)
        self.assertEqual(int(empty.sum()),25)
        self.assertTrue(np.all(g['epsinf'][empty]==3.75))
        self.assertTrue(np.all(g['ell'][empty]==.125))
        mass=np.bincount(g['support_id'],weights=g['prior_weight'],minlength=14)
        self.assertTrue(np.allclose(mass,g['support_prior']))
        # The shared material prior at fixed nonempty support has CC3 masses.
        chosen=np.flatnonzero(g['legal_supports'].reshape(14,4).any(axis=1))[0]
        selected=g['support_id']==chosen
        for node,expected in [(1.5,1/6),(3.75,4/6),(6.,1/6)]:
            actual=g['prior_weight'][selected&(g['epsinf']==node)].sum()/mass[chosen]
            self.assertAlmostEqual(actual,expected)
        # Pose CC5 masses are independent of all wall states.
        xmass=[g['prior_weight'][np.isclose(g['pose'][:,0],x)].sum() for x in np.unique(g['pose'][:,0])]
        self.assertTrue(np.allclose(xmass,[1/30,4/15,2/5,4/15,1/30]))

    def test_coarse_bank_has_no_posterior_claim(self):
        g=build_coarse_grid()
        self.assertEqual(g['state_count'],477)
        self.assertEqual(g['wall_group_count'],53)
        self.assertFalse(g['posterior_claims_allowed'])
        self.assertEqual(g['label'],'PROPOSAL_BANK_ONLY')
        self.assertAlmostEqual(g['prior_weight'].sum(),1.)

    def test_exact_unique_observations_and_stale_wall_approximation(self):
        g=tiny_grid()
        ll=np.log(np.array([1.,3.,1.,1.]))
        posterior=np.exp(update_logweights(g['prior_weight'],[('new',ll)]))
        expected=g['prior_weight']*np.exp(ll)
        expected/=expected.sum()
        self.assertTrue(np.allclose(posterior,expected))
        with self.assertRaises(ValueError):
            update_logweights(g['prior_weight'],[('same',ll),('same',ll)])
        stale=conditional_target_only(g['prior_weight'],ll,g['wall_group_id'])
        self.assertEqual(stale['label'],'STALE_WALL_MARGINAL')
        self.assertFalse(stale['exact_new_target_posterior'])
        self.assertTrue(np.allclose(stale['weights'],[.05,.15,.4,.4]))
        self.assertAlmostEqual(stale['weights'][:2].sum(),.2)
        self.assertGreater(posterior[:2].sum(),.2)

    def test_moments_conditionals_and_zero_occupancy(self):
        g=tiny_grid()
        s=summarize(g['prior_weight'],g)
        self.assertTrue(np.allclose(s['occupancy'],[.8,0,0,0]))
        self.assertEqual(s['joint_cov'].shape,(10,10))
        self.assertTrue(np.allclose(s['pose_mean'],[0,2]))
        self.assertAlmostEqual(s['pose_cov'][0,0],1.)
        self.assertAlmostEqual(s['full_wall_cov'][0,0],.16)
        self.assertAlmostEqual(s['joint_cov'][0,4],.016)
        self.assertAlmostEqual(s['event_px_ge_0'],.5)
        self.assertTrue(s['conditional_material'][0]['defined'])
        self.assertFalse(s['conditional_material'][1]['defined'])
        self.assertTrue(np.array_equal(s['conditional_material'][0]['epsinf']['values'],[2.]))
        self.assertTrue(np.allclose(s['conditional_material'][0]['joint']['mass'],[1.]))
        posecdf=s['coordinate_distributions'][8]
        self.assertTrue(np.allclose(posecdf['cdf'],[.5,1.]))
        self.assertTrue(np.array_equal(posecdf['quantiles'],[-1.,-1.,1.]))
        # Empty-wall placeholders never leak into contrast/loss physical states.
        empty=summarize([.5,.5,0,0],g)
        self.assertTrue(np.allclose(empty['full_wall_mean'],0))
        self.assertTrue(np.allclose(empty['full_wall_cov'],0))

    def test_metrics_zero_and_known_squared_loss_regret(self):
        g=tiny_grid()
        r=g['prior_weight']
        same=posterior_metrics(r,r,g)
        self.assertAlmostEqual(same['action_regret'],0.)
        self.assertTrue(np.allclose(same['coordinate_cdf_sup'],0))
        approximate=np.array([.2,0,.8,0])
        metrics=posterior_metrics(approximate,r,g)
        self.assertAlmostEqual(metrics['event_px_ge_0_error'],.5)
        self.assertAlmostEqual(metrics['coordinate_cdf_sup'][8],.5)
        self.assertAlmostEqual(metrics['action_regret'],1.)
        self.assertAlmostEqual(metrics['occupancy_max_error'],0.)
        self.assertAlmostEqual(metrics['mode_mass_max_error_ref_ge_005'],0.)

    def test_primary_target_risk_ignores_wall_material_differences(self):
        g=tiny_grid()
        reference=np.array([.1,.1,.4,.4])
        material_shift=np.array([.4,.4,.1,.1])
        metrics=posterior_metrics(material_shift,reference,g)
        self.assertAlmostEqual(metrics['action_regret'],0.)
        self.assertAlmostEqual(metrics['reference_bayes_risk'],1.)
        self.assertAlmostEqual(metrics['approximate_bayes_risk'],1.)
        self.assertEqual(metrics['action'].shape,(2,))
        self.assertGreater(metrics['occupancy_max_error'],.5)

    def test_matching_marginals_do_not_certify_joint_branches(self):
        g=tiny_grid()
        reference=np.array([.5,0,0,.5])
        swapped=np.array([0,.5,.5,0])
        metrics=posterior_metrics(swapped,reference,g)
        self.assertTrue(np.allclose(metrics['coordinate_cdf_sup'],0))
        self.assertAlmostEqual(metrics['mode_mass_max_error_ref_ge_005'],0.)
        self.assertAlmostEqual(metrics['action_regret'],0.)
        self.assertAlmostEqual(metrics['joint_total_variation'],1.)
        self.assertGreater(np.linalg.norm(metrics['wall_pose_cross_covariance_difference']),.5)
        self.assertFalse(metrics['joint_diagnostics_change_frozen_gate'])

    def test_prior_importance_and_tempered_smc_reproducibility(self):
        g=tiny_grid()
        ll=np.log(np.array([1.,3.,1.,2.]))
        exact=g['prior_weight']*np.exp(ll)
        exact/=exact.sum()
        for sampler in (importance_sample,tempered_smc):
            result=sampler(g,ll,n_particles=12000,seed=1807)
            again=sampler(g,ll,n_particles=12000,seed=1807)
            self.assertTrue(np.array_equal(result['indices'],again['indices']))
            self.assertTrue(np.array_equal(result['weights'],again['weights']))
            self.assertAlmostEqual(result['weights'].sum(),1.)
            self.assertFalse(result['global_posterior_certificate'])
            self.assertLess(np.abs(result['weights']-exact).max(),.03)
        smc=tempered_smc(g,ll,n_particles=64,seed=9,resample_ess_fraction=1.)
        self.assertEqual(smc['observation_likelihood_insertions'],1)
        self.assertEqual(smc['final_temperature'],1.)
        self.assertTrue(any(x['resampled'] for x in smc['history']))
        with self.assertRaises(ValueError):
            tempered_smc(g,ll,temperatures=[0,.7,.5,1])


if __name__=='__main__':
    unittest.main()

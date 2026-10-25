"""Substantive finite-dimensional checks; no physical solve."""
import unittest

import numpy as np

from wslam.information import (schur, dynamic_information,
    dynamic_information_full_rank, generalized_retention, conditional_innovation,
    local_gaussian_audit, target_structure, svd_audit)


class TestInformation(unittest.TestCase):
    def test_singular_hidden_target_counterexample(self):
        b = np.array([[1.,0,0,0],[0,0,1,0]])
        a = np.array([[1.],[-1.]])
        d = np.array([[-1.,-1,1,0],[0,-1,0,1]])
        c = np.linalg.pinv(b) @ a
        naive = c.T @ d.T @ np.linalg.solve(np.eye(2)+d @ np.linalg.pinv(b.T@b) @ d.T, d@c)
        self.assertAlmostEqual(naive.item(), 4/3)
        result = dynamic_information(a,b,d,np.eye(2))
        np.testing.assert_allclose(result['information'], 0, atol=1e-12)
        np.testing.assert_allclose(result['information'], schur(a,b,d.T@d)['information'], atol=1e-12)
        self.assertEqual(result['formula'], 'D2')
        with self.assertRaises(ValueError):
            dynamic_information_full_rank(a,b,d,np.eye(2))

    def test_sharing_and_finite_motion_prior_assistance(self):
        aa = [np.ones((1,1)), np.ones((1,1))]
        bb = [np.ones((1,1)), -np.ones((1,1))]
        self.assertAlmostEqual(target_structure(aa,bb,structure='shared')['information'].item(), 2)
        self.assertAlmostEqual(target_structure(aa,bb,structure='independent')['information'].item(), 0)
        motion = target_structure(aa,bb,structure='motion',d=np.array([[-1.,1.]]),q=np.array([[.3]]))
        self.assertTrue(motion['prior_assisted'])
        self.assertAlmostEqual(motion['data_only_information'].item(), 0)
        self.assertAlmostEqual(motion['information'].item(), 4/2.3)

    def test_retention_angles_and_metric_invariance(self):
        a = np.array([[2.,0,0],[0,3,0],[0,0,0]])
        b = np.array([[1.],[1.],[0.]])
        out = generalized_retention(a,b,material_metric=np.diag([4.,9.,2.]))
        self.assertEqual(out['wall_rank_audit']['rank'], 2)
        np.testing.assert_allclose(out['eigenvalues'], [0,1], atol=1e-12)
        np.testing.assert_allclose(np.sort(out['bare_angle_retention']),out['eigenvalues'],atol=1e-12)
        self.assertAlmostEqual(out['angle_sum'], out['retention_trace'])
        with_prior = generalized_retention(a,b,prior_precision=np.array([[2.]]))
        self.assertTrue(with_prior['prior_assisted'])
        self.assertGreater(with_prior['retention_trace'], with_prior['angle_sum'])

    def test_d3_matches_direct_joint_update(self):
        h = np.array([[4.,1.],[1.,2.]])
        a, b = np.array([[2.]]), np.array([[3.]])
        out = conditional_innovation(h,1,a,b)
        j = np.hstack([a,b])
        hn = h+j.T@j
        direct = hn[:1,:1]-hn[:1,1:]@np.linalg.solve(hn[1:,1:],hn[1:,:1])
        np.testing.assert_allclose(out['new_information'],direct)

    def test_gaussian_covariance_and_inexact_dual_interval(self):
        h = np.array([[3.,.4,.1],[.4,2.,.3],[.1,.3,1.5]])
        g, c = np.array([1.,-.5,.3]), np.array([[1.,2.,-1.],[0.,1.,1.]])
        exact = local_gaussian_audit(h,g,c,[0,2])
        approx = local_gaussian_audit(h,g,c,[0,2],approximate_duals=exact['exact_duals']+.03)
        lo, hi = approx['dual_interval']
        self.assertLessEqual(lo, approx['task_regret'])
        self.assertLessEqual(approx['task_regret'], hi)
        np.testing.assert_allclose(exact['task_covariance']-exact['pinned_task_covariance'],exact['omitted_task_covariance'])
        delta = c@(exact['mean']-exact['pinned_mean'])
        self.assertAlmostEqual(exact['task_regret'], delta@delta)
        self.assertEqual(exact['scope'], 'LOCAL_GAUSSIAN_EXACT_SMALL_SYSTEM')
        self.assertFalse(exact['true_posterior_certificate'])

    def test_d7_partial_reopen_can_increase_regret(self):
        h, g, c = np.eye(2), np.array([1.,-1.]), np.array([[1.,1.]])
        before = local_gaussian_audit(h,g,c,[0,1])
        after = local_gaussian_audit(h,g,c,[1])
        self.assertAlmostEqual(before['task_regret'], 0)
        self.assertAlmostEqual(after['task_regret'], 1)
        self.assertLess(after['objective_gap'], before['objective_gap'])

    def test_rank_zero_and_input_guards(self):
        out = svd_audit(np.zeros((2,3)))
        self.assertEqual(out['rank'], 0)
        self.assertEqual(out['null_basis'].shape, (3,3))
        with self.assertRaises(ValueError):
            schur(np.ones((1,1),complex),np.ones((1,1)))
        with self.assertRaises(ValueError):
            local_gaussian_audit(np.zeros((2,2)),np.zeros(2),np.eye(2),[0])

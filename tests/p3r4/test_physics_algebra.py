"""Small generic matrix algebra only: this file never imports a Maxwell backend."""
import unittest
import numpy as np
from scipy.stats import multivariate_normal
from p3r4.reuse import (FiniteMaterialCache, ReuseUnavailable, compose_jacobian,
    jacobian_evolution, broyden_update, secant_decomposition, line_events, ellipse_events,
    stability_from_saved_inverse, dual_corrected_prediction, same_uniform_decision)
from p3r4.geometry import (bracket_from_differentials, material_bracket_diagnostic,
    full_rank_covariance, gaussian_log_proposal, corrected_log_acceptance, soft_mobility)


class Algebra(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(2026102504)
        self.n, self.sources, self.receivers, self.latents = 9, 3, 4, 5
        complex_array = lambda shape: self.rng.normal(size=shape)+1j*self.rng.normal(size=shape)
        self.G = complex_array((9, 9))*.025
        self.S = complex_array((4, 9))
        self.incident = complex_array((9, 3))
        self.chi = complex_array((9,))*.1
        self.L0 = np.eye(9)-self.chi[:, None]*self.G
        self.R0 = np.linalg.inv(self.L0)  # SYNTHETIC ONLY; never in reuse module
        self.j0 = self.R0@(self.chi[:, None]*self.incident)
        self.E0 = self.incident+self.G@self.j0
        self.F0 = self.S@self.j0
        self.I = np.array([1, 3, 5])
        self.cache = FiniteMaterialCache(self.chi, self.j0, self.E0, self.F0,
            self.G, self.S, self.incident, self.I, self.R0[:, self.I],
            'synthetic_base', len(self.I), 'synthetic_fixture')

    def test_exact_finite_response_all_terms(self):
        chi1 = self.chi.copy()
        chi1[self.I] += .3+.1j
        result = self.cache.evaluate(chi1)
        R1 = np.linalg.inv(np.eye(9)-chi1[:, None]*self.G)
        j1 = R1@(chi1[:, None]*self.incident)
        E1 = self.incident+self.G@j1
        self.assertLess(np.linalg.norm(result['current']-j1), 1e-14)
        np.testing.assert_allclose(result['field'], E1, atol=1e-14)
        np.testing.assert_allclose(result['prediction'], self.S@j1, atol=1e-14)
        self.assertLess(max(v for v in self.cache.audit().values() if isinstance(v, float)), 1e-14)
        self.assertEqual(result['profile']['new_physical_RHS'], 0)
        D0 = self.rng.normal(size=(9, 5))+1j*self.rng.normal(size=(9, 5))
        D1 = D0+self.rng.normal(size=(9, 5))*.2
        T0, T1 = self.S@self.R0, self.S@R1
        delta = compose_jacobian(T1, E1, D1, .3)-compose_jacobian(T0, self.E0, D0, .3)
        terms = jacobian_evolution(T0, T1, self.E0, E1, D0, D1, .3)
        np.testing.assert_allclose(delta, terms['total'], atol=2e-13)
        for name in ('propagation', 'total_field', 'decoder_derivative'):
            self.assertGreater(np.linalg.norm(terms[name]), 1e-5)

    def test_cumulative_support_and_properties_fail_closed(self):
        chi1 = self.chi.copy()
        chi1[2] += .01
        with self.assertRaises(ReuseUnavailable):
            self.cache.evaluate(chi1, property_changed=True)
        chi1 = self.chi.copy()
        chi1[self.I] += .01
        with self.assertRaises(ReuseUnavailable):
            self.cache.evaluate(chi1, max_rank=2)
        # A second candidate is still compared to original chi, not first endpoint.
        self.cache.evaluate(chi1)
        chi2 = self.chi.copy()
        chi2[1] += .02
        np.testing.assert_allclose(self.cache.evaluate(chi2)['prediction'],
            self.S@np.linalg.solve(np.eye(9)-chi2[:, None]*self.G, chi2[:, None]*self.incident), atol=1e-14)
        same = self.cache.evaluate(self.chi)
        np.testing.assert_array_equal(same['prediction'], self.F0)

    def test_stale_inverse_columns_rejected(self):
        with self.assertRaises(ReuseUnavailable):
            FiniteMaterialCache(self.chi, self.j0, self.E0, self.F0,
                self.G, self.S, self.incident, self.I, np.zeros((9, 3)),
                'bad_base', 3, 'synthetic')

    def test_conditional_residual_bounds_and_same_uniform(self):
        chi1 = self.chi.copy()
        chi1[self.I] += .02+.01j
        L1 = np.eye(9)-chi1[:, None]*self.G
        rhs = chi1[:, None]*self.incident
        stability = stability_from_saved_inverse(self.L0, self.R0, L1)
        self.assertGreater(stability['beta'], 0)
        basis = np.linalg.qr(self.rng.normal(size=(9, 3)))[0]
        y = self.rng.normal(size=(4, 3))*.1
        result = dual_corrected_prediction(L1, rhs, self.S, basis,
            (self.S@self.R0).conj().T, stability, observation=y, sigma=.2)
        phi = np.linalg.norm(self.S@np.linalg.solve(L1, rhs)-y)**2/.2**2
        lo, hi = result['likelihood_interval']
        self.assertLessEqual(lo, phi)
        self.assertLessEqual(phi, hi)
        self.assertIn('NOT_INTERVAL_CERTIFIED', result['bound_label'])
        self.assertEqual(same_uniform_decision(None, [1, 2], 0, -.5)['action'], 'REFINE')
        self.assertEqual(same_uniform_decision([1, 1], [1, 1], 0, -.5)['action'], 'ACCEPT')
        self.assertEqual(same_uniform_decision([1, 1], [2, 2], 0, -.5)['action'], 'REJECT')

    def test_empty_support_and_input_ownership(self):
        cache = FiniteMaterialCache(self.chi, self.j0, self.E0, self.F0,
            self.G, self.S, self.incident, np.array([], int), np.zeros((9, 0)),
            'empty', 0, 'synthetic')
        np.testing.assert_array_equal(cache.evaluate(self.chi)['prediction'], self.F0)
        before = self.cache.chi0.copy()
        self.chi[0] += 2
        np.testing.assert_array_equal(before, self.cache.chi0)

    def test_secant_vector_accounting(self):
        J = self.rng.normal(size=(12, 5))
        h = self.rng.normal(size=5)
        df = self.rng.normal(size=12)
        np.testing.assert_allclose(broyden_update(J, h, df)@h, df, atol=1e-14)
        D = self.rng.normal(size=(9, 5))
        dchi = D@h+self.rng.normal(size=9)*.02
        response = self.S@(self.E0*dchi[:, None])
        out = secant_decomposition(response, compose_jacobian(self.S, self.E0, D, .2), h,
            self.S, self.E0, dchi, D, .2)
        np.testing.assert_allclose(out['linear']+out['decoder']+out['maxwell'], out['actual'])
        np.testing.assert_allclose(out['remainder'], out['decoder']+out['maxwell'])

    def test_event_roots(self):
        np.testing.assert_array_equal(line_events([0, 2], [1, 0]), [-1, 1])
        for theta in ellipse_events([2.], [.7], [.2]):
            self.assertLess(min(abs(2*np.cos(theta)+.7*np.sin(theta)+.2-s) for s in (-1, 1)), 2e-14)
        self.assertEqual(len(ellipse_events([0], [0], [1])), 0)

    def test_full_rank_reverse_proposal_density(self):
        Cx = full_rank_covariance(5, weak=self.rng.normal(size=(5, 2)), floor=.02)
        Cy = full_rank_covariance(5, brackets=self.rng.normal(size=(5, 1)), floor=.03)
        x, y = self.rng.normal(size=(2, 5))
        a = gaussian_log_proposal(y, x, Cx)
        self.assertAlmostEqual(a, multivariate_normal.logpdf(y, mean=x, cov=Cx), places=10)
        expected = min(0., -.3+.1+multivariate_normal.logpdf(x, mean=y, cov=Cy)-a)
        actual = corrected_log_acceptance(x, y, logtarget_x=-.1, logtarget_candidate=-.3,
            mean_at_x=x, covariance_at_x=Cx, mean_at_candidate=y, covariance_at_candidate=Cy)
        self.assertAlmostEqual(actual, expected, places=10)
        with self.assertRaises(np.linalg.LinAlgError):
            gaussian_log_proposal(y, x, np.diag([1, 1, 0, 0, 0]))
        self.assertGreater(np.linalg.eigvalsh(soft_mobility(self.rng.normal(size=(3, 5)))).min(), 0)

    def test_bracket_convention_and_material_metric(self):
        # Heisenberg X=(1,0,-b/2,0), Y=(0,1,a/2,0): bracket=e_c.
        X, Y = np.array([1., 0, -.2, 0]), np.array([0., 1., .3, 0])
        DX, DY = np.zeros((4, 4)), np.zeros((4, 4))
        DX[2, 1], DY[2, 0] = -.5, .5
        bracket = bracket_from_differentials(X, Y, DX, DY)
        np.testing.assert_array_equal(bracket, [0, 0, 1, 0])
        result = material_bracket_diagnostic(bracket, np.eye(4), np.eye(4)[:, :2])
        self.assertEqual(result['normal_fraction'], 1.)
        self.assertFalse(result['full_lie_rank_proven'])


if __name__ == '__main__':
    unittest.main()

"""Synthetic algebra only; no antenna/Green construction or physical job.

Every synthetic LU and RHS still uses the real reservation Book. These tests
are defined here but must not be invoked until the parent permits execution.
"""
import numpy as np
from research_w.accounting import Book
from research_w.vie import VIE


def synthetic(tmp_path):
    rng=np.random.default_rng(71)
    p=VIE.__new__(VIE)
    p.N=5;p.nx=5;p.ny=1;p.freq=1.;p.phase='repair';p.book=Book(tmp_path)
    matrix=.03*(rng.normal(size=(5,5))+1j*rng.normal(size=(5,5)))
    p.green_block=lambda rows,cols:matrix[np.ix_(np.asarray(rows),np.asarray(cols))]
    p.S=rng.normal(size=(12,5))+1j*rng.normal(size=(12,5))
    p.incident=rng.normal(size=(5,4))+1j*rng.normal(size=(5,4))
    p.direct=np.zeros((12,4),complex);p.setup_seconds=0.;p.residual_tolerance=1e-12
    return p


def test_jacobian_source_major_and_real_imag_fd(tmp_path):
    p=synthetic(tmp_path);eps=np.array([2+.04j,1.,3+.02j,1.,1.5+.03j])
    state=p.solve(eps,jacobian=True)
    rng=np.random.default_rng(9);direction=rng.normal(size=5)
    # Keep exact-zero reference cells inside legal material bounds on both sides.
    direction[[1,3]]=0;step=1e-6
    plus=p.solve(eps+step*direction);minus=p.solve(eps-step*direction)
    fd=(plus.data-minus.data).ravel()/(2*step)
    np.testing.assert_allclose(state.jacobian[:,:5]@direction,fd,rtol=2e-8,atol=2e-8)
    assert state.data.shape==(4,12)
    np.testing.assert_array_equal(state.jacobian[:,5:],1j*state.jacobian[:,:5])
    totals=p.book.totals()
    assert totals['factorizations']==3 and totals['rhs']==24


def test_inactive_cell_derivative_is_retained(tmp_path):
    p=synthetic(tmp_path);eps=np.array([2+.02j,1,1,1,1],complex)
    state=p.solve(eps,jacobian=True);step=1e-7
    shifted=eps.copy();shifted[3]+=step
    fd=(p.solve(shifted).data-state.data).ravel()/step
    np.testing.assert_allclose(state.jacobian[:,3],fd,rtol=2e-6,atol=2e-6)
    assert len(state.active)==1


def test_zero_contrast_has_no_native_solve_but_nonzero_jacobian(tmp_path):
    p=synthetic(tmp_path);state=p.solve(np.ones(5),jacobian=True)
    assert state.factor is None and len(state.active)==0
    assert p.book.totals()['rhs']==0 and p.book.totals()['factorizations']==0
    expected=np.einsum('rn,ns->srn',p.S,p.incident).reshape(48,5)
    np.testing.assert_allclose(state.jacobian[:,:5],expected)
    np.testing.assert_allclose(state.current,0)


def test_invalid_material_rejected_before_reservation(tmp_path):
    p=synthetic(tmp_path)
    import unittest
    for value in (0.9,9.1,2-.01j,2+.151j):
        with unittest.TestCase().assertRaises(ValueError):p.solve(np.full(5,value))
    assert p.book.totals()['factorizations']==0

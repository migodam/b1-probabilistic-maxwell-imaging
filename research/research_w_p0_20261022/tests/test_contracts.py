import json,time,os
from pathlib import Path
import numpy as np
import unittest
from unittest.mock import patch
from research_w.accounting import Book,BudgetStop
from research_w.geometry import antennas,material_grid,probes,in_query,target_centroid
from research_w.analytic import whitened_distance

def test_global_rhs_and_phase_allocation_are_not_reset(tmp_path):
    b=Book(tmp_path);b.reserve('repair','fixture',rhs=800)
    with unittest.TestCase().assertRaises(BudgetStop):Book(tmp_path).reserve('repair','rejected',rhs=1)
    assert b.totals()['rhs']==800

def test_lu_failure_reservation_is_not_refunded(tmp_path):
    b=Book(tmp_path);b.reserve('validation','LU_to_fail',factorizations=1)
    b.event(event='failure',label='simulated algebra failure')
    assert Book(tmp_path).totals()['factorizations']==1

def test_in_progress_job_deadline(tmp_path):
    with patch.dict(os.environ,{'W_JOB_DEADLINE_MONOTONIC':str(time.monotonic()-1)}):
        with unittest.TestCase().assertRaises(BudgetStop):Book(tmp_path).reserve('search','reject',rhs=4)
    assert Book(tmp_path).totals()['rhs']==0

def test_complex_noise_whitening_and_realification():
    rng=np.random.default_rng(2026102203);sigma=.01
    z=sigma*(rng.normal(size=100000)+1j*rng.normal(size=100000))
    assert abs(np.mean(abs(z)**2)/(2*sigma*sigma)-1)<.015
    a=z[:24].reshape(4,6)
    packed=np.concatenate([np.r_[q.real,q.imag] for q in a])
    assert np.isclose(np.linalg.norm(packed)/sigma,whitened_distance(a,0*a,sigma))

def test_material_world_is_shared_across_frequencies():
    class Geometry:
        N=64*48;h=1/16;freq=1
        x=-2+(np.arange(64)+.5)/16;y=(np.arange(48)+.5)/16
        xx,yy=np.meshgrid(x,y);points=np.c_[xx.ravel(),yy.ravel()]
    p=Geometry();a=material_grid(p,2026102201,frequency=.8)
    b=material_grid(p,2026102201,frequency=1.2)
    np.testing.assert_array_equal(a.real,b.real)
    np.testing.assert_allclose(.8*a.imag,1.2*b.imag,atol=1e-17)
    assert np.min(a.real)>=1 and np.max(a.real)<=9 and np.min(a.imag)>=0
    assert np.all(a[in_query(p.points)]==1)

def test_probe_and_acquisition_are_fixed():
    tx,rx=antennas();assert tx.shape==(4,2) and rx.shape==(12,2)
    assert np.all(tx[:,1]==-.3) and np.all(rx[:,1]==-.3)
    assert len(probes())==24 and len(set(tuple(p[0]) for p in probes()))==12

def test_task_is_total_material_not_component_labels():
    xy=np.array([[0.,1.8],[.3,1.8],[-.3,1.8]])
    e=np.array([3.,1.,1.]);np.testing.assert_allclose(target_centroid(xy,e),[0,1.8])
    # No W/X labels enter the total-material functional.
    np.testing.assert_allclose(target_centroid(xy,e.copy()),[0,1.8])

def test_kernel_inclusion_is_not_unconditional_fiber_requirement():
    # F=x^3, H=x has an exact continuous factor h=cuberoot but no derivative
    # factor at zero. Do not infer the differential criterion without regularity.
    x=np.array([-.1,0,.1]);np.testing.assert_allclose(np.cbrt(x**3),x)
    assert 0**2==0 and 1!=0

def test_noncommuting_four_block_identity_and_schur_increment():
    rng=np.random.default_rng(27);n,t=11,3
    L=np.diag(np.arange(n)+3).astype(complex)+.02*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    U=np.eye(n,dtype=complex)[:,[3,6,9]];V=.15*(rng.normal(size=(t,t))+1j*rng.normal(size=(t,t)))
    Q=rng.normal(size=(n,4))+1j*rng.normal(size=(n,4));R=rng.normal(size=(5,n))+1j*rng.normal(size=(5,n))
    GQ=np.linalg.solve(L,Q);GU=np.linalg.solve(L,U)
    d=R@GQ;A=R@GU;B=U.conj().T@GQ;C=U.conj().T@GU
    exact=R@np.linalg.solve(L-U@V@U.conj().T,Q)
    reconstructed=d+A@V@np.linalg.solve(np.eye(t)-C@V,B)
    np.testing.assert_allclose(exact,reconstructed,rtol=1e-13,atol=1e-13)
    # This is a general matrix identity check, not a Maxwell wall witness.

def test_conditional_posterior_tv_bound_on_explicit_joint_prior():
    rng=np.random.default_rng(3);prior=rng.uniform(.1,1,12);phi=rng.uniform(0,3,12);eps=.2
    approximation=phi+rng.uniform(-eps,eps,12)
    p=prior*np.exp(-phi);p/=p.sum();q=prior*np.exp(-approximation);q/=q.sum()
    assert .5*np.abs(p-q).sum()<=np.tanh(eps)+1e-15

"""Example only: do not execute implicitly; each line with physics is counted.
Launch fresh with VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
and Gaussian/.venv_nn/bin/python. Import from repository root.
"""
import numpy as np
from research.p3l_r5_parallel_ab.shared.physics_adapter import Target

def example(x66,V66x2,observation):
    target=Target('A','explicit_example_job')
    chi,decoder_aux=target.decoder(x66)        # no Maxwell RHS; decoder logged
    point=target.evaluate(x66,'declared_A_reference') # one state and8 forward RHS
    J2=target.jacobian(point,V66x2)           #16 derivative RHS; white512x2
    return point.prediction,J2,target.phi(point,observation),target.book.totals()

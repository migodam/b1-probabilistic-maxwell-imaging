import unittest
import numpy as np
from p3r4.dynamic_geometry import (Chart,GeometryInvalid,material_frame,projector_field,
    central_bracket,material_real,finite_material_normal)

class GeometryAlgebra(unittest.TestCase):
    def test_frame_coordinate_rotation_and_projector(self):
        rng=np.random.default_rng(20261025030);D=rng.normal(size=(11,5))+1j*rng.normal(size=(11,5));J=rng.normal(size=(17,5))
        frame=material_frame(J,D,2);R=np.linalg.qr(rng.normal(size=(5,5)))[0]
        transformed=material_frame(J@R,D@R,2)
        np.testing.assert_allclose(transformed['projector'],R.T@frame['projector']@R,atol=2e-13)
        np.testing.assert_allclose(frame['projector']@frame['projector'],frame['projector'],atol=2e-14)
        field=projector_field(frame,D,rng.normal(size=5))
        self.assertAlmostEqual(np.linalg.norm(material_real(D)@field),1.,places=13)

    def test_material_rank_fail_closed(self):
        with self.assertRaises(GeometryInvalid):material_frame(np.ones((3,5)),np.zeros((11,5)),2)

    def test_bracket_sign_and_normalization_derivative(self):
        # X=(1,0), Y=(0,x): [X,Y]=(0,1), not its negative.
        h=.002;X=lambda x:np.array([1.,0.]);Y=lambda x:np.array([0.,x[0]])
        x=np.array([.3,.4]);vx=X(x);vy=Y(x)
        bracket=central_bracket(X(x+h*vy),X(x-h*vy),Y(x+h*vx),Y(x-h*vx),h)
        np.testing.assert_allclose(bracket,[0.,1.],atol=1e-13)

    def test_probit_boundary_no_shift_and_finite_material_normal(self):
        q=np.r_[np.zeros(64),0.,.5]
        with self.assertRaises(GeometryInvalid):Chart(q,66)
        chart=Chart(q,64);np.testing.assert_array_equal(chart.native(chart.x0),q)
        D=np.eye(3,dtype=complex);frame={'Vw':np.eye(3)[:,:1]}
        value=finite_material_normal(np.array([0,0,1j]),D,frame)
        self.assertAlmostEqual(value['actual_normal_fraction'],1.)
        self.assertAlmostEqual(value['actual_material_RMS'],1/np.sqrt(3))

if __name__=='__main__':unittest.main()

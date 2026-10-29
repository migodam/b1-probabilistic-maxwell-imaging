"""Positive bilinear material prolongation on one fixed physical grid."""
import numpy as np
from scipy.sparse import coo_matrix


def interpolation_1d(fine,coarse):
    x=np.linspace(0,1,int(fine))
    t=x*(int(coarse)-1)
    lo=np.minimum(np.floor(t).astype(int),int(coarse)-2)
    a=t-lo
    B=np.zeros((int(fine),int(coarse)))
    B[np.arange(int(fine)),lo]=1-a
    B[np.arange(int(fine)),lo+1]=a
    return B


def material_basis(grid_shape,stage):
    ny,nx=map(int,grid_shape)
    cy,cx=(int(stage),int(stage)) if np.isscalar(stage) else tuple(map(int,stage))
    if (ny,nx)==(cy,cx):
        return np.eye(ny*nx)
    return np.kron(interpolation_1d(ny,cy),interpolation_1d(nx,cx))


def resample_material(coarse,old_shape,new_shape):
    # New coefficients sample the preceding material field with linear interpolation.
    return material_basis(new_shape,old_shape)@np.asarray(coarse)


def difference_operator(grid_shape,width=19.5):
    ny,nx=map(int,grid_shape)
    rows,cols,values=[],[],[]
    edge=0
    for y in range(ny):
        for x in range(nx):
            i=y*nx+x
            for j in ([i+1] if x+1<nx else [])+([i+nx] if y+1<ny else []):
                rows.extend([edge,edge]);cols.extend([i,j]);values.extend([-1/width,1/width]);edge+=1
    return coo_matrix((values,(rows,cols)),shape=(edge,ny*nx)).tocsr()


class Regularizer:
    def __init__(self,grid_shape,smooth=1e-4,tv=1e-4,epsilon=.01):
        self.D=difference_operator(grid_shape)
        self.smooth,self.tv,self.epsilon=float(smooth),float(tv),float(epsilon)
        self.edges=self.D.shape[0]

    def evaluate(self,m):
        d=self.D@m
        q=np.sqrt(d*d+self.epsilon**2)
        value=self.smooth*np.mean(d*d)+self.tv*np.mean(q-self.epsilon)
        gd=(2*self.smooth*d+self.tv*d/q)/self.edges
        g=np.asarray(self.D.T@gd).reshape(-1)
        weights=(2*self.smooth+self.tv/q)/self.edges
        return float(value),g,weights

"""Public material/target coordinates; no data, labels or truth-file loader."""
from dataclasses import dataclass
import numpy as np
from .geometry import ROI,validate_geometry,expand_support

def circle_vertices(center,radius,n):
    angle=2*np.pi*np.arange(n)/n
    return np.asarray(center)+radius*np.column_stack((np.cos(angle),np.sin(angle)))

def polygon_inside(points,center,radius,n):
    d=np.asarray(points)-np.asarray(center)
    angle=np.arctan2(d[:,1],d[:,0]);segment=2*np.pi/n
    delta=np.mod(angle,segment)-segment/2
    boundary=radius*np.cos(np.pi/n)/np.cos(delta)
    return np.linalg.norm(d,axis=1)<=boundary

def boundary_quadrature(center,radius,n,degree=4):
    v=circle_vertices(center,radius,n);edge=np.roll(v,-1,axis=0)-v
    length=np.linalg.norm(edge,axis=1)
    q,w=np.polynomial.legendre.leggauss(degree);t=(q+1)/2
    points=(v[:,None,:]+edge[:,None,:]*t[None,:,None]).reshape(-1,2)
    normal=np.column_stack((edge[:,1],-edge[:,0]))/length[:,None]
    normals=np.repeat(normal,degree,axis=0)
    ds=(length[:,None]*w[None,:]/2).ravel()
    radial=np.sum(normals*(points-np.asarray(center)),axis=1)/radius
    return points,normals,ds,radial

def region_indices(points,shape):
    points=np.asarray(points);ny,nx=shape
    ix=np.floor((points[:,0]-ROI[0])/(ROI[1]-ROI[0])*nx).astype(int)
    iy=np.floor((points[:,1]-ROI[2])/(ROI[3]-ROI[2])*ny).astype(int)
    valid=(points[:,0]>=ROI[0])&(points[:,0]<ROI[1])&(points[:,1]>=ROI[2])&(points[:,1]<ROI[3])
    return np.clip(iy,0,ny-1),np.clip(ix,0,nx-1),valid

@dataclass
class World:
    support: object
    epsinf: object
    ell: object
    pose: object
    radius: float=.12
    contrast: float=3.
    target_loss: float=.02
    def __post_init__(self):
        self.support=np.asarray(self.support,dtype=bool)
        self.pose=np.asarray(self.pose,dtype=float)
        check=validate_geometry(self.support,self.epsinf,self.ell,self.pose,self.radius)
        if not check['legal']:raise ValueError('Illegal candidate: '+str(check['reasons']))
        if not 0<=self.contrast<=8 or not 0<=self.target_loss<=.25:raise ValueError('Invalid intrinsic target')
    def occupied(self,points):
        iy,ix,inside=region_indices(points,self.support.shape)
        return inside&self.support[iy,ix]
    def wall_values(self,points,f):
        iy,ix,inside=region_indices(points,self.support.shape)
        mask=inside&self.support[iy,ix]
        e=np.broadcast_to(np.asarray(self.epsinf,float),self.support.shape)
        loss=np.broadcast_to(np.asarray(self.ell,float),self.support.shape)
        value=np.ones(len(points),complex)
        im=loss[iy,ix]/f
        # At the exact legal upper endpoint, avoid a division/multiplication
        # roundoff excursion above .25. One representational ULP, no jitter.
        too_high=f*im>.25
        im[too_high]=np.nextafter(im[too_high],0.)
        value[mask]=e[iy[mask],ix[mask]]+1j*im[mask]
        return value
    def coefficient(self,f,n,include_target=True):
        def value(points):
            a=self.wall_values(points,f)
            if include_target:
                target=polygon_inside(points,self.pose,self.radius,n)
                a[target]=1+self.contrast+1j*self.target_loss/f
            return a
        return value
    def material_directions(self,points,f,n,kind='shared',include_target=True):
        occupied=self.occupied(points)
        if kind=='shared':
            d=np.column_stack((occupied.astype(complex),1j*occupied/f))
        elif kind=='left_right':
            d=np.column_stack((occupied&(points[:,0]<0),occupied&(points[:,0]>=0))).astype(complex)
        else:raise ValueError(kind)
        target=polygon_inside(points,self.pose,self.radius,n) if include_target else np.zeros(len(points),bool)
        return np.column_stack((d,target.astype(complex),1j*target/f))
    def metadata(self):
        return dict(support=self.support.astype(int).tolist(),epsinf=np.asarray(self.epsinf).tolist(),ell=np.asarray(self.ell).tolist(),
                    target_pose=self.pose.tolist(),target_intrinsic=dict(radius=self.radius,contrast=self.contrast,loss=self.target_loss))

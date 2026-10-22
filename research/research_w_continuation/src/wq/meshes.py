"""Interface-conforming polygon meshes for qualification, not inverse input.

All boundaries sample the *same* continuous frozen worlds. Polygon tolerance is
refined independently from the field mesh. No truth geometry is exposed to a
future witness optimizer. The material coefficient is still evaluated through
the original continuous callable, not a truth label inferred from mesh region.
"""
import numpy as np

def _rdp(points,tol):
    if len(points)<3:return points
    a,b=points[0],points[-1];v=b-a
    if np.dot(v,v)==0:
        idx=np.argmax(np.linalg.norm(points-a,axis=1))
        return np.vstack((_rdp(points[:idx+1],tol)[:-1],_rdp(points[idx:],tol))) if idx not in (0,len(points)-1) else points[[0,-1]]
    t=np.clip((points-a)@v/np.dot(v,v),0,1);d=np.linalg.norm(points-a-t[:,None]*v,axis=1);i=np.argmax(d)
    if d[i]<=tol:return points[[0,-1]]
    return np.vstack((_rdp(points[:i+1],tol)[:-1],_rdp(points[i:],tol)))

def wall_loops(seed,tol=1e-4,resolution=1):
    import contourpy
    # Implicit equations exactly duplicate frozen wall_values support/cavity.
    rng=np.random.default_rng(seed);p=rng.uniform(0,2*np.pi,5)
    x=np.linspace(-1.65,1.65,1320*resolution+1);y=np.linspace(.3,1.45,576*resolution+1)
    xx,yy=np.meshgrid(x,y)
    low=.64+.09*xx+.055*np.sin(2.4*xx+p[0]);thick=.48+.085*np.cos(2.1*xx+p[1])
    angle=np.arctan2(yy-.92,xx+.20);rad=np.hypot((xx+.20)/.25,(yy-.92)/.15)
    phi=np.minimum.reduce([1.48+.08*np.sin(4*yy+p[2])-np.abs(xx),yy-low,low+thick-yy,rad-1-.18*np.sin(3*angle+p[3])])
    paths=contourpy.contour_generator(x=x,y=y,z=phi).lines(0.)
    loops=[]
    for a in paths:
        if np.linalg.norm(a[0]-a[-1])>1e-9:raise RuntimeError('Open material boundary')
        # Split closed path to avoid a degenerate RDP baseline.
        j=len(a)//2;b=np.vstack((_rdp(a[:j+1],tol)[:-1],_rdp(a[j:],tol)))
        loops.append(b[:-1])
    return loops

def circle_loop(center,radius,n):
    a=2*np.pi*np.arange(n)/n
    return np.asarray(center)+radius*np.c_[np.cos(a),np.sin(a)]

def interface_mesh(*,h,outer_box=(-3.5,3.5,-1.8,4.5),physical_box=(-2.5,2.5,-.8,3.5),loops=()):
    import triangle
    from skfem import MeshTri
    vertices=[];segments=[]
    def add(loop):
        start=len(vertices);vertices.extend(np.asarray(loop).tolist());n=len(loop)
        segments.extend((start+i,start+(i+1)%n) for i in range(n))
    def rect(box):
        a,b,c,d=box;return np.array([[a,c],[b,c],[b,d],[a,d]])
    add(rect(outer_box));add(rect(physical_box));add(rect((-2.,2.,0.,3.)))
    for loop in loops:add(loop)
    result=triangle.triangulate(dict(vertices=np.asarray(vertices),segments=np.asarray(segments,dtype=np.int32)),f'pq28a{h*h/2:.12g}Q')
    if 'triangles' not in result:raise RuntimeError('Meshing failed')
    return MeshTri(result['vertices'].T,result['triangles'].T)

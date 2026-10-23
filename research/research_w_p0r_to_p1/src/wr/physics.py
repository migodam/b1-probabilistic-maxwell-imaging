"""Shared full Maxwell evaluation; no oracle/world-generator imports.

Only a caller-supplied coefficient model and a public mesh are used. Each
frequency's field/current/prediction is saved before releasing its factor.
No factor strong references are stored in archives or cross-frequency lists.
"""
import gc,time,json
import numpy as np
from skfem import asm
from wq.fem import FEM,_forcing,PMLConfig
from wq.inverse_support import PixelGrid,pixel_mesh,pixel_material,volume_pixel_jacobian
from .accounting import ROOT,Book,write

def pack(a):
    a=np.asarray(a).ravel();return np.r_[a.real,a.imag]/.01
def packj(a): return np.vstack([a.real,a.imag])/.01

def coefficient_jacobian(p,s,derivative):
    """Exact discrete volume output chain rule for a differentiable decoder."""
    nr=len(p.rx); ns=s.field.shape[1]
    p.book.reserve(p.phase,'fem_general_receiver_transpose',rhs=nr,factorizations=0,freq=p.freq)
    rhs=np.empty((len(p._free_dofs),nr),complex,order='F')
    for r in range(nr):
        load=(p.k**2*(s.epsilon-1)*p._receiver_green[r]).reshape(p._material_shape)
        rhs[:,r]=asm(_forcing,p.material_basis,load=load)[p._free_dofs]
    v=s.factor.solve(rhs,trans='T')
    defect=np.linalg.norm(s._matrix.T@v-rhs,axis=0)/np.maximum(np.linalg.norm(rhs,axis=0),np.finfo(float).tiny)
    if defect.max()>p.residual_tolerance:raise RuntimeError('receiver-transpose residual')
    p._event('backend_solve',label='general_receiver_transpose',rhs=nr,transpose='T',residual=float(defect.max()))
    dofs=np.zeros((p.ndofs,nr),complex); dofs[p._free_dofs]=v
    recv=p._receiver_green.T+np.column_stack([np.asarray(p.material_basis.interpolate(dofs[:,r])).ravel() for r in range(nr)])
    d=derivative(s.material_points,p.freq)
    j=np.empty((ns*nr,d.shape[1]),complex)
    for source in range(ns):
        j[source*nr:(source+1)*nr]=(p.k**2*(s.material_weights*s.field[:,source])[:,None]*recv).T@d
    return j,float(defect.max())

class PixelModel:
    def __init__(self,grid,indices=None):
        self.grid=grid; self.scale=np.r_[np.full(grid.n_pixels,8.),np.full(grid.n_pixels,.15)]
        self.indices=np.arange(grid.n_parameters) if indices is None else np.asarray(indices,int)
        self.n=len(self.indices)
    def full(self,x):
        q=np.zeros(self.grid.n_parameters);q[self.indices]=np.asarray(x)*self.scale[self.indices];return q
    def coefficient(self,x,f):return pixel_material(self.full(x),f,self.grid)
    def jacobian(self,p,s,x):
        return volume_pixel_jacobian(p,s,self.grid)[:,self.indices]*self.scale[self.indices][None,:],0.
    def material_metadata(self,x):return dict(model='pixel',grid=[self.grid.ny,self.grid.nx],indices=self.indices,physical_parameters=self.full(x))

class Forward:
    def __init__(self,model,mesh,frequencies,book,folder,*,phase,identity,h=.18,order=3,pml_strength=4,max_dofs=60000,max_factor_bytes=1500000000):
        self.model=model;self.mesh=mesh;self.freq=list(frequencies);self.book=book;self.folder=ROOT/'results'/folder
        self.folder.mkdir(parents=True,exist_ok=True);self.phase=phase;self.identity=identity
        self.h=h;self.order=order;self.pm=PMLConfig(strength=pml_strength)
        self.max_dofs=max_dofs; self.max_factor_bytes=max_factor_bytes
        self.tx=np.c_[np.linspace(-1.4,1.4,4),np.full(4,-.3)]
        self.rx=np.c_[np.linspace(-1.5,1.5,12),np.full(12,-.3)]
        self.saved={};self.calls=0;self.last_jac=None
    def evaluate(self,x,jac=False):
        x=np.asarray(x,float);key=x.tobytes()
        if key in self.saved:
            prior=self.saved[key]
            if not jac or prior[1] is not None:return prior[0].copy(),None if prior[1] is None else prior[1].copy()
        call=self.calls;self.calls+=1;dest=self.folder/f'endpoint_{call:04d}';dest.mkdir(exist_ok=False)
        write(dest/'PARAMETERS.json',dict(x=x,identity=self.identity,jacobian_requested=jac,frequencies=self.freq,phase=self.phase,metadata=self.model.material_metadata(x)))
        data=[];js=[];logs=[]
        for idx,f in enumerate(self.freq):
            p=FEM(f,self.tx,self.rx,self.book,h=self.h,order=self.order,mesh=self.mesh,pml=self.pm,phase=self.phase,receiver_mode='volume',max_dofs=self.max_dofs,max_factor_bytes=self.max_factor_bytes)
            s=p.solve(self.model.coefficient(x,f));data.append(s.data.copy())
            # Save endpoint before optional derivative or any subsequent solve.
            np.savez_compressed(dest/f'f{idx}.npz',x=x,frequency=f,data=s.data,probe_data=s.probe_data,epsilon=s.epsilon,points=s.material_points,weights=s.material_weights,field=s.field,current=s.current,scattered_dofs=s.scattered_dofs,mesh_p=p.mesh.p,mesh_t=p.mesh.t,residuals=s.residuals)
            if jac:
                j,adjres=self.model.jacobian(p,s,x);js.append(j)
                np.save(dest/f'J_f{idx}.npy',j)
            else:adjres=None
            logs.append(dict(frequency=f,residual=s.residual,adjoint_residual=adjres,metadata=s.metadata,timings=s.timings))
            write(dest/'FREQUENCIES.json',logs)
            del p,s;gc.collect()
        d=np.stack(data);j=np.concatenate(js,axis=0) if jac else None
        np.savez_compressed(dest/'RESPONSE.npz',x=x,data=d,**({'jacobian':j} if j is not None else {}))
        write(dest/'COMPLETE.json',dict(identity=self.identity,frequency_count=len(self.freq),jacobian=jac,totals=self.book.totals()))
        # At most two response/J arrays; never retain LU/field references.
        self.saved[key]=(d.copy(),None if j is None else j.copy(),dest)
        if len(self.saved)>2:self.saved.pop(next(iter(self.saved)))
        return d,j
    def residual(self,y):
        def ev(x,jac=False):
            d,j=self.evaluate(x,jac=jac)
            return pack(d-y),None if j is None else packj(j)
        return ev

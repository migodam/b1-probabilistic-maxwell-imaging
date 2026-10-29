"""R3 common-mesh scalar TM full-wave model, exact moving polygon geometry.

Reuses read-only R2 P3 FEM/PML and analytic incident scattered-field forcing.
Target permittivity is 1+contrast+i*loss/f. Wall is eps+i*ell/f. Ellipse
semiaxes and angle are independent continuous coordinates. All six public
wall cells remain active derivative coordinates even at ambient material.
"""
from dataclasses import dataclass
import numpy as np
from r2 import physics as old
from r2.physics import S1_TX,RX16,TargetQuadrature,TargetAssembler,line_green
from bgga.fixed_mesh import TARGET_POLYGON_N,_ccw,_clip_convex,_area,_cross
from .accounting import fem_lock
NAMES_TARGET=('px','py','rx','ry','contrast','loss','angle')
def default_config():
    return dict(version='R3_CONTINUOUS_FULL_MACRO_3X2_V1',wall_x=[-1.5,-.5,.5,1.5],wall_y=[.5,.875,1.25],frequencies=[.85,1.,1.15],sigma=.01,h=.17,targetstep=.10,wall_refine=1,target_degree=8,phase='R3',pixel_nx=241,pixel_ny=181)
def parameter_spec(config=None):
    c=default_config() if config is None else config;n=(len(c['wall_x'])-1)*(len(c['wall_y'])-1)
    return dict(names=[f'wall_eps_{i}' for i in range(n)]+[f'wall_loss_{i}' for i in range(n)]+list(NAMES_TARGET),lower=[1.]*n+[0.]*n+[-.9,1.7,.07,.07,.5,0.,-np.pi/2],upper=[6.]*n+[.25]*n+[.9,2.5,.18,.18,8.,.25,np.pi/2],scale=[1.]*n+[.08]*n+[.5,.35,.06,.06,3.,.1,1.],initial=[3.]*n+[.06]*n+[0.,2.1,.125,.105,3.,.06,0.])
@dataclass(frozen=True,eq=False)
class PublicWall(old.WallGrid):
    def __post_init__(self):
        x=old._axis(self.x,old.WALL_BOX[:2],'wall x');y=old._axis(self.y,old.WALL_BOX[2:],'wall y');shape=(len(y)-1,len(x)-1)
        e=old._real_array(self.epsinf,shape,'epsinf');l=old._real_array(self.ell,shape,'ell')
        if np.any((e<1.)|(e>6.)) or np.any((l<0.)|(l>.25)):raise ValueError('R3 wall bounds [1,6],[0,.25]')
        for k,v in [('x',x),('y',y),('epsinf',e),('ell',l)]:object.__setattr__(self,k,old._readonly(v))
    @property
    def support(self):return np.ones(self.shape,bool)
@dataclass(frozen=True,eq=False)
class Ellipse(old.Target):
    # inherited radius=rx, aspect=ry/rx, orientation continuous
    def __post_init__(self):
        p=np.asarray(self.pose,float);rx=float(self.radius);ry=rx*float(self.aspect)
        if p.shape!=(2,) or not np.isfinite(p).all() or not (-1<=p[0]<=1 and 1.6<=p[1]<=2.6):raise ValueError('target public ROI')
        if not(.06<=rx<=.2 and .06<=ry<=.2):raise ValueError('R3 independent semiaxes [.06,.2]')
        if not(0<=self.contrast<=8 and 0<=self.loss<=.25) or not np.isfinite(self.orientation):raise ValueError('target material/angle')
        object.__setattr__(self,'pose',old._readonly(p));object.__setattr__(self,'shape','ellipse')
    @property
    def rx(self):return self.radius
    @property
    def ry(self):return self.radius*self.aspect
    def rotation(self):
        c,s=np.cos(self.orientation),np.sin(self.orientation);return np.array([[c,-s],[s,c]])
    def transform(self):return self.rotation()@np.diag([1.,self.aspect])
    def metadata(self):return dict(pose=self.pose.tolist(),rx=self.rx,ry=self.ry,contrast=self.contrast,epsilon_real=1+self.contrast,loss=self.loss,angle=self.orientation,shape='ellipse',polygon_vertices=TARGET_POLYGON_N)
class EllipseAssembler(TargetAssembler):
    """Affine polygon cut integral with arbitrary aspect/rotation."""
    def set_geometry(self,target):self.geometry=target
    def quadrature(self,pose,radius):
        g=self.geometry;center=g.pose;vertices=g.target_vertices();extent=np.sum(np.abs(g.radius*g.transform()),axis=1)
        near=np.asarray(sorted(self._tree.query_ball_point(center,max(g.rx,g.ry)+self._max_radius)),dtype=np.int64)
        if len(near):near=near[np.all(self._upper[near]>center-extent,axis=1)&np.all(self._lower[near]<center+extent,axis=1)]
        pieces=[];touched=[];count=0;inv=np.linalg.inv(radius*g.transform());incircle=np.cos(np.pi/TARGET_POLYGON_N)
        for local in near:
            tri=_ccw(self._triangles[local]);polygon=tri if np.all(np.linalg.norm((tri-center)@inv.T,axis=1)<=incircle) else _clip_convex(vertices,tri)
            if _area(polygon)<=0:continue
            element=int(self.roi_elements[local]);found=False
            for j in range(1,len(polygon)-1):
                a,b,c=polygon[0],polygon[j],polygon[j+1];det=float(_cross(b-a,c-a))
                if det<=0:continue
                count+=len(self.reference_weights)
                if count>self.max_points or count*32>self.max_work_bytes:raise MemoryError('R3 ellipse cut quadrature guard')
                pieces.append((a.copy(),b.copy(),c.copy(),det,element));found=True
            if found:touched.append(element)
        if not pieces:raise RuntimeError('empty ellipse cut')
        q=np.empty((count,2));w=np.empty(count);els=np.empty(count,dtype=np.int64);nq=len(self.reference_weights)
        for j,(a,b,c,det,e) in enumerate(pieces):
            sl=slice(j*nq,(j+1)*nq);q[sl]=a+self.reference_points[:,:1]*(b-a)+self.reference_points[:,1:]*(c-a);w[sl]=det*self.reference_weights;els[sl]=e
        area=.5*TARGET_POLYGON_N*g.rx*g.ry*np.sin(2*np.pi/TARGET_POLYGON_N)
        return TargetQuadrature(q,w,els,np.asarray(touched,dtype=np.int64),float(w.sum()),float(area),vertices,self.reference_points.copy(),self.reference_weights.copy(),self.degree,len(near))
class Eval(old.Eval):
    def _assembler(self,target):
        if not hasattr(self,'r3ass'):self.r3ass=EllipseAssembler(self.p.basis,self.p._free_dofs,**self._target_settings)
        self.r3ass.set_geometry(target);return self.r3ass
    def _target_names(self,names):
        if any(n not in NAMES_TARGET for n in names):raise ValueError('unknown target derivative')
        return tuple(names)
    def target_differentials(self,state,*,names=NAMES_TARGET):
        self._state(state);p=self.p;t=state.world;out={n:np.zeros_like(state.u) for n in names};q=state.target.quadrature
        if set(names)&{'contrast','loss'}:
            volume=np.zeros_like(state.u)
            for j in range(0,len(q.weights),self.ass.chunk_points):
                sl=slice(j,min(j+self.ass.chunk_points,len(q.weights)));probe=self.ass._known_element_probes(q.points[sl],q.elements[sl],p._free_dofs);total=line_green(q.points[sl],p.tx,p.k)+probe@state.u
                volume+=p.k**2*np.asarray(probe.T@(q.weights[sl,None]*total))
            if 'contrast' in out:out['contrast']=volume
            if 'loss' in out:out['loss']=1j/p.freq*volume
        shape=set(names)-{'contrast','loss'}
        if shape:
            v=t.target_vertices();edge=np.roll(v,-1,axis=0)-v;length=np.linalg.norm(edge,axis=1);nodes,ww=np.polynomial.legendre.leggauss(6);alpha=(nodes+1)/2
            points=(v[:,None,:]+edge[:,None,:]*alpha[None,:,None]).reshape(-1,2);normal=np.repeat(np.column_stack((edge[:,1],-edge[:,0]))/length[:,None],len(nodes),axis=0);ds=(length[:,None]*ww[None,:]/2).ravel()
            local=(points-t.pose)@t.rotation();rot=t.rotation();drx=local[:,0,None]/t.rx*rot[:,0];dry=local[:,1,None]/t.ry*rot[:,1];rel=points-t.pose;dang=np.column_stack((-rel[:,1],rel[:,0]))
            vel={'px':normal[:,0],'py':normal[:,1],'rx':np.einsum('qi,qi->q',normal,drx),'ry':np.einsum('qi,qi->q',normal,dry),'angle':np.einsum('qi,qi->q',normal,dang)}
            probe=self.probe(points);total=line_green(points,p.tx,p.k)+probe@state.u;scale=p.k**2*(t.contrast+1j*t.loss/p.freq)
            for n in shape:out[n]=scale*np.asarray(probe.T@((ds*vel[n])[:,None]*total))
        return out
@dataclass
class Prediction:
    data:np.ndarray
    jacobian:object=None
    diagnostics:object=None
class JointModel:
    def __init__(self,config,book):
        self.config={**default_config(),**config};self.book=book;self.spec=parameter_spec(self.config);self.evaluators=[];self.charts=[];self.audit_dir=None;self.audit_count=0;self.last_prediction=None;self.last_chart_z=None;self.last_chart_prediction=None;c=self.config
        book.check(new_job=True)
        with fem_lock(book.root):
            book.check(new_job=True)
            for f in c['frequencies']:
                self.evaluators.append(Eval(book,c.get('phase','R3'),f,S1_TX,c['wall_x'],c['wall_y'],wall_refine=c.get('wall_refine',1),h=c['h'],targetstep=c['targetstep'],target_degree=c.get('target_degree',8),max_factor_bytes=1500000000))
    def decode(self,z):
        z=np.asarray(z,float);n=(len(self.config['wall_x'])-1)*(len(self.config['wall_y'])-1)
        if z.shape!=(2*n+7,) or not np.isfinite(z).all():raise ValueError('R3 parameter shape')
        wall=PublicWall(self.config['wall_x'],self.config['wall_y'],z[:n].reshape(-1,len(self.config['wall_x'])-1),z[n:2*n].reshape(-1,len(self.config['wall_x'])-1));x,y,rx,ry,chi,loss,ang=z[2*n:]
        return wall,Ellipse([x,y],rx,chi,loss,'ellipse',ry/rx,ang)
    def initialize_charts(self,rank=16):
        from .charts import ForcingChart
        z=np.asarray(self.spec['initial']);z[-3: -1]=0.;wall,target=self.decode(z)
        with fem_lock(self.book.root):
            for ev in self.evaluators:
                st=ev.evaluate(wall,target,label='R3 public reference chart state');ch=ForcingChart(ev,st,rank=rank);self.charts.append(ch)
                if self.audit_dir is not None:
                    self.audit_dir.mkdir(parents=True,exist_ok=True);np.savez_compressed(self.audit_dir/f'REFERENCE_f{ev.p.freq:g}.npz',z=z,**ch.reference_audit)
    def audit_chart_step(self,z,direction,prediction,tolerance=.1):
        wall,target=self.decode(z);events=[];full=[];bad=False
        with fem_lock(self.book.root):
            for i,(ev,ch) in enumerate(zip(self.evaluators,self.charts)):
                st=ev.evaluate(wall,target,label='R3 adaptive base state');J=prediction.jacobian[i*64:(i+1)*64];event,F=ch.audit(st,direction,J,tolerance);events.append(event);full.append(F);bad|=F is not None
            if bad:
                for i,(ev,F) in enumerate(zip(self.evaluators,full)):
                    if F is None:
                        st=ev.evaluate(wall,target);A,B=ev.jacobians(st,target_names=NAMES_TARGET);full[i]=np.column_stack((A,B))
                fullJ=np.vstack(full)
                if self.audit_dir is not None:
                    self.audit_count+=1;np.savez_compressed(self.audit_dir/f'CURRENT_AUDIT_{self.audit_count:05d}.npz',z=z,data=prediction.data,exact_jacobian=fullJ,approximate_jacobian=prediction.jacobian,proposed_direction=direction,parameter_scales=np.asarray(self.spec['scale']))
                return events,Prediction(prediction.data,fullJ,prediction.diagnostics)
        return events,None
    def forward(self,z,jacobian=False,label='R3 full objective'):
        wall,target=self.decode(z);data=[];jac=[];diag=[]
        with fem_lock(self.book.root):
            for i,ev in enumerate(self.evaluators):
                st=ev.evaluate(wall,target,label=label);data.append(st.data.copy());diag.append(dict(residual=st.residual,timings=st.timings,ndofs=ev.p.ndofs))
                if jacobian:
                    if jacobian in ('FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART'):
                        jac.append(self.charts[i].jacobian(st))
                    else:
                        A,B=ev.jacobians(st,target_names=NAMES_TARGET);jac.append(np.column_stack((A,B)))
            prediction=Prediction(np.asarray(data),np.vstack(jac) if jacobian else None,diag)
            if jacobian in ('FIXED_CURRENT_CHART','ADAPTIVE_CURRENT_CHART'):
                self.last_chart_z=np.asarray(z).copy();self.last_chart_prediction=prediction
            return prediction
    def material_maps(self,z,x=None,y=None):
        c=self.config;x=np.linspace(-1.5,1.5,c['pixel_nx']) if x is None else np.asarray(x);y=np.linspace(.4,2.8,c['pixel_ny']) if y is None else np.asarray(y);xx,yy=np.meshgrid(x,y);points=np.column_stack((xx.ravel(),yy.ravel()));wall,t=self.decode(z);eps=wall.values(points,1.).real;loss=wall.values(points,1.).imag;loc=(points-t.pose)@t.rotation();inside=(loc[:,0]/t.rx)**2+(loc[:,1]/t.ry)**2<=1.;eps[inside]=1+t.contrast;loss[inside]=t.loss
        return dict(pixel_x=x,pixel_y=y,image_eps=eps.reshape(xx.shape),image_loss=loss.reshape(xx.shape),target_mask=inside.reshape(xx.shape))

"""Fake dense matrices only: no real Model/Book or Green kernels are constructed."""
import contextlib,json,sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.linalg import lu_factor
P=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(P/'src'))
from physics import Model
class FakeBook:
    def __init__(self):self.operations=[]
    @contextlib.contextmanager
    def op(self,kind,**kw):
        self.operations.append(dict(kind=kind,**kw));yield len(self.operations)
    def save(self,*args,**kw):pass
    def append(self,row):pass
rng=np.random.default_rng(901672)
n,ns,nr=6,2,2
G=.025*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
S=rng.normal(size=(nr,n))+1j*rng.normal(size=(nr,n))
E0=rng.normal(size=(n,ns))+1j*rng.normal(size=(n,ns))
m=np.linspace(.8,2.8,n)
L=np.eye(n)-m[:,None]*G
j=np.linalg.solve(L,m[:,None]*E0)
E=E0+G@j
state=dict(m=m,j=j,E=E)
fake=SimpleNamespace(n=n,ns=ns,nr=nr,state=state,factor=lu_factor(L),S=S,book=FakeBook())
fake.pack_jac=lambda J:Model.pack_jac(fake,J)
J=Model.jacobian(fake,state)
def response(x):
    current=np.linalg.solve(np.eye(n)-x[:,None]*G,x[:,None]*E0)
    return Model.pack(fake,S@current)
errors={}
for step in [1e-4,1e-5,1e-6]:
    fd=np.column_stack([(response(m+step*np.eye(n)[i])-response(m-step*np.eye(n)[i]))/(2*step) for i in range(n)])
    errors[str(step)]=float(np.linalg.norm(fd-J)/np.linalg.norm(J))
T_reference=S@np.linalg.inv(L)
J_reference=np.einsum('ri,is->sri',T_reference,E).reshape(ns*nr,n)
transpose_gap=float(np.linalg.norm(state['Jcomplex']-J_reference)/np.linalg.norm(J_reference))
label_y=np.array([[11+12j,21+22j],[13+14j,23+24j]])
packing=Model.pack(fake,label_y).tolist()
expected=[11,12,13,14,21,22,23,24]
wrong_H=S@np.linalg.inv(L).conj().T
transpose_vs_hermitian_gap=float(np.linalg.norm(wrong_H-T_reference)/np.linalg.norm(T_reference))
result=dict(scope='Algebraic fake matrices only; zero real Model/Book construction and zero new physics RHS',
    seed=901672,dimensions=dict(n=n,sources=ns,receivers=nr),
    central_difference_relative_errors=errors,explicit_inverse_reference_gap=transpose_gap,
    analytic_transpose_vs_wrong_hermitian_gap=transpose_vs_hermitian_gap,
    pack_observed=packing,pack_expected=expected,
    fakebook_operations=fake.book.operations,
    gate=bool(max(errors.values())<1e-7 and transpose_gap<1e-13 and packing==expected and len(fake.book.operations)==1 and fake.book.operations[0]['rhs']==nr))
(P/'audit/solver_review/ALGEBRAIC_QA.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
from scipy.linalg import null_space
from charts import CurrentChart,partition
# A separate fake system has four real responses and two free real materials.
n,ns,nr=6,2,1
G=.025*(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
S=rng.normal(size=(nr,n))+1j*rng.normal(size=(nr,n))
E0=rng.normal(size=(n,ns))+1j*rng.normal(size=(n,ns))
m=np.linspace(.8,2.8,n);L=np.eye(n)-m[:,None]*G
j=np.linalg.solve(L,m[:,None]*E0);E=E0+G@j;y=S@j
state=dict(m=m,j=j,E=E,y=y)
fake=SimpleNamespace(n=n,ns=ns,nr=nr,state=state,factor=lu_factor(L),S=S,E0=E0,G=G,N=null_space(S),book=FakeBook())
fake.pack_jac=lambda J:Model.pack_jac(fake,J)
J=Model.jacobian(fake,state);dep,free=partition(J)
chart=CurrentChart(fake,y,state,dep,free);z=chart.encode(j,m)
r0,D,(m0,j0,E00)=chart.residual_jac(z,m[free])
current_errors={}
for step in [1e-4,1e-5,1e-6]:
    eye=np.eye(len(z))
    fd=np.column_stack([(chart.residual_jac(z+step*v,m[free],False)[0]-chart.residual_jac(z-step*v,m[free],False)[0])/(2*step) for v in eye])
    current_errors[str(step)]=float(np.linalg.norm(fd-D)/np.linalg.norm(D))
tan=chart.material_tangent(z,m[free]);tan_ref=np.zeros((n,len(free)));tan_ref[free]=np.eye(len(free));tan_ref[dep]=-np.linalg.solve(J[:,dep],J[:,free])
current_result=dict(scope=result['scope'],central_difference_relative_errors=current_errors,
    center_residual=float(np.linalg.norm(r0)),receiver_compatibility_gap=float(np.linalg.norm(S@j0-y)),
    encode_decode_material_gap=float(np.linalg.norm(m0-m)),encode_decode_current_gap=float(np.linalg.norm(j0-j)),
    material_is_real=bool(np.isrealobj(m0)),material_tangent_relative_gap=float(np.linalg.norm(tan-tan_ref)/np.linalg.norm(tan_ref)),
    constitutive_equations=D.shape[0],real_variables=D.shape[1],free_material_coordinates=free.tolist(),
    gate=bool(max(current_errors.values())<1e-7 and np.linalg.norm(r0)<1e-13 and np.linalg.norm(tan-tan_ref)/np.linalg.norm(tan_ref)<1e-12 and np.isrealobj(m0)))
(P/'audit/solver_review/CURRENT_ALGEBRAIC_QA.json').write_text(json.dumps(current_result,indent=2)+'\n')
print(json.dumps(current_result,indent=2))

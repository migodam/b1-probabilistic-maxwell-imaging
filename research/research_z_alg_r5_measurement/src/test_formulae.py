from common import *
import certification as c
from experiment import quotient,solve_trf
from scipy.special import j0
book=Ledger('formula_tests');checks={};certificates=[]
for n in (4,6):
 g=base_geometry(n);pos,Q=dictionary('test',g);p=Physics(g,pos,book);m=np.linspace(2,8,n);old=np.eye(p.L)[:,:p.s];y,J=p.response(m,old,True)
 q=quotient(p,y);sys=q['raw']['system'];w=p.dual(m);cc=(w-sys.w0)@np.linalg.pinv(sys.N);x=np.r_[cc,cc.conj()];f,j=sys.evaluate(x,True)
 h=1e-6;fd=np.column_stack([(sys.evaluate(x+np.eye(n)[i]*h)-sys.evaluate(x-np.eye(n)[i]*h))/(2*h) for i in range(n)])
 checks[f'polynomial_{n}']=float(np.max(abs(f)));checks[f'poly_jac_{n}']=float(np.linalg.norm(j-fd)/np.linalg.norm(j))
 cert=c.certify(p.H,p.B[:,:p.s],p.R,y.reshape(1,-1),m);certificates.append(cert);checks[f'cert_{n}']=cert['status']=='CERTIFIED_LOCAL_REAL_ROOT'
 # Directed interval point and full-J must enclose independently evaluated values to float rounding tolerance.
 import mpmath as mp
 with c._interval_precision(128):
  X=[mp.iv.mpf([float(v-1e-7),float(v+1e-7)]) for v in m]
  F,JI,_,_=c._interval_forward_jacobian(X,p.H.tolist(),p.B[:,:p.s].tolist(),p.R.tolist(),y.tolist())
  checks[f'interval_J_{n}']=all(float(JI[i][k].a)<=J[i,k]<=float(JI[i][k].b) for i in range(n) for k in range(n))
 theta=np.arange(2048)*2*np.pi/2048;phase=np.exp(-1j*g['k']*(np.c_[np.cos(theta),np.sin(theta)]@pos.T));Qnum=phase.conj().T@phase/len(theta)
 checks[f'power_gram_{n}']=float(np.max(abs(Q-Qnum)))
 # Node: G=0 and one receiver component zero makes that material unobservable for every source.
 G=np.zeros((n,n),complex);S=np.ones((1,n),complex);S[0,0]=0;E=p.B[:,:p.s];yn=(S@(m[:,None]*E));mn=m.copy();mn[0]+=1
 checks[f'node_invariant_{n}']=float(np.linalg.norm(S@(mn[:,None]*E)-yn))
 cn=c.certify(G,E,S,yn,m);checks[f'node_unknown_{n}']=cn['status']=='UNKNOWN';certificates.append(cn)
 # Pole enclosure rejects singular A.
 cp=c.certify(np.eye(n,dtype=complex),E,S,yn,np.ones(n));checks[f'pole_unknown_{n}']=cp['status']=='UNKNOWN';certificates.append(cp)
passed=all(v if isinstance(v,bool) else v<1e-6 for v in checks.values())
save(ROOT/'results/formula_tests.json',dict(passed=passed,checks=checks,certificates=certificates,receipt=book.receipt()));print(checks);assert passed

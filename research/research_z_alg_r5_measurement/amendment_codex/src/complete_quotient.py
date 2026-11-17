"""Codex-authored exact-chart interval quotient coverage prototype.

Scope: original rounded discrete Maxwell constants, not rounded coefficients.
No GT is used. Bounds/roots are independently certified with interval arithmetic.
"""
import sys,time,json,itertools,math
from pathlib import Path
BASE=Path(__file__).resolve().parents[2];sys.path.insert(0,str(BASE/'src'))
from common import *
import certification as oldcert
import mpmath as mp
from scipy.linalg import qr
from src.general_quotient_solver import maximal_minors,multiplication_matrices,homogeneous_monomials
OUT=Path(__file__).resolve().parents[1]

def ivc(z):return mp.iv.mpc(float(complex(z).real),float(complex(z).imag))
def mid(z):return complex(float(z.real.mid),float(z.imag.mid))
def conj(z):return mp.iv.mpc(z.real,-z.imag)
def upper(x):return max(abs(x.a),abs(x.b))
def interval_json(z):
 if hasattr(z,'imag') and isinstance(z,mp.iv.mpc):return {'real':interval_json(z.real),'imag':interval_json(z.imag)}
 return {'lo':float(np.nextafter(float(z.a),-np.inf)),'hi':float(np.nextafter(float(z.b),np.inf))}
def addmat(a,b):return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]
def mm(a,b):return [[sum((a[i][k]*b[k][j] for k in range(len(b))),mp.iv.mpc(0)) for j in range(len(b[0]))] for i in range(len(a))]
def transpose(a):return [list(r) for r in zip(*a)]

def chart(H,B,R,t,data_radius=0.):
 n,s=B.shape;_,_,piv=qr(B.T,pivoting=True);I=list(map(int,piv[:s]));J=[i for i in range(n) if i not in I]
 BI=[[ivc(B[i,j]) for j in range(s)] for i in I];BJ=[[ivc(B[i,j]) for j in range(s)] for i in J]
 inv,pivots=oldcert._interval_inverse(BI)
 ti=[mp.iv.mpc(mp.iv.mpf([float(np.nextafter(z.real-data_radius,-np.inf)) if data_radius else float(z.real),float(np.nextafter(z.real+data_radius,np.inf)) if data_radius else float(z.real)]),mp.iv.mpf([float(np.nextafter(z.imag-data_radius,-np.inf)) if data_radius else float(z.imag),float(np.nextafter(z.imag+data_radius,np.inf)) if data_radius else float(z.imag)])) for z in t]
 wp=mm([ti],inv)[0];Npart=mm(BJ,inv)
 W=[[mp.iv.mpc(0) for _ in range(n)] for _ in range(s+1)]
 for j,i in enumerate(I):W[0][i]=wp[j]
 for l in range(s):
  for j,i in enumerate(I):W[l+1][i]=-Npart[l][j]
  W[l+1][J[l]]=mp.iv.mpc(1)
 V=mm(W,[[ivc(v) for v in row] for row in H])
 for i in range(n):V[0][i]+=ivc(R.reshape(-1)[i])
 coeff=[]
 for i in range(n):
  ci=[[W[a][i]*conj(V[b][i])-V[a][i]*conj(W[b][i]) for b in range(s+1)] for a in range(s+1)]
  scale=max(abs(mid(v)) for row in ci for v in row)
  if scale==0:raise ValueError('ZERO_EQUATION')
  coeff.append([[v/float(scale) for v in row] for row in ci])
 return dict(I=I,J=J,W=W,V=V,coeff=coeff)

def interval_minors(C):
 n=len(C);s=len(C[0])-1;width=s+1
 basis=[e for degree in range(s+1) for e in homogeneous_monomials(s,degree)];leading=homogeneous_monomials(s,s+1);mons=basis+leading;out=[]
 for rows in itertools.combinations(range(n),width):
  poly={}
  for perm in itertools.permutations(range(width)):
   sign=(-1)**sum(perm[i]>perm[j] for i in range(width) for j in range(i+1,width));term={(0,)*s:mp.iv.mpc(sign)}
   for row,col in zip(rows,perm):
    new={}
    for e,v in term.items():
     for axis in range(width):
      changed=list(e)
      if axis:changed[axis-1]+=1
      key=tuple(changed);new[key]=new.get(key,mp.iv.mpc(0))+v*C[row][axis][col]
    term=new
   for e,v in term.items():poly[e]=poly.get(e,mp.iv.mpc(0))+v
  out.append([poly.get(e,mp.iv.mpc(0)) for e in mons])
 return out,basis,leading

def invertibility(A):
 midpoint=np.array([[mid(z) for z in row] for row in A]);C=np.linalg.inv(midpoint);CA=mm([[ivc(z) for z in row] for row in C],A)
 norm=[]
 for i,row in enumerate(CA):norm.append(sum((upper((mp.iv.mpc(i==j)-z).real)+upper((mp.iv.mpc(i==j)-z).imag) for j,z in enumerate(row)),mp.iv.mpf(0)))
 bound=max(norm)
 return dict(certified=bool(bound<1),preconditioned_residual_norm_upper=float(bound.b),midpoint_condition=float(np.linalg.cond(midpoint)),dimension=len(A))

def evaluate(C,x):
 n=len(C);s=n//2;c=[mp.iv.mpc(1)]+x[:s];d=[mp.iv.mpc(1)]+x[s:]
 F=[sum((c[a]*d[b]*ci[a][b] for a in range(s+1) for b in range(s+1)),mp.iv.mpc(0)) for ci in C]
 J=[[sum((d[b]*ci[a][b] for b in range(s+1)),mp.iv.mpc(0)) for a in range(1,s+1)]+[sum((c[a]*ci[a][b] for a in range(s+1)),mp.iv.mpc(0)) for b in range(1,s+1)] for ci in C]
 return F,J

def evalfloat(C,x):
 s=len(x)//2;c=np.r_[1,x[:s]];d=np.r_[1,x[s:]];F=np.einsum('a,b,iab->i',c,d,C)
 J=np.column_stack([np.einsum('b,ib->i',d,C[:,a,:]) for a in range(1,s+1)]+[np.einsum('a,ia->i',c,C[:,:,b]) for b in range(1,s+1)])
 return F,J

def enumerate_midpoint(C):
 cf=np.array([[[mid(z) for z in row] for row in ci] for ci in C]);s=len(cf)//2
 minors,basis,leading=maximal_minors(cf);mats,diag=multiplication_matrices(minors,basis,leading,condition_limit=1e16)
 sep=[1,.371+.193j,-.227+.419j];linear=sum(c*m for c,m in zip(sep,mats));eig,vec=np.linalg.eig(linear.T);roots=[]
 for v in vec.T:
  c=np.array([np.vdot(v,M.T@v)/np.vdot(v,v) for M in mats]);A=np.einsum('a,iab->ib',np.r_[1,c],cf);d=np.linalg.lstsq(A[:,1:],-A[:,0],rcond=1e-14)[0];x=np.r_[c,d]
  for it in range(12):
   F,J=evalfloat(cf,x);dx=np.linalg.solve(J,-F);x+=dx
   if np.max(abs(dx))<1e-14:break
  physical=np.linalg.norm(x[s:]-x[:s].conj())<1e-7
  if physical:x=np.r_[x[:s],x[:s].conj()]
  roots.append((x,physical))
 return roots,dict(float_leading_condition=diag['leading_block_condition'],float_commutator=diag['maximum_relative_commutator'])

def root_certificate(C,x,physical,radii=(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,.1,1.)):
 start=time.process_time();n=len(x);xp=[ivc(z) for z in x];F0,J0=evaluate(C,xp);Jmid=np.array([[mid(z) for z in row] for row in J0]);pre=np.linalg.inv(Jmid);preiv=[[ivc(z) for z in row] for row in pre];CF=mm(preiv,[[f] for f in F0]);trials=[]
 for radius in radii:
  radius*=max(1.,np.max(abs(x)));X=[];domains=[]
  for z in x:
   re=[float(np.nextafter(z.real-radius,-np.inf)),float(np.nextafter(z.real+radius,np.inf))];im=[float(np.nextafter(z.imag-radius,-np.inf)),float(np.nextafter(z.imag+radius,np.inf))];X.append(mp.iv.mpc(mp.iv.mpf(re),mp.iv.mpf(im)));domains.append({'real':{'lo':re[0],'hi':re[1]},'imag':{'lo':im[0],'hi':im[1]}})
  _,JX=evaluate(C,X);CJ=mm(preiv,JX);M=[[mp.iv.mpc(i==j)-v for j,v in enumerate(row)] for i,row in enumerate(CJ)];K=[xp[i]-CF[i][0]+sum((M[i][j]*(X[j]-xp[j]) for j in range(n)),mp.iv.mpc(0)) for i in range(n)]
  norms=[sum((upper(v.real)+upper(v.imag) for v in row),mp.iv.mpf(0)) for row in M];norm=max(norms)
  inside=all(k.real.a>z.real.a and k.real.b<z.real.b and k.imag.a>z.imag.a and k.imag.b<z.imag.b for k,z in zip(K,X));good=bool(inside and norm<1)
  trials.append(dict(radius=radius,strict_inclusion=bool(inside),contraction_upper=float(norm.b)))
  if good:return dict(status='CERTIFIED_COMPLEX_ROOT',X=X,domain=domains,K=[interval_json(z) for z in K],conjugate_invariant_box=physical,trials=trials,CPU_seconds=time.process_time()-start)
 return dict(status='UNKNOWN',trials=trials,CPU_seconds=time.process_time()-start)

def nonzero(z):return bool(z.real.a>0 or z.real.b<0 or z.imag.a>0 or z.imag.b<0)
def disjoint(A,B):return any(a.real.b<b.real.a or b.real.b<a.real.a or a.imag.b<b.imag.a or b.imag.b<a.imag.a for a,b in zip(A,B))

def classify(ch,cr,H,R,box):
 if cr['status']!='CERTIFIED_COMPLEX_ROOT':return dict(status='UNKNOWN')
 X=cr['X'];s=len(X)//2;n=2*s
 # Exclude another d over the same c: full column rank is necessary for physical classification.
 cvec=[mp.iv.mpc(1)]+X[:s];C=ch['coeff']
 recovery=[[sum((cvec[a]*C[i][a][b] for a in range(s+1)),mp.iv.mpc(0)) for b in range(1,s+1)] for i in range(n)]
 midrec=np.array([[mid(z) for z in row] for row in recovery]);_,_,rp=qr(midrec.T,pivoting=True)
 rankproof=invertibility([recovery[int(i)] for i in rp[:s]])
 cr['recovery_rank_certificate']=rankproof
 if not rankproof['certified']:return dict(status='UNKNOWN_RECOVERY_RANK')
 if cr['conjugate_invariant_box'] and not all(X[s+j].real.a==X[j].real.a and X[s+j].real.b==X[j].real.b and X[s+j].imag.a==-X[j].imag.b and X[s+j].imag.b==-X[j].imag.a for j in range(s)):
  return dict(status='UNKNOWN_CONJUGATE_BOX_SYMMETRY')
 nonreal=any(nonzero(X[s+j]-conj(X[j])) for j in range(s))
 if nonreal:return dict(status='EXCLUDED_NONCONJUGATE')
 if not cr['conjugate_invariant_box']:return dict(status='UNKNOWN_CONJUGACY')
 w=mm([[mp.iv.mpc(1)]+X[:s]],ch['W'])[0];v=mm([[mp.iv.mpc(1)]+X[:s]],ch['V'])[0]
 if not all(nonzero(z) for z in v):return dict(status='UNKNOWN_NODE',v=[interval_json(z) for z in v])
 m=[a/b for a,b in zip(w,v)]
 if any(z.real.b<box[0] or z.real.a>box[1] for z in m):return dict(status='EXCLUDED_OUTSIDE_BOX',material=[interval_json(z.real) for z in m])
 if not all(z.real.a>box[0] and z.real.b<box[1] for z in m):return dict(status='UNKNOWN_BOX_BOUNDARY',material=[interval_json(z.real) for z in m])
 A=[[mp.iv.mpc(i==j)-m[i].real*ivc(H[i,j]) for j in range(n)] for i in range(n)];pole=invertibility(A)
 return dict(status='CERTIFIED_LEGAL_REAL_MATERIAL' if pole['certified'] else 'UNKNOWN_POLE',material=[interval_json(z.real) for z in m],pole=pole,w=[interval_json(z) for z in w],v=[interval_json(z) for z in v])

def run(case='T0',arm='clean',bits=128):
 mp.iv.prec=bits;start=time.process_time();wall=time.monotonic();book=Ledger('complete_quotient_'+case+'_'+arm);g=geometry(case);pos,Q=dictionary(case,g);p=Physics(g,pos,book);y=np.load(BASE/'inputs'/f'{case}.npz')['old_'+arm]
 ch=chart(p.H,p.B[:,:p.s],p.R,y);polys,basis,leading=interval_minors(ch['coeff']);lead=invertibility([r[len(basis):] for r in polys]);roots,diag=enumerate_midpoint(ch['coeff']);certs=[];classifications=[]
 for x,phys in roots:
  cr=root_certificate(ch['coeff'],x,phys);cls=classify(ch,cr,p.H,p.R,BOX);certs.append(cr);classifications.append(cls)
 allcert=all(c['status']=='CERTIFIED_COMPLEX_ROOT' for c in certs);distinct=all(disjoint(a['X'][:p.s],b['X'][:p.s]) for a,b in itertools.combinations(certs,2)) if allcert else False
 quotient_complete=bool(lead['certified'] and len(certs)==len(basis) and allcert and distinct)
 material_complete=quotient_complete and all(v['status'] in ('EXCLUDED_NONCONJUGATE','EXCLUDED_OUTSIDE_BOX','CERTIFIED_LEGAL_REAL_MATERIAL') for v in classifications)
 result=dict(case=case,arm=arm,bits=bits,scope='All regular real materials in [.5,30]^k for original rounded G/E/S and exact central binary64 old observation; no continuous PDE/Hankel error or noise feasible coverage.',leading=lead,chart_indices=ch['I'],quotient_upper_dimension=len(basis),distinct_c_boxes=distinct,quotient_complete=quotient_complete,material_box_complete=material_complete,legal_count=sum(v['status']=='CERTIFIED_LEGAL_REAL_MATERIAL' for v in classifications),classifications=classifications,certificates=[{k:v for k,v in cr.items() if k!='X'} for cr in certs],diagnostic=diag,CPU_seconds=time.process_time()-start,wall_seconds=time.monotonic()-wall,physics_receipt=book.receipt())
 save(OUT/'results'/f'complete_{case}_{arm}_{bits}.json',result);print({k:v for k,v in result.items() if k not in ('classifications','certificates','physics_receipt')},flush=True)
 return result
if __name__=='__main__':run(*(sys.argv[1:] or ['T0','clean']))

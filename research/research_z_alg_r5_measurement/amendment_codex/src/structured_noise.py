"""Codex: preserve data correlation before interval Krawczyk preconditioning.
Coefficient family C(t) is exact quadratic in the real/imag old data. Leading
maximal-minor coefficients are affine because their sole data-dependent column
is column zero. No noise reduction, no target change, no GT inputs.
"""
from complete_quotient import *

def tensor_zero(n,s):return [[[mp.iv.mpc(0) for b in range(s+1)] for a in range(s+1)] for i in range(n)]
def tensor_add(A,B):return [[[x+y for x,y in zip(ar,br)] for ar,br in zip(ai,bi)] for ai,bi in zip(A,B)]
def bilinear(W,V):
 s=len(W)-1;n=len(W[0]);return [[[W[a][i]*conj(V[b][i])-V[a][i]*conj(W[b][i]) for b in range(s+1)] for a in range(s+1)] for i in range(n)]
def polar(W,V,dW,dV):
 s=len(W)-1;n=len(W[0]);return [[[dW[a][i]*conj(V[b][i])+W[a][i]*conj(dV[b][i])-dV[a][i]*conj(W[b][i])-V[a][i]*conj(dW[b][i]) for b in range(s+1)] for a in range(s+1)] for i in range(n)]
def scale_tensor(A,sc):return [[[v/float(sc[i]) for v in row] for row in ai] for i,ai in enumerate(A)]
def family(H,B,R,t):
 ch=chart(H,B,R,t);W=ch['W'];V=ch['V'];n,s=B.shape;raw=bilinear(W,V);sc=[max(abs(mid(v)) for row in ai for v in row) for ai in raw];linear=[];DW=[];DV=[]
 inv,_=oldcert._interval_inverse([[ivc(B[i,j]) for j in range(s)] for i in ch['I']])
 for l in range(2*s):
  unit=[mp.iv.mpc(0) for _ in range(s)];unit[l//2]=ivc(1 if l%2==0 else 1j);row=mm([unit],inv)[0]
  dW=[[mp.iv.mpc(0) for _ in range(n)] for _ in range(s+1)]
  for j,i in enumerate(ch['I']):dW[0][i]=row[j]
  dV=mm(dW,[[ivc(z) for z in row] for row in H]);DW.append(dW);DV.append(dV);linear.append(scale_tensor(polar(W,V,dW,dV),sc))
 quadratic=[]
 for l in range(2*s):
  for k in range(l,2*s):
   if l==k:C=bilinear(DW[l],DV[l])
   else:C=polar(DW[l],DV[l],DW[k],DV[k])
   quadratic.append((l,k,scale_tensor(C,sc)))
 return ch,scale_tensor(raw,sc),linear,quadratic

def pre_tensor(P,C):
 n=len(C);s=len(C[0])-1
 return [[[sum((P[i][j]*C[j][a][b] for j in range(n)),mp.iv.mpc(0)) for b in range(s+1)] for a in range(s+1)] for i in range(n)]
def envelope(C,linear,quad,eta,P=None):
 if P is not None:
  C=pre_tensor(P,C);linear=[pre_tensor(P,d) for d in linear];quad=[(l,k,pre_tensor(P,d)) for l,k,d in quad]
 n=len(C);s=len(C[0])-1;delta=mp.iv.mpf([-eta,eta]);square=mp.iv.mpf([0,eta])*mp.iv.mpf([0,eta]);cross=delta*delta
 return [[[C[i][a][b]+sum((D[i][a][b]*delta for D in linear),mp.iv.mpc(0))+sum((D[i][a][b]*(square if l==k else cross) for l,k,D in quad),mp.iv.mpc(0)) for b in range(s+1)] for a in range(s+1)] for i in range(n)]

def affine_leading_certificate(C,L,eta):
 n=len(C);s=n//2;polys,basis,lead=interval_minors(C);T=[r[len(basis):] for r in polys];Tmid=np.array([[mid(v) for v in row] for row in T]);P=[[ivc(v) for v in row] for row in np.linalg.inv(Tmid)];PT=mm(P,T);der=[]
 for d in L:
  # Highest c degree: a>=1 only. Column b=0 is the sole data-dependent column.
  F=[[[mp.iv.mpc(0) if a==0 else (d[i][a][b] if b==0 else C[i][a][b]) for b in range(s+1)] for a in range(s+1)] for i in range(n)]
  pm,_,_=interval_minors(F);dT=[r[len(basis):] for r in pm];der.append(mm(P,dT))
 delta=mp.iv.mpf([-eta,eta]);M=[[mp.iv.mpc(i==j)-PT[i][j]-sum((D[i][j]*delta for D in der),mp.iv.mpc(0)) for j in range(len(T))] for i in range(len(T))]
 bound=max(sum((upper(v.real)+upper(v.imag) for v in row),mp.iv.mpf(0)) for row in M)
 absbound=np.array([[float((upper(v.real)+upper(v.imag)).b) for v in row] for row in M]);ev,vec=np.linalg.eig(absbound);weights=np.abs(vec[:,np.argmax(ev.real)].real)+1e-15
 weighted=max(sum(((upper(M[i][j].real)+upper(M[i][j].imag))*mp.iv.mpf(float(weights[j])) for j in range(len(T))),mp.iv.mpf(0))/mp.iv.mpf(float(weights[i])) for i in range(len(T)))
 return dict(certified=bool(min(bound,weighted)<1),preconditioned_residual_norm_upper=float(bound.b),weighted_residual_norm_upper=float(weighted.b),dimension=len(T),explanation='Exact affine leading block with data dependence combined before interval bounding')

def preconditioned_parametric_root(C,L,Q,x,physical,eta):
 n=len(x);s=n//2;xp=[ivc(z) for z in x];_,J0=evaluate(C,xp);P=[[ivc(v) for v in row] for row in np.linalg.inv(np.array([[mid(z) for z in row] for row in J0]))]
 PC=pre_tensor(P,C);PL=[pre_tensor(P,d) for d in L];PQ=[(l,k,pre_tensor(P,d)) for l,k,d in Q]
 Fc,Jc=evaluate(PC,xp);Lvals=[evaluate(d,xp) for d in PL];Qvals=[(l,k,evaluate(d,xp)[0]) for l,k,d in PQ]
 delta=mp.iv.mpf([-eta,eta]);square=mp.iv.mpf([0,eta])*mp.iv.mpf([0,eta]);cross=delta*delta
 # Evaluate and precondition exact Taylor coefficients BEFORE taking data ranges.
 F0=[Fc[i]+sum((f[i]*delta for f,j in Lvals),mp.iv.mpc(0))+sum((f[i]*(square if l==k else cross) for l,k,f in Qvals),mp.iv.mpc(0)) for i in range(n)]
 displacement=max(float(max(upper(z.real),upper(z.imag)).b) for z in F0)
 radii=[max(1e-10,displacement)*v for v in (1.05,1.25,1.5,2,3,4,6,8,12,16)];trials=[]
 for radius in radii:
  X=[];domains=[]
  for z in x:
   re=[float(np.nextafter(z.real-radius,-np.inf)),float(np.nextafter(z.real+radius,np.inf))];im=[float(np.nextafter(z.imag-radius,-np.inf)),float(np.nextafter(z.imag+radius,np.inf))];X.append(mp.iv.mpc(mp.iv.mpf(re),mp.iv.mpf(im)));domains.append({'real':{'lo':re[0],'hi':re[1]},'imag':{'lo':im[0],'hi':im[1]}})
  h=[X[j]-xp[j] for j in range(n)];JX=[]
  for i in range(n):
   row=[]
   for j in range(n):
    # The original F is quadratic jointly in x and real/imag data, hence exact affine Jacobian.
    shift=sum((PC[i][j+1][b+1]*h[s+b] for b in range(s)),mp.iv.mpc(0)) if j<s else sum((PC[i][a+1][j-s+1]*h[a] for a in range(s)),mp.iv.mpc(0))
    row.append(Jc[i][j]+sum((jl[i][j]*delta for fl,jl in Lvals),mp.iv.mpc(0))+shift)
   JX.append(row)
  M=[[mp.iv.mpc(i==j)-v for j,v in enumerate(row)] for i,row in enumerate(JX)];K=[xp[i]-F0[i]+sum((M[i][j]*h[j] for j in range(n)),mp.iv.mpc(0)) for i in range(n)];norm=max(sum((upper(v.real)+upper(v.imag) for v in row),mp.iv.mpf(0)) for row in M);inside=all(k.real.a>z.real.a and k.real.b<z.real.b and k.imag.a>z.imag.a and k.imag.b<z.imag.b for k,z in zip(K,X));trials.append(dict(radius=radius,strict_inclusion=bool(inside),contraction_upper=float(norm.b)))
  if inside and norm<1:return dict(status='CERTIFIED_COMPLEX_ROOT',X=X,domain=domains,K=[interval_json(z) for z in K],conjugate_invariant_box=physical,trials=trials)
 return dict(status='UNKNOWN',trials=trials)

def run_structured(case='T0',arm='noisy',eta=1.2e-4):
 mp.iv.prec=128;begin=time.process_time();book=Ledger('structured_noise');g=geometry(case);pos,Qpower=dictionary(case,g);p=Physics(g,pos,book);y=np.load(BASE/'inputs'/f'{case}.npz')['old_'+arm];ch,C,L,Q=family(p.H,p.B[:,:p.s],p.R,y);lead=affine_leading_certificate(C,L,eta);roots,_=enumerate_midpoint(C);certs=[preconditioned_parametric_root(C,L,Q,x,physical,eta) for x,physical in roots];allcert=all(c['status']=='CERTIFIED_COMPLEX_ROOT' for c in certs);distinct=all(disjoint(a['X'][:p.s],b['X'][:p.s]) for a,b in itertools.combinations(certs,2)) if allcert else False
 chunc=chart(p.H,p.B[:,:p.s],p.R,y,data_radius=eta);classes=[classify(chunc,cr,p.H,p.R,BOX) for cr in certs]
 out=dict(case=case,arm=arm,eta=eta,leading=lead,certified_root_count=sum(c['status']=='CERTIFIED_COMPLEX_ROOT' for c in certs),total_roots=len(certs),distinct_c_boxes=distinct,noise_feasible_root_coverage=bool(lead['certified'] and allcert and distinct),classifications=classes,certificates=[{k:v for k,v in c.items() if k!='X'} for c in certs],CPU_seconds=time.process_time()-begin,physics_receipt=book.receipt(),scope='Original rounded discrete model, uniform old-data component box; model-error assumption still requires independent validity.')
 save(OUT/'results'/f'structured_noise_{case}_{arm}.json',out);print({k:v for k,v in out.items() if k not in ('classifications','certificates','physics_receipt')},flush=True);return out
if __name__=='__main__':run_structured(*(sys.argv[1:] or ['T0','noisy']))

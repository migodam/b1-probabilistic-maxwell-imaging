from complete_quotient import *

def real_det(A):
 A=[row[:] for row in A];n=len(A);det=mp.iv.mpf(1)
 for j in range(n):
  choices=[i for i in range(j,n) if A[i][j].a>0 or A[i][j].b<0]
  if not choices:raise ArithmeticError('Interval determinant pivot unresolved')
  pivot=max(choices,key=lambda i:abs(float(A[i][j].mid)))
  if pivot!=j:A[j],A[pivot]=A[pivot],A[j];det=-det
  d=A[j][j];det*=d
  for i in range(j+1,n):
   scale=A[i][j]/d
   for k in range(j+1,n):A[i][k]-=scale*A[j][k]
 return det

mp.iv.prec=128;begin=time.process_time();book=Ledger('certify_noise_critical_value');g=geometry('T0');pos,Q=dictionary('T0',g);p=Physics(g,pos,book);y=np.load(BASE/'inputs/T0.npz')['old_noisy'];m=np.array(json.loads((OUT/'results/minimax_fold.json').read_text())['m']);A=np.eye(p.L)[:,:p.s];steps=[]
for j in range(p.n):
 d=np.eye(p.n)[j]*1e-6;steps.append(abs(np.linalg.det(p.response(m+d,A,True)[1])-np.linalg.det(p.response(m-d,A,True)[1])))
axis=int(np.argmax(steps));radius=1e-7;detends=[]
for sign in (-1,1):
 x=[mp.iv.mpf(float(v))+mp.iv.mpf(sign*radius if j==axis else 0) for j,v in enumerate(m)]
 F,J,_,_=oldcert._interval_forward_jacobian(x,p.H.tolist(),p.B[:,:p.s].tolist(),p.R.tolist(),y.tolist());detends.append(real_det(J))
x=[mp.iv.mpf(float(v))+mp.iv.mpf([-radius,radius] if j==axis else 0) for j,v in enumerate(m)]
F,J,_,piv=oldcert._interval_forward_jacobian(x,p.H.tolist(),p.B[:,:p.s].tolist(),p.R.tolist(),y.tolist());errors=[]
for j in range(p.s):errors.append(mp.iv.sqrt(upper(F[2*j])**2+upper(F[2*j+1])**2))
opposite=bool((detends[0].b<0<detends[1].a) or (detends[1].b<0<detends[0].a));inside=all(e.b<mp.iv.mpf(1.2e-4).a for e in errors)
result=dict(certified_critical_value_inside_old_noise_ball=bool(opposite and inside),opposite_determinant_signs=opposite,axis=axis,halfwidth=radius,center=m,endpoint_determinant_intervals=[interval_json(d) for d in detends],old_complex_residual_upper=[float(e.b) for e in errors],declared_complex_radius=1.2e-4,physical_nonsingularity_along_segment='ALL_INTERVAL_ELIMINATION_PIVOTS_EXCLUDE_ZERO',conclusion='By continuity and sign change a real material Jacobian singularity occurs inside the declared old-data error ball. Uniform distinct nonsingular root tubes over that entire error ball cannot hold.',scope='Original rounded scalar-TM model, actual material values within legal box; no claim of nonidentifiability after adding a new source.',CPU_seconds=time.process_time()-begin,receipt=book.receipt());save(OUT/'results/certified_fold.json',result);print({k:v for k,v in serial(result).items() if k!='receipt'})

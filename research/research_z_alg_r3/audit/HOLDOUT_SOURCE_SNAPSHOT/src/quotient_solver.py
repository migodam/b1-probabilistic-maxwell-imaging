"""Six-dimensional quotient algebra for four bilinear Maxwell equations.

Mathematical prototype by the principal investigator.  No private GT access.
All construction is on observed data and fixed physical operators.  Failure of
regular quotient conditions is explicit, never silently reported as coverage.
"""
import itertools,time
import numpy as np
try:
 from .branch_enumerator import BilinearDualSystem
except ImportError:
 from branch_enumerator import BilinearDualSystem
BASIS=((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
CUBICS=((3,0),(2,1),(1,2),(0,3))
MONOMIALS=BASIS+CUBICS

class QuotientFailure(ValueError):
 def __init__(self,message,diagnostic=None):
  super().__init__(message);self.diagnostic=diagnostic or {}

def _multiply(p,q):
 r={}
 for a,x in p.items():
  for b,y in q.items():
   e=(a[0]+b[0],a[1]+b[1]);r[e]=r.get(e,0)+x*y
 return r

def maximal_minors(coeff):
 """Four exact-in-form maximal minors of A(c)_{ib}=sum_a C_iab c_a."""
 polynomials=[]
 for omit in range(4):
  rows=[i for i in range(4) if i!=omit];out={}
  for perm in itertools.permutations(range(3)):
   parity=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
   p={(0,0):complex(parity)}
   for row,col in zip(rows,perm):
    p=_multiply(p,{(0,0):coeff[row,0,col],(1,0):coeff[row,1,col],(0,1):coeff[row,2,col]})
   for e,z in p.items():out[e]=out.get(e,0)+z
  polynomials.append([out.get(e,0) for e in MONOMIALS])
 return np.asarray(polynomials,complex)

def multiplication_matrices(minors):
 scales=np.max(np.abs(minors),axis=1)
 if np.any(scales==0):raise QuotientFailure('DEGENERATE_ZERO_MAXIMAL_MINOR')
 p=minors/scales[:,None];T=p[:,6:];R=p[:,:6]
 condition=float(np.linalg.cond(T))
 if not np.isfinite(condition) or condition>1e12:raise QuotientFailure('UNRESOLVED_CUBIC_LEADING_BLOCK',dict(leading_block_condition=condition))
 reduction=np.linalg.solve(T,-R)
 matrices=[]
 for axis in range(2):
  M=np.zeros((6,6),complex)
  for col,b in enumerate(BASIS):
   exponent=list(b);exponent[axis]+=1;exponent=tuple(exponent)
   if exponent in BASIS:M[BASIS.index(exponent),col]=1
   else:M[:,col]=reduction[CUBICS.index(exponent)]
  matrices.append(M)
 Mx,My=matrices
 commutator=float(np.linalg.norm(Mx@My-My@Mx)/max(np.linalg.norm(Mx)*np.linalg.norm(My),1e-300))
 return Mx,My,dict(leading_block_condition=condition,commutator_relative=commutator,reduction=reduction,minors=minors)

def enumerate_quotient(kernel,target,separating=.371+.193j,tolerance=2e-10,max_polish=5):
 start=time.perf_counter();system=BilinearDualSystem(kernel,target);paths=[];roots=[]
 counts=dict(maximal_minors=4,polynomial_multiplications=72,quotient_linear_solves=1,quotient_eigenproblems=1,
             quotient_recovery_least_squares=0,quotient_polish_linear_solves=0,quotient_condition_SVDs=0)
 with kernel.model.operation('QUOTIENT_ALGEBRA',maximal_minors=4,polynomial_multiplications=72,
                             quotient_linear_solves=1,quotient_eigenproblems=1):
  minors=maximal_minors(system.coefficients);Mx,My,diagnostic=multiplication_matrices(minors)
  counts['quotient_condition_SVDs']+=1
  if diagnostic['commutator_relative']>1e-7:raise QuotientFailure('UNRESOLVED_NONCOMMUTING_MULTIPLICATION',diagnostic)
  eig,vectors=np.linalg.eig((Mx+separating*My).T)
  counts['quotient_condition_SVDs']+=1
  gaps=[abs(eig[i]-eig[j])/(1+abs(eig[i])+abs(eig[j])) for i in range(6) for j in range(i)]
  diagnostic.update(eigenvector_condition=float(np.linalg.cond(vectors)),multiplication_x=Mx,multiplication_y=My,separating=separating)
  diagnostic['minimum_relative_separating_eigengap']=float(min(gaps))
  if diagnostic['eigenvector_condition']>1e12 or min(gaps)<=1e-10:
   raise QuotientFailure('UNRESOLVED_NONSEPARATING_OR_ILL_CONDITIONED_EIGENPROBLEM',diagnostic)
  for i,v in enumerate(vectors.T):
   c=np.array([np.vdot(v,Mx.T@v),np.vdot(v,My.T@v)])/np.vdot(v,v)
   A=np.einsum('a,iab->ib',np.r_[1,c],system.coefficients)
   d,_,rank,recovery_s=np.linalg.lstsq(A[:,1:],-A[:,0],rcond=1e-12)
   recovery_condition=float(recovery_s[0]/max(recovery_s[-1],1e-300))
   counts['quotient_recovery_least_squares']+=1
   x=np.r_[c,d];polish=0
   for k in range(max_polish):
    f,J=system.evaluate(x,jac=True)
    if np.linalg.norm(f,np.inf)<=tolerance:break
    try:
     counts['quotient_polish_linear_solves']+=1
     step=np.linalg.solve(J,-f)
    except np.linalg.LinAlgError:break
    x+=step;polish+=1
   f,J=system.evaluate(x,jac=True);residual=float(np.linalg.norm(f,np.inf));cond=float(np.linalg.cond(J))
   counts['quotient_condition_SVDs']+=1
   status='RESOLVED_FINITE_NUMERICAL' if rank==2 and recovery_condition<1e12 and residual<=tolerance and cond<1e12 else 'UNRESOLVED_QUOTIENT_ENDPOINT'
   path=dict(path=i,status=status,x=x,endpoint_residual=residual,endpoint_condition=cond,polish_steps=polish,recovery_rank=int(rank),recovery_condition=recovery_condition,separating_eigenvalue=eig[i])
   paths.append(path)
   if status=='RESOLVED_FINITE_NUMERICAL':roots.append(x.copy())
 distinct=[]
 for root in roots:
  if all(np.linalg.norm(root-r)>1e-6*(1+np.linalg.norm(root)) for r in distinct):distinct.append(root)
 elapsed=time.perf_counter()-start
 counts.update(polynomial_evaluations=system.evaluations,polynomial_J_evaluations=system.jacobian_evaluations)
 with kernel.model.operation('QUOTIENT_ENDPOINT_SUMMARY',**{k:v for k,v in counts.items() if k not in ('maximal_minors','polynomial_multiplications','quotient_linear_solves','quotient_eigenproblems')}):pass
 return dict(method='BILINEAR_QUOTIENT_SIX_EIGENROOTS',paths=paths,roots=np.asarray(roots).reshape(-1,4),resolved_paths=len(roots),unresolved_paths=6-len(roots),distinct_finite_endpoints=len(distinct),all_paths_numerically_resolved=len(roots)==len(distinct)==6,completeness='NOT_RIGOROUSLY_CERTIFIED',certificate_flag='NOT_RIGOROUSLY_CERTIFIED',system=system,diagnostic=diagnostic,wall_seconds=elapsed,counters=counts)

from complete_quotient import *
from scipy.optimize import minimize

def run_fold():
 begin=time.process_time();book=Ledger('nearest_fold_two_start_diagnostic');g=geometry('T0');pos,Q=dictionary('T0',g);p=Physics(g,pos,book);y=np.load(BASE/'inputs/T0.npz')['old_noisy'];old=np.eye(p.L)[:,:p.s];line=json.loads((OUT/'results/fold_diagnostic.json').read_text());raw=json.loads((BASE/'results/T0/noisy/Q.json').read_text());a,b=map(np.array,raw['accepted_roots']);starts0=[np.array(line['critical_material']),(a+b)/2];rows=[]
 def obj(m):
  yp,J=p.response(m,old,True);r=pack(yp-y)/1e-4;return float(r@r),2*J.T@r/1e-4
 def con(m):return float(np.linalg.det(p.response(m,old,True)[1])/1e-17)
 for m in starts0:
  res=minimize(obj,m,jac=True,constraints=[{'type':'eq','fun':con}],bounds=[BOX]*p.n,method='SLSQP',options={'maxiter':400,'ftol':1e-12})
  yp,J=p.response(res.x,old,True);rows.append(dict(m=res.x,objective=res.fun,success=res.success,message=res.message,nit=res.nit,complex_residual=abs(yp-y),determinant=np.linalg.det(J),sigma_min=np.linalg.svd(J,compute_uv=False)[-1]))
 best=min(rows,key=lambda r:r['objective']);save(OUT/'results/nearest_fold.json',dict(rows=rows,best=best,CPU_seconds=time.process_time()-begin,receipt=book.receipt(),scope='Post-result diagnostic, no GT input, no new target, no proof until interval bracketing'))
 print(serial(best),flush=True)
if __name__=='__main__':run_fold()

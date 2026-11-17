from complete_quotient import *
from scipy.optimize import minimize
begin=time.process_time();book=Ledger('minimax_fold_diagnostic');g=geometry('T0');pos,Q=dictionary('T0',g);p=Physics(g,pos,book);y=np.load(BASE/'inputs/T0.npz')['old_noisy'];old=np.eye(p.L)[:,:p.s];m=np.array(json.loads((OUT/'results/nearest_fold.json').read_text())['best']['m']);x0=np.r_[m,max(abs(p.response(m,old)-y))/1e-4]
def constraint(x):
 yp,J=p.response(x[:-1],old,True);r=abs(yp-y)/1e-4;return x[-1]**2-r**2
def detcon(x):return float(np.linalg.det(p.response(x[:-1],old,True)[1])/1e-17)
res=minimize(lambda x:x[-1],x0,jac=lambda x:np.r_[np.zeros(p.n),1.],bounds=[BOX]*p.n+[(0.,100.)],constraints=[{'type':'eq','fun':detcon},{'type':'ineq','fun':constraint}],method='SLSQP',options={'maxiter':400,'ftol':1e-12})
yc,J=p.response(res.x[:-1],old,True);row=dict(m=res.x[:-1],complex_residual=yc-y,max_complex_error=max(abs(yc-y)),success=res.success,message=res.message,nit=res.nit,determinant=np.linalg.det(J),sigma_min=np.linalg.svd(J,compute_uv=False)[-1],CPU_seconds=time.process_time()-begin,receipt=book.receipt());save(OUT/'results/minimax_fold.json',row);print({k:v for k,v in serial(row).items() if k!='receipt'})

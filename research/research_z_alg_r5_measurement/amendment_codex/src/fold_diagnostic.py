from complete_quotient import *
from scipy.optimize import brentq

def diagnose():
 begin=time.process_time();book=Ledger('fold_diagnostic');g=geometry('T0');pos,Q=dictionary('T0',g);p=Physics(g,pos,book);y=np.load(BASE/'inputs/T0.npz')['old_noisy'];raw=json.loads((BASE/'results/T0/noisy/Q.json').read_text());m0,m1=map(np.array,raw['accepted_roots']);direction=m1-m0;A=np.eye(p.L)[:,:p.s];rows=[]
 for t in np.linspace(0,1,201):
  m=m0+t*direction;yp,J=p.response(m,A,True);rows.append(dict(t=t,residual_max=float(max(abs(yp-y))),jacobian_determinant=float(np.linalg.det(J)),sigma_min=float(np.linalg.svd(J,compute_uv=False)[-1])))
 def determinant(t):return np.linalg.det(p.response(m0+t*direction,A,True)[1])
 critical=brentq(determinant,0,1,xtol=1e-14) if rows[0]['jacobian_determinant']*rows[-1]['jacobian_determinant']<0 else None
 out=dict(sampled_path=rows,maximum_sampled_old_residual=max(r['residual_max'] for r in rows),endpoint_determinants=[rows[0]['jacobian_determinant'],rows[-1]['jacobian_determinant']],critical_line_parameter=critical,critical_material=(m0+critical*direction) if critical is not None else None,scope='Sampled diagnostic, not yet interval path proof',CPU_seconds=time.process_time()-begin,receipt=book.receipt())
 save(OUT/'results/fold_diagnostic.json',out);print({k:v for k,v in serial(out).items() if k not in ('sampled_path','receipt')})
if __name__=='__main__':diagnose()

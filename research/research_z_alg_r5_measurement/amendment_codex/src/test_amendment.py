"""Independent checks of exact chart, polynomial derivatives and source KKT."""
from complete_quotient import *
from structured_noise import family,bilinear,scale_tensor
from constrained_source import constrained_match

def array(C):return np.array([[[mid(v) for v in r] for r in ci] for ci in C])
def main():
 begin=time.process_time();mp.iv.prec=128;book=Ledger('amendment_formula_tests');g=geometry('T0');pos,Q=dictionary('T0',g);p=Physics(g,pos,book);B=p.B[:,:3];y=np.load(BASE/'inputs/T0.npz')['old_noisy'];ch,C,L,T=family(p.H,B,p.R,y);rng=np.random.default_rng(618);checks={}
 W=np.array([[mid(z) for z in row] for row in ch['W']]);checks['chart_nullspace_error']=float(np.max(abs(W[1:]@B)));checks['chart_particular_error']=float(np.max(abs(W[0]@B-y)))
 # Arbitrary real material, independently forward propagated, must satisfy chart and conjugacy polynomials.
 m=np.array([2.,4.,6.,8.,10.,12.]);ym=p.response(m,np.eye(p.L)[:,:3],independent=True);cm=chart(p.H,B,p.R,ym);w=p.dual(m);x=np.r_[w[cm['J']],w[cm['J']].conj()];F,J=evalfloat(array(cm['coeff']),x);checks['full_forward_polynomial_error']=float(np.max(abs(F)))
 d=rng.normal(size=6)+1j*rng.normal(size=6);h=1e-6;num=(evalfloat(array(C),x+h*d)[0]-evalfloat(array(C),x-h*d)[0])/(2*h);checks['polynomial_directional_derivative_relative_error']=float(np.linalg.norm(num-evalfloat(array(C),x)[1]@d)/np.linalg.norm(num))
 # Data Taylor family against independently rebuilt unscaled W,V at perturbed data (fixed original scaling).
 dt=rng.normal(size=6)*1e-5;yp=y+dt[::2]+1j*dt[1::2];cp=chart(p.H,B,p.R,yp);raw=bilinear(ch['W'],ch['V']);sc=[max(abs(mid(v)) for row in ci for v in row) for ci in raw];direct=array(scale_tensor(bilinear(cp['W'],cp['V']),sc));pred=array(C)+sum(v*array(t) for v,t in zip(dt,L))+sum(dt[l]*dt[k]*array(t) for l,k,t in T);checks['exact_data_quadratic_error']=float(np.max(abs(pred-direct)))
 # Verify affine leading block construction with actual maximal minors, fixed scale.
 def leading(c):pm,ba,le=maximal_minors(c);return np.asarray(pm)[:,len(ba):]
 # Historical maximal_minors returns dense coefficient matrix.
 t0=leading(array(C));t1=leading(direct);aff=t0.copy()
 for q,t in enumerate(L):
  dd=np.zeros_like(array(C));dd[:,1:,1:]=array(C)[:,1:,1:];dd[:,1:,0]=array(t)[:,1:,0];aff+=dt[q]*leading(dd)
 checks['leading_affine_error']=float(np.max(abs(t1-aff)))
 # KKT dual upper bound: for lambda>=0, |d a|² <= d (Q+lambda I)^-1 d* (1+lambda cap²).
 direction=rng.normal(size=p.L)+1j*rng.normal(size=p.L);sol=constrained_match(direction,Q);a=sol['a'];lam=sol['lagrange_ratio'];bound=np.real(direction@np.linalg.solve(Q+lam*np.eye(p.L),direction.conj()))*(1+16*lam);checks['source_KKT_relative_dual_gap']=float(abs(bound-abs(direction@a)**2)/bound);checks['source_power_error']=float(abs(np.vdot(a,Q@a)-1));checks['source_feed_excess']=float(max(0,np.linalg.norm(a)-4))
 limits={'chart_nullspace_error':1e-12,'chart_particular_error':1e-12,'full_forward_polynomial_error':1e-10,'polynomial_directional_derivative_relative_error':1e-7,'exact_data_quadratic_error':1e-12,'leading_affine_error':1e-12,'source_KKT_relative_dual_gap':1e-8,'source_power_error':1e-10,'source_feed_excess':1e-8};ok=all(checks[k]<v for k,v in limits.items());save(OUT/'results/amendment_formula_tests.json',dict(passed=ok,checks=checks,limits=limits,CPU_seconds=time.process_time()-begin,receipt=book.receipt()));print(checks);assert ok
if __name__=='__main__':main()

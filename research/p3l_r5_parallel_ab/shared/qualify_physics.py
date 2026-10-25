"""Counted qualification: three forwards, rank2 JVP, receiver J, adjoint."""
import time
import numpy as np
from .physics_adapter import Target,OUT,write

def run():
    if (OUT/'audit/PHYSICS_QUALIFICATION.json').exists():raise FileExistsError('Qualification already retained; no implicit paid replay')
    target=Target('A','qualification');rng=np.random.default_rng(20261026031)
    x=rng.normal(scale=.3,size=66);V=np.linalg.qr(rng.normal(size=(66,2)))[0]
    before=target.book.totals();start=time.perf_counter()
    p=target.evaluate(x,'qualification_base',retain_reference=True)
    J2=target.jacobian(p,V,purpose='qualification_rank2_jvp')
    J=target.jacobian(p,purpose='qualification_receiver32')
    y=p.prediction+target.sigma*(.1+.2j)
    gradient=target.gradient(p,y,purpose='qualification_adjoint8')
    eps=1e-5
    plus=target.evaluate(x+eps*V[:,0],'qualification_fd_plus')
    minus=target.evaluate(x-eps*V[:,0],'qualification_fd_minus')
    fd=target.white(plus.prediction-minus.prediction)/(2*eps)
    err_fd=float(np.linalg.norm(fd-J2[:,0])/np.linalg.norm(J2[:,0]))
    err_j=float(np.linalg.norm(J@V-J2)/np.linalg.norm(J2))
    g2=J.T@target.white(p.prediction-y)
    err_g=float(np.linalg.norm(gradient-g2)/max(1.,np.linalg.norm(g2)))
    after=target.book.totals()
    result=dict(status='PASS' if max(err_fd,err_j,err_g)<1e-6 else 'FAIL',finite_difference_relative_error=err_fd,direct_vs_receiver_J_relative_error=err_j,adjoint_vs_J_gradient_relative_error=err_g,before=before,after=after,delta={k:after[k]-before[k] for k in ('states','forward_rhs','derivative_rhs','all_rhs','cpu_seconds','wall_seconds')},benchmark_wall_seconds=time.perf_counter()-start,backward_residual=p.state.residual,numpy=np.__version__,seed=20261026031,scope='original discrete TM24 backend and P1 decoder derivative qualification; no posterior or full-continuum claim')
    write(OUT/'audit/PHYSICS_QUALIFICATION.json',result)
    assert result['status']=='PASS',result
if __name__=='__main__':run()

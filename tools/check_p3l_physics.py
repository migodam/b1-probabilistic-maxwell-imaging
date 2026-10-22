from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from p3l.core import OUT,Book,Target,write,to_probit,from_probit,probit_derivative

def main():
    b=Book('A','initial_physics_contract',max_wall=120);status='COMPLETE';error=None
    try:
        T=Target(b,keep_fields=True);y=T.observation('DEV_000');q=np.r_[np.zeros(64),.5,.5]
        p=T.evaluate(q,'fixed_prior_center_derivative_contract');J=T.jacobian(p);g=T.phi_gradient(p,y)
        expected=J.T@T.white(p.state.data-y);rng=np.random.default_rng(2026102209);v=rng.normal(size=66);v/=np.linalg.norm(v);h=1e-5
        pp=T.evaluate(q+h*v,'centered_FD_plus');pm=T.evaluate(q-h*v,'centered_FD_minus')
        fd=(T.phi(pp,y)-T.phi(pm,y))/(2*h);gg=g@v
        y2=y+(.1+.2j)*T.sigma;g2=T.phi_gradient(p,y2)
        g2expected=J.T@T.white(p.state.data-y2)
        x=to_probit(q);xp=x+h*v;xm=x-h*v
        ppx=T.evaluate(from_probit(xp),'probit_FD_plus');pmx=T.evaluate(from_probit(xm),'probit_FD_minus')
        fdx=(T.phi(ppx,y)+.5*xp@xp-T.phi(pmx,y)-.5*xm@xm)/(2*h)
        gx=(g*probit_derivative(x)+x)@v
        result=dict(source_gradient_vs_receiver_J_relative=float(np.linalg.norm(g-expected)/np.linalg.norm(expected)),source_gradient_direction_FD_relative=float(abs(fd-gg)/max(1.,abs(gg))),observation_cache_invalidation_error=float(np.linalg.norm(g2-g2expected)/np.linalg.norm(g2expected)),probit_logdensity_gradient_FD_relative=float(abs(fdx-gx)/max(1.,abs(gx))),backward_residual=p.state.residual,gradient_incremental_RHS=8,full_J_incremental_RHS=32,truth_used=False)
        result['PASS']=all(result[k]<1e-5 for k in ['source_gradient_vs_receiver_J_relative','source_gradient_direction_FD_relative','observation_cache_invalidation_error','probit_logdensity_gradient_FD_relative'])
        write(OUT/'AUDIT/PHYSICS_PROBABILITY_CONTRACT.json',result);print(json.dumps(result))
        if not result['PASS']:raise RuntimeError('initial physical derivative contract failed')
    except BaseException as e:status='ERROR';error=repr(e);raise
    finally:b.close(status,error)

if __name__=='__main__':main()

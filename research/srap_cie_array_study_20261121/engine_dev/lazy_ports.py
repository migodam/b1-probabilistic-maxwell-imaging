"""Separate development experiment comparing eager and on-demand port Schur solves."""
from collections import Counter
import json
import platform
import time
from pathlib import Path
from check_interventions import read_public, model_for, HISTORICAL, HERE, engine, legacy, np
from scipy import linalg as la
from scipy.optimize import brentq


class TimedEager(engine.BallProjector):
    def solve(self, *args, **kwargs):
        t = time.perf_counter()
        x = super().solve(*args, **kwargs)
        self.costs['inner_H_solve_s'] += time.perf_counter()-t
        return x


class OnDemand(TimedEager):
    """Graph-first projection. Active balls use per-source H-metric port Krylov
    spaces and one secular equation for the total discrepancy over all sources.
    The full-space shifted linear-system residual is checked before accepting.
    """
    def __init__(self,model,Y,delta,recycled):
        self.model,self.E,self.sources = model,model.E,model.ns
        self.delta,self.costs = delta,model.clock.counts
        self.Y = Y if model.ports.valid_mask is None else np.where(model.ports.valid_mask,Y,0.)
        self.cg_tolerance = 2e-10
        self.warm = None
        self._coarse = None
        gb = self.B(recycled)
        vals,rot = la.eigh(gb.conj().T@gb,check_finite=False)
        basis = recycled@rot
        self._coarse = basis, (vals/(1+vals))[:,None]*basis.conj().T
        self.last_A_info = {}

    def T(self,x):
        return 6*self.model.ports.receive(x)

    def TT(self,r):
        return 6*self.model.ports.adjoint(r)

    def affine_ball(self,x,q):
        xu = self.solve(x+self.BH(q-self.E),start=self.warm)
        self.warm=xu.copy()
        r = self.T(xu)-self.Y
        raw = np.linalg.norm(r)
        if raw <= self.delta:
            self.last_A_info = dict(multiplier=0.,krylov_dimension=0,raw_residual=float(raw))
            return xu,self.E+self.B(xu)
        rhs_norm = np.linalg.norm(r,axis=0)
        v = r/np.maximum(rhs_norm,1e-300)
        vectors, images, solutions = [],[],[]
        for degree in range(1,r.shape[0]+1):
            self.model.clock.check()
            z = self.solve(self.TT(v))
            w = self.T(z)
            vectors.append(v)
            images.append(w)
            solutions.append(z)
            # Batch columns are physically distinct port maps. Galerkin bases
            # are therefore separate. Only the Lagrange multiplier is shared.
            U = np.stack(vectors,axis=2).transpose(1,0,2)
            KU = np.stack(images,axis=2).transpose(1,0,2)
            if degree == 1 or degree % 4 == 0 or degree == r.shape[0]:
                projected = U.conj().transpose(0,2,1)@KU
                ev,S = np.linalg.eigh(.5*(projected+projected.conj().transpose(0,2,1)))
                ev = np.maximum(ev,0.)
                rhs = (U.conj().transpose(0,2,1) @ r.T[:,:,None])[:,:,0]
                b = (S.conj().transpose(0,2,1)@rhs[:,:,None])[:,:,0]
                power = abs(b)**2
                def excess(mu):
                    return np.sum(power/(1+mu*ev)**2)-self.delta**2
                hi=1./float(np.max(ev))
                for _ in range(80):
                    if excess(hi) < 0:
                        break
                    hi *= 4.
                if excess(hi) < 0:
                    mu = brentq(excess,0.,hi,xtol=np.nextafter(0.,1.),rtol=8*np.finfo(float).eps)
                    coeff = (S@(b/(ev+1/mu))[:,:,None])[:,:,0]
                    port_dual = (U@coeff[:,:,None])[:,:,0].T
                    rr = r - (KU@coeff[:,:,None])[:,:,0].T
                    linerr = np.linalg.norm(rr-port_dual/mu)/raw
                    if linerr <= 2e-9:
                        Z = np.stack(solutions,axis=2).transpose(1,0,2)
                        xf = xu-(Z @ coeff[:,:,None])[:,:,0].T
                        self.last_A_info=dict(multiplier=float(mu),krylov_dimension=degree,
                                              shifted_schur_residual=float(linerr),raw_residual=float(raw))
                        return xf,self.E+self.B(xf)
            # Fully reorthogonalized block of independent Arnoldi recurrences.
            wt = w.T[:,:,None].copy()
            for _ in range(2):
                wt -= U @ (U.conj().transpose(0,2,1)@wt)
            w = wt[:,:,0].T
            n = np.linalg.norm(w,axis=0)
            v = np.divide(w,n,out=np.zeros_like(w),where=n>1e-30)
        raise ArithmeticError('On-demand active projection missed the full Schur residual tolerance')


def run():
    name='fresnel_twin_8GHz_half180_tx'
    pub,case = read_public(name,HISTORICAL)
    shared=np.load(HERE/'checkpoints'/f'{name}__shared.npz')
    pmodel=model_for(pub,case)
    rt=time.perf_counter()
    weak,strong=engine.route_basis(pmodel)
    route_s=time.perf_counter()-rt
    recycled=np.column_stack((strong,weak))
    Y=pub['Y']
    delta=float(pub['delta'])
    baselines={}
    records={}
    for label,typ in [('eager',TimedEager),('on_demand',OnDemand)]:
        model=model_for(pub,case)
        t0=time.perf_counter()
        P=typ(model,Y,delta,recycled)
        setup_s=time.perf_counter()-t0
        setup_cost=dict(model.clock.counts)
        records[label]=dict(setup_s=setup_s,setup_cost=setup_cost,inputs={})
        if label=='eager':
            cx,cq,_=P.material(shared['before_J']/6,shared['before_Q'])
        for input_tag,x,q in [('first_AP',cx,cq),('active_air',np.zeros_like(cx),model.E)]:
            t=time.perf_counter();c0=time.process_time()
            before=Counter(model.clock.counts)
            xf,qf=P.affine_ball(x,q)
            secs,cpu=time.perf_counter()-t,time.process_time()-c0
            cc=dict(Counter(model.clock.counts)-before)
            mu=P.last_A_info['multiplier']
            r=6*model.ports.receive(xf)-Y
            rhs=x+P.BH(q-model.E)
            grad=P.H(xf)-rhs+mu*6*model.ports.adjoint(r)
            denom=np.linalg.norm(rhs)+mu*np.linalg.norm(6*model.ports.adjoint(Y))
            rec=dict(wall_s=secs,cpu_s=cpu,extra_cost=cc,
                     data_to_radius=float(np.linalg.norm(r)/delta),
                     normalized_stationarity=float(np.linalg.norm(grad)/max(denom,1e-300)),
                     **P.last_A_info)
            if label=='eager':
                baselines[input_tag]=(xf,qf,mu)
            else:
                xref,qref,muref=baselines[input_tag]
                rec.update(x_relative=engine.ratio(xf-xref,xref),q_relative=engine.ratio(qf-qref,qref),
                           multiplier_relative=float(abs(mu-muref)/max(abs(muref),1.)))
            records[label]['inputs'][input_tag]=rec
            print(label,input_tag,'setup',setup_s,'projection',secs,'delta',rec['data_to_radius'],
                  'stationarity',rec['normalized_stationarity'],flush=True)
            (HERE/'lazy_ports.json').write_text(json.dumps(dict(case=name,host=platform.platform(),
                  shared_route_s=route_s,records=records),indent=2,default=float))
    return records


if __name__=='__main__':
    run()

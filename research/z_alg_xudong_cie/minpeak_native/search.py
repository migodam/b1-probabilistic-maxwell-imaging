"""Public-data box search on the exact native 96² FFT forward model."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
from time import perf_counter, process_time

import numpy as np
from scipy.fft import fft2, ifft2
from scipy.optimize import minimize
from scipy.sparse.linalg import LinearOperator, gmres
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
STUDY = HERE.parents[1]/'z_alg_optical_inertia'
sys.path.insert(0, str(STUDY/'submission/curve'))
from curve_inputs import full_ring_inputs

SIDE = 96
N = SIDE**2
RTOL = 7.e-7


def clocks():
    return perf_counter(), process_time()


def elapsed(start):
    return dict(wall_seconds=perf_counter()-start[0], cpu_seconds=process_time()-start[1])


def norm_columns(z):
    return np.sqrt(np.sum(z.real*z.real+z.imag*z.imag, axis=0))


class Model:
    def __init__(self, spectrum, sensor, incident, measurement):
        self.spectrum, self.s = spectrum, sensor
        self.b, self.y = incident, measurement
        self.scale = np.linalg.norm(measurement)**2
        self.count = Counter()

    def green(self, z, adjoint=False):
        ncols = z.shape[1]
        self.count['fft_columns'] += ncols
        z0 = fft2(z.reshape(SIDE,SIDE,ncols), s=self.spectrum.shape, axes=(0,1))
        s = self.spectrum.conj() if adjoint else self.spectrum
        return ifft2(z0*s[:,:,None], axes=(0,1))[:SIDE,:SIDE].reshape(N,ncols)

    @staticmethod
    def dot(a,b):
        return np.sum(a.conj()*b, axis=0)

    def linear_solve(self, material, rhs, adjoint, seed):
        """Vectorized independent BiCGSTAB recurrences; GMRES guards difficult columns."""
        if not np.any(material):
            return rhs.copy()
        self.count['adjoint_rhs' if adjoint else 'state_rhs'] += rhs.shape[1]
        m = material[:,None]
        op = (lambda z: z-self.green(m*z,True)) if adjoint else (lambda z: z-m*self.green(z))
        sol = np.zeros_like(rhs) if seed is None else seed.copy()
        r = rhs.copy() if seed is None else rhs-op(sol)
        tolerance = RTOL*np.maximum(norm_columns(rhs), 1e-20)
        initial = r.copy()
        ns = rhs.shape[1]
        one = np.ones(ns,complex)
        alpha, omega, rho = one.copy(), one.copy(), one.copy()
        p, v = np.zeros_like(r), np.zeros_like(r)
        to_solve = norm_columns(r) > tolerance
        for _ in range(85):
            ids = np.flatnonzero(to_solve)
            if not len(ids):
                break
            rn = self.dot(initial[:,ids],r[:,ids])
            ok = (np.abs(rn)>1e-55)&(np.abs(omega[ids])>1e-30)
            to_solve[ids[~ok]] = False
            ids, rn = ids[ok], rn[ok]
            if len(ids)==0:
                break
            beta = rn/rho[ids] * alpha[ids]/omega[ids]
            rho[ids] = rn
            p[:,ids] = r[:,ids] + beta*(p[:,ids]-omega[ids]*v[:,ids])
            v[:,ids] = op(p[:,ids])
            bottom = self.dot(initial[:,ids],v[:,ids])
            ok = np.abs(bottom)>1e-55
            to_solve[ids[~ok]] = False
            ids, bottom = ids[ok],bottom[ok]
            if not len(ids):
                break
            alpha[ids] = rho[ids]/bottom
            sm = r[:,ids]-alpha[ids]*v[:,ids]
            sol[:,ids] += alpha[ids]*p[:,ids]
            ok = norm_columns(sm)>tolerance[ids]
            to_solve[ids[~ok]] = False
            ids,sm = ids[ok],sm[:,ok]
            if not len(ids):
                break
            t = op(sm)
            tnorm = np.maximum(self.dot(t,t).real,1e-55)
            omega[ids] = self.dot(t,sm)/tnorm
            sol[:,ids] += omega[ids]*sm
            r[:,ids] = sm-omega[ids]*t
            to_solve[ids[norm_columns(r[:,ids])<=tolerance[ids]]] = False
            self.count['bicgstab_iterations'] += 1
        for k in np.flatnonzero(norm_columns(rhs-op(sol))>tolerance):
            self.count['gmres_fallback_rhs'] += 1
            action=LinearOperator((N,N),matvec=lambda z:op(z[:,None])[:,0],dtype=complex)
            start = sol[:,k] if np.all(np.isfinite(sol[:,k])) else None
            sol[:,k], info = gmres(action,rhs[:,k],x0=start,rtol=RTOL,atol=0., restart=50,maxiter=12)
            if info:
                self.count['failed_state_solves'] += 1
                raise ArithmeticError(f'GMRES failed at right hand side {k}, info={info}.')
        return sol

    def full(self, material, previous=None, gradient=True):
        self.count['native_full_wave_calls'] += 1
        seed = previous[0] if previous is not None else None
        j = self.linear_solve(material,material[:,None]*self.b,False,seed)
        ports = self.s @ j
        residual = ports-self.y
        fr = np.linalg.norm(residual)
        if not gradient:
            return fr, ports, j
        incident = self.b+self.green(j)
        ob = self.s.conj().T @ residual
        seed = previous[1] if previous is not None else None
        adj = self.linear_solve(material,ob,True,seed)
        g = np.sum((incident*adj.conj()).real,axis=1)/self.scale
        return .5*fr*fr/self.scale, g, fr, (j,adj), ports


def groups(side):
    rr,cc = np.indices((SIDE,SIDE))
    return ((rr*side//SIDE)*side+cc*side//SIDE).ravel()


def coarsen(material, index, n):
    return np.bincount(index,weights=material,minlength=n)/np.bincount(index,minlength=n)


def born_initial(model,index,active,upper):
    n = int(index.max())+1
    field = np.einsum('ri,is->rsi',model.s,model.b).reshape(-1,N)
    a = np.stack([np.bincount(index,weights=f.real,minlength=n)+
                  1j*np.bincount(index,weights=f.imag,minlength=n) for f in field])[:,active]
    target = model.y.ravel()
    hess = (a.conj().T@a).real
    b = (a.conj().T@target).real
    lips = np.linalg.norm(hess,ord=2)
    ridge = .001*lips
    m = w = np.zeros(len(active))
    momentum = 1.
    for _ in range(120):
        m0 = np.clip(w-(hess@w-b+ridge*w)/(lips+ridge),0,upper)
        next_momentum = .5*(1+np.sqrt(1+4*momentum*momentum))
        w = m0+(momentum-1)/next_momentum*(m0-m)
        m,momentum = m0,next_momentum
    return m


class BudgetEnd(Exception):
    pass


class Feasible(Exception):
    pass


def box_fit(model, material, active, side, cap, target, max_calls, max_cpu, keep_fit):
    index, size = groups(side), side*side
    start = clocks()
    call0 = model.count['native_full_wave_calls']
    base = np.zeros(size)
    best = {'residual':float('inf')}
    previous = None

    def loss(x):
        nonlocal previous
        if model.count['native_full_wave_calls']-call0>=max_calls or process_time()-start[1]>=max_cpu:
            raise BudgetEnd()
        base[active] = x
        m = base[index]
        try:
            f,g,fr,state,ports = model.full(m,previous)
        except ArithmeticError:
            previous = None
            if 'm' not in best:
                raise
            delta = x-best['x']
            return .5 + .01*np.dot(delta,delta), .02*delta
        previous = state
        if fr < best['residual']:
            best.update(residual=float(fr),m=m.copy(),x=x.copy(),ports=ports.copy(),j=state[0].copy())
        if keep_fit and fr < target:
            raise Feasible()
        g0 = np.bincount(index,weights=g,minlength=size)[active]
        return f,g0

    initial = np.clip(coarsen(material,index,size)[active],0,cap)
    try:
        result = minimize(loss, initial, jac=True, method='L-BFGS-B',bounds=[(0,cap)]*len(active),
                          options={'maxiter':max_calls,'maxfun':max_calls,'gtol':1.e-10,
                                   'ftol':2.e-12,'maxls':25,'maxcor':10})
        status = str(result.message)
    except BudgetEnd:
        status = 'stage budget'
    except Feasible:
        status = 'reached numerical interior of noise ball'
    return best, {'level':side,'upper_cap':cap,'stop':status,'active_variables':len(active),
                  'native_calls':model.count['native_full_wave_calls']-call0,
                  'best_residual':best['residual'],'best_peak':float(best['m'].max()),
                  'positive_cells':int(np.count_nonzero(best['m'])),**elapsed(start)}


def bounded_support(material, side, max_cells):
    index = groups(side)
    intensities = coarsen(material,index,side*side)
    pad = np.pad(intensities.reshape(side,side),1)
    halo = np.maximum.reduce([pad[dy:dy+side,dx:dx+side] for dy in range(3) for dx in range(3)])
    key = intensities+.15*halo.ravel()
    k = max_cells//((SIDE//side)**2)
    return np.sort(np.argsort(-key,kind='stable')[:k])


def run_one(spectrum,sensor,incident,obs,args):
    start = clocks()
    model = Model(spectrum,sensor,incident,obs.y)
    history = []
    initial_cap = args.initial_cap
    first_side = 24
    idx = groups(first_side)
    active = np.arange(first_side**2)
    initial = born_initial(model,idx,active,initial_cap)[idx]
    target = obs.padded_radius * .998
    state, row = box_fit(model,initial,active,first_side,initial_cap,target,
                         args.stage_calls, args.stage_cpu, True)
    history.append(row)
    print(obs.id,row,flush=True)
    best = state
    eligible = state if state['residual']<target else None
    for side in (48,96):
        active = bounded_support(state['m'],side,args.max_cells)
        state,row = box_fit(model,state['m'],active,side,initial_cap,target,
                            args.stage_calls,args.stage_cpu,True)
        history.append(row)
        print(obs.id,row,flush=True)
        if state['residual'] < best['residual']:
            best = state
        if state['residual']<target and (eligible is None or (np.count_nonzero(eligible['m'])>args.max_cells or state['m'].max()<=eligible['m'].max())):
            eligible=state
    if eligible is not None:
        active = bounded_support(eligible['m'],96,args.max_cells)
        left, right = .5, float(eligible['m'].max())
        for _ in range(args.peaks):
            cap = left + .67*(right-left)
            trial,row=box_fit(model,eligible['m'],active,96,cap,target,args.peak_calls,
                              args.stage_cpu,True)
            history.append(row)
            print(obs.id,row,flush=True)
            if trial['residual']<best['residual']:
                best=trial
            if trial['residual']<target:
                eligible=trial
                right=float(eligible['m'].max())
            else:
                left=cap
    selected=eligible if eligible is not None else best
    arrays = {'m_selected':selected['m'],'Y_selected_numeric':selected['ports'],
              'j_selected_numeric':selected['j'], 'm_best_residual':best['m'],
              'Y_best_residual_numeric':best['ports'], 'Y':obs.y}
    np.savez_compressed(args.out/(obs.id+'.npz'),**arrays)
    receipt = {'observation':obs.id,'radius':obs.padded_radius,'native_side':SIDE,
               'noise_pct':obs.noise_pct,'object':obs.object_id,
               'selected_numeric_residual':selected['residual'],
               'selected_peak':float(selected['m'].max()),'selected_positive_cells':int(np.count_nonzero(selected['m'])),
               'eligible_for_independent_validation':eligible is not None,
               'best_numeric_residual':best['residual'],
               'stages':history,'counts':dict(model.count),**elapsed(start)}
    (args.out/(obs.id+'.json')).write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
    print('DONE',obs.id,receipt['selected_peak'],receipt['selected_numeric_residual'],receipt['wall_seconds'],flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=HERE/'runs')
    parser.add_argument('--object',type=str,action='append')
    parser.add_argument('--noise',type=float,action='append')
    parser.add_argument('--initial-cap',type=float,default=3.8)
    parser.add_argument('--max-cells',type=int,default=2600)
    parser.add_argument('--stage-calls',type=int,default=70)
    parser.add_argument('--stage-cpu',type=float,default=28)
    parser.add_argument('--peak-calls',type=int,default=52)
    parser.add_argument('--peaks',type=int,default=5)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    s,p,b,obs=full_ring_inputs()
    for observation in obs:
        if observation.noise_pct not in (args.noise or [1.,2.]):
            continue
        if args.object and observation.object_id not in args.object:
            continue
        run_one(s,p,b,observation,args)


if __name__=='__main__':
    with threadpool_limits(1):
        main()

"""Paid exact dense rational responses on complete box-feasible rays."""
import numpy as np
from scipy.linalg import lu_solve, schur, solve_triangular
from scipy.optimize import minimize_scalar


def feasible_interval(m, p, lower=.5, upper=20.):
    active = np.abs(p) > 1e-14
    a = (lower-m[active])/p[active]
    b = (upper-m[active])/p[active]
    return float(np.max(np.minimum(a,b))), float(np.min(np.maximum(a,b)))


class RationalRay:
    def __init__(self, model, state, direction):
        model._matching(state)
        self.model, self.state = model, state
        self.p = np.asarray(direction, float)
        n = model.n
        with model.operation('RAY_DENSE_KERNEL', rhs=n+model.ns,
                             ray_kernel_RHS=n, ray_illumination_RHS=model.ns,
                             dense_ray_kernel_builds=1, ray_kernel_dimension=n):
            kernel = lu_solve(model.factor, self.p[:,None]*model.G, check_finite=False)
            illumination = lu_solve(model.factor, self.p[:,None]*state['Q'], check_finite=False)
        with model.operation('RAY_COMPLEX_SCHUR', dense_complex_schur=1, schur_dimension=n):
            self.T, Z = schur(kernel, output='complex', check_finite=False)
            self.left = model.S@Z
            self.right = Z.conj().T@illumination
        self.identity = np.eye(n,dtype=complex)

    def response(self, t):
        with self.model.operation('RAY_TRANSFER_EVALUATION', ray_transfer_evaluations=1,
                                  triangular_system_dimension=self.model.n,
                                  triangular_small_transfer_RHS=self.model.ns):
            w = solve_triangular(self.identity-float(t)*self.T, self.right,
                                 lower=False, check_finite=False)
            return self.state['y']+float(t)*(self.left@w)


def scan_ray(model, state, direction, target, scale2, regularizer, config, label, record=None):
    """All scan nodes are recorded; every separated grid valley is refined."""
    p = np.asarray(direction,float)
    if np.sqrt(np.mean(p*p)) < 1e-10:
        return [], dict(direction=label,status='ZERO_DIRECTION',nodes=[],minima=[])
    lo,hi = feasible_interval(state['m'],p,config['lower'],config['upper'])
    if record is None:record={}
    record.update(direction=label,status='SCANNING',feasible_interval=[lo,hi],nodes=[],minima=[],
                  anchor_material=state['m'].copy(),material_direction=p.copy(),
                  physical_dimension=model.n,approximation='none; full dense complex Schur')
    nodes=record['nodes']
    ray = RationalRay(model,state,p)
    def objective(t, refined=False):
        prediction=ray.response(t)
        m=state['m']+float(t)*p
        residual=float(np.linalg.norm(prediction-target)/np.sqrt(scale2))
        reg=regularizer.evaluate(m)[0]
        value=.5*residual*residual+reg
        nodes.append(dict(t=float(t),raw_score=value,relative_residual=residual,
                          regularization=reg,refined=bool(refined)))
        return value
    grid=np.linspace(lo,hi,int(config['ray_points']))
    values=np.asarray([objective(t) for t in grid])
    valleys=[i for i in range(len(grid)) if
             (i==0 or values[i]<=values[i-1]) and
             (i==len(grid)-1 or values[i]<=values[i+1])]
    proposals=[]
    for i in valleys:
        t=float(grid[i]); value=float(values[i])
        if 0<i<len(grid)-1:
            result=minimize_scalar(lambda x:objective(x,True),bounds=(grid[i-1],grid[i+1]),
                                   method='bounded',options={'maxiter':int(config['ray_refine_iterations']),
                                                            'xatol':float(config['ray_refine_tolerance'])})
            if result.fun<value:t,value=float(result.x),float(result.fun)
        m=np.clip(state['m']+t*p,config['lower'],config['upper'])
        proposals.append(dict(m=m,raw_score=value,t=t,ray=label,grid_valley_index=i,
                              feasible_interval=[lo,hi],corrected=False))
    record.update(status='SCANNED',minima=[{k:v for k,v in q.items() if k!='m'} for q in proposals])
    return proposals,record

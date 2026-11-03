"""Common mature bounded SciPy TRF and analytic-adjoint L-BFGS-B controls."""
from collections import Counter
import time
import numpy as np
from scipy.optimize import least_squares, minimize


def corrector(model,Y,m0,config=None):
    cfg = dict(config or {})
    box = cfg.get('box',cfg.get('material_box',(.02,20.)))
    # Common box projection is part of handing a bounded QP seed to TRF;
    # adding its h to the anchor can round an active coordinate outside by ulps.
    m0 = np.clip(np.asarray(m0, float), box[0], box[1])
    method = cfg.get('method','trf')
    maximum = cfg.get('max_nfev',cfg.get('corrector_nfev',120))
    start = time.perf_counter()
    initial_cost = Counter(model.cost)
    scale = np.linalg.norm(Y)
    history, archive = [], []
    incumbent = [np.asarray(m0,float).copy(),np.inf]
    deadline = start + cfg.get('wall_seconds',np.inf)

    class EndBudget(Exception):
        pass

    def record(m, state=None):
        state = model.state(m) if state is None else state
        rr = state['Y']-Y
        data = .5*np.vdot(rr,rr).real/(scale*scale)
        value = float(data+model.reg.evaluate(m)[0])
        history.append(dict(wall=time.perf_counter()-start,physical_score=value,
                            relative_residual=float(np.linalg.norm(rr)/scale),
                            LU=int(model.cost['LU']-initial_cost['LU'])))
        if value < incumbent[1]:
            incumbent[:] = [np.asarray(m).copy(),value]
            archive.append(dict(m=np.asarray(m).copy(),physical_score=value,
                relative_residual=float(np.linalg.norm(rr)/scale),
                wall=time.perf_counter()-start,cost=dict(Counter(model.cost)-initial_cost)))
        if time.perf_counter()>=deadline:
            raise EndBudget()
        return state

    reciprocal = method == 'reciprocal'
    def to_material(x):
        return 1/np.asarray(x) if reciprocal else np.asarray(x)

    def fun(x):
        m = to_material(x)
        state = record(m)
        prior,_ = model.regularizer_residual(m)
        return np.r_[model.pack(state['Y']-Y)/scale,prior]

    def jac(x):
        m = to_material(x)
        _,jp = model.regularizer_residual(m)
        jj = np.vstack((model.jacobian(m)/scale,jp))
        return jj*(-m*m)[None,:] if reciprocal else jj

    stopped = 'CONVERGED_OR_MAX_NFEV'
    try:
        if method in ('trf','native','reciprocal'):
            x0 = 1/np.asarray(m0) if reciprocal else np.asarray(m0)
            bounds = (1/box[1],1/box[0]) if reciprocal else box
            result = least_squares(fun,x0,jac=jac,bounds=bounds,method='trf',x_scale='jac',
                                   max_nfev=maximum,ftol=1e-10,xtol=1e-10,gtol=1e-8)
            record(to_material(result.x))
            stopped = str(result.message)
        else:
            beta = cfg.get('beta',.2) if method in ('feedback','globaldual') else None
            if method=='globaldual':
                beta=0.
            reference = model.state(np.full(model.n,.1)) if beta is not None else None
            def objective(m):
                value,g = model.data_value_gradient(m,Y,beta,reference)
                record(m)
                return value,g
            result = minimize(objective,np.asarray(m0),jac=True,method='L-BFGS-B',
                              bounds=[box]*model.n,options=dict(maxiter=maximum,maxfun=4*maximum,
                                                               ftol=1e-12,gtol=1e-8,maxls=30))
            record(result.x)
            stopped=str(result.message)
    except EndBudget:
        stopped='WALL_BUDGET'
    m = incumbent[0]
    return dict(m=m,candidates=archive,history=history,
                wall=time.perf_counter()-start,cost=dict(Counter(model.cost)-initial_cost),stop=stopped)


run_corrector = corrector

"""Dense full-wave projected Gauss–Newton on the original real material grid."""
from time import perf_counter, process_time

import numpy as np
from scipy.linalg import cho_factor, cho_solve, lu_factor, lu_solve


class Wave:
    def __init__(self, model):
        self.G,self.S,self.B = (model[x] for x in ("G","S","E0"))
        self.eye = np.eye(len(self.G),dtype=complex)
        self.lu_count = self.rhs = 0

    def state(self, material):
        matrix = self.eye-material[:,None]*self.G
        factor = lu_factor(matrix,check_finite=False)
        current = lu_solve(factor,material[:,None]*self.B,check_finite=False)
        self.lu_count += 1
        self.rhs += self.B.shape[1]
        return self.S@current,self.B+self.G@current,factor

    def transfer(self, factor):
        transfer = lu_solve(factor,self.S.T,trans=1,check_finite=False).T
        self.rhs += self.S.shape[0]
        return transfer


class Prior:
    def __init__(self, shape):
        grid = np.arange(np.prod(shape)).reshape(shape)
        self.left = np.r_[grid[:-1,:].ravel(),grid[:,:-1].ravel()]
        self.right = np.r_[grid[1:,:].ravel(),grid[:,1:].ravel()]
        self.tv = self.smooth = 1e-4

    def evaluate(self, material):
        delta = material[self.left]-material[self.right]
        q = np.hypot(delta,.01)
        value = np.mean(self.smooth*delta**2+self.tv*(q-.01))
        weight = (2*self.smooth+self.tv/q)/len(delta)
        return value,weight*delta,weight

    def gradient(self, edge, pixels):
        return np.bincount(self.left,edge,minlength=pixels)-np.bincount(self.right,edge,minlength=pixels)

    def add_curvature(self, H, weight):
        n = len(H)
        H[np.diag_indices(n)] += (np.bincount(self.left,weight,minlength=n) +
                                  np.bincount(self.right,weight,minlength=n))
        H[self.left,self.right] -= weight
        H[self.right,self.left] -= weight


def realify(x):
    return np.concatenate((x.real,x.imag),axis=0)


def voxel_gn(model, case, rhs_budget, iterations):
    wave,prior = Wave(model),Prior(model["shape"])
    y, bounds = case["Y"],model["bounds"]
    scale = np.linalg.norm(y)
    m = np.clip(np.ones(len(wave.G)),*bounds)
    wall,cpu = perf_counter(),process_time()
    prediction,field,factor = wave.state(m)
    residual = realify((prediction-y).T.ravel())/scale
    reg,gap,weight = prior.evaluate(m)
    objective = .5*residual@residual+reg
    damping = 1e-3
    stop,trace = "iteration_cap",[]
    for iteration in range(iterations):
        if wave.rhs+wave.S.shape[0]+wave.B.shape[1]>rhs_budget:
            stop = "rhs_cap"
            break
        transfer = wave.transfer(factor)
        complex_jac = np.einsum("ri,is->sri",transfer,field).reshape(-1,len(m))
        jac = realify(complex_jac)/scale
        H = jac.T@jac
        gradient = jac.T@residual
        gradient += prior.gradient(gap,len(m))
        prior.add_curvature(H,weight)
        diag = np.diag_indices(len(H))
        unit_damping = max(float(np.trace(H)/len(H)),1e-12)
        old_objective = objective
        accepted = False
        count_before = wave.lu_count
        for _ in range(4):
            H[diag] += damping*unit_damping
            step = cho_solve(cho_factor(H,check_finite=False),-gradient,check_finite=False)
            H[diag] -= damping*unit_damping
            alpha0 = min(1.,4./max(np.max(abs(step)),1e-30))
            for shrink in range(7):
                if wave.rhs+wave.B.shape[1]>rhs_budget:
                    break
                alpha = alpha0*.5**shrink
                trial = np.clip(m+alpha*step,*bounds)
                difference = trial-m
                if np.linalg.norm(difference,np.inf)<1e-10:
                    continue
                trial_y,trial_field,trial_factor = wave.state(trial)
                trial_residual = realify((trial_y-y).T.ravel())/scale
                trial_reg,trial_gap,trial_weight = prior.evaluate(trial)
                trial_objective = .5*trial_residual@trial_residual+trial_reg
                if trial_objective < objective:
                    trace.append(dict(iteration=iteration+1,objective=float(trial_objective),
                                      residual_rel=float(np.linalg.norm(trial_residual)),rhs=wave.rhs,
                                      wall_s=perf_counter()-wall,alpha=float(alpha),
                                      max_step=float(np.max(abs(difference))),damping=damping,
                                      relative_retraction_change=float(
                                          np.linalg.norm(difference-alpha*step)/np.linalg.norm(alpha*step)),
                                      touched_box_fraction=float(np.mean(
                                          (trial<=bounds[0]+1e-12)|(trial>=bounds[1]-1e-12))),
                                      rejected_trial_solves=wave.lu_count-count_before-1))
                    m,residual,objective = trial,trial_residual,trial_objective
                    prediction,field,factor = trial_y,trial_field,trial_factor
                    gap,weight = trial_gap,trial_weight
                    damping = max(1e-10,damping*.5)
                    accepted = True
                    break
            if accepted:
                break
            damping *= 10
        if not accepted:
            stop = "no_decrease" if wave.rhs+wave.B.shape[1]<=rhs_budget else "rhs_cap"
            break
        if old_objective-objective<1e-11:
            stop = "gain_tolerance"
            break
    return m,dict(stop=stop,trace=trace,
                  residual_over_sigma_rms=float(np.linalg.norm(prediction-y)/(case["sigma"]*np.sqrt(y.size))),
                  online_cost=dict(wall=perf_counter()-wall,cpu=process_time()-cpu,
                                   full_complex_rhs=wave.rhs,full_factorizations=wave.lu_count))

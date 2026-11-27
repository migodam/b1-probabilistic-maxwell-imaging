"""Bounded, full-wave real-material Gauss–Newton on an intersection of data balls."""
from collections import Counter

from ports import elapsed, now, project, realify
import numpy as np
from scipy.linalg import cho_factor, cho_solve

UPPER = 4.5


def mask_family(family, roi, patch, low, high, floor):
    lb = np.zeros(roi.size)
    ub = np.full(roi.size, UPPER)
    mean_min, mean_max = 0., UPPER
    if family == "zero":
        ub[roi] = 0
    elif family == "patch":
        lb[patch] = floor
    elif family == "mean_low":
        mean_max = low
    elif family == "mean_high":
        mean_min = high
    elif family == "unrestricted":
        pass
    elif isinstance(family, float):
        mean_min = mean_max = family
    else:
        raise ValueError(family)
    return lb, ub, roi, mean_min, mean_max


def box_project(x, bounds):
    lb,ub,roi,a,b=bounds
    y=np.clip(x,lb,ub)
    if a==0 and b==UPPER:
        return y
    t=float(np.clip(y[roi].mean(),a,b))
    if abs(t-y[roi].mean()) < 1e-13:
        return y
    values=x[roi]
    lower,upper=lb[roi],ub[roi]
    lo=float(np.min(lower-values))
    hi=float(np.max(upper-values))
    for _ in range(42):
        m=(lo+hi)/2
        if np.clip(values+m,lower,upper).mean()<t:
            lo=m
        else:
            hi=m
    # Round toward the feasible side for one-sided constraints.
    y[roi]=np.clip(values+(hi if t==a else lo),lower,upper)
    return y


class Observation:
    def __init__(self, tx, rx, y, radius):
        self.tx=np.asarray(tx)
        self.rx=np.asarray(rx)
        self.y=y
        self.radius=radius

    def project(self, matrix):
        return project(matrix,self.tx,self.rx)

    def residual(self, matrix):
        return (self.project(matrix)-self.y).ravel()/self.radius

    def j_block(self, jac):
        return realify(self.project(jac).reshape(self.y.size,-1))/self.radius


def describe(st):
    return dict(residual_over_each_radius=st["ratios"], objective=st["objective"],
                original_state_defect=st["state_defect"])


class Attack:
    def __init__(self,engine,observations,target,projection=box_project):
        self.engine=engine
        self.observations=observations
        self.target=target
        self.project=projection

    def state(self,x):
        st=self.engine.state(x)
        blocks=[o.residual(st["delta"]) for o in self.observations]
        st["ratios"]=[float(np.linalg.norm(r)) for r in blocks]
        # Use the actual per-layout balls. The increasing weight on the worst
        # residual helps a fit find the intersection after unequal starts.
        norm=np.asarray(st["ratios"])
        scale=np.maximum(norm,1.)
        st["weights"]=(scale/scale.max())**.8
        st["objective"]=float(np.linalg.norm(norm,ord=2.8)**2)
        st["residual"]=np.concatenate([w*realify(r) for r,w in zip(blocks,st["weights"])])
        return st

    def jacobian(self,st):
        jac=self.engine.jacobian(st)
        return np.vstack([w*o.j_block(jac) for o,w in zip(self.observations,st["weights"])])

    def fit(self,initial,bounds,cap):
        start=now()
        before=Counter(self.engine.counters)
        x=self.project(initial,bounds)
        st=self.state(x)
        first=describe(st)
        lb,ub,roi,lo_mean,hi_mean=bounds
        trace=[]
        trust=.85
        damping=3e-4
        accepted=0
        stop="major_iteration_budget"
        for it in range(cap):
            if max(st["ratios"])<=self.target:
                stop="joint_discrepancy_target"
                break
            jac=self.jacobian(st)
            residual=st["residual"]
            grad=jac.T@residual
            pg=x-self.project(x-grad,bounds)
            pg_max=float(np.max(abs(pg)))
            entry=dict(iteration=it+1,projected_gradient_inf=pg_max,damping=damping,trust=trust,trials=[])
            trace.append(entry)
            if pg_max<1e-9:
                stop="projected_stationarity"
                break
            free=(ub>lb)&~((x<=lb+1e-12)&(grad>0))&~((x>=ub-1e-12)&(grad<0))
            jfree=jac[:,free]
            mean_now=float(x[roi].mean())
            active_mean=(lo_mean>0 and mean_now<=lo_mean+1e-10) or (hi_mean<UPPER and mean_now>=hi_mean-1e-10)
            roi_in_free=roi[free]
            if active_mean and roi_in_free.any():
                jfree[:,roi_in_free]-=jfree[:,roi_in_free].mean(axis=1,keepdims=True)
            gram=jfree@jfree.T
            gram.flat[::len(gram)+1]+=max(damping*np.trace(gram)/len(gram),1e-13)
            direction=np.zeros(len(x))
            direction[free]=-jfree.T@cho_solve(cho_factor(gram,check_finite=False),residual,check_finite=False)
            direction*=min(1.,trust/max(np.max(abs(direction)),1e-300))
            for alpha in (1., .5, .2, .08):
                candidate=self.project(x+alpha*direction,bounds)
                trial=self.state(candidate)
                entry["trials"].append(dict(alpha=alpha,**describe(trial)))
                if trial["objective"] < st["objective"]*(1-1e-11):
                    x,st=candidate,trial
                    accepted+=1
                    trust=min(1.75,trust*1.2)
                    damping=max(3e-7,damping*.7)
                    break
            else:
                direction=-pg*min(1.,trust/pg_max)
                candidate=self.project(x+direction,bounds)
                trial=self.state(candidate)
                entry["trials"].append(dict(alpha="projected_gradient",**describe(trial)))
                if trial["objective"] < st["objective"]*(1-1e-11):
                    x,st=candidate,trial
                    accepted+=1
                else:
                    trust*=.6
                    damping=min(5.,damping*5)
                    if trust<.006:
                        stop="trust_region_stalled"
                        break
        if max(st["ratios"])<=self.target:
            stop="joint_discrepancy_target"
        if np.any(x<lb) or np.any(x>ub) or not lo_mean-1e-10<=x[roi].mean()<=hi_mean+1e-10:
            raise ArithmeticError("Fit escaped registered material family")
        return st,dict(initial=first,final=describe(st),accepted=accepted,stop=stop,
                       trace=trace,counts=dict(self.engine.counters-before),**elapsed(start))


def born_initial(engine, observation):
    j=np.einsum("ri,ij->rji",engine.c,engine.b)
    jac=observation.j_block(j)
    target=realify(observation.y.ravel())/observation.radius
    center=np.full(engine.n,.3)
    gram=jac@jac.T
    gram.flat[::len(gram)+1]+=0.015*np.trace(gram)/len(gram)
    lift=jac.T@cho_solve(cho_factor(gram,check_finite=False),target-jac@center,check_finite=False)
    return np.clip(center+lift,0.,UPPER)

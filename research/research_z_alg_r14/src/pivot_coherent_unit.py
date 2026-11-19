"""Search an old-compatible known-unit-gain pair for coherent physical port feeds."""
import argparse
import json
from time import process_time

from flint import acb_mat, arb
import numpy as np
from scipy.optimize import minimize
from scipy.special import j0
from scipy.stats import qmc

from continuous_band import Physical, lower, norm, pnorm, pt, up
from data14 import Fit, FrequencyTM, ROOT, cn_radius, load_cases, ri
from pivot_menu_witnesses import OFFSETS

SELECTED = ('fresh6-16', 'fresh6-10', 'fresh6-02', 'fresh6-05', 'fresh6-15', 'fresh6-03')
NODES = (-.02,-.01,0.,.01,.02)
CAP = .745


def rms(a,b):
    return float(np.linalg.norm(a-b) / np.sqrt(len(a)))


class Search:
    def __init__(self,case,nodes):
        self.case = case
        self.model=FrequencyTM(case.geometry)
        self.nodes=nodes
        self.mid=nodes.index(0.)
        self.noise=cn_radius(case.sigma,1)
        self.phy=Physical(case)
        self.strict_ops={}
        p=np.asarray(case.geometry['sources'])
        srcdist=np.linalg.norm(p[:,None]-p[None,:],axis=-1)
        self.whiten=[]
        for f in nodes:
            p=j0(case.geometry['k']*(1+f)*srcdist)
            lam, u=np.linalg.eigh(p)
            self.whiten.append((u * np.sqrt(3/lam)) @ u.T)
        self.whiten=np.array(self.whiten)

    def numerical(self,x):
        f,J=zip(*(self.model.response(x,d) for d in self.nodes))
        F, J=np.asarray(f),np.asarray(J)
        p=np.einsum('fij,fj->fi',self.whiten,F) / (2*self.noise)
        j=np.einsum('fij,fjk->fik',self.whiten,J) / (2*self.noise)
        r=(F[self.mid]-self.case.old)/self.case.additive_radius
        v=J[self.mid]/self.case.additive_radius
        return r,v,p,j

    def constraints(self,a,b):
        r,v,p,j=self.numerical(a)
        s,w,q,k=self.numerical(b)
        d=p-q
        value=np.r_[.997-np.vdot(r,r).real,.997-np.vdot(s,s).real,CAP**2-np.sum(abs(d)**2,axis=1)]
        grad=np.zeros((len(value),12))
        grad[0,:6]=-2*np.real(r.conj()@v)
        grad[1,6:]=-2*np.real(s.conj()@w)
        grad[2:,:6]=-2*np.real(np.einsum('fi,fij->fj',d.conj(),j))
        grad[2:,6:]=2*np.real(np.einsum('fi,fij->fj',d.conj(),k))
        return value,grad

    def transform(self,x):
        _,v,_,j=self.numerical(x)
        A=np.vstack((ri(v),ri(j.reshape(-1,6))/np.sqrt(len(self.nodes))))
        _,s,vh=np.linalg.svd(A,full_matrices=False)
        return vh.T / np.maximum(s,.002)

    def old_cost(self,x):
        return float(np.linalg.norm((self.model.response(x)[0]-self.case.old)/self.case.additive_radius)**2)

    def line(self,m,d):
        lo,hi=self.case.box
        dirs=[d.copy(),-d.copy()]
        ends=[]
        for z in dirs:
            active=((z>0)&(m>hi-1e-8))|((z<0)&(m<lo+1e-8))
            z[active]=0.
            span=[(hi-x)/v if v>0. else (lo-x)/v for x,v in zip(m,z) if v]
            ends.append(min(4.,*span) if span else 0.)
        t0,t1=0.,min(x for x in ends if x>0.) if any(x>0. for x in ends) else 0.
        for k in range(13):
            t=t1 if k==0 else (t0+t1)/2
            a,b=[np.clip(m+min(t,e)*v,lo,hi) for v,e in zip(dirs,ends)]
            if min(self.constraints(a,b)[0])>0.:
                t0=t
                if k==0:break
            else:
                t1=t
        return tuple(np.clip(m+min(t0,e)*v,lo,hi) for e,v in zip(ends,dirs))

    def strict(self,m):
        old=[]; gaps={}; limit=pt(2*self.noise)
        for f in sorted(set(OFFSETS)|set(NODES)):
            if f not in self.strict_ops:
                G,B,S,P=self.phy.kernels(pt(f),arb(0))
                if any(not acb_mat([[P[i,j] for j in range(n)] for i in range(n)]).det().real>0
                       for n in range(1,P.nrows()+1)):
                    raise ValueError('Positive definiteness of the port Gram was not enclosed')
                self.strict_ops[f]=((G,B,S,P),P.inv())
            kernels,Pi=self.strict_ops[f]
            vals=[self.phy.response(m0,kernels,False)[0] for m0 in m]
            gap=arb(3).sqrt()*pnorm(vals[0]-vals[1],Pi)
            gaps[f]=up(gap/limit)
            if f==0.:
                obs=acb_mat([[complex(z) for z in self.case.old]])
                old=[up((norm(u-obs)/pt(self.case.additive_radius))**2) for u in vals]
        rm=sum(((pt(a)-pt(b))**2 for a,b in zip(*m)),arb(0))/len(m[0])
        best=max(gaps.values()); off=max(gaps,key=gaps.get)
        return dict(materials=np.asarray(m).tolist(),gains_zeta_ri=[[0.,0.],[0.,0.]],rms_lower=lower(rm.sqrt()),
                    old_score_upper=old,one_coherent_noise_radius=self.noise,
                    peak_gap_ratio_upper=best,peak_offset=off,per_frequency_gap_ratio_upper=gaps,
                    strictly_inside_box=bool(np.min(m)>=self.case.box[0] and np.max(m)<=self.case.box[1]),
                    certified_finite_075=bool(max(old)<1 and best<.75 and np.min(m)>=self.case.box[0]
                                               and np.max(m)<=self.case.box[1]))

    def opt(self,pts,maxiter):
        tic=process_time();s0=self.model.factorizations
        aim=self.case.box
        base=np.r_[pts[0],pts[1]]
        a,b=[self.transform(x) for x in pts]
        T=np.block([[a,np.zeros_like(a)],[np.zeros_like(b),b]])
        values=[]; cache=None
        scale=max(.05,rms(*pts))**2
        def objective(y):
            p=base+T@y
            d=p[:6]-p[6:]
            return -np.dot(d,d)/(6*scale),-(T.T@np.r_[d,-d])/(3*scale)
        def c(y):
            nonlocal cache
            if cache is not None and np.array_equal(y,cache[0]):return cache[1:]
            x=base+T@y
            v,J=self.constraints(x[:6],x[6:])
            vals=np.r_[v,x-aim[0],aim[1]-x]
            jac=np.vstack((J@T,T,-T))
            if min(vals)>-1e-8:values.append((rms(x[:6],x[6:]),x.reshape(2,6)))
            cache=y.copy(),vals,jac
            return vals,jac
        res=minimize(objective,np.zeros(12),jac=True,method='SLSQP',constraints=[
            dict(type='ineq',fun=lambda y:c(y)[0],jac=lambda y:c(y)[1])],options=dict(maxiter=maxiter,ftol=1e-10))
        v,_=c(res.x)
        return values,dict(status=int(res.status),iterations=int(res.nit),minimum_constraint=float(min(v)),
                           feasible_iterates=len(values),best_rms=max((x[0] for x in values),default=None),
                           cpu_s=process_time()-tic,forward_derivative_calls=self.model.factorizations-s0)


def run(case, archive, gains, choices, maxiter):
    r=Search(case,NODES)
    tic=process_time(); attempts=[]; tried=[]
    known=gains['cases'][case.case_id]['modes']['known_unit']['best']
    incoming=[(x['witness']['materials'], name) for name,x in archive['cases'][case.case_id]['actions'].items()]
    if known:
        incoming.insert(0,(known['materials'],'known_unit'))
    fit=Fit(r.model,[0.],[case.old])
    points=[]; pairs=[]
    for m,name in incoming:
        m=np.array(m)
        points.extend(m)
        pairs.append((m,name+'_old_pair'))
        pairs.append((np.array([fit.solve(x,120)[0] for x in m]) ,name+'_known_gain_refit'))
    pairs.sort(key=lambda pair:-rms(*pair[0]))
    neutral=[np.full(6,12.),*list(case.box[0]+(case.box[1]-case.box[0])*qmc.Sobol(6,scramble=True,seed=29118).random_base2(2))]
    oldstarts=[np.array(x) for x in choices[case.case_id]['old_starts']]
    neutral.extend(oldstarts)
    points.extend(fit.solve(x,120)[0] for x in neutral)
    points.extend(x for m,name in pairs for x in m if name.endswith('refit'))
    distinct=[]
    for p in points:
        score=r.old_cost(p)
        if score<.997 and all(rms(p,x)>1e-4 for x in distinct):distinct.append(p)
    for p in distinct:
        for j in [-1,-2]:
            a,b=r.line(p,r.transform(p)[:,j]); d=rms(a,b)
            if d>1e-6:pairs.append((np.array([a,b]),f'known_anchor_weak_{j}'))
    best=None
    ordered=sorted(pairs,key=lambda p:-rms(*p[0]))
    lines=sorted([x for x in ordered if 'weak' in x[1]],key=lambda p:-rms(*p[0]))
    seeds=[]
    for group in (lines[:3],ordered[:5]):
        for pts,name in group:
            if all(np.linalg.norm(pts-p)>1e-5 and np.linalg.norm(pts[::-1]-p)>1e-5 for p,_ in seeds):seeds.append((pts,name))
    for pts,name in seeds:
        result,work=r.opt(pts,maxiter)
        attempts.append(dict(source=name,initial_rms=rms(*pts),**work))
        for d,m in sorted([(rms(*pts),pts),*result],key=lambda v:-v[0])[:4]:
            if d<=0 or (best is not None and d<=best['rms_lower']):continue
            v=r.strict(m)
            tried.append(dict(source=name,**v))
            if v['certified_finite_075']:
                best=v;break
        if best is not None and best['rms_lower']>1.:break
    return dict(case_id=case.case_id,best=best,attempts=attempts,verified_trials=tried,
                eligible_anchors=len(distinct),anchor_scores=[r.old_cost(p) for p in distinct],
                full_forward_derivative_calls=r.model.factorizations,arb_solves=r.phy.solves,
                arb_kernel_calls=r.phy.kernel_calls,cpu_s=process_time()-tic)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',action='append')
    p.add_argument('--maxiter',type=int,default=110)
    a=p.parse_args(); r=ROOT/'results'; known={x.case_id:x for x in load_cases()}
    archive=json.loads((r/'pivot_diameter_witnesses.json').read_text())
    gains=json.loads((r/'pivot_menu_gain_controls.json').read_text())
    choices={x['case_id']:x for x in json.loads((r/'measure_choices.json').read_text())['cases']}
    result={}
    for name in a.case or SELECTED:
        row=run(known[name],archive,gains,choices,a.maxiter)
        result[name]=row
        (r/'pivot_coherent_unit.json').write_text(json.dumps(dict(cases=result,target_gap_ratio=.75,
               algorithm_nodes=NODES,verified_offsets=sorted(set(NODES)|set(OFFSETS))),indent=2,allow_nan=False)+'\n')
        print(name,row['eligible_anchors'],row['cpu_s'],None if row['best'] is None else
              (row['best']['rms_lower'],row['best']['peak_gap_ratio_upper']),flush=True)
        if row['best'] is not None and row['best']['rms_lower']>1:break


if __name__=='__main__':
    main()

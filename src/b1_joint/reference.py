"""Coherent mixed-measure radial posterior and direct truncated CDF integrals.

All summaries use one continuous coefficient prior and one frozen response
density. Inactive geometry coordinates are integrated analytically. The same
prior still has no atom at a complete image because complex epsilon is free.
"""
from __future__ import annotations
import heapq
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import logsumexp
from scipy.optimize import brentq
from .material import BOUNDS,WIDTH,RMAX,SCALES,event_intervals,geometry,event_ids
from .likelihood import log_likelihood,context_loglik

def unit_rule(n):
    x,w=leggauss(n);return (x+1)/2,w/2

class Reference:
    def __init__(self,response):
        self.response=response;self.cache={}
    def cell(self,event,order,known_eps=None,cut=None,box_override=None):
        u,w=unit_rule(order);ci,si=divmod(event,3); intervals=event_intervals()
        boxes=[intervals[0][ci],intervals[1][si],tuple(BOUNDS[2]),tuple(BOUNDS[3])]
        if box_override is not None:boxes=box_override
        if known_eps is not None:boxes[2]=(known_eps[0],known_eps[0]);boxes[3]=(known_eps[1],known_eps[1])
        coef=[];weights=[];physical=[]
        for j,box in enumerate(boxes):
            if j<2 and box is None:
                if cut is not None and cut[0]==j and cut[1]<0:return None
                values=np.array([0.]);ww=np.array([.5]);pp=values
            elif j>=2 and known_eps is not None:
                values=np.array([box[0]]);ww=np.array([1.]);pp=values
                if cut is not None and cut[0]==j and cut[1]<box[0]:return None
            else:
                lo,hi=box
                if cut is not None and cut[0]==j:
                    threshold=cut[1]
                    if j<2:
                        if threshold<=0:return None
                        if threshold<RMAX[j]:lo=max(lo,-(1-(threshold/SCALES[j])**2)**-2)
                    else:hi=min(hi,threshold)
                if hi<=lo:return None
                values=lo+(hi-lo)*u;ww=w*(hi-lo)/WIDTH[j]
                pp=SCALES[j]*np.sqrt(1-1/np.sqrt(-values)) if j<2 else values
            coef.append(values);weights.append(ww);physical.append(pp)
        grids=np.meshgrid(*coef,indexing='ij');x=np.stack([a.ravel() for a in grids],axis=1)
        gw=np.meshgrid(*weights,indexing='ij');prior=np.prod(np.stack(gw),axis=0).ravel()
        gphys=np.meshgrid(*physical,indexing='ij');v=np.stack([a.ravel() for a in gphys],axis=1)
        f1=self.response.mesh(physical,1.6);f2=self.response.mesh(physical,2.4)
        return dict(x=x,physical=v,prior=prior,F1=f1,F2=f2,event=event,boxes=boxes)
    def bank(self,order,known_eps=None):
        key=(order,tuple(known_eps) if known_eps is not None else None)
        if key not in self.cache:
            cells=[self.cell(j,order,known_eps) for j in range(12)]
            self.cache[key]=dict(x=np.concatenate([c['x'] for c in cells]),physical=np.concatenate([c['physical'] for c in cells]),prior=np.concatenate([c['prior'] for c in cells]),F1=np.concatenate([c['F1'] for c in cells]),F2=np.concatenate([c['F2'] for c in cells]),events=np.concatenate([np.full(len(c['x']),j) for j,c in enumerate(cells)]),cells=cells)
        return self.cache[key]
    def solve(self,y,sigma,order=13,condition='C',known_eps=None,context=None):
        b=self.bank(order,known_eps);F=b['F1'] if condition!='D' else np.concatenate((b['F1'],b['F2']),axis=1)
        ll=log_likelihood(F,y,sigma,0. if condition in ('A','B') else .15)
        if context is not None:ll+=context_loglik(b['x'],context)
        logw=ll+np.log(b['prior']);logz=logsumexp(logw);mass=np.exp(logw-logz)
        v=b['physical'];mu=mass@v;cov=((v-mu).T*mass)@(v-mu)
        logp=np.array([logsumexp(logw[b['events']==j])-logz for j in range(12)])
        p=np.exp(logp)
        cond=[]
        for j in range(12):
            ids=b['events']==j;w=np.exp(logw[ids]-logsumexp(logw[ids]));a=v[ids];m=w@a;cv=((a-m).T*w)@(a-m)
            cond.append(dict(event_id=j,mean=m,SD=np.sqrt(np.maximum(np.diag(cv),0.)),covariance=cv))
        return dict(event_mass=p,log_event_mass=logp,mean=mu,SD=np.sqrt(np.maximum(np.diag(cov),0.)),covariance=cov,conditional=cond,log_evidence=float(logz),mass=mass,bank=b,loglik=ll,order=order,condition=condition,known_eps=known_eps,context=context,y=np.array(y),sigma=np.array(sigma))
    def cdf(self,result,axis,cut,order=None,event=None):
        if order is None:order=result['order']
        eids=range(12) if event is None else [event]
        logs=[]
        for j in eids:
            c=self.cell(j,order,result['known_eps'],(axis,cut))
            if c is None:continue
            F=c['F1'] if result['condition']!='D' else np.concatenate((c['F1'],c['F2']),axis=1)
            ll=log_likelihood(F,result['y'],result['sigma'],0. if result['condition'] in ('A','B') else .15)
            if result['context'] is not None:ll+=context_loglik(c['x'],result['context'])
            logs.append(logsumexp(ll+np.log(c['prior'])))
        if not logs:return 0.
        den=result['log_evidence']+(result['log_event_mass'][event] if event is not None else 0.)
        return float(np.exp(logsumexp(logs)-den))
    def continuous_summary(self,result,grid_size=17,conditional_events=()):
        bounds=np.array([[0.,RMAX[0]],[0.,RMAX[1]],[2.2,2.8],[.04,.16]])
        outputs=[]
        for event in [None]+list(conditional_events):
            curves=[];quant=[];grids=[];memo={}
            for j,(lo,hi) in enumerate(bounds):
                def cdf(t):
                    key=(j,float(t))
                    if key not in memo:memo[key]=self.cdf(result,j,t,event=event)
                    return memo[key]
                ts=np.linspace(lo,hi,grid_size);cs=np.array([cdf(t) for t in ts]);grids.append(ts);curves.append(cs)
                if result['known_eps'] is not None and j>=2:
                    q=np.full(3,result['known_eps'][j-2])
                else:
                    atom=cdf(lo)
                    q=np.array([lo if atom>=p else brentq(lambda t:cdf(t)-p,lo,hi,xtol=2e-9) for p in [.05,.5,.95]])
                quant.append(q)
            quant=np.array(quant)
            source=result if event is None else result['conditional'][event]
            mean=np.array(source['mean']);cov=np.array(source['covariance']);width=bounds[:,1]-bounds[:,0]
            lower=.05*(np.maximum(mean-quant[:,0],0)**2+np.maximum(quant[:,2]-mean,0)**2)
            outputs.append(dict(event_id=event,mean=mean,SD=np.array(source['SD']),quantiles_05_50_95=quant,CDF_grid=grids,CDF=curves,coherence_pass=bool(np.all(np.diag(cov)+1e-12*width**2>=lower)),variance_lower_bound=lower,measure='same continuous coefficient prior, geometry pushforward atoms; direct truncated integrals'))
        return outputs
    def particles(self,result,n,rng):
        cdf=np.cumsum(result['mass']);cdf[-1]=1.
        # Randomized stratification is a quadrature resampling scheme, not iid.
        ids=np.searchsorted(cdf,(np.arange(n)+rng.random(n))/n)
        x=result['bank']['x'][ids].copy()
        inactive=x[:,:2]>=-1
        x[:,:2]=np.where(inactive,rng.uniform(-1,1,size=(n,2)),x[:,:2])
        return x
    def adaptive(self,result,max_cells=512,relative=2e-6,cut=None,event=None):
        """Independent 3/7-point partitioned cubature, with a posterior scale.

        Error estimates are heuristic and reported; never a uniform bound.
        The adaptive partition differs from tensor quadrature nodes.
        """
        shift=float(np.max(result['loglik']));heap=[];value=np.zeros((12,21));error=np.zeros_like(value);counter=0
        def evaluate(j,box,n):
            c=self.cell(j,n,result['known_eps'],cut,box)
            if c is None:return np.zeros(21)
            F=c['F1'] if result['condition']!='D' else np.concatenate((c['F1'],c['F2']),axis=1)
            ll=log_likelihood(F,result['y'],result['sigma'],0. if result['condition'] in ('A','B') else .15)
            if result['context'] is not None:ll+=context_loglik(c['x'],result['context'])
            w=c['prior']*np.exp(ll-shift);v=c['physical'];normalized=(v-np.array([0,0,2.2,.04]))/np.array([*RMAX,.6,.12])
            feats=np.column_stack((np.ones(len(v)),normalized,(normalized[:,:,None]*normalized[:,None,:]).reshape(len(v),16)))
            return w@feats
        def make(j,box):
            nonlocal counter
            lo=evaluate(j,box,3);hi=evaluate(j,box,7);err=abs(hi-lo);counter+=1
            return (-float(err.max()),counter,j,box,hi,err)
        intervals=event_intervals()
        for j in range(12) if event is None else [event]:
            ci,si=divmod(j,3);box=[intervals[0][ci],intervals[1][si],tuple(BOUNDS[2]),tuple(BOUNDS[3])]
            if result['known_eps'] is not None:
                box[2]=(result['known_eps'][0],result['known_eps'][0]);box[3]=(result['known_eps'][1],result['known_eps'][1])
            q=make(j,box);heapq.heappush(heap,q);value[j]+=q[4];error[j]+=q[5]
        def done():return float(error.sum(axis=0).max())<=1e-12+relative*float(value[:,0].sum())
        while not done() and len(heap)<max_cells:
            _,_,j,box,hi,err=heapq.heappop(heap);value[j]-=hi;error[j]-=err
            spans=[0 if b is None else (b[1]-b[0])/WIDTH[i] for i,b in enumerate(box)]
            axis=int(np.argmax(spans));a,b=box[axis];mid=(a+b)/2
            for ab in [(a,mid),(mid,b)]:
                part=list(box);part[axis]=ab;q=make(j,part);heapq.heappush(heap,q);value[j]+=q[4];error[j]+=q[5]
        z=value[:,0].sum();p=value[:,0]/z;feat=value.sum(axis=0)/z;scale=np.array([*RMAX,.6,.12]);offset=np.array([0,0,2.2,.04]);mu=offset+feat[1:5]*scale;cov=(feat[5:].reshape(4,4)-np.outer(feat[1:5],feat[1:5]))*np.outer(scale,scale)
        return dict(converged=bool(done()),cells=len(heap),event_mass=p,mean=mu,SD=np.sqrt(np.maximum(np.diag(cov),0)),log_evidence=float(np.log(z)+shift),estimated_relative_error=float(error.sum(axis=0).max()/z),certified=False)

    def adaptive_physical(self,result,max_cells=2048,relative=2e-6,prior_measure='coefficient_uniform'):
        """Independent smooth-coordinate 5/9 cubature, SAME prior measure.

        The induced active geometry density is r/a²(1-(r/a)²)^(-3).
        Its integral is .5, while the inactive atom also has mass .5.
        Transforming coordinates removes the coefficient-root cusp without
        altering the physical prior. Leaves define a computable continuous
        defensive importance proposal, not a discrete posterior prior.
        """
        intervals=event_intervals();scale=np.array([*RMAX,.6,.12]);offset=np.array([0,0,2.2,.04]);shift=float(max(result['loglik']))
        heap=[];value=np.zeros((12,21));error=np.zeros_like(value);counter=0
        def initial(j):
            c,s=divmod(j,3);boxes=[]
            for i,coef in enumerate([intervals[0][c],intervals[1][s]]):
                boxes.append(None if coef is None else tuple(SCALES[i]*np.sqrt(1-1/np.sqrt(-np.array(coef[::-1])))))
            return boxes+[tuple(BOUNDS[2]),tuple(BOUNDS[3])]
        def evaluate(box,n):
            u,w=unit_rule(n);axes=[];weights=[]
            for j,b in enumerate(box):
                if b is None:axes.append(np.array([0.]));weights.append(np.array([.5]));continue
                lo,hi=b;v=lo+(hi-lo)*u;ww=w*(hi-lo)
                if j<2:
                    if prior_measure=='coefficient_uniform':ww*=v/SCALES[j]**2*(1-(v/SCALES[j])**2)**-3
                    elif prior_measure=='radius_uniform_with_atoms':ww*=.5/RMAX[j]
                    else:raise ValueError('undeclared integration prior')
                else:ww/=WIDTH[j]
                axes.append(v);weights.append(ww)
            mesh=np.meshgrid(*axes,indexing='ij');v=np.stack([a.ravel() for a in mesh],axis=1)
            wg=np.meshgrid(*weights,indexing='ij');pw=np.prod(np.stack(wg),axis=0).ravel()
            f=self.response.mesh(axes,1.6)
            if result['condition']=='D':f=np.concatenate((f,self.response.mesh(axes,2.4)),axis=1)
            ll=log_likelihood(f,result['y'],result['sigma'],0. if result['condition'] in ('A','B') else .15)
            if result['context'] is not None:
                volume=v[:,0]**3+np.where(v[:,1]>0,(.65+v[:,1])**3-(.65-v[:,1])**3,0.)
                ll+=-.5*((volume-result['context'])/.01)**2-np.log(.01*np.sqrt(2*np.pi))
            mass=pw*np.exp(ll-shift);z=(v-offset)/scale;features=np.column_stack((np.ones(len(v)),z,(z[:,:,None]*z[:,None,:]).reshape(len(v),16)))
            return mass@features
        def make(j,box):
            nonlocal counter
            low=evaluate(box,5);high=evaluate(box,9);err=abs(high-low);counter+=1
            return (-float(err.max()),counter,j,box,high,err)
        for j in range(12):
            box=initial(j);parts=[box]
            if result['context'] is not None and box[1] is not None:
                # Start near the declared sensor ridge. Splitting the integration
                # domain changes neither the prior nor its likelihood/target.
                lo,hi=box[1];cuts=[lo,hi]
                radii=[0.,RMAX[0]] if box[0] is None else [*box[0],sum(box[0])/2]
                for radius in radii:
                    for noise in [-4,-2,0,2,4]:
                        h=(result['context']-radius**3+noise*.01)/(6*.65**2)
                        if lo<h<hi:cuts.append(h)
                cuts=sorted(set(cuts));parts=[]
                for a,b in zip(cuts[:-1],cuts[1:]):
                    p=list(box);p[1]=(a,b);parts.append(p)
            for part in parts:
                q=make(j,part);heapq.heappush(heap,q);value[j]+=q[4];error[j]+=q[5]
        def done():return error.sum(axis=0).max()<=1e-13+relative*value[:,0].sum()
        while not done() and len(heap)<max_cells:
            _,_,j,box,hi,err=heapq.heappop(heap);value[j]-=hi;error[j]-=err
            spans=[0 if b is None else (b[1]-b[0])/scale[i] for i,b in enumerate(box)]
            axis=int(np.argmax(spans));lo,up=box[axis];mid=(lo+up)/2
            for interval in [(lo,mid),(mid,up)]:
                new=list(box);new[axis]=interval;q=make(j,new);heapq.heappush(heap,q);value[j]+=q[4];error[j]+=q[5]
        norm=value[:,0].sum();stats=value.sum(axis=0)/norm;mu=offset+scale*stats[1:5];cov=(stats[5:].reshape(4,4)-np.outer(stats[1:5],stats[1:5]))*np.outer(scale,scale)
        leaves=[dict(event_id=q[2],box=q[3],scaled_mass=float(q[4][0])) for q in heap]
        return dict(converged=bool(done()),cells=len(heap),event_mass=value[:,0]/norm,mean=mu,SD=np.sqrt(np.maximum(np.diag(cov),0)),log_evidence=float(np.log(norm)+shift),estimated_relative_error=float(error.sum(axis=0).max()/norm),leaves=leaves,coordinate='physical radii, with declared prior density and atoms',prior_measure=prior_measure,certified=False)

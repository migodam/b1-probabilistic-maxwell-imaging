"""One bounded roundoff repair: interval normalization of finite log potentials.

Scope is deliberately narrow: input float64 log likelihoods and prior masses
are interpreted as exact finite numerical-model constants. It does NOT bound
Maxwell discretization, the rounding used to form the Gaussian potential, or
an unseen forward model. Decimal.exp is correctly rounded; its next neighbors
enclose that operation, and positive arithmetic uses directed rounding.
"""
from decimal import Decimal,Context,ROUND_FLOOR,ROUND_CEILING,ROUND_HALF_EVEN
import numpy as np

def _contexts(precision):
    kw=dict(prec=precision,Emin=-9999999,Emax=9999999)
    return (Context(rounding=ROUND_FLOOR,**kw),Context(rounding=ROUND_CEILING,**kw),
            Context(rounding=ROUND_HALF_EVEN,**kw))

def _mass_intervals(prior,loglike,known,precision):
    low,up,near=_contexts(precision);lo=[];hi=[]
    for p,l,k in zip(prior,loglike,known):
        a=Decimal.from_float(float(p))
        if not k:lo.append(Decimal(0));hi.append(a);continue
        if not np.isfinite(l) or l>0:raise ValueError('finite nonpositive log potentials required')
        e=near.exp(Decimal.from_float(float(l)))
        emin=max(Decimal(0),near.next_minus(e));emax=min(Decimal(1),near.next_plus(e))
        lo.append(low.multiply(a,emin));hi.append(up.multiply(a,emax))
    return lo,hi,(low,up,near)

def _prefix(values,context):
    ans=[Decimal(0)]
    for v in values:ans.append(context.add(ans[-1],v))
    return ans

def _outward(value,upper):
    f=float(value)
    if (upper and Decimal.from_float(f)<value) or (not upper and Decimal.from_float(f)>value):
        f=np.nextafter(f,np.inf if upper else -np.inf)
    return min(1.,max(0.,float(f)))

def normalized_intervals(prior,loglike,known,precision=60):
    prior=np.asarray(prior,float);loglike=np.asarray(loglike,float);known=np.asarray(known,bool)
    if prior.shape!=loglike.shape or known.shape!=prior.shape or np.any(prior<=0):raise ValueError('positive matching finite-prior vectors required')
    lo,hi,(low,up,near)=_mass_intervals(prior,loglike,known,precision)
    lp,hp=_prefix(lo,low),_prefix(hi,up)
    ls,hs=_prefix(lo[::-1],low)[::-1],_prefix(hi[::-1],up)[::-1]
    outlo=[];outhi=[]
    for i,(a,b) in enumerate(zip(lo,hi)):
        other_hi=up.add(hp[i],hs[i+1]);other_lo=low.add(lp[i],ls[i+1])
        if b==0:lower=upper=Decimal(0)
        elif other_hi==0:lower=upper=Decimal(1)
        else:
            den=up.add(a,other_hi);lower=low.divide(a,den) if a else Decimal(0)
            upper=up.divide(b,low.add(b,other_lo))
        outlo.append(_outward(lower,False));outhi.append(_outward(upper,True))
    un=sum((Decimal.from_float(float(p)) for p,k in zip(prior,known) if not k),Decimal(0))
    # sum() above uses the ambient context; redo the small positive sum directed.
    un=_prefix([Decimal.from_float(float(p)) for p,k in zip(prior,known) if not k],up)[-1]
    if un:
        denom=low.add(un,lp[-1]);unknown=up.divide(un,denom)
    else:unknown=Decimal(0)
    return dict(lower=np.array(outlo),upper=np.array(outhi),unknown_mass_upper=_outward(unknown,True),
        precision=precision,label='DIRECTED_NORMALIZATION_CONDITIONAL_ON_F64_LOG_POTENTIALS',
        forward_and_potential_roundoff_certified=False,continuum_certified=False)

def independent_reference(prior,loglike,precision=80):
    _,_,ctx=_contexts(precision)
    values=[ctx.multiply(Decimal.from_float(float(p)),ctx.exp(Decimal.from_float(float(l)))) for p,l in zip(prior,loglike)]
    normal=_prefix(values,ctx)[-1]
    return [ctx.divide(v,normal) for v in values]

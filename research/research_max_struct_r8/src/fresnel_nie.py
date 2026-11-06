"""Source-specific measured-data independent real-material NIE/LSIE adaptation.

Implements NIE2016 Eq3 and Eq14 LS modified contrast, sourcewise PR+ and
native projected Fourier coordinates. This is a new independent implementation,
not an author-code/mature-software reproduction. Real nonnegative material
box and finite iteration ceiling are declared adaptations. No extra prior.
"""
import time
from collections import Counter
import numpy as np
from scipy.linalg import svd
from research.research_max_struct_r8.src.fresnel_cie_fft import NativeFourierAPCS,SourcePorts
from research.research_max_struct_r8.src.paper_physics import EqualAreaGreen


def run_paper_tsom(model,Y,m0,stages,lsie=False,rank=15,maxiter=500):
    t,c=time.perf_counter(),time.process_time();counts=Counter();E,S=model.E,SourcePorts(model.Ss);ns=E.shape[1]
    geom=model.paper_geometry;green=EqualAreaGreen(geom['coords'],model.shape,geom['radius'],geom['k'],counts)
    V=[];JD=np.zeros_like(E)
    for source in range(ns):
        U,s,Vh=svd(S.Ss[source],full_matrices=False,check_finite=False);v=Vh[:rank].conj().T;V.append(v)
        JD[:,source]=v@((U[:,:rank].conj().T@Y[:,source])/s[:rank]);counts['receiver_SVD']+=1
    GJD=green.apply(JD);SJD=S@JD
    wd=1/np.maximum(np.sum(abs(Y)**2,axis=0),1e-300);ws=1/np.sum(abs(E)**2,axis=0)
    old={};m=np.array(m0).copy();history=[];states=[]
    for stage_index,stage in enumerate(stages):
        beta=None if lsie else float(stage['beta_material']);basis=NativeFourierAPCS(model.shape,stage['mf'],V,counts)
        A=np.array([old.get(tuple(f),np.zeros(ns,complex)) for f in basis.frequencies]);og=od=None
        for it in range(maxiter):
            factor=np.ones(len(m)) if beta is None else beta/(1+beta*m)
            def H(a):
                f=basis.apply(a);return factor[:,None]*(f-m[:,None]*green.apply(f))
            def HH(r):
                f=factor[:,None]*r;return basis.adjoint(f-green.adjoint(m[:,None]*f))
            f=basis.apply(A);J=JD+f;Q=E+GJD+green.apply(f);D=S@J-Y;R=factor[:,None]*(J-m[:,None]*Q)
            g=basis.adjoint(S.conj().T@(D*wd))+HH(R*ws)
            pr=0 if og is None else np.maximum(0,np.sum(np.real(g.conj()*(g-og)),axis=0)/np.maximum(np.sum(abs(og)**2,axis=0),1e-300))
            d=-g if od is None else -g+od*pr;bad=np.sum(np.real(g.conj()*d),axis=0)>=0;d[:,bad]=-g[:,bad]
            dJ=basis.apply(d);dD=S@dJ;dR=H(d)
            den=wd*np.sum(abs(dD)**2,axis=0)+ws*np.sum(abs(dR)**2,axis=0)
            alpha=-np.sum(np.real(g.conj()*d),axis=0)/np.maximum(den,1e-300);change=d*alpha
            delta=float(np.sqrt(np.mean(np.sum(abs(change)**2,axis=0)/np.maximum(np.sum(abs(A)**2,axis=0),1e-300))))
            A+=change;f=basis.apply(A);J=JD+f;Q=E+GJD+green.apply(f)
            T=J if beta is None else beta*J;B=Q if beta is None else Q+beta*J
            r=np.sum(ws*np.real(B.conj()*T),axis=1)/np.maximum(np.sum(ws*abs(B)**2,axis=1),1e-300)
            if beta is None:m=np.clip(r,0,20)
            else:
                r=np.clip(r,0,beta*20/(1+beta*20));m=r/(beta*(1-r))
            og,od=g,d;counts['PR_steps']+=1
            history.append(dict(stage=stage_index,iteration=it+1,delta2D=delta))
            if it>=1 and delta<float(stage.get('paper_delta',.001)):break
        state=model.state(m);states.append(dict(m=m.copy(),Y_prediction=state['Y'].copy(),stage=stage_index))
        old={tuple(f):a.copy() for f,a in zip(basis.frequencies,A)}
    return dict(m=m,history=history,stage_states=states,receipt=dict(status='COMPLETED',wall_seconds=time.perf_counter()-t,cpu_seconds=time.process_time()-c,counts=dict(counts),implementation='independent real-material paper NIE/LSIE adaptation'))

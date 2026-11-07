"""Offline squared-error geometry; no inverse solver, data selection or GT feedback."""
import importlib.util
from pathlib import Path
import numpy as np
from scipy.linalg import svd, lu_factor, lu_solve

REF=Path(__file__).resolve().parents[2]/'research_gn_green_error_20261106/src/analyze.py'
spec=importlib.util.spec_from_file_location('historical_green_definitions',REF)
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
energy=ref.energy
BANDS=ref.BANDS

def orth(x,tol=1e-10):
    u,s,_=svd(x,full_matrices=False,check_finite=False)
    return u[:,s>tol]

def physical_current(G,E,m):
    f=lu_factor(np.eye(len(m),dtype=complex)-m[:,None]*G,check_finite=False)
    return lu_solve(f,m[:,None]*E,check_finite=False)

def fourier_columns(shape,mf,components=1):
    axes=[np.unique(np.r_[np.arange(mf),np.arange(n-mf,n)]) for n in shape]
    ids=np.array(np.meshgrid(*axes,indexing='ij')).reshape(len(shape),-1).T
    work=np.zeros(tuple(shape)+(len(ids),),complex)
    work[tuple(ids.T)+(np.arange(len(ids)),)]=1
    F=np.fft.ifftn(work,axes=tuple(range(len(shape))),norm='ortho').reshape(-1,len(ids))
    return np.kron(F,np.eye(components))

class Geometry:
    def __init__(self,G,S,joint=True):
        self.G,self.S=G,S
        self.frames={}
        self.meta=[]
        for family,gi,go in [('complex',G,S),('real',np.vstack([G.real,G.imag]),np.vstack([S.real,S.imag]))]:
            a=ref.spectrum(gi,family+'_G_in');b=ref.spectrum(go,family+'_G_ob')
            for tau in (.03,.1,.3):
                aa,bb=ref.split(a,tau),ref.split(b,tau)
                z=ref.joint_spaces(aa,bb) if joint else dict(inn=aa,ob=bb)
                if joint:
                    vo,vi=bb['bases'][0],aa['bases'][0]
                    z['legacy_second']=orth(vi-vo@(vo.conj().T@vi))
                self.frames[family,tau]=z
                self.meta.append(dict(family=family,tau=tau,in_ranks=[len(i) for i in aa['ids']],
                    ob_ranks=[len(i) for i in bb['ids']],groups=z.get('meta'),joint_computed=joint))
        self.receiver_u,self.receiver_s,self.receiver_vh=svd(S,full_matrices=False,check_finite=False)
        self.affines={}

    def affine(self,Y,rank,shape,mf,components=1):
        key=(rank,tuple(shape),mf,components)
        if key not in self.affines:
            vd=self.receiver_vh[:rank].conj().T
            F=fourier_columns(shape,mf,components)
            U=F-vd@(vd.conj().T@F)
            q=orth(U)
            self.affines[key]=(vd,U,q)
        vd,U,q=self.affines[key]
        JD=vd@((self.receiver_u[:,:rank].conj().T@Y)/self.receiver_s[:rank,None])
        return JD,U,q

    def decompose(self,x,family='complex'):
        x=np.asarray(x);x=x[:,None] if x.ndim==1 else x
        total=energy(x);single=[];two=[];rem=[]
        for tau in (.03,.1,.3):
            z=self.frames[family,tau]
            for name,op in [('G_in',z['inn']),('G_ob',z['ob'])]:
                coeff=op['V'].conj().T@x
                modes=np.sum(abs(coeff)**2,axis=1);out=float(np.sum(op['s']**2*modes))
                for band,ids in zip(BANDS,op['ids']):
                    e=float(np.sum(modes[ids]));oe=float(np.sum(op['s'][ids]**2*modes[ids]))
                    single.append(dict(tau=tau,operator=name,band=band,rank=len(ids),total_energy=total,
                        component_energy=e,fraction=e/total if total else 0,output_energy=oe,output_total_energy=out))
            if 'caps' not in z:continue
            pi=ref.components(z['inn'],x);po=ref.components(z['ob'],x)
            for order in ('in_then_ob','ob_then_in'):
                parts=[]
                for first in range(3):
                    pp=ref.components(z['ob'],pi[first]) if order=='in_then_ob' else ref.components(z['inn'],po[first])
                    for second,p in enumerate(pp):
                        a,b=(first,second) if order=='in_then_ob' else (second,first)
                        parts.append(p)
                        row=dict(tau=tau,definition=order,in_band=BANDS[a],ob_band=BANDS[b],total_energy=total,
                            component_energy=energy(p),fraction=energy(p)/total if total else 0)
                        if tau==.1:
                            row.update(in_output_energy=energy(self.G@p),ob_output_energy=energy(self.S@p))
                        two.append(row)
                if tau==.1:
                    for name,A in [('G_in',self.G),('G_ob',self.S)]:
                        isolated=sum(energy(A@p) for p in parts);actual=energy(A@x)
                        rem.append(dict(tau=tau,definition=order,operator=name,total_energy=total,
                            propagated_actual_energy=actual,propagated_isolated_sum=isolated,
                            signed_cross_energy=actual-isolated))
            residual=x.astype(complex if family=='complex' else float,copy=True)
            for (a,b),q in z['caps'].items():
                p=q@(q.conj().T@x);residual-=p
                two.append(dict(tau=tau,definition='numerical_intersection',in_band=BANDS[a],ob_band=BANDS[b],
                    rank=q.shape[1],total_energy=total,component_energy=energy(p),fraction=energy(p)/total if total else 0))
            rem.append(dict(tau=tau,definition='intersection_remainder',total_energy=total,component_energy=energy(residual),
                fraction=energy(residual)/total if total else 0))
            vo=z['ob']['bases'][0];q=z['legacy_second']
            ev,ef=energy(vo.conj().T@x),energy(q.conj().T@x)
            rem.append(dict(tau=tau,definition='legacy_twofold_chart_not_FFT_TSom',total_energy=total,
                receiver_strong_energy=ev,second_fold_energy=ef,remainder_energy=total-ev-ef,
                receiver_strong_rank=vo.shape[1],second_fold_rank=q.shape[1]))
            for (a,b),q in z['projected'].items():
                ee=energy(q.conj().T@x)
                two.append(dict(tau=tau,definition='overlapping_projected_spans_not_intersections',in_band=BANDS[a],ob_band=BANDS[b],
                    rank=q.shape[1],total_energy=total,component_energy=ee,fraction=ee/total if total else 0))
        return single,two,rem

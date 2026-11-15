"""Sourcewise/global PR-CG and scalar-isotropic real A+ adaptation.
One deterministic exact-quadratic current line step per joint cycle; low-frequency
FFT basis and DPCS fixed per stage. Complex material uses legal 2D optimization.
"""
import time
import numpy as np
from scipy.optimize import minimize
from research.research_max_struct_r8.src.vector_cie3d import VectorCIE3D

class VectorPR3D(VectorCIE3D):
    def current_pr(self,z,Y,metric,beta,mf,global_pr):
        x=self.decode(z);d=np.repeat(x,3);K=np.eye(self.M)-d[:,None]*self.G;rhs=d[:,None]*self.inc
        if metric!='LSIE_TSOM':
            b=(beta-self.s) if metric=='A_TSOM' else beta;t=beta/(1+b*d);K=t[:,None]*K;rhs=t[:,None]*rhs
            if metric=='A_TSOM':
                aa=d/(1+b*d);K=K+aa[:,None]*(self.Hop@K);rhs=rhs+aa[:,None]*(self.Hop@rhs)
        B=self.bases[mf];SB=self.S@B;KB=K@B;Y2=np.sum(abs(Y)**2,axis=0)
        if self.prstate is None or self.prstate['mf']!=mf:
            gam=np.zeros((B.shape[1],self.P),complex);prevw=None;prevdir=None
        else:gam=self.prstate['gam'];prevw=self.prstate['w'];prevdir=self.prstate['direction']
        D=self.S@self.Jd+SB@gam-Y;R=K@self.Jd+KB@gam-rhs
        w=SB.conj().T@(D/Y2)+self.rho*KB.conj().T@(R/self.E2)
        axis=None if global_pr else 0
        if prevw is None:direction=-w
        else:
            pr=np.maximum(0,np.real(np.sum(np.conj(w)*(w-prevw),axis=axis))/np.maximum(np.sum(abs(prevw)**2,axis=axis),1e-300));direction=-w+prevdir*pr
        slope=np.real(np.sum(np.conj(w)*direction,axis=axis));bad=slope>=0
        if global_pr:
            if bad:direction=-w
        else:direction[:,bad]=-w[:,bad]
        slope=np.real(np.sum(np.conj(w)*direction,axis=axis));ds=SB@direction;dk=KB@direction
        den=np.sum(np.sum(abs(ds)**2,axis=0)/Y2+self.rho*np.sum(abs(dk)**2,axis=0)/self.E2) if global_pr else np.sum(abs(ds)**2,axis=0)/Y2+self.rho*np.sum(abs(dk)**2,axis=0)/self.E2
        alpha=-slope/np.maximum(den,1e-300);gam+=direction*alpha;self.prstate=dict(gam=gam,w=w,direction=direction,mf=mf);self.counts['current_PR_steps']=self.counts.get('current_PR_steps',0)+1
        return self.Jd+B@gam
    def aplus(self,x,J,beta):
        assert not self.lossy
        x=x.copy();b=beta-self.s;Q=self.inc+self.G@J;v=-beta*(Q+b*J);u0=beta*J;z=x/(1+b*x);u=u0+np.repeat(z,3)[:,None]*v;R=u+np.repeat(z,3)[:,None]*(self.Hop@u);calls=0;updates=0
        for sweep in range(2):
            for i in range(self.N):
                rows=np.flatnonzero(np.any(abs(self.H[:,3*i:3*i+3])>0,axis=1));rows=np.unique(np.r_[rows,np.arange(3*i,3*i+3)])
                D=np.repeat(z,3)[rows,None]*(self.H[np.ix_(rows,np.arange(3*i,3*i+3))]@v[3*i:3*i+3]);own=np.flatnonzero(np.isin(rows,np.arange(3*i,3*i+3)));D[own]=v[3*i:3*i+3]+self.H[3*i:3*i+3]@u
                C=R[rows]-z[i]*D;q=self.rho/self.P*np.sum(abs(D)**2/self.E2);hh=self.rho/self.P*np.sum(np.conj(C)*D/self.E2)
                aa=q*b.real+np.real(hh*np.conj(b)**2);bb=q+2*np.real(hh*np.conj(b));cc=hh.real;roots=np.roots([aa,bb,cc]) if abs(aa)>1e-30 else np.array([-cc/bb]) if abs(bb)>1e-30 else np.array([])
                candidates=[x[i],0.,4.]+[float(t.real) for t in roots if abs(t.imag)<1e-9 and 0<=t.real<=4]
                e=self.edges;neighbors=np.r_[e[e[:,0]==i,1],e[e[:,1]==i,0]];a=self.lam*len(neighbors)/len(e)+self.l2/self.N;c=-self.lam*np.sum(x[neighbors])/len(e)
                def loss(t):zz=t/(1+b*t);return q*abs(zz)**2+2*np.real(hh*zz)+a*t*t+2*c*t
                new=min(candidates,key=loss);zz=new/(1+b*new);R[rows]+=(zz-z[i])*D;u[3*i:3*i+3]=u0[3*i:3*i+3]+zz*v[3*i:3*i+3];updates+=int(new!=x[i]);x[i]=new;z[i]=zz;calls+=len(candidates)
        self.counts['aplus_arc_candidates']=self.counts.get('aplus_arc_candidates',0)+calls;self.counts['aplus_coordinate_changes']=self.counts.get('aplus_coordinate_changes',0)+updates
        return x
    def optimize_pr(self,Y,kind,outer=48,inner=12,corrector=300):
        self.reset();self.prstate=None;t=time.perf_counter();cpu=time.process_time();z=self.encode(np.full(self.N,.25+(.05j if self.lossy else 0)));initial=z.copy();best=z.copy();bestphi=self.physical(z,Y,False)[0];logs=[];curve=[];trials=[];bounds=[(0,4)]*self.N+([(0,1.2)]*self.N if self.lossy else [])
        def callback(x):
            nonlocal best,bestphi
            phi=self.physical(x,Y,False)[0];trials.append(x.copy());curve.append(dict(phi=phi,physical_LU=self.counts['physical_LU'],cpu_s=time.process_time()-cpu))
            if phi<bestphi:bestphi=phi;best=x.copy()
        prep=self.prepare(Y);metric='LSIE_TSOM' if kind=='LSIE_sourcewise' else 'A_TSOM' if kind in ['A_global','Aplus_global'] else 'CIE_TSOM'
        for it in range(outer):
            mf=1 if it<outer//2 else 2;beta=6 if it<outer//2 else 1;bc=3 if kind=='Dual_sourcewise' and it>=outer//2 else beta;J=self.current_pr(z,Y,metric,bc,mf,kind.endswith('global'))
            if kind=='Aplus_global':
                before=self.material(z,J,metric,beta)[0];candidate=self.aplus(z,J,beta);after=self.material(candidate,J,metric,beta)[0];accepted=bool(after<=before+1e-13);z=candidate if accepted else z;log=dict(nit=2,nfev=2,message='real scalar-isotropic arc stationary proposals + full quadratic prior acceptance; complex unavailable',success=accepted,before=before,after=after)
            else:
                res=minimize(lambda q:self.material(q,J,metric,beta),z,jac=True,bounds=bounds,method='L-BFGS-B',options=dict(maxiter=inner,ftol=1e-12,gtol=1e-8,maxls=25));z=res.x;log=dict(nit=res.nit,nfev=res.nfev,message=str(res.message),success=bool(res.success))
            callback(z);logs.append(dict(stage='joint',cycle=it,mf=mf,beta_material=beta,beta_current=bc,**log))
        raw=best.copy();pre=self.physical(raw,Y,False);res=minimize(lambda q:self.physical(q,Y),raw,jac=True,bounds=bounds,method='L-BFGS-B',callback=callback,options=dict(maxiter=corrector,ftol=1e-13,gtol=1e-9,maxls=30));callback(res.x);f,df,J,E=self.physical(best,Y,False);logs.append(dict(stage='common_long_corrector',nit=res.nit,nfev=res.nfev,message=str(res.message),success=bool(res.success)))
        return dict(initial=self.decode(initial),raw=self.decode(raw),prediction=self.decode(best),J=J,E=E,trials=np.array([self.decode(q) for q in trials])),dict(method=kind,wall_s=time.perf_counter()-t,cpu_s=time.process_time()-cpu,counts=self.counts.copy(),physical_phi=f,physical_data_residual=np.sqrt(df),pre_corrector_phi=pre[0],pre_corrector_residual=np.sqrt(pre[1]),optimizer=logs,curve=curve,subspace=prep,truth_access='NONE during optimization/selection')

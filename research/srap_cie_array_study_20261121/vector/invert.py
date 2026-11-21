"""Public-data-only cold trajectory, deadline-aware audits and physical correctors."""
import argparse,traceback
from scipy.optimize import minimize
from projection import *
class Fourier:
    def __init__(self,mod,V,mf):
        self.n=mod.n;self.V=V
        f=np.rint(np.fft.fftfreq(self.n)*self.n);a=(f>=-mf)&(f<mf)
        self.ids=torch.as_tensor(np.flatnonzero((a[:,None,None]&a[None,:,None]&a[None,None,:]).ravel()),device=DEV)
        self.FV=self.fft(V);self.small=None
    def fft(self,x):
        n=self.n;return torch.fft.fftn(x.reshape(n,n,n,3,-1),dim=(0,1,2),norm='ortho').reshape(n**3,3,-1)[self.ids].reshape(len(self.ids)*3,-1)
    def ortho(self,x):return x-self.V@(self.V.mH@x)
    def apply(self,c):
        n=self.n;grid=torch.zeros((n**3,3,c.shape[1]),device=DEV,dtype=CD);grid[self.ids]=c.reshape(len(self.ids),3,-1)
        return self.ortho(torch.fft.ifftn(grid.reshape(n,n,n,3,-1),dim=(0,1,2),norm='ortho').reshape(n**3*3,-1))
    def adj(self,x):return self.fft(self.ortho(x))
    def fit(self,x):
        b=self.adj(x);D=self.FV.to(torch.complex128)
        if self.small is None:
            g=torch.eye(D.shape[1],dtype=D.dtype,device=DEV)-D.mH@D
            self.small=torch.linalg.pinv(g,hermitian=True,rtol=1e-7)
        return b+(D@(self.small@(D.mH@b.to(D.dtype)))).to(CD)
def readout(J,Q,w,beta):
    n=J.shape[0]//3;V=(Q+beta*J).reshape(n,3,-1);T=(beta*J).reshape(n,3,-1)
    r=((V.conj()*T).real*w).sum((1,2))/((abs(V)**2*w).sum((1,2)).clamp_min(1e-30))
    r=r.clamp(0,4*beta/(1+4*beta));return r/(beta*(1-r))
def basic(mod,Y):
    gram=mod.S.to(torch.complex128)@mod.S.mH.to(torch.complex128);ev,U=torch.linalg.eigh(gram)
    keep=ev>=.01*ev[-1];s=ev[keep].sqrt();uy=U[:,keep];vd=mod.S.mH.to(s.dtype if not s.is_floating_point() else torch.complex128)@uy/s
    V=vd.to(CD);jd=V@((uy.mH@Y.to(uy.dtype)/s[:,None]).to(CD));return V,jd
def optimize(geom,idx,arm):
    out=ROOT/'predictions'/f'{geom}_{idx:02d}_{arm}';out.mkdir(parents=True,exist_ok=True)
    if (out/'SEALED.json').is_file():return
    torch.cuda.synchronize();torch.cuda.reset_peak_memory_stats();clock=time.perf_counter();cpu=time.process_time()
    data=np.load(ROOT/'data'/f'{geom}_{idx:02d}.npz');Y=t(data['Y']);delta=float(data['delta'])
    mod=Model(12,geom);m=torch.full((mod.N,),.25,device=DEV)
    records=[];mats=[];checkpoints=[];stages=[];proj=None;cyc=0
    def elapsed():torch.cuda.synchronize();return time.perf_counter()-clock
    def audit(m,step,stage):
        phi,rel=mod.physical(m,Y);done=elapsed()
        records.append(dict(time=done,phi=phi,data_rel=rel,step=step,stage=stage,cycles=cyc,counts=dict(mod.count)));mats.append(npy(m).copy())
    audit(m,0,-1)
    V,jd=basic(mod,Y);wd=1/normcols(Y);ws=1/normcols(mod.E)
    J=jd.clone();Q=mod.E+mod.G(J);D=mod.S@J-Y;mx=torch.zeros_like(J);mq=torch.zeros_like(J);step=0;oldbasis=None
    for stage,(mf,bm,bj) in enumerate(zip([2,3,4],[6.,3.,1.],[6.,4.,3.])):
        basis=Fourier(mod,V,mf);g=torch.zeros((len(basis.ids)*3,16),device=DEV,dtype=CD)
        if oldbasis is not None:
            pos=torch.searchsorted(basis.ids,oldbasis.ids)
            g.reshape(-1,3,16)[pos]=oldg.reshape(-1,3,16)
        J=jd+basis.apply(g);Q=mod.E+mod.G(J);D=mod.S@J-Y
        oldgrad=None
        for it in range(500):
            d=m.repeat_interleave(3)[:,None];p=bj/(1+bj*d);v=p*p*(J-d*Q)*ws
            grad=basis.adj(mod.S.mH@(D*wd)+v-mod.G(d*v,True));direct=-grad
            if oldgrad is not None:
                beta=((grad.conj()*(grad-oldgrad)).real.sum(0)/normcols(oldgrad).clamp_min(1e-30)).clamp_min(0)
                direct=direct+olddir*beta;bad=(grad.conj()*direct).real.sum(0)>=0;direct=torch.where(bad,-grad,direct)
            u=basis.apply(direct);gu=mod.G(u);su=mod.S@u;hu=p*(u-d*gu)
            alpha=-(grad.conj()*direct).real.sum(0)/(normcols(su)*wd+normcols(hu)*ws).clamp_min(1e-30)
            update=direct*alpha;change=(normcols(update)/normcols(g).clamp_min(1e-30)).mean().sqrt().item()
            g+=update;J+=u*alpha;Q+=gu*alpha;D+=su*alpha;m=readout(J,Q,ws,bm)
            oldgrad=grad;olddir=direct;step+=1;mod.count['PR_steps']+=1
            end=(it==499 or (it>=19 and change<.001));event=(step%16==0 or end)
            if event and arm=='SRAP':
                if proj is None:proj=Projector(mod,V,Y,delta)
                gx=J/6;cx,cq,_=material(gx+mx,Q+mq);ax,aq=proj.affine(cx,cq);mx=proj.route(ax-gx);mq=proj.route(aq-Q)
                sx,_,m=material(gx+mx,Q+mq);g=basis.fit(6*sx-jd);J=jd+basis.apply(g);Q=mod.E+mod.G(J);D=mod.S@J-Y
                oldgrad=None;cyc+=1
            if event:audit(m,step,stage)
            if elapsed()>=180 or end:break
        stages.append(dict(stage=stage,steps=it+1,finished=bool(end),coef_change=change))
        oldbasis=basis;oldg=g
        if elapsed()>=180:break
    raw_s=elapsed();raw_cpu=time.process_time()-cpu;totalcounts=dict(mod.count);corrlog=[]
    selected={};outm={}
    for limit in (60,180):
        eligible=[i for i,r in enumerate(records) if r['time']<=limit]
        prefixes=[i for i,r in enumerate(records) if r['time']<=.8*limit]
        ibest=min(eligible,key=lambda i:records[i]['phi']);iseed=min(prefixes,key=lambda i:records[i]['phi'])
        key=f'raw_{limit}';selected[key]=records[ibest];outm[key]=mats[ibest]
        allowance=.2*limit;start=time.perf_counter();cost0=Counter(mod.count);base=mats[iseed];best=base.copy();bestphi=records[iseed]['phi'];best_time=0.;good=0
        trials=[];last=[None,None];over=False
        class Deadline(Exception):pass
        def objective(a):
            nonlocal over
            if time.perf_counter()-start>=allowance:over=True;raise Deadline
            value,grad=mod.physical(t(a).real,Y,True);duration=time.perf_counter()-start
            trials.append(dict(time=duration,phi=value))
            if duration>allowance:over=True;raise Deadline
            last[:]=[a.copy(),value];return value,grad
        def accepted(a):
            nonlocal good,best,bestphi,best_time
            if last[0] is None or not np.array_equal(last[0],a):raise ArithmeticError('callback cache mismatch')
            duration=time.perf_counter()-start
            if duration>allowance:raise Deadline
            good+=1
            if last[1]<bestphi:best=a.copy();bestphi=last[1];best_time=duration
        try:opt=minimize(objective,base.astype(float),method='L-BFGS-B',jac=True,bounds=[(0,4)]*mod.N,callback=accepted,options=dict(maxiter=20,maxls=12,ftol=1e-12,gtol=1e-8));reason=str(opt.message)
        except Deadline:reason='deadline'
        correct_time=time.perf_counter()-start
        key=f'light_{limit}';selected[key]=dict(seed=records[iseed],phi=bestphi,correction_s=best_time,accepted=good);outm[key]=best
        corrlog.append(dict(budget=limit,seconds=correct_time,cap=allowance,num_accepted=good,time_exceeded=over,reason=reason,trials=trials,counts=dict(mod.count-cost0)))
    np.savez_compressed(out/'materials.npz',**outm)
    write(out/'SEALED.json',dict(case=idx,array=geom,method=arm,raw_s=raw_s,raw_cpu=raw_cpu,raw_counts=totalcounts,total_s=elapsed(),cpu=time.process_time()-cpu,stages=stages,rank=V.shape[1],cycle_count=cyc,route=proj.info if proj else {},projection_log=proj.logs if proj else [],raw_audits=records,selected=selected,correctors=corrlog,gpu_peak_scope='per_arm_reset',peak_gpu_bytes=torch.cuda.max_memory_allocated()))
    print(geom,idx,arm,step,cyc,round(raw_s,2),[x['num_accepted'] for x in corrlog],flush=True)
if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('--first',type=int,default=0);pa.add_argument('--last',type=int,default=12);a=pa.parse_args()
    for idx in range(a.first,a.last):
        for geom in ('plane_grid','hemi'):
            for arm in ('Dual','SRAP'):
                optimize(geom,idx,arm)

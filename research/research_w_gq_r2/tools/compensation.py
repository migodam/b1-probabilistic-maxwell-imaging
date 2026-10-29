"""Full Maxwell bidirectional information and bounded finite compensation, not online imaging."""
from pathlib import Path
from dataclasses import replace
import argparse,json,gc,time,numpy as np
from scipy.optimize import lsq_linear
from r2.accounting import BudgetBook,plain
from r2.physics import Eval,WallGrid,Target,S1_TX
from r2.worlds import white,whiteJ
ROOT=Path(__file__).resolve().parents[1]
def base():
    x=np.linspace(-1.5,1.5,25);y=np.linspace(.5,1.25,7)
    xx,yy=np.meshgrid((x[:-1]+x[1:])/2,(y[:-1]+y[1:])/2)
    return WallGrid(x,y,3.1+.15*np.cos(2*xx)+.1*(yy-.8),.09+.01*np.cos(xx)),Target([.2,2.],.12,3.,.05)
def information(A,B):
    u,s,vh=np.linalg.svd(A,full_matrices=False);eps=np.finfo(float).eps*max(A.shape)*s[0];rk=int(np.sum(s>eps));U=u[:,:rk]
    rem=B-U@(U.T@B);Kw=A.T@A-A.T@B@np.linalg.solve(B.T@B+np.eye(B.shape[1]),B.T@A)
    precision=np.eye(A.shape[1]+B.shape[1])+np.column_stack((A,B)).T@np.column_stack((A,B));cov=np.linalg.inv(precision)
    return dict(wall_singular=s.tolist(),algebraic_float_rank=rk,rows=A.shape[0],columns=A.shape[1],rank_threshold=float(eps),exact_rank='NOT_CERTIFIED_FROM_FLOAT_MATRIX',data_target_Keffective=rem.T@rem,prior_target_covariance=cov[-B.shape[1]:,-B.shape[1]:],prior_wall_variance_retained=float(np.trace(cov[:A.shape[1],:A.shape[1]])/A.shape[1]),wall_target_cross=cov[:A.shape[1],-B.shape[1]:],proper_prior='untruncated standardized N(0,I) tangent diagnostic, not posterior over boxed nonlinear wall',no_uniform_Jerror_bound=True)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--validate',action='store_true');a=ap.parse_args()
    b=BudgetBook();out=ROOT/'results/compensation';out.mkdir(exist_ok=True);w0,t0=base();n=w0.epsinf.size
    fs=[.85,1.,1.15];scale=np.r_[np.full(n,.5),np.full(n,.05)];ts=np.array([.2,.2,2.,.08,.04])
    evs=[Eval(b,'E2',f,S1_TX,w0.x,w0.y,h=.15 if f>1 else .17) for f in fs]
    if a.validate:
        rows=json.loads((out/'CURVES.json').read_text())
        eligible=[r for r in rows if r['distance_A0']<=1.]
        if not eligible:
            (out/'VALIDATION.json').write_text(json.dumps(dict(status='NO_WITNESS_FOUND')));return
        r=max(eligible,key=lambda r:r['dx']);d=np.load(out/r['endpoint']);w=WallGrid(w0.x,w0.y,d['eps'],d['ell']);t=replace(t0,pose=t0.pose+[r['dx'],0])
        res=[]
        for f in [1.,1.15]:
            e=Eval(b,'E2',f,S1_TX,w0.x,w0.y,wall_refine=2,h=.145,targetstep=.1,max_dofs=120000)
            s0=e.evaluate(w0,t0,label='finite compensation independent base');z0=s0.data.copy();s=e.evaluate(w,t,label='finite compensation independent endpoint');z1=s.data.copy()
            coarse=np.load(out/r['endpoint']);i=fs.index(f);ref0=np.load(out/'BASE.npz')['data'][i];ref1=coarse['nested'][i]
            res.append(dict(f=f,base_error=float(np.linalg.norm(z0-ref0)/.01),endpoint_error=float(np.linalg.norm(z1-ref1)/.01),pair_A0=float(np.linalg.norm((z1-z0)[:2])/.01),pair_all4=float(np.linalg.norm(z1-z0)/.01),dx=r['dx'],branch=r['branch']))
            np.savez_compressed(out/f'VALIDATED_f{f}.npz',base=z0,endpoint=z1,eps=w.epsinf,ell=w.ell)
        (out/'VALIDATION.json').write_text(json.dumps(dict(rows=res,scope='selected finite pair only, empirical independent refinement not exact continuum equality'),indent=2));return
    states=[ev.evaluate(w0,t0,label='bidirectional base') for ev in evs];data=np.asarray([s.data for s in states]);AJ=[];BJ=[]
    for ev,s in zip(evs,states):
        A,B=ev.jacobians(s);AJ.append(A);BJ.append(B)
    np.savez_compressed(out/'BASE.npz',data=data,A=np.asarray(AJ),B=np.asarray(BJ),eps=w0.epsinf,ell=w0.ell,x=w0.x,y=w0.y)
    infos={}
    for name,ac,bc in [('A0',AJ[1].reshape(4,16,-1)[:2].reshape(32,-1),BJ[1].reshape(4,16,-1)[:2].reshape(32,-1)),('A1',AJ[1],BJ[1]),('A2',np.concatenate(AJ),np.concatenate(BJ))]:
        infos[name]=information(whiteJ(ac)*scale[None,:],whiteJ(bc)*ts[None,:])
    # Independent physical directional FD; one observed direction does not bound all columns.
    ev=evs[1];rng=np.random.default_rng(28102922);v=rng.normal(size=n);v/=np.linalg.norm(v)
    fd=[]
    for h in [1e-4,2e-5]:
        wp=WallGrid(w0.x,w0.y,w0.epsinf+h*v.reshape(w0.shape),w0.ell);wm=WallGrid(w0.x,w0.y,w0.epsinf-h*v.reshape(w0.shape),w0.ell)
        dp=ev.evaluate(wp,t0,label='wall FDplus').data.copy();dm=ev.evaluate(wm,t0,label='wall FDminus').data.copy();g=(dp-dm)/(2*h);analytic=(AJ[1][:,:n]@v).reshape(4,16)
        fd.append(dict(kind='wall',h=h,relative=float(np.linalg.norm(g-analytic)/np.linalg.norm(analytic)),white_error=float(np.linalg.norm(g-analytic)/.01)))
        dp=ev.evaluate(w0,replace(t0,pose=t0.pose+[h,0]),label='target FDplus').data.copy();dm=ev.evaluate(w0,replace(t0,pose=t0.pose-[h,0]),label='target FDminus').data.copy();g=(dp-dm)/(2*h);analytic=BJ[1][:,0].reshape(4,16)
        fd.append(dict(kind='target_x',h=h,relative=float(np.linalg.norm(g-analytic)/np.linalg.norm(analytic)),white_error=float(np.linalg.norm(g-analytic)/.01)))
    (out/'INFORMATION.json').write_text(json.dumps(dict(info=infos,FD=fd,parameter_scales=dict(wall_eps=.5,wall_loss=.05,target=ts.tolist()),scope='local information plus proper tangent prior; no exact physical rank certificate'),default=plain,indent=2))
    A=whiteJ(AJ[1].reshape(4,16,-1)[:2].reshape(32,-1))*scale[None,:]
    U,sv,Vh=np.linalg.svd(A,full_matrices=True);rk=np.sum(sv>np.finfo(float).eps*max(A.shape)*sv[0]);null=Vh[int(rk):]
    direction=np.zeros(2*n)
    if len(null):direction=null[0]*.03
    curves=[]
    for branch in [0,1]:
        wall=w0
        for dx in [.005,.02,.05]:
            target=replace(t0,pose=t0.pose+[dx,0]);history=[]
            if branch and dx==.005:
                wall=WallGrid(w0.x,w0.y,w0.epsinf+(direction[:n]*scale[:n]).reshape(w0.shape),w0.ell+(direction[n:]*scale[n:]).reshape(w0.shape))
            for it in range(6):
                s=ev.evaluate(wall,target,label=f'compensation b{branch} x{dx} it{it}');r=white((s.data-data[1])[:2]);norm=float(np.linalg.norm(r));history.append(dict(iteration=it,d=norm,eps_min=float(wall.epsinf.min()),eps_max=float(wall.epsinf.max()),loss_min=float(wall.ell.min()),loss_max=float(wall.ell.max())))
                if norm<1e-6:break
                Ac=ev.wall_jacobian(s).reshape(4,16,-1)[:2].reshape(32,-1);Ar=whiteJ(Ac)*scale[None,:]
                cur=np.r_[wall.epsinf.ravel(),wall.ell.ravel()];lo=(np.r_[np.full(n,1.5),np.zeros(n)]-cur)/scale;hi=(np.r_[np.full(n,6.),np.full(n,.25)]-cur)/scale
                mat=np.vstack((Ar,1e-4*np.eye(2*n)));rhs=np.r_[-r,np.zeros(2*n)]
                sol=lsq_linear(mat,rhs,bounds=(lo,hi),tol=1e-10,max_iter=100,lsmr_tol=1e-12)
                accepted=False
                for alpha in [1.,.5]:
                    z=cur+alpha*scale*sol.x
                    trial=WallGrid(w0.x,w0.y,z[:n].reshape(w0.shape),z[n:].reshape(w0.shape))
                    ss=ev.evaluate(trial,target,label='bounded full-objective corrector trial')
                    nr=float(np.linalg.norm(ss.data[:2]-data[1,:2])/.01)
                    history[-1].setdefault('trials',[]).append(dict(alpha=alpha,d=nr))
                    if nr<norm:wall=trial;accepted=True;break
                if not accepted:break
            nested=np.asarray([e.evaluate(wall,target,label='finite compensation endpoint nested').data.copy() for e in evs]);d0=float(np.linalg.norm((nested-data)[1,:2])/.01);d1=float(np.linalg.norm((nested-data)[1])/.01);d2=float(np.linalg.norm(nested-data)/.01)
            endpoint=f'endpoint_b{branch}_dx{dx}.npz';np.savez_compressed(out/endpoint,eps=wall.epsinf,ell=wall.ell,nested=nested)
            row=dict(branch=branch,dx=dx,distance_A0=d0,distance_A1=d1,distance_A2=d2,material_RMS=float(np.sqrt(np.mean((wall.epsinf-w0.epsinf)**2))),history=history,endpoint=endpoint,stop='local six-corrector cap or small residual; not global optimum',legal=True)
            curves.append(row);(out/'CURVES.json').write_text(json.dumps(curves,default=plain,indent=2));print(json.dumps(row,default=plain),flush=True)
    print(json.dumps(dict(curves=curves,physics=b.snapshot()),default=plain))
if __name__=='__main__':main()

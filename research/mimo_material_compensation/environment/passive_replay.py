"""Enclose a complete passive feedback disk and independently replay its wave operators."""
import argparse,json,sys
from itertools import combinations
from pathlib import Path
from time import perf_counter,process_time
import numpy as np
from flint import acb,arb,ctx
from threadpoolctl import threadpool_limits
from scipy.linalg import lu_factor,lu_solve
from probe import HERE, PRIOR, Circuit, R, at_load
sys.path.insert(0,str(PRIOR/'diagnostics'))
from verified_wave import VerifiedWave,rd,fro,opnorm,up,arb_upper


def lower(x):
    return float(np.nextafter(x-abs(x)*rd(3),-np.inf))


def stored_response(w,v,x):
    """Positive contrast confines the bare current exactly to the supplied active union."""
    ix=w.ids; n=len(ix); peak=float(x.max())
    yy,xx=divmod(ix,w.side)
    g=v.lag[(yy[:,None]-yy)%v.p,(xx[:,None]-xx)%v.p]
    ge=up(n*v.lagerr)
    a=np.eye(n)-x[:,None]*g
    ae=up(peak*ge+rd(2)*(np.sqrt(n)+peak*fro(g)))
    lu=lu_factor(a,check_finite=False)
    q=lu_solve(lu,np.eye(n),check_finite=False)
    qn=opnorm(q)
    defect=up(fro(np.eye(n)-a@q)+rd(n)*(np.sqrt(n)+fro(a)*fro(q))+ae*qn)
    if defect>=1:raise ArithmeticError('bare inverse')
    ainv=up(qn/(1-defect))
    f=x[:,None]*w.b
    j=lu_solve(lu,f,check_finite=False)
    residual=up(fro(f-a@j)+rd(n)*(fro(f)+fro(a)*fro(j))+ae*fro(j)+rd(2)*peak*fro(w.b))
    je=up(ainv*residual)
    z=w.z0-1j*(w.c@j)
    ze=up(opnorm(w.c)*je+rd(n)*(fro(w.z0)+fro(w.c)*fro(j)))
    mat=z+R*np.eye(25);me=up(ze+rd(2)*(fro(z)+5*R))
    inv=np.linalg.solve(mat,np.eye(25));ni=opnorm(inv)
    invdef=up(fro(np.eye(25)-mat@inv)+rd(25)*(5+fro(mat)*fro(inv))+me*ni)
    if invdef>=1:raise ArithmeticError('port inverse')
    s=np.eye(25)-2*R*inv
    serr=up(2*R*ni*invdef/(1-invdef)+rd(3)*(5+2*R*fro(inv)))
    return s,serr,dict(operator_lag_entry_error=v.lagerr,material_solve_enclosure=je,
                  inverse_defect=defect,circuit_inverse_defect=invdef,scattering_operator_error=serr,
                  inverse_bound=ainv,minimum_material=float(min(x)),maximum_material=float(max(x)))


def singular_upper(a):
    u,s,vh=np.linalg.svd(a,full_matrices=False)
    n=len(s)
    un=up(np.sqrt(1+fro(u.conj().T@u-np.eye(n))+rd(n)*(fro(u)**2+np.sqrt(n))))
    vn=up(np.sqrt(1+fro(vh@vh.conj().T-np.eye(n))+rd(n)*(fro(vh)**2+np.sqrt(n))))
    err=up(fro(a-(u*s)@vh)+rd(n+4)*(fro(a)+fro(u)*fro(s)*fro(vh)))
    return up(un*vn*float(s[0])+err)


def load_cover(s,e,k,h,kr):
    """All loads whose distance to the rounded center is at most h + kr."""
    b=s[:16,24];c=s[24,:16];ss=s[24,24]; b0=fro(b); c0=fro(c)
    p=up(abs(ss)*(1+rd(2))+e)
    if p>=1: raise ArithmeticError('pole not excluded in closed feedback disk')
    center_gap=lower(abs(1-k*ss)-rd(6)*(1+abs(ss))-e-kr*p)
    gap=lower(center_gap-h*p)
    if gap<=0: raise ArithmeticError('feedback inverse')
    # Changes in A,b,c,s and unit-circle node rounding.
    center_error=up(e+((b0+c0+e)*e)/center_gap+(b0*c0)*e/center_gap**2+
                    (b0+e)*(c0+e)*kr/center_gap**2+
                    rd(24)*(1+fro(s)+(1+b0)*(1+c0)/center_gap**2))
    arc=up(h*(b0+e)*(c0+e)/gap**2)
    return at_load(s,k),center_error,arc


def sweep(s1,e1,s0,e0,r,steps):
    den=2*r; full=-1.; comp=[];unit=[];deltas=[]
    edge=np.zeros((16,16)); largest_e=0.
    al=up(arb_upper(arb.pi())/steps)
    for i in range(steps):
        kk=acb(0,2*i*arb.pi()/steps).exp();k=complex(kk)
        kr=arb_upper(abs(kk-acb(float(k.real),float(k.imag))))
        a,ea,qa=load_cover(s1,e1,k,al,kr);b,eb,qb=load_cover(s0,e0,k,al,kr)
        diff=a-b
        err=up(ea+eb+qa+qb+rd(3)*(fro(a)+fro(b)))
        full=max(full,up(singular_upper(diff)+err))
        edge=np.maximum(edge,(abs(diff)*(1+rd(3))+err)*(1+rd(4)))
        # Concrete test at the rounded nearby unit-modulus boundary node.
        u,s,vh=np.linalg.svd(diff,full_matrices=False)
        uu=u[:,0];vv=vh[0,:].conj()
        score=abs(np.vdot(uu,diff@vv))
        comp.append(max(0.,lower((score-rd(64)*fro(diff)*fro(uu)*fro(vv))/up(fro(uu)*fro(vv))-ea-eb)))
        largest_e=max(largest_e,err)
    cuts=np.array([(0,*t) for t in combinations(range(1,16),7)])
    rx=np.array([[i for i in range(16) if i not in row] for row in cuts])
    # Max both orientations; exact comparison does not assume floating reciprocity.
    def cover(c,rr):
        val=np.linalg.norm(edge[rr[:,:,None],c[:,None,:]],axis=(1,2))
        return up(max(val)*(1+rd(128))/np.sqrt(8))
    maxcut=max(cover(cuts,rx),cover(rx,cuts))
    qa=up(np.linalg.norm(edge[8:16,:8])*(1+rd(128))/np.sqrt(8))
    qb=up(np.linalg.norm(edge[1::2,::2])*(1+rd(128))/np.sqrt(8))
    # The max-cut envelope also allows changing the termination at each shot.
    u1,ue1,_=load_cover(s1,e1,1.,0.,0.);u0,ue0,_=load_cover(s0,e0,1.,0.,0.)
    raw=u1-u0
    rawerr=up(ue1+ue0+rd(3)*(fro(u1)+fro(u0)))
    def fixed_lower(x):
        f=np.linalg.norm(x); guarded=lower(f/(1+rd(160))-rawerr)
        return max(0.,lower(guarded/(np.sqrt(8)*den)))
    return dict(open_A_lower_over_2r=fixed_lower(raw[8:16,:8]),open_B_lower_over_2r=fixed_lower(raw[1::2,::2]),arbitrary_wave_upper_over_2r=up(full/den),any_uniform_eight_cut_upper_over_2r=up(maxcut/den),
                concrete_wave_lower_over_2r=lower(max(comp)/den),A_upper_over_2r=up(qa/den),B_upper_over_2r=up(qb/den),
                max_operator_error_plus_arc=largest_e,n_arcs=steps,poles_abs_upper=[up(abs(s1[24,24])+e1),up(abs(s0[24,24])+e0)])


def run(side,t,file,steps,bands):
    p=np.load(HERE/file);ix=p['ix'];group=p['group']
    xp,xm=(p[k] if k in p else p[f'{k}_{t}'] for k in ['plus','minus']);r=float(p['r'])
    pref=perf_counter(),process_time()
    record=dict(core=t,side=side,file=file,noise_radius=r,core_cells=int(sum(group==0)),core_area_m2=sum(group==0)*(2/side)**2,outputs={})
    outside=group!=0;delta=xp[outside]-xm[outside]
    record['exterior_budget']=dict(L1_m2=float(sum(abs(delta))*(2/side)**2),L2_m=float(np.linalg.norm(delta)*2/side),Linf=float(max(abs(delta))))
    record['roi_exactly_flat']=bool(np.all(xp[group==0]==t) and np.all(xm[group==0]==0))
    assert record['roi_exactly_flat'] and min(*xp,*xm)>=0 and max(*xp,*xm)<=4.5
    for band in bands:
        w=Circuit(side,band,ix);n=side*side
        v=VerifiedWave(dict(spectrum=w.symbol,grid_shape=[side,side],S=w.c,E0=w.b,low_rank_left=np.zeros((n,1)),low_rank_right=np.zeros((1,n))))
        s1,e1,r1=stored_response(w,v,xp);s0,e0,r0=stored_response(w,v,xm)
        empty=np.eye(25)-2*R*np.linalg.inv(w.z0+R*np.eye(25))
        z1=at_load(s1,1);z0=at_load(s0,1);air=at_load(empty,1)
        born=[];rels=[]
        mref=np.linalg.inv(w.z0[:24,:24]+R*np.eye(24))
        for x,obs in [(xp,z1),(xm,z0)]:
            lin=(-2j*R*mref@(w.c[:24]@(x[:,None]*w.b[:,:24]))@mref)[8:16,:8]
            actual=(obs-air)[8:16,:8]
            born.append(float(np.linalg.norm(actual-lin)/np.linalg.norm(actual)))
            rels.append(float(r/(np.linalg.norm(actual)/np.sqrt(8))))
        out=sweep(s1,e1,s0,e0,r,steps)
        out.update(born_A_relative_diagnostic=born,noise_as_fraction_of_each_scattered_A=rels,plus=r1,minus=r0)
        record['outputs'][str(band)]=out
        np.savez_compressed(HERE/f'verified_responses_{side}_{t}_{Path(file).stem}_{band}.npz',s1=s1,s0=s0,e1=e1,e0=e0,probe_incident=w.b[:,-1],probe_sensor=w.c[-1],Z_empty=w.z0)
    record['cost']=dict(wall_s=perf_counter()-pref[0],cpu_s=process_time()-pref[1],bare_lu=2*len(bands),port_solutions=2*len(bands))
    fn=HERE/f'cover_{side}_{t}_{Path(file).stem}.json'
    fn.write_text(json.dumps(record,indent=2,default=lambda x: x.item())+'\n')
    print(side,t,file,{band:(v['any_uniform_eight_cut_upper_over_2r'],v['arbitrary_wave_upper_over_2r'],v['concrete_wave_lower_over_2r']) for band,v in record['outputs'].items()},flush=True)


if __name__=='__main__':
    from pathlib import Path
    arg=argparse.ArgumentParser();arg.add_argument('side',type=int);arg.add_argument('t',type=float);arg.add_argument('file');arg.add_argument('--steps',type=int,default=1024);arg.add_argument('--bands',type=int,nargs='+',default=[400,460]);a=arg.parse_args()
    with threadpool_limits(1):run(a.side,a.t,a.file,a.steps,a.bands)

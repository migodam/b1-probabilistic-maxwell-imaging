"""Post-study trial of a one-cell optical multiplier on two public material boxes.
Reads only the stored prior, model and old readings. Computes the first old-data
saturated indefinite form of the prescribed family for each class and its future
support. Never changes the previous acquisition or reads the future fields.
"""
import json
import sys
from pathlib import Path
from time import process_time, perf_counter
import numpy as np
from scipy.linalg import eigh
from flint import acb, arb, ctx
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
C=HERE.parent
sys.path.insert(0,str(C))
from model import OLD, ANGLES, geometry, noise_radius, radius_squared, read_json, energy, eye, diag, star, ball, minor, hjoin, vjoin, inertia, mid, norm
from design import Border,envelope,source_pair,split_query,certify_pair

def herm(a):return (a+a.conj().T)/2

def parts(g,l,u,z,mu,cell,ratio,exact):
    n,m=g.b.shape
    p=len(g.s)
    if exact:
        j=hjoin(hjoin(eye(n),ball(np.zeros((n,m)))),ball(np.zeros((n,p))))
        e=hjoin(hjoin(g.ex_g,g.ex_b),ball(np.zeros((n,p))))
        r=hjoin(hjoin(g.ex_s,ball(np.zeros((p,m)))),-eye(p))
        lo=list(map(lambda x:arb(float(x)),l)); hi=list(map(lambda x:arb(float(x)),u))
        c=diag([acb(-mu*(x+y)/2,-mu*ratio/2 if i==cell else 0) for i,(x,y) in enumerate(zip(lo,hi))])
        a=star(j)*c*e
        form=star(r)*r+mu*star(j)*j+a+star(a)+star(e)*diag([mu*x*y for x,y in zip(lo,hi)])*e
        ids,ports=range(n),range(n,n+m+p)
        k=minor(form,ids,ids);f=minor(form,ids,ports)
        return k,minor(form,ports,ports)-star(f)*k.solve(f)
    j=np.column_stack((np.eye(n),np.zeros((n,m+p)))); e=np.column_stack((g.g,g.b,np.zeros((n,p)))); r=np.column_stack((g.s,np.zeros((p,m)),-np.eye(p)))
    v=-(l+u)*mu/2+np.zeros(n,complex);v[cell]-=.5j*mu*ratio
    a=j.conj().T@(v[:,None]*e)
    form=herm(r.conj().T@r+mu*j.conj().T@j+a+a.conj().T+e.conj().T@((mu*l*u)[:,None]*e))
    k=form[:n,:n]; f=form[:n,n:]
    return k,herm(form[n:,n:]-f.conj().T@np.linalg.solve(k,f))

def find(g,l,u,z,nr):
    zn=mid(z); m=g.b.shape[1]; tried=0; stat={}
    for mu in (128,8,1):
        for cell in range(len(l)):
            for ratio in (2,4,8,16,32,64):
                tried+=1
                k,q=parts(g,l,u,z,mu,cell,ratio,False)
                nu=int(np.count_nonzero(eigh(k,eigvals_only=True)<0))
                if not nu:continue
                po=herm(zn.conj().T@q@zn-float(radius_squared(len(OLD),nr))*np.eye(len(OLD)))
                ko=int(np.count_nonzero(eigh(po,eigvals_only=True)>0))
                if ko!=nu:continue
                b=q@zn;qc=herm(q-b@np.linalg.solve(po,b.conj().T))
                if np.min(eigh(qc[m:,m:],eigvals_only=True))<=0:continue
                K,Q=parts(g,l,u,z,mu,cell,ratio,True)
                sK,_,_=inertia(K)
                if not sK['full_count_verified'] or sK['negative']==0:continue
                E,rec=envelope(Border(mu,sK['negative'],Q,{}),g,z,nr)
                if E is not None:
                    return E,dict(trials=tried,mu=mu,pixel=cell,optical_ratio=ratio,K=sK,old=rec,**stat)
    return None,dict(trials=tried,strict_bounded_indefinite_form=False)

def common_clean_reading(first,second,index):
    # A point on the line of centres can prove overlap of the two outer sets.
    p=first.center.nrows()
    delta=minor(second.center,range(p),[index])-minor(first.center,range(p),[index])
    costs=[(star(delta)*e.receive*delta)[0,0].real for e in (first,second)]
    radii=[e.gram[index,index].real for e in (first,second)]
    if any(v.lower()<=0 for v in costs+radii): return dict(confirmed=False)
    reach=[float((r/c).sqrt().mid()) for r,c in zip(radii,costs)]
    a=arb(float(reach[0]/sum(reach)))
    gap=[r-c*z**2 for r,c,z in zip(radii,costs,(a,1-a))]
    return dict(confirmed=all(v.lower()>0 for v in gap),centre_fraction=float(a),
                slack_lower=[float(np.nextafter(float(v.lower()),-np.inf)) for v in gap])

def main():
    ctx.prec=192; ctx.threads=1
    reg=read_json(C/'registration.json'); t=perf_counter();cput=process_time()
    output=dict(scope='post-study old-data only; original occupied compact TM classes; pre-existing positive-state control also separates each pair',proof='interval Hankel/Richmond discrete model; exact power, joint error majorant',rows=[])
    with np.load(C/'public_model.npz') as pub:
        for scene,info in reg['scenes'].items():
            g=geometry(info['layout']);nr=noise_radius(reg,scene);power,_=energy('added',1)
            with np.load(C/'old'/f'{scene}__added.npz') as f:z=source_pair(OLD,power.sqrt(),f['Y'])
            opts=[];rows=[]
            for l,u in zip(pub[scene+'__lo'],pub[scene+'__hi']):
                e,rec=find(g,l,u,z,nr);opts.append(e);rows.append(rec)
            edges=[]
            if all(e is not None for e in opts):
                for i in range(len(ANGLES)):
                    if i in OLD:continue
                    _,w=split_query(opts,(i,),arb(1),nr)
                    if w is not None:
                        c=certify_pair(opts,w,(i,),arb(1),nr)
                        edges.append(dict(beam=i,angle_degrees=ANGLES[i],gap=c['noisy_separation_gap_lower'],strict=c['strict_edge'],w_real=w.real.tolist(),w_imag=w.imag.tolist(),clean_overlap=common_clean_reading(*opts,i)))
            output['rows'].append(dict(scene=scene,classes=rows,edges=edges,strict_edges=sum(x['strict'] for x in edges)))
            output['wall_seconds']=perf_counter()-t;output['cpu_seconds']=process_time()-cput
            (HERE/'results.json').write_text(json.dumps(output,indent=2))
            print(scene,[r['trials'] for r in rows],output['rows'][-1]['strict_edges'],flush=True)
if __name__=='__main__':
    with threadpool_limits(limits=1):main()

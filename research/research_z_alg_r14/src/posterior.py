"""Sequential set inversion after a real simulated read. No truth in the solver.

The old and new residual numerators at distinct frequencies are concatenated
WITHOUT multiplying their denominators; the common complex gain is retained
when forming every directional error bound.
"""
import argparse
from collections import deque
import json
from time import perf_counter, process_time

from flint import arb
import numpy as np

from data14 import ACTIONS, Action, load_cases
from frequency_cover import group_cover
from material_controls import (GRID, HI, K, LO, MASS, RESULTS, ceilroot, complex_ball,
                              int_directions, material_ranges, quantize, restrict, wave_table)


class Joint:
    def __init__(self, case, action, new, numerator=False):
        g = json.dumps(case.geometry,sort_keys=True)
        tabs=[wave_table(g,delta) for delta in (0.,*action.deltas)]
        data=[[complex_ball(x) for x in block] for block in (case.old,*new)]
        residuals=[]; den=[]; numerators=[]
        for v in range(1<<K):
            rp=[];rq=[];rn=[]
            for tab,D in zip(tabs,data):
                q,N=tab[v]
                raw=[n-y*q for n,y in zip(N,D)]
                rp.extend(z.real for z in raw);rp.extend(z.imag for z in raw)
                rn.extend(z.real for z in N);rn.extend(z.imag for z in N)
                rq.extend((q.real,q.imag))
            residuals.append(rp);den.append(rq);numerators.append(rn)
        self.numerator=numerator
        self.nout=len(residuals[0]);self.s=len(case.old)
        pc,ep=quantize([a+b for a,b in zip(residuals,numerators)] if numerator else residuals)
        qc,eq=quantize(den)
        self.table=np.column_stack([pc,qc])
        gamma=arb(min(case.gain_radius,case.hard_gain_radius or 1.))+arb(2)**-54
        s=arb(MASS)*arb(2)**(eq-ep)/(1-gamma)
        def ceil(x):
            return int(x.upper().ceil().unique_fmpz())
        def near(x):
            v=int((x.mid()+arb(.5)).floor().unique_fmpz())
            assert abs(x-v)<1
            return v
        eps=arb(case.additive_radius)+arb(2)**-54
        r=arb(case.new_noise_radius(action))+arb(2)**-54
        self.noise=[ceil(s*eps)]+[ceil(s*r/arb(float(e)).sqrt()) for e in action.energies]
        if numerator:
            scale=arb(MASS)*arb(2)**(eq-ep)
            self.noise=[ceil(scale*eps)]+[ceil(scale*r/arb(float(e)).sqrt()) for e in action.energies]
            self.gamma=ceil(arb(MASS)*(arb(case.gain_radius)+arb(2)**-54))
            self.hard=ceil(arb(MASS)*gamma)
            return
        self.gains=[]
        for D in data:
            z=[x*s*gamma for x in D]
            self.gains.append((np.array([near(x.real) for x in z],dtype=object),
                               np.array([near(x.imag) for x in z],dtype=object)))

    def test(self, lo, hi):
        z=restrict(self.table,lo,hi)
        p=z[:,:self.nout]
        q=z[:,self.nout*(1+int(self.numerator)):].reshape((len(z),-1,2))
        lower=q.min(axis=0)-(1+K);upper=q.max(axis=0)+(1+K)
        mid=(lower+upper)//2
        size=np.maximum(mid-lower,upper-mid)
        maxima=np.maximum(abs(lower),abs(upper))
        weights=int_directions(p)
        proj=p@weights.T
        coeff=weights.reshape((len(weights),-1,2*self.s))
        l1=np.abs(coeff).sum(axis=2)
        l2=(coeff**2).sum(axis=2)
        if self.numerator:
            n=z[:,self.nout:2*self.nout]
            # imag of sum (wr-i wi) (Nr+i Ni): real multipliers [-wi,wr].
            wi=np.concatenate((-coeff[:,:,self.s:],coeff[:,:,:self.s]),axis=2).reshape(weights.shape)
            dx=np.abs(n@weights.T).max(axis=0)
            dy=np.abs(n@wi.T).max(axis=0)
            mag_error=(1+K)*l1.sum(axis=1)
            old_noise=None;future=np.zeros(len(weights),dtype=object)
            for k,xy in enumerate(maxima):
                qmax=ceilroot(sum(int(x)**2 for x in xy))
                val=(qmax*self.noise[k])**2*l2[:,k].astype(object)
                if k: future+=val
                else: old_noise=val
            radius=[]
            for a,b,m,t1,t2 in zip(dx.tolist(),dy.tolist(),mag_error.tolist(),old_noise,future):
                cal=ceilroot((a+m)**2+(b+m)**2)
                joint=ceilroot(t1+(self.gamma*cal)**2)
                hard=ceilroot(t1)+self.hard*cal
                radius.append((min(joint,hard)+ceilroot(t2)+MASS-1)//MASS+m)
            radii=np.array(radius,dtype=np.int64)
            if np.any(proj.min(axis=0)>radii) or np.any(proj.max(axis=0)<-radii):
                return None
            return material_ranges(lo,hi,proj,radii)
        real=np.zeros(len(weights),dtype=object);imag=real.copy()
        rema=real.copy();future2=real.copy();old_radius=None
        for k,(dr,di) in enumerate(self.gains):
            ar=coeff[:,k,:self.s].astype(object);ai=coeff[:,k,self.s:].astype(object)
            vr=ar@dr+ai@di;vi=ar@di-ai@dr
            qr,qi=map(int,mid[k]);xr,xi=map(int,size[k])
            real+=qr*vr-qi*vi;imag+=qr*vi+qi*vr
            # Each coefficient of the scaled calibration vector is enclosed by +/-1.
            # L1 controls the disk enclosing complex products and quantization.
            l=l1[:,k].astype(object)
            qmax=ceilroot(sum(int(x)**2 for x in maxima[k]))
            qrad=ceilroot(xr*xr+xi*xi)
            rema+=qrad*(np.abs(vr)+np.abs(vi)) + (2*qmax)*l
            noise=(self.noise[k]*qmax)**2*l2[:,k].astype(object)
            if k:
                future2+=noise
            else:
                old_radius=noise
        rtot=[]
        for a,b,c,d,e,l in zip(real,imag,rema,old_radius,future2,l1.sum(axis=1)):
            t=ceilroot(a*a+b*b)+c+ceilroot(d)+ceilroot(e)
            rtot.append((t+MASS-1)//MASS + (1+K)*int(l))
        radii=np.array(rtot,dtype=np.int64)
        if np.any(proj.min(axis=0)>radii) or np.any(proj.max(axis=0)<-radii):
            return None
        return material_ranges(lo,hi,proj,radii)


def work(obj, seeds, target, nodes):
    queue=deque((x[0],x[1]) for x in seeds)
    leaves=[];seen=discarded=contracted=0
    while queue and seen<nodes:
        lo,hi=queue.popleft()
        seen+=1
        bound=obj.test(lo,hi)
        if bound is None:
            discarded+=1
            continue
        l,h=bound
        size=h-l
        if max(size)*(HI-LO)/GRID<=target:
            leaves.append(bound)
        elif size.sum()<.78*(hi-lo).sum():
            queue.append(bound);contracted+=1
        else:
            ax=np.argmax(size);mid=(l[ax]+h[ax])//2
            left=h.copy();right=l.copy()
            left[ax]=right[ax]=mid
            queue.extend(((l,left),(right,h)))
    unfinished=len(queue)
    leaves.extend(queue)
    box=np.array(leaves,dtype=np.int64).reshape(-1,2,K)
    return box,dict(visited=seen,discarded=discarded,contracted=contracted,unfinished=unfinished)


def hull_diam(box):
    if not len(box):
        return 0.
    grouped=group_cover(box, min(96,len(box)))
    i,j=np.triu_indices(len(grouped))
    d=np.maximum(grouped[i,1]-grouped[j,0],grouped[j,1]-grouped[i,0])
    return float(np.nextafter((HI-LO)/GRID*np.sqrt((d*d).sum(axis=1).max()/K)*(1.+1.e-14),np.inf))


def run(case,reads,maximum=4500, numerator=False, seed=None):
    if seed is None:
        seed=np.load(RESULTS/"material_covers_seeded"/(case.case_id+".npz"))["boxes"]
    first=Action("first_half",(-.02,),(.5,))
    both=ACTIONS["split_0.02"]
    summary={};saved={}
    for stage,action,blocks in (("first_half",first,reads[:1]),("both_halves",both,reads)):
        wall,cpu=perf_counter(),process_time()
        joint=Joint(case,action,blocks,numerator=numerator)
        start=group_cover(seed,96)
        out,stats=work(joint,start,.085,maximum)
        diameter=hull_diam(out)
        summary[stage]=dict(**stats,retained=len(out),D_upper=diameter,
                reached_half=bool(len(out) and diameter<=.5),new_energy=action.energy,
                cpu_s=process_time()-cpu,wall_s=perf_counter()-wall)
        saved[stage]=out
    return summary,saved


def main():
    p=argparse.ArgumentParser();p.add_argument("--case",action="append")
    p.add_argument("--nodes",type=int,default=4500);p.add_argument("--tag",default="")
    p.add_argument("--numerator",action="store_true")
    a=p.parse_args();folder=RESULTS/("posterior"+a.tag);folder.mkdir(exist_ok=True)
    with np.load(RESULTS/"measure_new_reads.npz") as src:
        for case in load_cases():
            if a.case and case.case_id not in a.case:
                continue
            obs=src[case.case_id+"/split_0.02"]
            info,arrays=run(case,obs,a.nodes,numerator=a.numerator)
            (folder/(case.case_id+".json")).write_text(json.dumps(info,indent=2)+"\n")
            np.savez_compressed(folder/(case.case_id+".npz"),**arrays)
            print(case.case_id,[(k, round(v["D_upper"],4),v["visited"],round(v["wall_s"],2)) for k,v in info.items()],flush=True)
if __name__=="__main__":
    main()
